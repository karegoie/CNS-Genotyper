#!/bin/python3

import datetime
import os
import click
from Bio import SeqIO
import gzip

# type hinting for lower version of python
from typing import List

from src.reference import Reference
from src.line_set import Line_Set
from src.indel_counter_for_genotype import InDel_Counter_for_Genotype
from src.log_writer import write_main_log, write_main_html_log, write_sub_log, write_main_csv_log, get_main_log_name, write_raw_data_log
import src.globals as glv

DATA_ADDRESS = "./data/"
REF_ADDRESS = "./ref/"
GUIDE_RNA_SET_ADDRESS = "./ref/guide_RNA_set.fasta"
REF_SET_ADDRESS = "./ref/reference_seq_set.fasta"

# with spans, 3 second for 1000 lines: 20000 for a minute, 1,000,000 ~ 1,500,000 for an hour
# > total 800 nt of ref, 150 nt for a line: 40,000,000 for a second.


def get_best_line_set(read: SeqIO.SeqRecord, reference_list: list):
    if len(reference_list) == 0:
        return None
    best_line_set = Line_Set(read_raw=read, reference=reference_list[0])

    if len(reference_list) > 1:
        for reference in reference_list[1:]:
            test_line_set = Line_Set(read_raw=read, reference=reference)

            if test_line_set.score > best_line_set.score:
                best_line_set = test_line_set
    return best_line_set


def get_reference_list_from_file():
    # Get reference and guide RNA sequences,
    # and match them as a 'Reference' class

    reference_list = []

    ref_raw_iter = SeqIO.parse(REF_SET_ADDRESS, "fasta")
    ref_raw_list = []
    for item in ref_raw_iter:
        ref_raw_list.append(item)

    g_rna_seq_iter = SeqIO.parse(GUIDE_RNA_SET_ADDRESS, "fasta")
    g_rna_raw_list = []
    for item in g_rna_seq_iter:
        g_rna_raw_list.append(item)

    for i, ref_raw in enumerate(ref_raw_list):
        if i >= len(g_rna_raw_list):
            i = len(g_rna_raw_list) - 1
        reference = Reference(ref_raw=ref_raw, guide_rna_raw=g_rna_raw_list[i])
        reference_list.append(reference)

    return reference_list


def get_file_data_file_list():
    data_file_list = [file_name for file_name in os.listdir(DATA_ADDRESS)
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_name))
                    and file_name[-6:] in ('.fastq', 'stq.gz')]

    for ignore_text in glv.READ_IGNORE:
        if len(ignore_text) > 0:
            data_file_list = [address for address in data_file_list if address.find(ignore_text) < 0]

    data_file_list.sort(key=lambda f: int(''.join(filter(str.isdigit, f)) + '0'))
    return data_file_list


def get_total_number_of_reads(data_file_list: List[str]):
    total_reads_count = 0
    reads_count_list = []
    for i, file_name in enumerate(data_file_list):
        reads_count = 0
        print(f"\r({i + 1}/{len(data_file_list)}) reading {file_name}", end="")

        if file_name[-5:] == str("file.fastq.gz")[-5:]:
            read_raw_iter = SeqIO.parse(gzip.open(str(os.path.join(DATA_ADDRESS, file_name)), "rt"), "fastq")
        else:
            read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")

        for _ in read_raw_iter:
            reads_count += 1
            total_reads_count += 1
        reads_count_list.append(reads_count)
    print(f"\rTotal reads :{total_reads_count} for {len(data_file_list)} files                     ")
    print()
    return total_reads_count, reads_count_list


def key_for_sorting_err(line_set: Line_Set):
    if line_set.indel_type == 'err':
        return 1
        # return len(line_set)
    return 0


def test_all_input_files():

    is_data_exist = False
    is_guide_rna_exist = False
    is_reference_exist = False

    if not os.path.exists(DATA_ADDRESS):
        os.makedirs(DATA_ADDRESS)
    if not os.path.exists(REF_ADDRESS):
        os.makedirs(REF_ADDRESS)

    if len(get_file_data_file_list()) > 0:
        is_data_exist = True

    if os.path.exists(GUIDE_RNA_SET_ADDRESS):
        try:
            start = True
            seq_iter = SeqIO.parse(GUIDE_RNA_SET_ADDRESS, 'fasta')
            for seq in seq_iter:
                for a in str(seq.seq):
                    if a not in 'ATGCatgc':
                        start = False
            is_guide_rna_exist = start
        except (TypeError, ValueError, PermissionError):
            pass
    else:
        with open(GUIDE_RNA_SET_ADDRESS, 'w') as file:
            file.write(">Guide_RNA_SPACER\n"
                       "ATTATAGGAAGAAAGGGGAA # insert the spacer sequence of your guide RNA without PAM sequence here\n")

    if os.path.exists(REF_SET_ADDRESS):
        try:
            start = True
            seq_iter = SeqIO.parse(REF_SET_ADDRESS, 'fasta')
            for seq in seq_iter:
                for a in str(seq.seq):
                    if a not in 'ATGCatgc':
                        start = False
            is_reference_exist = start
        except (TypeError, ValueError, PermissionError):
            pass
    else:
        with open(REF_SET_ADDRESS, 'w') as file:
            file.write(">Rererence_Sequence_around_sequencing_part\n"
                       "ATTATAGGAAGAAAGGGGAA # insert the reference sequence here "
                       "(length > 300 nt, margin > 50 nt is recommended)\n")

    if not is_data_exist:
        print("ERROR: No Data(NGS result files: fastq, fastq.gz) found in the ./data folder")
    if not is_reference_exist:
        print("ERROR: Something is wrong with ./ref/reference_seq_set.fasta file!")
    if not is_guide_rna_exist:
        print("ERROR: Something is wrong with ./ref/guide_RNA_set.fasta file!")
    if is_guide_rna_exist and is_data_exist and is_reference_exist:
        return True
    print("Please check the files and try again!")
    print("Terminating...")
    return False


@click.command()
@click.option('-x', '--read_ignore', default=['R2', 'Undetermined'], multiple=True,
              help=glv.EXPLANATION_MAP['read_ignore'])
#
@click.option('-e', '--err_ratio_max', default=0.03,
              help=glv.EXPLANATION_MAP['err_ratio_max'])
@click.option('-p', '--err_padding_for_seq', default=1,
              help=glv.EXPLANATION_MAP['err_padding_for_seq'])
@click.option('-d', '--cut_pos_from_pam', default=-3,
              help=glv.EXPLANATION_MAP['cut_pos_from_pam'])
@click.option('-r', '--cut_pos_radius', default=5,
              help=glv.EXPLANATION_MAP['cut_pos_radius'])
@click.option('-s', '--phred_meaningful_score_min', default=30,
              help=glv.EXPLANATION_MAP['phred_meaningful_score_min'])
#
@click.option('--pam_distance_max', default=5,
              help=glv.EXPLANATION_MAP['pam_distance_max'])
@click.option('--score_match', default=2,
              help=glv.EXPLANATION_MAP['score_match'])
@click.option('--score_mismatch', default=-1,
              help=glv.EXPLANATION_MAP['score_mismatch'])
@click.option('--score_gap_open', default=-30,
              help=glv.EXPLANATION_MAP['score_gap_open'])
@click.option('--score_gap_extend', default=-4,
              help=glv.EXPLANATION_MAP['score_gap_extend'])
#
@click.option('-t', '--task_title', default="task " + str(datetime.datetime.now())[5:-10],
              help=glv.EXPLANATION_MAP['task_title'])
@click.option('-o', '--open_xlsx_auto', default=False, is_flag=True,
              help=glv.EXPLANATION_MAP['open_xlsx_auto'])
@click.option('--debug', default=False, is_flag=True,
              help=glv.EXPLANATION_MAP['debug'])
def main(read_ignore, err_ratio_max, err_padding_for_seq, cut_pos_from_pam, cut_pos_radius,
         phred_meaningful_score_min, pam_distance_max,
         score_match, score_mismatch, score_gap_open, score_gap_extend,
         task_title, open_xlsx_auto, debug):
    # set global variables
    glv.READ_IGNORE = read_ignore

    glv.ERR_RATIO_MAX = err_ratio_max
    glv.ERR_PADDING_FOR_SEQ = err_padding_for_seq
    glv.CUT_POS_FROM_PAM = cut_pos_from_pam
    glv.CUT_POS_RADIUS = cut_pos_radius

    glv.PHRED_MEANINGFUL_MIN = phred_meaningful_score_min
    glv.PAM_DISTANCE_MAX = pam_distance_max

    glv.MAT = score_match
    glv.MIS = score_mismatch
    glv.GAP_OPEN = score_gap_open
    glv.GAP_EXTEND = score_gap_extend

    glv.TASK_TITLE = task_title
    glv.OPEN_XLSX_AUTO = open_xlsx_auto
    glv.DEBUG = debug

    is_good_to_go = test_all_input_files()
    if not is_good_to_go:
        return -1

    # Get sorted file address list from a folder, list[str]
    data_file_list = get_file_data_file_list()
    print(f"File list: files with '.fastq.gz' or '.fastq' in {DATA_ADDRESS} only")
    print(f"File list: ignoring keyword {glv.READ_IGNORE} in the name")
    print(f"File list: {data_file_list}")
    print()

    # get a list[Reference]
    reference_list = get_reference_list_from_file()

    # for total result, making list[list[InDel_Counter]]
    all_indel_counter_list_list = []

    debug_data = {}

    # # # for counting expected time left,
    # # # count total number of reads,
    # # # count total number of finished number of reads,
    # # # and check the time of initiation
    '''This function will make a text print: opening large file takes some time'''
    total_reads_count, reads_count_list = get_total_number_of_reads(data_file_list=data_file_list)
    finish_reads_count = 0
    start_time = datetime.datetime.now()
    start_time_for_file_before = datetime.datetime.now()
    start_time_for_file = datetime.datetime.now()

    # to show that the program is running by someone else in the folder:
    # this will make a no_name file of 'using' the folder.
    # file = open(os.path.join(DATA_ADDRESS, ".program_is_running_here"), 'w')
    # file.close()

    for file_no, file_name in enumerate(data_file_list):
        '''
        For each file, This happens:
            get a list of 'reads' from file (NGS-fastq only),
            
            for each read:
                for each reference:
                    (make line_set from 1 read and 1 reference)
                    make possible aligns with each reference sequences,
                    align the phred score line,
                    get the position of guide_RNA and PAM sequence,
                    get the indel type by checking around the PAM starting point,
                    get the score(mismatch ratio without main indel) and check if it is error seq or not.
                
                pick the best aligned line_set by score, 
                and add it to the list.
                
            for each line_set in list:
                for each indel_counter in list:
                    if line_set.ref == indel_counter.ref:
                        indel_counter.count(line_set) < count the each indel type
        
        '''
        # # # for counting expected time left
        start_time_for_file_before = start_time_for_file
        start_time_for_file = datetime.datetime.now()

        # build list[InDel_Counter_For_Ref]
        indel_counter_list = []
        for reference in reference_list:
            indel_counter = InDel_Counter_for_Genotype(reference=reference)
            indel_counter.set_file_name(file_name=file_name)
            indel_counter_list.append(indel_counter)

        if file_name[-5:] == str("file.fastq.gz")[-5:]:
            read_raw_iter = SeqIO.parse(gzip.open(str(os.path.join(DATA_ADDRESS, file_name)), "rt"), "fastq")
        elif file_name[-5:] == str(".fastq")[-5:]:
            read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
        else:
            print(f"({file_no + 1}/{len(data_file_list)}) {file_name} is not readable: is it .fastq or .fastq.gz ?")
            continue

        # # # move the file to the 'used' folder, as quick as possible
        # os.replace(os.path.join(DATA_ADDRESS, file_name), os.path.join(USED_DATA_ADDRESS, file_name))
        # not working...

        # build list[Bio.SeqRecord]
        read_raw_list = [read_raw for read_raw in read_raw_iter]

        # build list[Line_Set]
        line_set_list = []
        for i, read_raw in enumerate(read_raw_list):

            # # # for showing expected time left
            finish_reads_count += 1

            if (i % 100) == 0:
                now_time = datetime.datetime.now()
                delta_time = now_time - start_time
                delta_count = finish_reads_count
                if file_no > 0:
                    delta_time = now_time - start_time_for_file_before
                    delta_count = i + reads_count_list[file_no-1]

                print(f"\r({file_no + 1}/{len(data_file_list)}) "
                      f"for {file_name}: {((i+1)/len(read_raw_list)):.3f} / "
                      f"remaining: {(delta_time/delta_count)*(total_reads_count-finish_reads_count)} "
                      f"(for this file: {(delta_time/delta_count)*(len(read_raw_list)-(i+1))}) "
                      f"(length: {len(read_raw_list)})                              ", end="")

            best_line_set = get_best_line_set(read_raw, reference_list)
            best_line_set.set_file_name(file_name=file_name)

            line_set_list.append(best_line_set)

        # # # for showing expected time left / while log writing
        print(f"\r({file_no + 1}/{len(data_file_list)}) for {file_name}: Complete / "
              f"Writing log files (length: {len(line_set_list)})                              ", end="")

        # count the number of each indel type,
        # also setting the best line_set for each indel type happens here
        for line_set in line_set_list:
            for indel_counter in indel_counter_list:
                if indel_counter.ref_name == line_set.ref_name:
                    indel_counter.count(line_set)

        # # save the indel type count to each line set
        # # just for the log file
        for line_set in line_set_list:
            for indel_counter in indel_counter_list:
                if indel_counter.ref_name == line_set.ref_name:
                    line_set.set_indel_same_type_count(indel_counter.count_map)

        #
        # add the file result to the total result
        all_indel_counter_list_list.append(indel_counter_list)

        # # Sorting Line Set List!
        # # err for the last / biggest indel type first / higher score / higher phred score / longer one first
        # # just for the log file
        indel_counter_map = {}
        for indel_counter in indel_counter_list:
            indel_counter_map[indel_counter.ref_name] = indel_counter
        line_set_list.sort(key=lambda l: l.phred_score, reverse=True)
        line_set_list.sort(key=lambda l: l.score, reverse=True)
        line_set_list.sort(key=lambda l: indel_counter_map[l.ref_name].count_map[l.indel_type], reverse=True)
        line_set_list.sort(key=lambda l: key_for_sorting_err(l))

        # Writing sub log
        for indel_counter in indel_counter_list:
            write_sub_log(line_set_list=[l for l in line_set_list if l.ref_name == indel_counter.ref_name],
                          indel_counter=indel_counter, file_name=file_name)

        # # for debug data
        if glv.DEBUG:
            pos_phred_score = [0] * 25
            pos_error_count = [0] * 25

            for line_set in line_set_list:
                if len(line_set) < 20:
                    print(line_set)
                    continue
                for i in range(-10, 10):
                    pos_phred_score[i] += (ord(line_set.phred_line[i]) - glv.PHRED_ENCODING)

                    if line_set.match_line[i] != '|':
                        pos_error_count[i] += 1

            debug_data[file_name] = {
                'length': len(line_set_list),
                'pos_phred_score': pos_phred_score,
                'pos_error_count': pos_error_count
            }

        # # for showing time used
        end_time_for_file = datetime.datetime.now()
        print(f"\r({file_no + 1}/{len(data_file_list)}) for {file_name}: Complete / Log written / "
              f"{end_time_for_file - start_time} ({end_time_for_file - start_time_for_file} for this file) is passed "
              f"(length: {len(line_set_list)})       ")

        # end.

    # Writing total log
    write_main_log(indel_counter_list_list=all_indel_counter_list_list, total_length=total_reads_count)
    write_main_html_log(indel_counter_list_list=all_indel_counter_list_list, total_length=total_reads_count)
    write_main_csv_log(indel_counter_list_list=all_indel_counter_list_list, ref_set_list=reference_list)
    write_raw_data_log(indel_counter_list_list=all_indel_counter_list_list, debug_data=debug_data)

    # # for showing time used
    print(f"Work Completed! (total time: {datetime.datetime.now() - start_time})")

    # # delete the file of 'program is running here'
    # os.remove(os.path.join(DATA_ADDRESS, ".program_is_running_here"))

    #
    if glv.OPEN_XLSX_AUTO:
        os.system(f"start EXCEL.EXE {get_main_log_name('xlsx')}")


if __name__ == '__main__':
    print("CNS-Genotyper ver. "+glv.VERSION)
    print()
    main()

    # this is how 'click' works...
    print("Things never works below here... << error message")




