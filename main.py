#!/bin/python3

import datetime
import os
import click
from Bio import SeqIO
import gzip
import configparser
import gc

# type hinting for lower version of python
from typing import List

from src.reference import Reference
from src.aligned import Aligned_Key, Aligned_Read, Aligned
from src.genotyper_for_reference import Genotyper_For_Reference
from src.log_writer import write_main_log, write_main_html_log, write_sub_log, write_main_csv_log, get_main_log_name, \
    write_raw_data_log
import src.globals as glv
# TODO: use asnycio module for faster processing

DATA_ADDRESS = "./data/"
REF_ADDRESS = "./ref/"
GUIDE_RNA_SET_ADDRESS = "./ref/guide_RNA_set.txt"
REF_SET_ADDRESS = "./ref/reference_seq_set.txt"
CONFIG_ADDRESS = "./config.txt"


# with spans, 3 second for 1000 lines: 20000 for a minute, 1,000,000 ~ 1,500,000 for an hour
# > total 800 nt of ref, 150 nt for a line: 40,000,000 for a second.

# Now it works like, 2,000,000 reads for only 10 minutes.


def get_best_aligned_key(seq_key: str, reference_list: list):
    if len(reference_list) == 0:
        return None
    best_aligned_key = None

    for reference in reference_list:
        test_aligned_key = Aligned_Key(key_seq=seq_key, reference=reference, strand='+')
        if best_aligned_key is None or best_aligned_key.score < test_aligned_key.score:
            best_aligned_key = test_aligned_key

        test_aligned_key = Aligned_Key(key_seq=glv.get_opposite_strand(seq_key), reference=reference, strand='-')
        if best_aligned_key is None or best_aligned_key.score < test_aligned_key.score:
            best_aligned_key = test_aligned_key

    return best_aligned_key


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
    print(f"File list: files with '.fastq.gz' or '.fastq' in {DATA_ADDRESS} only\n"
          f"File list: ignoring keyword {glv.READ_IGNORE} in the name\n"
          f"File list: {data_file_list}\n")
    return data_file_list


# def get_total_number_of_reads(data_file_list: List[str]):
#     total_reads_count = 0
#     reads_count_list = []
#     for i, file_name in enumerate(data_file_list):
#         reads_count = 0
#         print(f"\r({i + 1}/{len(data_file_list)}) reading {file_name}", end="")
#
#         if file_name[-5:] == str("file.fastq.gz")[-5:]:
#             read_raw_iter = SeqIO.parse(gzip.open(str(os.path.join(DATA_ADDRESS, file_name)), "rt"), "fastq")
#         else:
#             read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
#
#         for _ in read_raw_iter:
#             reads_count += 1
#             total_reads_count += 1
#         reads_count_list.append(reads_count)
#     print(f"\rTotal reads :{total_reads_count} for {len(data_file_list)} files                     ")
#     print()
#     return total_reads_count, reads_count_list


def key_for_sorting_err(aligned_read: Aligned_Read):
    if aligned_read.indel.indel_type == 'err':
        return 1
        # return len(aligned_read)
    return 0


def test_all_input_files():
    is_data_exist = False
    is_guide_rna_exist = False
    is_reference_exist = False

    if not os.path.exists(DATA_ADDRESS):
        os.makedirs(DATA_ADDRESS)
    if not os.path.exists(REF_ADDRESS):
        os.makedirs(REF_ADDRESS)

    data_file_list = get_file_data_file_list()
    if len(data_file_list) > 0:
        is_data_exist = True

    if os.path.exists(GUIDE_RNA_SET_ADDRESS):
        try:
            start = True
            seq_iter = SeqIO.parse(GUIDE_RNA_SET_ADDRESS, 'fasta')
            for seq in seq_iter:
                for a in str(seq.seq):
                    if a not in 'ATGCatgcN':
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
        print("ERROR: Something is wrong with ./ref/reference_seq_set.txt file!")
    if not is_guide_rna_exist:
        print("ERROR: Something is wrong with ./ref/guide_RNA_set.txt file!")
    if is_guide_rna_exist and is_data_exist and is_reference_exist:
        return True
    print("Please check the files and try again!")
    print("Terminating...")
    return False


def get_seq_key(seq: str):
    return seq[glv.ERR_PADDING_FOR_SEQ:-glv.ERR_PADDING_FOR_SEQ]


def get_key_list_of_all_seq(data_file_list: list):
    hashmap_seq_key = dict()
    total_reads = 0

    for file_no, file_name in enumerate(data_file_list):
        print(f"\r({file_no + 1}/{len(data_file_list)}) reading {file_name}           ", end="")

        if file_name[-5:] == str("file.fastq.gz")[-5:]:
            read_raw_iter = SeqIO.parse(gzip.open(str(os.path.join(DATA_ADDRESS, file_name)), "rt"), "fastq")
        elif file_name[-5:] == str(".fastq")[-5:]:
            read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
        else:
            print(f"({file_no + 1}/{len(data_file_list)}) {file_name} is not readable: is it .fastq or .fastq.gz ?")
            continue

        for read_raw in read_raw_iter:
            total_reads += 1
            seq = str(read_raw.seq)
            if len(seq) < glv.ERR_PADDING_FOR_SEQ * 2 + 10:
                print("\n", seq, "the length looks very wrong...")
            seq_key = get_seq_key(seq)
            hashmap_seq_key[seq_key] = None

    key_list = list(hashmap_seq_key.keys())

    print(f"Total reads: {total_reads}\n"
          f"Total unique reads: {len(key_list)}")
    return key_list, total_reads


def get_aligned_hashmap(seq_key_list, reference_list):
    print("\n")
    hashmap = dict()
    start_time = datetime.datetime.now()
    for i, seq_key in enumerate(seq_key_list):
        hashmap[seq_key] = get_best_aligned_key(seq_key, reference_list)
        # print(seq_key, hashmap[seq_key])
        if i % 100 == 0:
            print(f"\rUnique Sequencing aligning: {(i + 1) / len(seq_key_list):.03f} / "
                  f"remaining: {(datetime.datetime.now() - start_time) * (len(seq_key_list) - (i + 1)) / (i + 1)} "
                  f"({datetime.datetime.now() - start_time} is passed) "
                  f"(length: {len(seq_key_list)})                        ",
                  end="")

    print()
    return hashmap


def set_global_variables(read_ignore, err_ratio_max, err_padding_for_seq, cut_pos_from_pam, cut_pos_radius,
                         phred_meaningful_score_min, pam_distance_max, score_match, score_mismatch, score_gap_open,
                         score_gap_extend, task_title, open_xlsx_auto, debug):
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
    return


def get_genotyper_list_for_reference_list(reference_list: list, file_name: str, aligned_read_list: list):
    genotyper_list = []
    for reference in reference_list:
        genotyper = Genotyper_For_Reference(reference=reference, file_name=file_name)
        genotyper_list.append(genotyper)

    # count the number of each indel type,
    # also setting the best aligned_read for each indel type happens here
    for aligned_read in aligned_read_list:
        for genotyper in genotyper_list:
            if genotyper.ref_name == aligned_read.ref_name:
                genotyper.count(aligned_read)
    return genotyper_list


def get_aligned_read_list_for_file(file_no, file_name, aligned_hashmap, total_file_length,
                                   total_reads_count, finish_reads_count, start_time_3):
    if file_name[-5:] == str("file.fastq.gz")[-5:]:
        read_raw_iter = SeqIO.parse(gzip.open(str(os.path.join(DATA_ADDRESS, file_name)), "rt"), "fastq")
    elif file_name[-5:] == str(".fastq")[-5:]:
        read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
    else:
        return None

    aligned_read_list = []
    read_raw_list = [read_raw for read_raw in read_raw_iter]
    for i, read_raw in enumerate(read_raw_list):

        seq = str(read_raw.seq)
        seq_key = get_seq_key(seq)
        aligned_key = aligned_hashmap[seq_key]
        if aligned_key.strand == '-':
            aligned_read = Aligned_Read(aligned_key, read_raw.reverse_complement(), file_name)
        else:
            aligned_read = Aligned_Read(aligned_key, read_raw, file_name)

        # print(aligned_read)
        aligned_read_list.append(aligned_read)

        # # # # for showing expected time left
        finish_reads_count += 1

        if (i % 1000) == 0:
            now_time = datetime.datetime.now()
            delta_time = now_time - start_time_3
            delta_count = finish_reads_count

            print(f"\r({file_no}/{total_file_length}) "
                  f"for {file_name}: {((i + 1) / len(read_raw_list)):.3f} / "
                  f"remaining: {(delta_time / delta_count) * (total_reads_count - finish_reads_count)} "
                  f"(for this file: {(delta_time / delta_count) * (len(read_raw_list) - (i + 1))}) "
                  f"(length: {len(read_raw_list)})                              ", end="")

    # # # for showing expected time left / while log writing
    print(f"\r({file_no}/{total_file_length}) for {file_name}: Complete / "
          f"Writing log files (length: {len(aligned_read_list)})                              ", end="")

    return aligned_read_list


def get_sorted_aligned_read_list(aligned_read_list, genotyper_list):
    genotyper_map = {}
    for genotyper in genotyper_list:
        genotyper_map[genotyper.ref_name] = genotyper
    aligned_read_list.sort(key=lambda l: l.score, reverse=True)
    aligned_read_list.sort(key=lambda l: genotyper_map[l.ref_name].count_map[l.indel.indel_type], reverse=True)
    aligned_read_list.sort(key=lambda l: key_for_sorting_err(l))

    return aligned_read_list


@click.command()
@click.option('-x', '--read_ignore', default=['R2', 'Undetermined'], multiple=True,
              help=glv.EXPLANATION_MAP['read_ignore'])
#
@click.option('-e', '--err_ratio_max', default=glv.ERR_RATIO_MAX,
              help=glv.EXPLANATION_MAP['err_ratio_max'])
@click.option('-p', '--err_padding_for_seq', default=glv.ERR_PADDING_FOR_SEQ,
              help=glv.EXPLANATION_MAP['err_padding_for_seq'])
@click.option('-d', '--cut_pos_from_pam', default=glv.CUT_POS_FROM_PAM,
              help=glv.EXPLANATION_MAP['cut_pos_from_pam'])
@click.option('-r', '--cut_pos_radius', default=glv.CUT_POS_RADIUS,
              help=glv.EXPLANATION_MAP['cut_pos_radius'])
@click.option('-s', '--phred_meaningful_score_min', default=glv.PHRED_MEANINGFUL_MIN,
              help=glv.EXPLANATION_MAP['phred_meaningful_score_min'])
#
@click.option('--pam_distance_max', default=glv.PAM_DISTANCE_MAX,
              help=glv.EXPLANATION_MAP['pam_distance_max'])
@click.option('--score_match', default=glv.MAT,
              help=glv.EXPLANATION_MAP['score_match'])
@click.option('--score_mismatch', default=glv.MIS,
              help=glv.EXPLANATION_MAP['score_mismatch'])
@click.option('--score_gap_open', default=glv.GAP_OPEN,
              help=glv.EXPLANATION_MAP['score_gap_open'])
@click.option('--score_gap_extend', default=glv.GAP_EXTEND,
              help=glv.EXPLANATION_MAP['score_gap_extend'])
#
@click.option('-t', '--task_title', default=glv.TASK_TITLE,
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
    set_global_variables(read_ignore, err_ratio_max, err_padding_for_seq, cut_pos_from_pam, cut_pos_radius,
                         phred_meaningful_score_min, pam_distance_max,
                         score_match, score_mismatch, score_gap_open, score_gap_extend,
                         task_title, open_xlsx_auto, debug)

    is_good_to_go = test_all_input_files()
    if not is_good_to_go:
        # input("Press any key to finish...")
        return -1

    data_file_list = get_file_data_file_list()
    reference_list = get_reference_list_from_file()

    start_time_0 = start_time_1 = datetime.datetime.now()
    seq_key_list, total_reads_count = get_key_list_of_all_seq(data_file_list)
    finish_reads_count = 0
    print(f"\n\tUnique sequence reading: {datetime.datetime.now() - start_time_1}")

    start_time_2 = datetime.datetime.now()
    aligned_hashmap = get_aligned_hashmap(seq_key_list, reference_list)
    print(f"\n\tUnique sequence aligning: {datetime.datetime.now() - start_time_2}")

    # for total result, making list[list[genotyper]]
    start_time_3 = datetime.datetime.now()
    all_genotyper_list_list = []
    debug_data = dict()

    for file_no, file_name in enumerate(data_file_list):
        start_time_for_file = datetime.datetime.now()

        # build list[Aligned_Read]
        aligned_read_list = get_aligned_read_list_for_file(file_no+1, file_name, aligned_hashmap, len(data_file_list),
                                                           total_reads_count, finish_reads_count, start_time_3)
        # build list[genotyper_For_Ref] with count_result
        genotyper_list = get_genotyper_list_for_reference_list(reference_list, file_name, aligned_read_list)

        # log order: err for the last / biggest indel type first / higher score / higher phred score / longer one first
        sorted_aligned_read_list = get_sorted_aligned_read_list(aligned_read_list, genotyper_list)

        # add the file result to the total result
        all_genotyper_list_list.append(genotyper_list)

        # Writing sub log
        for genotyper in genotyper_list:
            write_sub_log(aligned_read_list=[l for l in sorted_aligned_read_list if l.ref_name == genotyper.ref_name],
                          genotyper=genotyper, file_name=file_name)
        # # for showing time used
        finish_reads_count += len(sorted_aligned_read_list)

        end_time_for_file = datetime.datetime.now()
        print(f"\r({file_no + 1}/{len(data_file_list)}) for {file_name}: Complete / Log written / "
              f"{end_time_for_file - start_time_0} ({end_time_for_file - start_time_for_file} for this file) is passed "
              f"(length: {len(sorted_aligned_read_list)})       ")

        # end.
    print(f"\n\tAll read matching: {datetime.datetime.now() - start_time_3}")

    # Writing total log
    write_main_log(genotyper_list_list=all_genotyper_list_list, total_length=total_reads_count)
    write_main_html_log(genotyper_list_list=all_genotyper_list_list, total_length=total_reads_count)
    write_main_csv_log(genotyper_list_list=all_genotyper_list_list, ref_set_list=reference_list)
    write_raw_data_log(genotyper_list_list=all_genotyper_list_list, debug_data=debug_data)

    # # for showing time used
    print(f"Work Completed!\n\t(total time: {datetime.datetime.now() - start_time_0})")

    if glv.OPEN_XLSX_AUTO:
        os.system(f"start EXCEL.EXE {get_main_log_name('xlsx')}")
    # input("Press enter to finish the program...")


if __name__ == '__main__':
    print("CNS-Genotyper ver. " + glv.VERSION)
    print()
    main()

    # this is how 'click' works...
    print("Things never works below here... << error message")

## testing pull request
