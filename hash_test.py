
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

from main import get_file_data_file_list, get_reference_list_from_file, DATA_ADDRESS, get_total_number_of_reads


if __name__ == '__main__':

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

    # Hash-Map'uh
    read_seq_dict = {}

    for file_no, file_name in enumerate(data_file_list):

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

        # build list[Bio.SeqRecord]
        read_raw_list = [read_raw for read_raw in read_raw_iter]

        for read_raw in read_raw_list:
            seq = str(read_raw.seq)[glv.ERR_PADDING_FOR_SEQ:-glv.ERR_PADDING_FOR_SEQ]
            if seq not in read_seq_dict.keys():
                read_seq_dict[seq] = [read_raw]
            else:
                read_seq_dict[seq] += [read_raw]

        print(len(read_seq_dict.keys()), len(read_raw_list), round(len(read_seq_dict.keys()) / len(read_raw_list), 2))
        print(datetime.datetime.now() - start_time_for_file, len(read_raw_list))

    print(datetime.datetime.now() - start_time, total_reads_count)
