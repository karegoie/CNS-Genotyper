import datetime
import csv
import os
import re
from xlsxwriter.workbook import Workbook

# type hinting for lower version of python
from typing import List, Dict

from src.aligned import Aligned_Read
from src.reference import Reference
from src.genotyper_for_reference import Genotyper_For_Reference
import src.globals as glv

# Variables for uh... zero division
Z = 0.000000001


def get_main_log_name(extension: str, is_folder=False, is_sub=False):
    valid_task_title = "".join([c for c in glv.TASK_TITLE if c not in "\/:*?<>| -"])

    if not os.path.exists(f"./Result_{valid_task_title}/"):
        os.makedirs(f"./Result_{valid_task_title}/")
    if not os.path.exists(f"./Result_{valid_task_title}/log/"):
        os.makedirs(f"./Result_{valid_task_title}/log/")

    if is_folder:
        return f"./Result_{valid_task_title}/"
    if is_sub:
        return f"./Result_{valid_task_title}/log/{valid_task_title}_result.{extension}"
    return f"./Result_{valid_task_title}/{valid_task_title}_result.{extension}"


def get_sub_log_address(file_name: str, ref_name: str, extension: str):
    valid_task_title = "".join([c for c in glv.TASK_TITLE if c not in "\/:*?<>| -"])
    if file_name[-5:] == str(".fastq")[-5:]:
        valid_tested_file_name = "".join([c for c in file_name[:-6] if c not in "\/:*?<>| -"])
    else:
        valid_tested_file_name = "".join([c for c in file_name[:-9] if c not in "\/:*?<>| -"])
    valid_ref_name = "".join([c for c in ref_name if c not in "\/:*?<>| -"])

    address = get_main_log_name('txt', is_folder=True) + 'log/' + valid_task_title + "--" + valid_tested_file_name + "--" \
              + valid_ref_name + "." + extension

    return address


def write_sub_log(aligned_read_list: List[Aligned_Read], genotyper: Genotyper_For_Reference, file_name: str):
    file_log = open(get_sub_log_address(file_name, genotyper.ref_name, 'txt'), 'w')

    file_log.write(f""
                   f"# <CNS-Genotyper {glv.VERSION} Side Log for {file_name}>\n"
                   f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n"
                   f"# \n"
                   f"# Task Title: {glv.TASK_TITLE}\n"
                   f"# {file_name} as a data / {genotyper.ref_name} as a reference sequence\n"
                   f"\n"
                   f"{glv.get_text_of_global_variables()}"
                   f"\n"
                   f"{_showing_selected_area_to_text(guide_rna_seq=genotyper.guide_rna_seq)}"
                   f"\n"
                   f"\n")

    file_log.write(genotyper.get_abstract_text())
    file_log.write("\n"
                   "\n"
                   "----------------------\n")

    # file_log.write(genotyper.get_examples_text())
    # file_log.write("\n"
    #                "\n"
    #                "----------------------\n")

    file_log.write("[Raw Data]\n"
                   "\n"
                   "\n")

    for aligned_read in aligned_read_list:
        if aligned_read.ref_name == genotyper.ref_name:
            file_log.write(str(aligned_read))
            file_log.write("\n"
                           "\n")

    file_log.close()


def write_main_log(genotyper_list_list: List[List[Genotyper_For_Reference]], total_length: int):
    main_log_name = get_main_log_name("txt", is_sub=True)
    file_log = open(main_log_name, "w")

    file_log.write(f"# <CNS-Genotyper {glv.VERSION} Main Log>\n"
                   f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n"
                   f"# \n"
                   f"# Task Title: {glv.TASK_TITLE}\n"
                   f"\n"
                   f"{glv.get_text_of_global_variables()}\n"
                   f"\n"
                   f"total length: {total_length}\n"
                   f"\n"
                   f"\n")

    for genotyper_list in genotyper_list_list:
        file_log.write(f"<{genotyper_list[0].file_name}>\n"
                       f"\n")
        for genotyper in genotyper_list:
            file_log.write(genotyper.get_abstract_text())
            file_log.write("\n"
                           "\n"
                           "\n")

    file_log.close()


def write_main_html_log(genotyper_list_list: List[List[Genotyper_For_Reference]], total_length: int):
    main_log_name = get_main_log_name("html")
    file_log = open(main_log_name, "w")

    file_log.write(f"<!DOCTYPE html>\n"
                   "<html lang=\"en\">\n"
                   "<head>\n"
                   "    <meta charset=\"UTF-8\">\n"
                   "    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
                   f"    <title>Result for {glv.TASK_TITLE}</title>\n"
                   "\n"
                   "    <style>\n"
                   "        body {font-family:Consolas;}\n"
                   "        footer {text-align: center; background-color: #EBEFF3; padding: 10px;}\n"
                   "        .file_name {color: blue; font-weight: bold; display: inline;}\n"
                   "        .ref_name {color: red; font-weight: bold; display: inline;}\n"
                   "        .important {color: green; font-weight: bold; display: inline;}\n"
                   "        .not_important {color: gray; font-weight: bold; display: inline;}\n"
                   "        .point {color: red; font-weight: bold; display: inline;}\n"
                   "        .file {margin: 3%}\n"
                   "        .heading_info {margin: 3%}\n"
                   "\n"
                   "    </style>\n"
                   "</head>\n"
                   "<body>\n"
                   "    \n")

    html_global_variables = glv.get_text_of_global_variables().replace('\n', '<br>')

    file_log.write(f"    <div class=heading_info><h2> [CNS-Genotyper {glv.VERSION} Main Log] </h2>\n"
                   f"    <h2> Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()}) </h2>\n"
                   f"    <h2> Task Title: {glv.TASK_TITLE}</h2>\n"
                   f"    <br>\n"
                   f"    {html_global_variables}<br>\n"
                   f"    <br>\n"
                   f"    total reads: {total_length} reads<br>\n"
                   f"    </div>\n")

    for genotyper_list in genotyper_list_list:
        file_log.write(
            f"<div class=file id={genotyper_list[0].file_name}><h2 class=file_name>[{genotyper_list[0].file_name}]</h2>\n")
        for genotyper in genotyper_list:
            example_text = genotyper.get_abstract_text(is_html=True).replace("\n", "<br>").replace('  ', '&nbsp;&nbsp;')
            file_log.write(f"<div class=ref id={genotyper.ref_name}>")

            file_log.write(example_text)
            file_log.write("\n"
                           "<br></div>\n")
        file_log.write("</div>\n")

    file_log.write("<footer>"
                   "<b3><div class=point>C</div>as9 Mutant Genotyper with "
                   "<div class=point>N</div>ext Generation <div class=point>S</div>equencing Data</b3><br>"
                   "<b3>CNS-Genotyper by Octo Moon</b3><br>\n"
                   "<b3>Contact: octahedron00@gmail.com "
                   "or <a href=\"https://github.com/octahedron00\">GitHub</a></b3>\n"
                   "</footer>\n"
                   "</body>\n"
                   "</html>")
    file_log.close()


def write_main_csv_log(genotyper_list_list: List[List[Genotyper_For_Reference]], ref_set_list: List[Reference]):
    csv_log_name = get_main_log_name("csv")
    xlsx_log_name = get_main_log_name("xlsx")
    file_csv = open(csv_log_name, 'w', newline="")
    file_csv_writer = csv.writer(file_csv)

    file_csv_writer.writerow([f"<CNS-Genotyper {glv.VERSION} Main Log>"])
    file_csv_writer.writerow(
        [f"Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})"])
    file_csv_writer.writerow(["Task Title", glv.TASK_TITLE])
    file_csv_writer.writerows(glv.get_row_of_global_variables())
    file_csv_writer.writerow([])

    file_csv_writer.writerow(["References:"])
    file_csv_writer.writerow(["Name", "seq", "Guide RNA name", "Guide RNA seq"])
    for ref_set in ref_set_list:
        file_csv_writer.writerow([ref_set.ref_name, ref_set.ref_seq, ref_set.guide_rna_name, ref_set.guide_rna_seq])

    file_csv_writer.writerow([])
    file_csv_writer.writerow(["Results:"])
    file_csv_writer.writerow(["file", "reference", "warning", "genotype", "_of",
                              "# total", "# error", "# not err", "# of genotype"])

    for genotyper_list in genotyper_list_list:
        total_lines_for_file = 0

        for genotyper in genotyper_list:
            total_lines_for_file += len(genotyper)

        for i, genotyper in enumerate(genotyper_list):
            row = [""]
            genotype = genotyper.get_genotype()
            sorted_count_map_list = genotyper.get_sorted_count_map_list()
            len_without_err = genotyper.get_len(with_err=False) + Z
            len_total = len(genotyper) + Z
            if i == 0:
                row[0] = genotyper.file_name
            if i == 1:
                row[0] = f"{total_lines_for_file} lines"
            row += [genotyper.ref_name,
                    genotype.warning.strip().replace("\n", ", "),
                    genotype.allele_set_text,
                    str(genotype).splitlines()[0].strip(),
                    len(genotyper),
                    f"err {genotyper.count_map['err']}({round(genotyper.count_map['err'] / len_total, 3)})",
                    genotyper.get_len(with_err=False)]

            for key, value in sorted_count_map_list:
                if key == 'err':
                    continue
                row.append(f"{key} {value}({round(value / len_without_err, 3)})")
            file_csv_writer.writerow(row)
        file_csv_writer.writerow([])

    file_csv.close()

    #
    # Write Excel file: CSV is not working well in Excel software
    workbook = Workbook(xlsx_log_name)
    worksheet = workbook.add_worksheet()

    italic_format = workbook.add_format({'italic': True, 'font_color': 'silver'})
    with open(csv_log_name, 'rt', encoding='utf8') as f:
        reader = csv.reader(f)
        for r, row in enumerate(reader):
            cell_format = workbook.add_format({})
            if len(row) > 2 and len(str(row[2])) > 8:
                if str(row[2])[:5] in (str("Reads not enough")[:5], str("Error only")[:5]):
                    cell_format = italic_format

            for c, col in enumerate(row):
                worksheet.write(r, c, col)
                if c > 0:
                    worksheet.write(r, c, col, cell_format)
        worksheet.set_column_pixels(0, 0, 255 / 1.4275)
        worksheet.set_column_pixels(2, 2, 255 / 1.4275)
        worksheet.set_column_pixels(4, 4, 660 / 1.4275)
        worksheet.set_column_pixels(6, 6, 127 / 1.4275)
        worksheet.set_column_pixels(8, 9, 150 / 1.4275)
        worksheet.set_column_pixels(10, 80, 136 / 1.4275)
    workbook.close()

    if os.path.exists(csv_log_name):
        os.remove(csv_log_name)


def _showing_selected_area_to_text(guide_rna_seq: str):
    # a little hardcoded: not for change...
    pos_line = ""
    ref_line = guide_rna_seq
    selected_area_line = ""

    # pos_rna_end = position rna ending
    pos_rna_end = len(ref_line)
    std_pos = len(ref_line)
    cut_pos = len(ref_line) + glv.CUT_POS_FROM_PAM

    ref_line += "NGG_______"

    for i, a in enumerate(ref_line):
        if i < pos_rna_end:
            pos_line += '>'
        elif std_pos <= i < std_pos + 3:
            pos_line += '<'
        else:
            pos_line += ' '

        if (cut_pos - glv.CUT_POS_RADIUS) <= i < cut_pos:
            selected_area_line += '['
        elif cut_pos <= i < (cut_pos + glv.CUT_POS_RADIUS):
            selected_area_line += ']'
        else:
            selected_area_line += ' '

    pos_line += '    '
    ref_line += ' or '
    selected_area_line += '    '

    new_starting_point = len(ref_line)

    ref_line += guide_rna_seq

    pos_rna_end = len(ref_line)
    std_pos = len(ref_line) + 1
    cut_pos = len(ref_line) + 1 + glv.CUT_POS_FROM_PAM

    ref_line += "_NGG_______"

    for i, a in enumerate(ref_line):
        if i < new_starting_point:
            continue

        if i < pos_rna_end:
            pos_line += '>'
        elif std_pos <= i < std_pos + 3:
            pos_line += '<'
        else:
            pos_line += ' '

        if (cut_pos - glv.CUT_POS_RADIUS) <= i < cut_pos:
            selected_area_line += '['
        elif cut_pos <= i < (cut_pos + glv.CUT_POS_RADIUS):
            selected_area_line += ']'
        else:
            selected_area_line += ' '

    return f"Selected area for estimated indel position is: \n" \
           f"{pos_line}\n" \
           f"{ref_line}\n" \
           f"{selected_area_line}\n"


def write_raw_data_log(genotyper_list_list: List[List[Genotyper_For_Reference]], debug_data: Dict[str, Dict]):
    debug_log_name = get_main_log_name("raw", is_sub=True)

    with open(debug_log_name, 'w') as file_log:
        if glv.DEBUG:
            for file_name, data_dict in debug_data.items():
                file_log.write(f"\n{file_name}\n")
                for key, value in data_dict.items():
                    file_log.write(f"\t{key}: {value}\n")
            file_log.write("\n\n")

        file_log.write(f"indel_type\tinsertion\tdeletion\tindel_length\tindel_position\t# count\n")

        for genotyper_list in genotyper_list_list:
            for genotyper in genotyper_list:
                file_log.write(f"\n"
                               f"{genotyper.file_name}--{genotyper.ref_name}\n")

                for key, value in genotyper.get_sorted_count_map_list():
                    indel_i, indel_d, indel_length, indel_pos = get_indel_type_text_unpacked(key)
                    file_log.write(f"{key}\t{indel_i}\t{indel_d}\t{indel_length}\t{indel_pos}\t{value}\n")


def get_indel_type_text_unpacked(indel_type: str):
    indel_i = 0
    indel_d = 0
    indel_length = 0
    indel_pos = 0

    if indel_type in ('err', 'WT'):
        return indel_i, indel_d, indel_length, indel_pos

    set_of_letter_between_number = ('I', 'D', ':', ';', ',', 'F')
    # Including phred score to show

    indel_type_lines = indel_type
    for letter in set_of_letter_between_number:
        indel_type_lines = indel_type_lines.replace(letter, '\n' + letter + '\n')

    str_set = indel_type_lines.splitlines(keepends=False)

    a = 0
    for s in str_set:
        if s == 'I':
            indel_i = a
            a = 0
        elif s == 'D':
            indel_d = a
            a = 0
        elif re.match('^-?\d+$', s):
            a = int(s)

    indel_pos = a
    indel_length = max(indel_i, indel_d)

    return indel_i, indel_d, indel_length, indel_pos
