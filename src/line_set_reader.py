

DIV_LETTER = ('@',)
SEQ_LETTER = ('A', 'T', 'G', 'C', '-')
MATCH_LETTER = ('|', 'x', '+', '\\', '/', '.')


# variables / functions naming rule
# text_ = english message
# line_ = sequence or string to perform with
# is_ = boolean

# etc. = int / char / etc.

# The sequence must have at least 5 NT span from the guide RNA sequence


def get_reformat_line_set_list(match_result_list: list):
    line_set_list = []

    # Making multiple-line seq to one long line
    # p_1: data / p_2: match / p_3: ref_seq position,
    p11, p12, p13, p21, p22, p23 = -1, -1, -1, -1, -1, -1
    for i, line in enumerate(match_result_list):
        if line[0] in DIV_LETTER:
            p11, p12, p13, p21, p22, p23 = -1, -1, -1, -1, -1, -1
        if line[0] in SEQ_LETTER:
            if p12 < 0:
                p11 = i
            elif p12 >= 0 and p13 < 0:
                p13 = i
            elif p13 >= 0 and p21 < 0:
                p21 = i
            elif p22 >= 0 and p23 < 0:
                p23 = i
                # print(f"{p21} to {p11}, {p22} to {p12}, {p23} to {p13}")
                match_result_list[p11] = match_result_list[p11][:-1] + match_result_list[p21]
                match_result_list[p12] = match_result_list[p12][:-1] + match_result_list[p22]
                match_result_list[p13] = match_result_list[p13][:-1] + match_result_list[p23]
                match_result_list[p21] = match_result_list[p22] = match_result_list[p23] = '.'
                p21 = p22 = p23 = -1

        if line[0] in MATCH_LETTER:
            if p12 < 0:
                p12 = i
            else:
                p22 = i

    # Making a list
    p11 = p12 = p13 = -1
    for i, line in enumerate(match_result_list):
        if line[0] in DIV_LETTER:
            line_set_list.append([match_result_list[p11], match_result_list[p12], match_result_list[p13]])
            p11 = p12 = p13 = -1
        if line[0] in SEQ_LETTER:
            if p11 < 0:
                p11 = i
            elif p12 >= 0 and p13 < 0:
                p13 = i
        if line[0] in MATCH_LETTER:
            p12 = i

    return line_set_list

