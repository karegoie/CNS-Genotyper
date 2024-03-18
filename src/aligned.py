import re

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO
from Bio.SeqIO import SeqRecord

from src.reference import Reference

# Variables for uh... making the results:
# Can be changed by command line interface; All in globals.py
import src.globals as glv

# Variables for letter recognition
MATCH_LETTER = ('|', 'x', '+', '\\', '/', '.')
MATCH_ERR_LETTER = MATCH_LETTER[1:]

# Variables for uh... zero division
Z = 0.000000001


#
# def get_guide_rna_seq_position(ref_line: str, guide_rna_seq: str):
#     pos_rna_start = -1
#     pos_rna_end = -1
#     if len(guide_rna_seq) > 0:
#         for i in range(len(ref_line) - len(guide_rna_seq)):
#             k, is_good = 0, False
#             pos_rna_start = i
#             for j, n in enumerate(guide_rna_seq):
#                 while ref_line[i + j + k] == '-' and (i+j+k) < len(ref_line):
#                     k += 1
#                 if (i+j+k) >= len(ref_line) or ref_line[i + j + k] != n:
#                     break
#                 if j == (len(guide_rna_seq) - 1):
#                     is_good = True
#                     pos_rna_end = (i + j + k + 1)
#             if is_good:
#                 break
#             else:
#                 pos_rna_start, pos_rna_end = -1, -1
#     return pos_rna_start, pos_rna_end


class Aligned:

    pos_line = ""
    ref_line = ""
    match_line = ""
    read_line = ""
    phred_line = ""

    read_name = ""
    ref_name = ""
    file_name = "(file_name not defined)"

    reference = None
    indel = None
    strand = "?"

    score = 0
    int_score = 0

    def __len__(self):
        return len(self.read_line)

    def __str__(self):
        return f"@\tref: {self.ref_name}\t /read: {self.read_name}\t /file: {self.file_name}\n" \
               f"@\tscore: {self.int_score}\t /indel: {self.indel}\n" \
               f"position: {self.pos_line.replace('(', '>').replace(')', '>').replace('|', ' ')}\n" \
               f"ref_line: {self.ref_line}\n" \
               f"match   : {self.match_line}\n" \
               f"readline: {self.read_line}\n" \
               f"phred   : {self.phred_line}\n"

    # # only for the beauty of the sub-log
    def get_str_simple(self):
        return f"r read: {self.read_name}\t /score: {self.int_score}\n" \
               f"| indel: {self.indel}\n" \
               f"| position: {self.pos_line.replace('(', '>').replace(')', '>').replace('|', ' ')}\n" \
               f"| ref_line: {self.ref_line}\n" \
               f"| match   : {self.match_line}\n" \
               f"| readline: {self.read_line}\n" \
               f"L phred   : {self.phred_line}\n"

class Indel_Type:
    indel_type = ""
    indel_reason = ""
    indel_length = 0
    indel_type_pos = 0

    def __init__(self, align_result: Aligned):

        match_line = align_result.match_line
        read_line = align_result.read_line
        ref_line = align_result.ref_line
        pos_line = align_result.pos_line

        std_pos = pos_line.find('<')
        cut_pos = pos_line.find('(')

        p1 = p2 = -1
        indel_i, indel_d, indel_length = 0, 0, 0

        indel_type = "WT"
        indel_reason = "(WT)"
        indel_pos = -9999

        for i, m in enumerate(match_line + str("X" * len(match_line))):
            if i < glv.ERR_PADDING_FOR_SEQ:
                continue
            if i >= (len(match_line) - glv.ERR_PADDING_FOR_SEQ) or match_line[i] == '$':
                break

            if m in MATCH_ERR_LETTER:
                p2 = i + 1
                indel_length += 1
                if p1 < 0:
                    p1 = i
                if m == 'x':
                    indel_d += 1
                if m == '+':
                    indel_i += 1
                if m == '.':
                    indel_d += 1
                    indel_i += 1

            elif p2 >= 0:
                if ((cut_pos - glv.CUT_POS_RADIUS) < p2) and (p1 < (cut_pos + glv.CUT_POS_RADIUS)):
                    indel_type = _get_indel_shape_text(indel_i, indel_d, p2 - std_pos, read_line[p1:p2])
                    indel_pos = p2 - std_pos
                    indel_reason = "Indel position confirmed by guide RNA sequence"
                    break
                p1 = p2 = -1
                indel_i = indel_d = indel_length = 0
                if cut_pos + glv.CUT_POS_RADIUS <= i:
                    break
        if indel_type == "WT":
            indel_pos = -9999
            indel_length = 0

        if indel_type[:4] == '1I1D':
            indel_length = 1

        self.indel_type = indel_type
        self.indel_pos = indel_pos
        self.indel_length = indel_length
        self.indel_reason = indel_reason

    def __str__(self):
        return f"{self.indel_type} ({self.indel_reason})"

    def test_err(self, score: float, align_result: Aligned, reference: Reference):
        pos_line = align_result.pos_line.replace('(', '>').replace(')', '>').replace(' ', '').replace('-', '')

        if pos_line.find('><') > 0:
            if (1 - score) < glv.ERR_RATIO_MAX:
                return

        self.indel_type = 'err'
        self.indel_pos = 0
        self.indel_length = 0
        self.indel_reason = 'Too many mismatch and error'


def _get_indel_shape_text(indel_i: int, indel_d: int, pos: int, insert_part: str = '', phred_part: str = ''):

    if indel_i == indel_d == 1 and glv.DEBUG:
        return f"1I1D{pos}{phred_part}"

    insert_part = insert_part.replace('-', '')
    if len(insert_part) > 2:
        insert_part = insert_part[:2]

    if indel_i == 0:
        if indel_d == 0:
            return 'WT'
        else:
            return f"{indel_d}D{pos}"
    else:
        if indel_d == 0:
            return f"{indel_i}I{pos}{insert_part}"
        else:
            return f"{indel_i}I{indel_d}D{pos}{insert_part}"


def get_score(match_line, indel_length):

    error = 0
    for a in match_line:
        if a not in ('|', '$'):
            error += 1

    return 1 - ((error - indel_length) / (len(match_line) - indel_length))


class Aligned_Key(Aligned):

    def __init__(self, key_seq: str, reference: Reference, strand: str):

        align_result = get_align_result_for_key(key_seq, reference)

        indel = Indel_Type(align_result=align_result)

        score = get_score(match_line=align_result.match_line, indel_length=indel.indel_length)
        int_score = int((score-1)*1000)

        indel.test_err(score, align_result, reference)

        # set variables below
        self.pos_line = align_result.pos_line
        self.ref_line = align_result.ref_line
        self.match_line = align_result.match_line
        self.read_line = align_result.read_line
        self.phred_line = ""

        self.read_name = "(aligned_key)"
        self.ref_name = reference.ref_name
        self.file_name = "(aligned_key)"

        self.reference = reference
        self.indel = indel
        self.strand = strand

        self.score = score
        self.int_score = int_score


def get_aligned_read_line(read_seq, match_line):
    read_line = ""

    # list: for the faster calculation
    read_line_list = []
    i = 0
    for m in match_line:
        if m == 'x':
            read_line_list.append('-')
        else:
            read_line_list.append(read_seq[i])
            i += 1

    read_line = "".join(read_line_list)
    return read_line


def get_aligned_phred_line(quality_list: list, match_line):
    phred_line = ""
    if len(quality_list) < 1:
        return phred_line

    for m in match_line:
        if m == 'x':
            phred_line += ' '
        else:
            phred_line += str(chr(quality_list[0] + glv.PHRED_ENCODING))
            quality_list = quality_list[1:]

    return phred_line


class Aligned_Read(Aligned):

    def __init__(self, aligned_key: Aligned_Key, read_raw: SeqRecord, file_name: str = "?"):

        read_line = get_aligned_read_line(str(read_raw.seq), aligned_key.match_line)

        quality_list = read_raw.letter_annotations['phred_quality']
        phred_line = get_aligned_phred_line(quality_list, aligned_key.match_line)

        # set variables below

        self.pos_line = aligned_key.pos_line
        self.ref_line = aligned_key.ref_line
        self.match_line = aligned_key.match_line
        self.read_line = read_line
        self.phred_line = phred_line

        self.read_name = read_raw.name
        self.ref_name = aligned_key.ref_name
        self.file_name = file_name

        self.reference = aligned_key.reference
        self.indel = aligned_key.indel

        self.score = aligned_key.score
        self.int_score = aligned_key.int_score


def get_pos_line_untrimmed(ref_line_untrimmed: str, ref_pos_line: str):

    is_in_focus = False

    pos_line_untrimmed = ""
    pos_line_untrimmed_list = []
    i = 0
    for s in ref_line_untrimmed:
        if s in ('-', 'X'):
            if is_in_focus:
                pos_line_untrimmed_list.append('-')
            else:
                pos_line_untrimmed_list.append(' ')
        else:
            if len(ref_pos_line) > i+2:
                if ref_pos_line[i] == '>':
                    is_in_focus = True
                if ref_pos_line[i] == '<' and ref_pos_line[i+1] == ' ':
                    is_in_focus = False
            pos_line_untrimmed_list.append(ref_pos_line[i])
            i += 1

    pos_line_untrimmed = "".join(pos_line_untrimmed_list)
    return pos_line_untrimmed


def get_alignment(ref_seq: str, read_seq: str):
    if len(read_seq) == 0:
        read_seq = 'A'

    aligner = PairwiseAligner()
    matrix = substitution_matrices.Array(data=glv.get_align_matrix_for_subsequence_positioning())
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = glv.GAP_OPEN
    aligner.extend_gap_score = glv.GAP_EXTEND
    alignments = aligner.align(ref_seq, read_seq)
    return alignments[0]


def get_longer_ref_line(ref_line: str, ref_line_untrimmed: str, padding: int):

    p = ref_line_untrimmed.find(ref_line)

    return ref_line_untrimmed[p-padding:p+len(ref_line)+padding]


def get_align_result_for_key(key_seq: str, reference: Reference):

    ref_seq = "X"*(glv.ERR_PADDING_FOR_SEQ+1) + reference.ref_seq + "X"*(glv.ERR_PADDING_FOR_SEQ+1)
    read_seq = key_seq.replace("N", "-")

    # align with the pre-set global variables for alignments.
    # get the untrimmed alignment result, based on the reference sequence
    alignment = get_alignment(ref_seq=ref_seq, read_seq=read_seq)
    ref_line_untrimmed = alignment[0]
    read_line_untrimmed = alignment[1]
    pos_line_untrimmed = get_pos_line_untrimmed(ref_line_untrimmed, reference.ref_pos_line)

    # trim the sequence, and build the match_line
    ref_line = match_line = read_line = pos_line = ""
    ref_line_list = []
    match_line_list = []
    read_line_list = []
    pos_line_list = []
    count_base = 0

    for i, c in enumerate(read_line_untrimmed):
        if count_base >= len(read_seq):
            break
        if len(read_line_list) == 0 and c == "-":
            continue
        read_line_list.append(c)
        ref_line_list.append(ref_line_untrimmed[i])
        pos_line_list.append(pos_line_untrimmed[i])
        count_base += 1
        if read_line_untrimmed[i] == '-':
            match_line_list.append('x')
            count_base -= 1
        elif ref_line_untrimmed[i] == '-':
            match_line_list.append('+')
        elif read_line_untrimmed[i] == ref_line_untrimmed[i]:
            match_line_list.append('|')
        else:
            match_line_list.append('.')

    ref_line = "".join(ref_line_list)
    match_line = "".join(match_line_list)
    read_line = "".join(read_line_list)
    pos_line = "".join(pos_line_list)

    align_result = Aligned()
    align_result.ref_line = get_longer_ref_line(ref_line, ref_line_untrimmed, glv.ERR_PADDING_FOR_SEQ)
    align_result.match_line = '$'*glv.ERR_PADDING_FOR_SEQ + match_line + '$'*glv.ERR_PADDING_FOR_SEQ
    align_result.pos_line = ' '*glv.ERR_PADDING_FOR_SEQ + pos_line + ' '*glv.ERR_PADDING_FOR_SEQ
    align_result.read_line = 'N'*glv.ERR_PADDING_FOR_SEQ + read_line + 'N'*glv.ERR_PADDING_FOR_SEQ

    return align_result

#
# class Line_Set:
#     '''
#     class Line_Set
#     a set of line(sequence and string) to show the alignment and indel type,
#         with a validated indel type.
#
#     The class has:
#         name of reference and read(id with number from sequencer), guide_RNA
#             to calculate / show / debug,
#
#         ref(erence)_line
#         read_line
#         match_line
#             which are aligned & trimmed to show the mismatch and indel.
#         pos(ition)_line
#             to show the position of guide_RNA with '>',
#             and the possible PAM sequence with '<'.
#         phred_line
#             aligned to show the signal strength(score) of each NT of read,
#
#         cut_pos(ition)
#             a calculated point... starting position to look forward.
#         std_pos(ition)
#             a starting point of PAM sequence: starting position to look forward.
#         indel_type
#         indel_reason
#             which are set by calculation,
#         indel_length
#             for the calculation of the score, like...
#             2I2D means two mismatch, but the indel_length is 2
#
#         score
#             calculated by the perfect match rate, without the main indel mismatch.
#         int_score
#             (score-1) * 1000 : <more big 'minus' value means more error>
#         phred_score
#             average phred score value, 'the quality' of the read.
#             all score is used to find the best line_set,
#             to check the best reference for a Read to match, and
#             to get the best line set for an InDel_Type
#
#
#     '''
#
#     reference = None
#     pos_line = ""
#     ref_line = ""
#     match_line = ""
#     read_line = ""
#     phred_line = ""
#
#     guide_rna_seq = ""
#
#     read_name = ""
#     ref_name = ""
#     guide_rna_name = ""
#     file_name = "(file_name not defined)"
#
#     indel_type = ""
#     indel_length = 0
#     indel_reason = ""
#     indel_type_pos = 0
#     cut_pos = 0
#     std_pos = 0
#     rna_pos = 0
#
#     # # only for the beauty of the log
#     indel_same_type_count = 0
#
#     score = 0
#     int_score = 0
#     phred_score = 0
#     #
#     # def __init__(self, read_raw: SeqIO.SeqRecord, reference: Reference):
#     #
#     #     self.reference = reference
#     #     self.ref_name = reference.ref_name
#     #     self.guide_rna_seq = reference.guide_rna_seq
#     #     self.guide_rna_name = reference.guide_rna_name
#     #
#     #     self.read_name = read_raw.name
#     #
#     #     # try aligning, and get ref, match, seq line here
#     #     self._set_align_line_set(ref_seq=reference.ref_seq, read_seq=str(read_raw.seq).upper())
#     #
#     #     # get aligned phred line here
#     #     self._set_aligned_phred_line_and_score(read_raw=read_raw)
#     #
#     #     # get pos line and cut position (endpoint of guide_RNA / startpoint of PAM)
#     #     self._set_indel_pos_line(guide_rna_seq=reference.guide_rna_seq)
#     #
#     #     # get indel type from the position... only if guide RNA sequence matches.
#     #     self._set_main_indel(cut_pos=self.cut_pos)
#     #
#     #     # check whether the type is error or not... and set align score with ratio of mismatch without main indel
#     #     self._set_score_and_check_err()
#
#     def __len__(self):
#         return len(self.ref_line)
#
#     def __str__(self):
#         return f"@\tref: {self.ref_name}\t /read: {self.read_name}\t /file: {self.file_name}\n" \
#                f"@\tscore: {self.int_score}\t /indel: {self.indel_type}({self.indel_same_type_count})\t /reason: \"{self.indel_reason}\"\n" \
#                f"position: {self.pos_line}\n" \
#                f"ref_line: {self.ref_line}\n" \
#                f"match   : {self.match_line}\n" \
#                f"readline: {self.read_line}\n" \
#                f"phred   : {self.phred_line}\n"
#
#     # # only for the beauty of the sub-log
#     def get_str_simple(self):
#         return f"r read: {self.read_name}\t /score: {self.int_score}\n" \
#                f"| indel: {self.indel_type}\t /reason: \"{self.indel_reason}\"\n" \
#                f"| position: {self.pos_line}\n" \
#                f"| ref_line: {self.ref_line}\n" \
#                f"| match   : {self.match_line}\n" \
#                f"| readline: {self.read_line}\n" \
#                f"L phred   : {self.phred_line}\n"

    # def set_file_name(self, file_name: str):
    #     self.file_name = file_name

    # def _set_align_line_set(self, ref_seq, read_seq):
    #     # 'X' means that the read must be between the reference
    #     # = Aligning the subsequence in the sequence
    #     ref_seq = "X" + ref_seq + "X"
    #
    #     read_seq = read_seq.replace("N", "-")
    #
    #     # align with the pre-set global variables for alignments.
    #     aligner = PairwiseAligner()
    #     matrix = substitution_matrices.Array(data=glv.get_align_matrix_for_subsequence_positioning())
    #     aligner.substitution_matrix = matrix
    #     aligner.open_gap_score = glv.GAP_OPEN
    #     aligner.extend_gap_score = glv.GAP_EXTEND
    #     alignments = aligner.align(ref_seq, read_seq)
    #
    #     # get the untrimmed alignment result, based on the reference sequence
    #     ref_line_untrimmed = alignments[0][0]
    #     read_line_untrimmed = alignments[0][1]
    #
    #     # trim the sequence, and build the match_line
    #     ref_line = match_line = read_line = ""
    #     count_base = 0
    #
    #     match_list = []
    #
    #     for i, c in enumerate(read_line_untrimmed):
    #         if count_base >= len(read_seq):
    #             break
    #         if len(read_line) == 0 and c == "-":
    #             continue
    #         read_line += c
    #         ref_line += ref_line_untrimmed[i]
    #         count_base += 1
    #         if read_line_untrimmed[i] == '-':
    #             match_list.append('x')
    #             count_base -= 1
    #         elif ref_line_untrimmed[i] == '-':
    #             match_list.append('+')
    #         elif read_line_untrimmed[i] == ref_line_untrimmed[i]:
    #             match_list.append('|')
    #         else:
    #             match_list.append('.')
    #
    #     match_line = "".join(match_list)
    #
    #     # set the lines
    #     self.ref_line = ref_line
    #     self.match_line = match_line
    #     self.read_line = read_line
    #
    # def _set_aligned_phred_line_and_score(self, read_raw: SeqIO.SeqRecord):
    #     # get quality score list[int] from SeqRecord: fastq files
    #     quality_list = read_raw.letter_annotations['phred_quality']
    #     if len(quality_list) < 1:
    #         return
    #
    #     # Build unaligned phred line from the quality score list, use 33-encoding
    #     phred_line = ""
    #     phred_score_sum = 0
    #     for a in quality_list:
    #         phred_line += str(chr(a + glv.PHRED_ENCODING))
    #         phred_score_sum += a
    #
    #     # Align the phred line
    #     aligned_phred_line = ""
    #     for c in self.read_line:
    #         if c == '-':
    #             aligned_phred_line += " "
    #         else:
    #             aligned_phred_line += phred_line[0]
    #             phred_line = phred_line[1:]
    #
    #     self.phred_score = phred_score_sum / (len(phred_line) + Z)
    #     self.phred_line = aligned_phred_line
    #
    # def _set_indel_pos_line(self, guide_rna_seq: str):
    #     # get the guide_rna_seq alignment, from the function outside
    #     pos_rna_start, pos_rna_end = get_guide_rna_seq_position(self.ref_line, self.guide_rna_seq)
    #
    #     # if guide RNA is not aligned: the sequence will be considered as an error,
    #     # since no indel is close the cut_pos(ition)
    #     if pos_rna_start < 0 or pos_rna_end < 0:
    #         self.cut_pos = -1000
    #         self.std_pos = -1000
    #         self.rna_pos = -1000
    #         return
    #
    #     # buffer: for finding the PAM sequence
    #     ref_line_buffer = self.ref_line + "X" * (glv.PAM_DISTANCE_MAX * 2)
    #
    #     # build the pos(ition)_line
    #     # add ' ' to the starting point of guide RNA
    #     pos_line = " " * pos_rna_start
    #     # add '>' at the position of guide RNA
    #     for i in range(pos_rna_start, pos_rna_end):
    #         if ref_line_buffer[i] == '-':
    #             pos_line += "-"
    #         else:
    #             pos_line += ">"
    #
    #     # find the first possible PAM sequence by 'NGG' (range: PAM_DISTANCE_MAX)
    #     k0 = 0
    #     for i in range(glv.PAM_DISTANCE_MAX):
    #         k1 = 1
    #         k2 = 1
    #         pos_line_add = ''
    #         # find first NT
    #         while ref_line_buffer[pos_rna_end + i + k0] == '-':
    #             k0 += 1
    #             pos_line += '-'
    #         pos_line_add += '<'
    #         # find second NT
    #         while ref_line_buffer[pos_rna_end + i + k0 + k1] == '-':
    #             k1 += 1
    #             pos_line_add += '-'
    #         pos_line_add += '<'
    #         # find third NT
    #         while ref_line_buffer[pos_rna_end + i + k0 + k1 + k2] == '-':
    #             k2 += 1
    #             pos_line_add += '-'
    #         pos_line_add += '<'
    #
    #         # check whether the second and the third NT is all 'G' -> check NGG PAM sequence
    #         if ref_line_buffer[pos_rna_end + i + k0 + k1] == ref_line_buffer[pos_rna_end + i + k0 + k1 + k2] == 'G':
    #             pos_line += pos_line_add
    #             break
    #         else:
    #             pos_line += ' '
    #     # cut or add more ' ' to make the both length match
    #     if len(pos_line) < len(self.ref_line):
    #         pos_line += " " * (len(self.ref_line) - len(pos_line))
    #     else:
    #         pos_line = pos_line[:len(self.ref_line)]
    #
    #     # So, where was the PAM?
    #     pos_pam = -2000
    #     for i, c in enumerate(pos_line):
    #         if c == '<':
    #             pos_pam = i
    #             break
    #
    #     # set the lines
    #     self.pos_line = pos_line
    #     # cut position = pam position + pre-set delta value(default: -3) - insertion
    #     self.cut_pos = pos_pam + glv.CUT_POS_FROM_PAM - k0
    #     # standard position = naming indels:
    #     # Each sequence of PAM will have the position number 1, 2, 3.
    #     self.std_pos = pos_pam
    #     self.rna_pos = pos_rna_end
    #
    # def _set_main_indel(self, cut_pos: int):
    #     '''
    #     set main indel, type and reason.
    #
    #     How it works:
    #         by using the end position of 'RNA', uh...
    #
    #     :param cut_pos:
    #     '''
    #
    #     read_line = self.read_line
    #     ref_line = self.ref_line
    #     match_line = self.match_line
    #     phred_line = self.phred_line
    #     std_pos = self.std_pos
    #     cut_pos = self.cut_pos
    #
    #     indel = _InDel(read_line=read_line, ref_line=ref_line, match_line=match_line, phred_line=phred_line, std_pos=std_pos, cut_pos=cut_pos)
    #
    #     self.indel_type_pos = indel.indel_type_pos
    #     self.indel_type = indel.indel_type
    #     self.indel_length = indel.indel_length
    #     self.indel_reason = indel.indel_reason
    #
    # def _set_score_and_check_err(self):
    #     if len(self) < 5:
    #         return 0
    #     a = 0
    #     for c in self.match_line[glv.ERR_PADDING_FOR_SEQ:-glv.ERR_PADDING_FOR_SEQ]:
    #         if c == '|':
    #             a += 1
    #     score = a / (len(self) - (2 * glv.ERR_PADDING_FOR_SEQ) - self.indel_length)
    #
    #     self.score = score
    #
    #     if score < (1 - glv.ERR_RATIO_MAX):
    #         self.indel_type = 'err'
    #         self.indel_reason = 'Too many mismatch and error'
    #     self.int_score = int((score - 1) * 1000)
    #
    # def set_indel_same_type_count(self, counter_dict: dict):
    #     self.indel_same_type_count = counter_dict[self.indel_type]

