import datetime

VERSION = '1.5.0'

# Set of Global Variables for editing things:
global EXPLANATION_MAP
EXPLANATION_MAP = dict()
EXPLANATION_MAP[''] = ""

# 2 variables for user experience
global TASK_TITLE, OPEN_XLSX_AUTO, READ_IGNORE, DEBUG

TASK_TITLE = "task " + str(datetime.datetime.now())[5:-10]
OPEN_XLSX_AUTO = False
READ_IGNORE = ['R2', 'Undetermined']
DEBUG = False
EXPLANATION_MAP['task_title'] = "Title for this task."
EXPLANATION_MAP['open_xlsx_auto'] = "Open the excel log file automatically if finished"
EXPLANATION_MAP['read_ignore'] = "ignoring some files by its name containing: " \
                                 "insert R2 to ignore all Read 2(default), " \
                                 "Undetermined to ignore all extra data(default), " \
                                 "or '' to ignore nothing."
EXPLANATION_MAP['debug'] = "Debug mode."


# 1 variables for formatting the data : will never change
# phred encoding for output! (for encoding, automatically set in the Bio.SeqIO)
global PHRED_ENCODING
PHRED_ENCODING = 33

# 9 variables for changing results
global PAM_DISTANCE_MAX, PHRED_MEANINGFUL_MIN, ERR_RATIO_MAX, ERR_PADDING_FOR_SEQ, CUT_POS_FROM_PAM, CUT_POS_RADIUS

ERR_RATIO_MAX = 0.03
ERR_PADDING_FOR_SEQ = 1
CUT_POS_FROM_PAM = -3
CUT_POS_RADIUS = 5
EXPLANATION_MAP['err_ratio_max'] = "The threshold of mismatch ratio in aligned line set, without the main indel."
EXPLANATION_MAP['err_padding_for_seq'] = "The mismatch in this length from both end will not be considered as an error."
EXPLANATION_MAP['cut_pos_from_pam'] = "Set estimated 'cut' position from the starting point of PAM sequence."
EXPLANATION_MAP['cut_pos_radius'] = "The indel within this radius from the estimated 'cut' pos will be considered."

PAM_DISTANCE_MAX = 5
PHRED_MEANINGFUL_MIN = 30
EXPLANATION_MAP['pam_distance_max'] = "The max distance of PAM sequence from guide RNA sequence"
EXPLANATION_MAP['phred_meaningful_score_min'] = "One mismatch will be shown as '1I1D', " \
                                          "only if the NT's phred score is higher than this; " \
                                          "For ignoring all 'one mismatch', make it higher than 100"

global MAT, MIS, GAP_OPEN, GAP_EXTEND
MAT = 2
MIS = -1
GAP_OPEN = -30
GAP_EXTEND = -4
EXPLANATION_MAP['score_match'] = "Score for align: for Match"
EXPLANATION_MAP['score_mismatch'] = "Score for align: for Mismatch"
EXPLANATION_MAP['score_gap_open'] = "Score for align: for Gap Open"
EXPLANATION_MAP['score_gap_extend'] = "Score for align: for Gap Extension; " \
                                    "low penalty for 'gaps' will make fake indels more often"


# for generating align matrix
def get_align_matrix_for_subsequence_positioning():
    # 'X' = for both end of main sequence, meaning the subsequence must be between the sequence.
    matrix = {
        ('A', 'A'): MAT, ('A', 'T'): MIS, ('A', 'G'): MIS, ('A', 'C'): MIS, ('A', '-'): MIS, ('A', 'X'): -1000,
        ('T', 'A'): MIS, ('T', 'T'): MAT, ('T', 'G'): MIS, ('T', 'C'): MIS, ('T', '-'): MIS, ('T', 'X'): -1000,
        ('G', 'A'): MIS, ('G', 'T'): MIS, ('G', 'G'): MAT, ('G', 'C'): MIS, ('G', '-'): MIS, ('G', 'X'): -1000,
        ('C', 'A'): MIS, ('C', 'T'): MIS, ('C', 'G'): MIS, ('C', 'C'): MAT, ('C', '-'): MIS, ('C', 'X'): -1000,
        ('-', 'A'): MIS, ('-', 'T'): MIS, ('-', 'G'): MIS, ('-', 'C'): MIS, ('-', '-'): MAT, ('-', 'X'): -1000,
        ('X', 'A'): -1000, ('X', 'T'): -1000, ('X', 'G'): -1000, ('X', 'C'): -1000, ('X', '-'): -1000, ('X', 'X'): MAT,
    }
    return matrix


# for txt log
def get_text_of_global_variables():
    txt = f"" \
          f"ERR_RATIO_MAX: {ERR_RATIO_MAX}\n" \
          f"ERR_PADDING_FOR_SEQ: {ERR_PADDING_FOR_SEQ}\n" \
          f"CUT_POS_FROM_PAM: {CUT_POS_FROM_PAM}\n" \
          f"CUT_POS_RADIUS: {CUT_POS_RADIUS}\n" \
          f"PHRED_MEANINGFUL_SCORE_MIN: {PHRED_MEANINGFUL_MIN}\n" \
          f"PAM_DISTANCE_MAX: {PAM_DISTANCE_MAX}\n" \
          f"SCORE_MATCH: {MAT} / SCORE_MISMATCH: {MIS}\n" \
          f"SCORE_GAP_OPEN: {GAP_OPEN} / SCORE_GAP_EXTEND: {GAP_EXTEND}\n"
    return txt


# for csv, xlsx log
def get_row_of_global_variables():
    row = [["ERR_RATIO_MAX", ERR_RATIO_MAX, "ERR_PADDING_FOR_SEQ", ERR_PADDING_FOR_SEQ],
           ["CUT_POS_FROM_PAM", CUT_POS_FROM_PAM, "CUT_POS_RADIUS", CUT_POS_RADIUS],
           ["PHRED_MEANINGFUL_SCORE_MIN", PHRED_MEANINGFUL_MIN, "PAM_DISTANCE_MAX", PAM_DISTANCE_MAX],
           ["SCORE_MATCH", MAT, "SCORE_MISMATCH", MIS],
           ["SCORE_GAP_OPEN,", GAP_OPEN, "SCORE_GAP_EXTEND", GAP_EXTEND]]
    return row


def get_opposite_strand(seq: str):
    opposite = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'
    }

    seq = seq.upper()
    new_seq = ""
    for s in seq[::-1]:
        if s not in opposite.keys():
            continue
        new_seq += opposite[s]

    return new_seq


# for cas9
PAM_IS_AFTER = True
PAM_SEQ_REGEX_AS_DOWNSTREAM = r"[ATGC]GG"
