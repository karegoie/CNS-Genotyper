from Bio import SeqIO
import re

import src.globals as glv


def get_pos_line_of_ref_seq(ref_seq, guide_rna_seq):

    pos_line = ""

    a = ref_seq.find(guide_rna_seq)
    if a < 0:
        return None

    if glv.PAM_IS_AFTER:
        p = a + len(guide_rna_seq)

        pam_match = re.search(glv.PAM_SEQ_REGEX_AS_DOWNSTREAM, ref_seq[p:])
        if pam_match is None:
            return None
        pam_start = pam_match.start()
        pam_end = pam_match.end()

        if pam_start > glv.PAM_DISTANCE_MAX:
            return None

        pos_line += " " * a
        pos_line += ">" * len(guide_rna_seq)
        pos_line += " " * pam_start
        pos_line += "<" * (pam_end - pam_start)
        pos_line += " " * (len(ref_seq) - len(pos_line))

        # set cut position
        std_pos = pos_line.find('<')
        cut_pos = std_pos + glv.CUT_POS_FROM_PAM
        pos_line = pos_line[:cut_pos-1] + ')(' + pos_line[cut_pos+1:]
        pos_line = pos_line[:pos_line.find(">") - 5] + '|' + pos_line[pos_line.find(">") - 4:]
        pos_line = pos_line[:std_pos + pam_end + 5] + '|' + pos_line[std_pos + pam_end + 6:]
        print(pos_line)

    else:
        # not finished yet... for Cpf1
        opp_seq = glv.get_opposite_strand(ref_seq[:a])

        pam_match = re.search(glv.PAM_SEQ_REGEX_AS_DOWNSTREAM, opp_seq)
        if pam_match is None:
            return None
        pam_start = pam_match.start()
        pam_end = pam_match.end()

        if pam_start > glv.PAM_DISTANCE_MAX:
            return None

        pos_line += " " * (a - pam_end)
        pos_line += ">" * (pam_end - pam_start)
        pos_line += " " * pam_start
        pos_line += "<" * len(guide_rna_seq)
        pos_line += " " * (len(ref_seq) - len(pos_line))

    return pos_line


class Reference:
    '''
    class Reference
    a set of reference for alignment and indel validation.

    The class has:
        a reference sequence to align,
        a guide_RNA (spacer) sequence to validate the main indels,
        a name for both sequence to show.
    '''

    ref_seq = ""
    ref_name = ""
    guide_rna_seq = ""
    guide_rna_name = ""
    ref_pos_line = ""

    def __init__(self, ref_raw: SeqIO.SeqRecord, guide_rna_raw: SeqIO.SeqRecord):
        ref_seq = str(ref_raw.seq).upper()
        guide_rna_seq = str(guide_rna_raw.seq).upper()

        ref_pos_line = get_pos_line_of_ref_seq(ref_seq, guide_rna_seq)
        if ref_pos_line is None:
            ref_seq = glv.get_opposite_strand(ref_seq)
            ref_pos_line = get_pos_line_of_ref_seq(ref_seq, guide_rna_seq)
        if ref_pos_line is None:
            raise ValueError

        #
        # set variables below
        self.ref_seq = ref_seq
        self.guide_rna_seq = guide_rna_seq
        self.ref_pos_line = ref_pos_line

        self.ref_name = str(ref_raw.name)
        self.guide_rna_name = str(guide_rna_raw.name)

        print(ref_seq)
        print(ref_pos_line)

    def __len__(self):
        return len(self.ref_seq)

    def __str__(self):
        ref_seq_for_print = self.ref_seq
        if len(ref_seq_for_print) < 30:
            ref_seq_for_print = ref_seq_for_print[:20] + "..." + ref_seq_for_print[-5:]

        return f"<Class Reference>\n" \
               f"ref_seq: {ref_seq_for_print}\n" \
               f"ref_name: {self.ref_name}\n" \
               f"guide_rna_seq: {self.guide_rna_seq}\n" \
               f"guide_rna_name: {self.guide_rna_name}\n"
