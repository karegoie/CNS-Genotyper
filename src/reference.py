from Bio import SeqIO


class Reference:
    '''
    class Reference
    a set of reference for alignment and indel validation.

    The class has:
        a reference sequence to align,
        a guide_RNA (spacer) sequence to validate the main indels,
        a name for both sequence to show.
    '''

    ref_raw = None
    ref_seq = ""
    ref_name = ""
    guide_rna_raw = None
    guide_rna_seq = ""
    guide_rna_name = ""

    def __init__(self, ref_raw: SeqIO.SeqRecord, guide_rna_raw: SeqIO.SeqRecord):
        self.ref_raw = ref_raw
        self.ref_seq = str(ref_raw.seq).upper()
        self.ref_name = ref_raw.name
        self.guide_rna_raw = guide_rna_raw
        self.guide_rna_seq = str(guide_rna_raw.seq).upper()
        self.guide_rna_name = guide_rna_raw.name

    def __len__(self):
        return len(self.ref_seq)

    def __str__(self):
        ref_seq_for_print = self.ref_seq
        if len(ref_seq_for_print) < 30:
            ref_seq_for_print = ref_seq_for_print[:20] + "..." + ref_seq_for_print[-5:]

        return f"<Class Reference>\n" \
               f"ref_raw: {self.ref_raw}\n" \
               f"ref_seq: {ref_seq_for_print}\n" \
               f"ref_name: {self.ref_name}\n" \
               f"guide_rna_raw: {self.guide_rna_raw}\n" \
               f"guide_rna_seq: {self.guide_rna_seq}\n" \
               f"guide_rna_name: {self.guide_rna_name}\n"
