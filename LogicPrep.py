from Bio import SeqIO

import Valid
class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def get_r1_r2_pairs(self, path, file_form, fastq_pairs):
        valid = Valid.Validation()

        r1_dict = {}
        r2_dict = {}
        for fastq_r1_r2 in fastq_pairs:
            tmp_r1 = list(SeqIO.parse(path + fastq_r1_r2[0] + "." + file_form, file_form))
            tmp_r2 = list(SeqIO.parse(path + fastq_r1_r2[1] + "." + file_form, file_form))

            if not valid.equal_lens(tmp_r1, tmp_r2, fastq_r1_r2):
                continue

            r1_dict[fastq_r1_r2[0] + "^" + fastq_r1_r2[1]] = [str(tmp_r1[j].seq).upper() for j in range(len(tmp_r1))]
            r2_dict[fastq_r1_r2[0] + "^" + fastq_r1_r2[1]] = [str(tmp_r2[k].seq).upper() for k in range(len(tmp_r2))]

        return r1_dict, r2_dict