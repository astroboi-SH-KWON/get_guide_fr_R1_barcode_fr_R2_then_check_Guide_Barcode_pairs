from Bio.Seq import Seq

class Logics:
    def __init__(self, init):
        self.bef_guide_R1 = init[0].upper()
        self.len_guide_R1 = init[1]
        self.rp_seq_R2 = init[2].upper()
        self.blank_size_R2 = init[3]
        self.len_brcd_R2 = init[4]

    def find_guide_from_r1(self, r1_seq):
        len_bef_guide_r1 = len(self.bef_guide_R1)
        pos_bef_guide = r1_seq.find(self.bef_guide_R1)
        if pos_bef_guide == -1:
            return ""
        strt_idx = pos_bef_guide + len_bef_guide_r1
        return r1_seq[strt_idx: strt_idx + self.len_guide_R1]

    def find_brcd_from_r2(self, r2_seq):
        len_rp_seq = len(self.rp_seq_R2)
        pos_rp_seq = r2_seq.find(self.rp_seq_R2)
        if pos_rp_seq == -1:
            return "", ""
        strt_idx = pos_rp_seq + len_rp_seq + self.blank_size_R2
        result_seq = r2_seq[strt_idx: strt_idx + self.len_brcd_R2]
        return str(Seq(result_seq).reverse_complement()), result_seq

    def check_PAM_idx_by_guide_brcd(self, guide_seq, brcd_seq, idx_g_b_list):
        g_PAM_arr = []
        b_PAM_idx = -1
        for idx_g_b in range(len(idx_g_b_list)):
            g_seq = idx_g_b_list[idx_g_b][1]
            b_seq = idx_g_b_list[idx_g_b][2]
            if g_seq == guide_seq:
                g_PAM_arr.append(idx_g_b)
            if b_seq == brcd_seq:
                b_PAM_idx = idx_g_b  # barcode is unique

        return g_PAM_arr, b_PAM_idx

    def analyze_(self, r1_dict, r2_dict, idx_g_b_list):
        err_dict = {}
        result_dict = {}
        for fastq_key, fastq_r1_list in r1_dict.items():
            result_dict.update({fastq_key: []})
            err_dict.update({fastq_key: []})
            fastq_r2_list = r2_dict[fastq_key]
            for fastq_idx in range(len(fastq_r1_list)):
                r1_seq = fastq_r1_list[fastq_idx]
                r2_seq = fastq_r2_list[fastq_idx]
                g_seq = self.find_guide_from_r1(r1_seq)
                rvrsed_brcd_seq, brcd_seq = self.find_brcd_from_r2(r2_seq)
                tmp_arr = [r1_seq, r2_seq, g_seq, rvrsed_brcd_seq, brcd_seq]

                if g_seq == "":
                    tmp_arr.append("no_guide")
                    err_dict[fastq_key].append(tmp_arr)
                    continue
                if rvrsed_brcd_seq == "":
                    tmp_arr.append("no_rp_seq")
                    err_dict[fastq_key].append(tmp_arr)
                    continue

                g_PAM_arr, b_PAM_idx = self.check_PAM_idx_by_guide_brcd(g_seq, rvrsed_brcd_seq, idx_g_b_list)

                if b_PAM_idx == -1:
                    tmp_arr.append("no_barcode_in_Guide_barcode.xlsx")
                    err_dict[fastq_key].append(tmp_arr)
                    continue

                if len(g_PAM_arr) == 0:
                    tmp_arr.append("no_guide_in_Guide_barcode.xlsx")
                    err_dict[fastq_key].append(tmp_arr)
                    continue

                if b_PAM_idx in g_PAM_arr:
                    tmp_arr.append(idx_g_b_list[b_PAM_idx][0])
                    tmp_arr.append("O")
                    result_dict[fastq_key].append(tmp_arr)
                else:
                    tmp_arr.append(idx_g_b_list[b_PAM_idx][0])
                    tmp_arr.append("X")
                    for idx_g_b in g_PAM_arr:
                        tmp_arr.append(idx_g_b_list[idx_g_b][0])
                    result_dict[fastq_key].append(tmp_arr)

        return result_dict, err_dict




