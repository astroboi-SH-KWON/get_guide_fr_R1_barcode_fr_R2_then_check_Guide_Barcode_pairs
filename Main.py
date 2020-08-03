import time
import os


import Util
import LogicPrep
import Logic
############### start to set env ###############
# WORK_DIR = os.getcwd() + "/"
WORK_DIR = "D:/000_WORK/KimNahye/20200803/WORK_DIR/"
PROJECT_NAME = WORK_DIR.split("/")[-2]
INPUT = "input/"
OUTPUT = "output/"
FASTQ = "NGS_READ/"
Guide_barcode = "Guide_barcode.txt"
FASTQ_pair = "FASTQ_list.txt"

BEFORE_GUIDE_R1 = "CCAGCAGGTCCCATGGTGTAATGGTtAGCACTCTGGACTTTGAATCCAGCGaTCCGAGTTCAAATCTCGGTGGGACCT"
LEN_GUIDE_R1 = 20
RP_SEQ_R2 = "CAGAAGACGGCATACGA"
BLANK_BP_R2 = 33
LEN_BRCD_R2 = 15

INIT = [BEFORE_GUIDE_R1, LEN_GUIDE_R1, RP_SEQ_R2, BLANK_BP_R2, LEN_BRCD_R2]
############### end setting env ################

def main():
    util = Util.Utils()
    logic_prep = LogicPrep.LogicPreps()
    logic = Logic.Logics(INIT)

    idx_g_b_list = util.csv_to_list_ignr_header(WORK_DIR + INPUT + Guide_barcode, "\t")

    fastq_pairs = util.csv_to_list_ignr_header(WORK_DIR + INPUT + FASTQ_pair, "\t")
    fastq_r1_dict, fastq_r2_dict = logic_prep.get_r1_r2_pairs(WORK_DIR + FASTQ, "fastq", fastq_pairs)

    result_dict, err_dict = logic.analyze_(fastq_r1_dict, fastq_r2_dict, idx_g_b_list)

    util.make_excel(WORK_DIR + OUTPUT + "result", result_dict)

if __name__ == '__main__':
    start_time = time.perf_counter()
    print("start [ " + PROJECT_NAME + " ]>>>>>>>>>>>>>>>>>>")
    main()
    print("::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))