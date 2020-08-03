

class Validation:
    def __init__(self):
        pass

    def equal_lens(self, tmp_r1, tmp_r2, fastq_r1_r2):
        if len(tmp_r1) == len(tmp_r2):
            return True
        print("lens are different : ", fastq_r1_r2)
        return False
