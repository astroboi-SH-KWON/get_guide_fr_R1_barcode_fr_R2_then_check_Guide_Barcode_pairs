import openpyxl

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def csv_to_list_ignr_header(self, path, deli_str=","):
        result_list = []
        with open(path, "r") as f:
            header = f.readline()
            print(header)
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == '':
                    break

                result_list.append(tmp_line.split(deli_str))
        return result_list

    def make_row(self, sheet, row, data_arr, col=1):
        for idx in range(len(data_arr)):
            sheet.cell(row=row, column=(col + idx), value=data_arr[idx])

    def make_excel(self, path, data_dict):
        for key, data_list in data_dict.items():
            workbook = openpyxl.Workbook()
            sheet = workbook.active
            row = 1
            header_arr = ['r1_seq', 'r2_seq', 'guide_seq', 'rvrsed_brcd_seq', 'brcd_seq', 'PAM_index_by_barcode',
                          'flag', 'PAM_index_by_guide']
            self.make_row(sheet, row, header_arr)

            for data_arr in data_list:
                row += 1
                self.make_row(sheet, row, data_arr)

            workbook.save(filename=path + "_" + str(key) + self.ext_xlsx)

