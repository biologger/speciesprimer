#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
BASEPATH = "/home/bio/primerdesign/Chlamydia_pneumoniae/Pangenome/results/primer"
p3_output = os.path.join(BASEPATH, "primer3_output")

settings = 1
primerdatasets = []
primerdata = []
with open(p3_output) as f:
    for line in f:
        line = line.strip()
        if "P3_SETTINGS_FILE_END=":
            settings = 0
        if settings == 0:
            if not line == "=":
                primerdata.append(line)
#            if 'SEQUENCE_ID=' in line:
#                seqid = line.split("=")[1]
#            if "SEQUENCE_TEMPLATE=" in line:
#                temp = line.split("=")[1]
#            if 'PRIMER_NUM_RETURNED=' in line:
#                primernum = int(line.split('=')[1])
                
        if line == "=":
            primerdatasets.append(primerdata)
            primerdata = []
            
    print(primerdatasets)
#    print(primerdata)
            


def parse_Primer3_output(self, path_to_file):

    def parseSeqId(key, value):
        if key.startswith("SEQUENCE_ID"):
            if "group" in value:
                value = "g" + value.split("group_")[1]
            seq_id = value
            self.p3dict.update({seq_id: {"Primer_pairs": None}})
            p3list.append(seq_id)

    def parseTemplate(key, value):
        if key.startswith("SEQUENCE_TEMPLATE"):
            self.p3dict[p3list[-1]].update(
                {"template_seq": value})


    def countPrimer(key, value):
        if key.startswith("PRIMER_PAIR_NUM_RETURNED"):
            lookup = p3list[-1]
            self.p3dict[lookup]["Primer_pairs"] = int(value)
            for i in range(0, int(value)):
                self.p3dict[p3list[-1]].update({"Primer_pair_"+str(i): {}})

    def parseRightPrimer(key, value):
        pp = "Primer_pair_"
        if key.startswith("PRIMER_RIGHT"):
            i = key.split("_")[2]
            if key.endswith("_PENALTY"):
                primer_rpen = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_R_penalty": float(primer_rpen)})
            if key.endswith("_SEQUENCE"):
                primer_rseq = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_R_sequence": primer_rseq})
            if key.endswith("_TM"):
                right_TM = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_R_TM": float(right_TM)})

    def parseLeftPrimer(key, value):
        pp = "Primer_pair_"
        if key.startswith("PRIMER_LEFT"):
            i = key.split("_")[2]
            if key.endswith("_PENALTY"):
                primer_lpen = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_L_penalty": float(primer_lpen)})
            if key.endswith("_SEQUENCE"):
                primer_lseq = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_L_sequence": primer_lseq})
            if key.endswith("_TM"):
                left_TM = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_L_TM": float(left_TM)})

    def parseInternalProbe(key, value):
        pp = "Primer_pair_"
        if key.startswith("PRIMER_INTERNAL"):
            i = key.split("_")[2]
            if key.endswith("_PENALTY"):
                primer_ipen = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_I_penalty": float(primer_ipen)})
            if key.endswith("_SEQUENCE"):
                primer_iseq = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_I_sequence": primer_iseq})
            if key.endswith("_TM"):
                int_TM = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_I_TM": float(int_TM)})

    def parsePrimerPair(key, value):
        pp = "Primer_pair_"
        if key.startswith("PRIMER_PAIR"):
            i = key.split("_")[2]
            if (key.endswith('_PENALTY') and "WT" not in key):
                primer_ppen = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"primer_P_penalty": float(primer_ppen)})
            if (key.endswith('_PRODUCT_SIZE') and "WT" not in key):
                prod_size = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"product_size": int(prod_size)})
            if (key.endswith("_PRODUCT_TM") and "WT" not in key):
                prod_TM = value
                self.p3dict[p3list[-1]][pp + str(i)].update(
                    {"product_TM": float(prod_TM)})

    p3list = []
    info = "Run: parse_Primer3_output(" + self.target + ")"
    print(info)
    G.logger(info)
    problem = []
    with open(path_to_file, "r") as p:
        for line in p:
            if "PRIMER_ERROR=Cannot open " in line:
                problem.append(line)
            key = line.split("=")[0]
            value = line.strip().split("=")[1]
            if not ("MIN" in key or "MAX" in key or "OPT" in key):
                parseSeqId(key, value)
                parseTemplate(key, value)
                countPrimer(key, value)
                parseRightPrimer(key, value)
                parseLeftPrimer(key, value)
                parseInternalProbe(key, value)
                parsePrimerPair(key, value)

    if len(problem) > 0:
        os.remove(path_to_file)
        for item in problem:
            errors.append([self.target, item])
        msg = "Error during primer design "
        print(msg)
        G.logger(msg)
        G.logger(problem)
        return 1
    return 0