#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:17:48 2019

@author: ags-bioinfo
"""

# mfe3parser

def MFEprimer_template(self, primerinfo):
    num_amp = None
    result = []
    [nameF, seqF, nameR, seqR, temp_seq, amp_seq] = primerinfo
    pp_name = "_".join(nameF.split("_")[0:-1])
    with tempfile.NamedTemporaryFile(
        mode='w+', dir=self.primer_qc_dir, prefix="primer",
        suffix=".fa", delete=False
    ) as primefile:
        primefile.write(
            ">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")
    db = "template.sequences"
    cmd = [
        "mfeprimer3.1", "spec", "-i", primefile.name, "--misMatch",
        str(self.mismatches), "-s", "50", "-d", db ]
    cmd = " ".join(cmd)
    result = G.read_shelloutput(
        cmd, printcmd=False, logcmd=False, printoption=False)
    os.unlink(primefile.name)
    datalist = []
    for index, item in enumerate(result):
        try:
            if "Descriptions of " in item:
                num_amp = int(item.split(" ")[3])
            elif "Amp " in item:
                ampsize = result[index+1].split(" ")[2]
            elif ">Amp_" in item:
                amplist = item.split(" ")
                amp_c = int(amplist[0].split("_")[1])
                fp = amplist[1]
                rp = amplist[3]
                target = amplist[5]
                seq = []
                for x in range(index, len(result)):
                    if not (
                        "Parameters" in result[x+1] or "Amp " in result[x+1]):
                        seq.append(result[x+1])
                    else:
                        ampseq = "".join(seq)
                        break
                datalist.append([amp_c, fp, rp, target, ampsize, ampseq])
        except IndexError:
            self.MFEprimer_template(primerinfo)

    if datalist == []:
        if num_amp == 0:
            return [[pp_name], [["no amplicon"]]]
        else:
            return self.MFEprimer_template(primerinfo)
    elif num_amp == 1:
        return [primerinfo, datalist]
    else:
        return [[pp_name], datalist]
    
    
def MFEprimer_template(self, primerinfo):
    num_amp = None
    result = []
    [nameF, seqF, nameR, seqR, temp_seq, amp_seq] = primerinfo
    pp_name = "_".join(nameF.split("_")[0:-1])
    with tempfile.NamedTemporaryFile(
        mode='w+', dir=self.primer_qc_dir, prefix="primer",
        suffix=".fa", delete=False
    ) as primefile:
        primefile.write(
            ">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")
    db = "template.sequences"
    cmd = [
        "mfeprimer3.1", "spec", "-i", primefile.name, "--misMatch",
        str(self.mismatches), "-s", "50", "-d", db ]
    cmd = " ".join(cmd)
    result = G.read_shelloutput(
        cmd, printcmd=False, logcmd=False, printoption=False)
    os.unlink(primefile.name)
    datalist = []
    for index, item in enumerate(result):
        try:
            if "Descriptions of " in item:
                num_amp = int(item.split(" ")[3])
            elif "Amp " in item:
                ampsize = result[index+1].split(" ")[2]
            elif ">Amp_" in item:
                amplist = item.split(" ")
                amp_c = int(amplist[0].split("_")[1])
                fp = amplist[1]
                rp = amplist[3]
                target = amplist[5]
                seq = []
                for x in range(index, len(result)):
                    if not (
                        "Parameters" in result[x+1] or "Amp " in result[x+1]):
                        seq.append(result[x+1])
                    else:
                        ampseq = "".join(seq)
                        break
                datalist.append([amp_c, fp, rp, target, ampsize, ampseq])
        except IndexError:
            self.MFEprimer_template(primerinfo)

    if datalist == []:
        if num_amp == 0:
            return [[pp_name], [["no amplicon"]]]
        else:
            return self.MFEprimer_template(primerinfo)
    elif num_amp == 1:
        return [primerinfo, datalist]
    else:
        return [[pp_name], datalist]