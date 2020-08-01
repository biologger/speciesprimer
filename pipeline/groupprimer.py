#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import csv
import time
import shutil
import itertools
from itertools import islice
import multiprocessing
from datetime import timedelta
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
from basicfunctions import GeneralFunctions as G
from basicfunctions import ParallelFunctions as P
from basicfunctions import HelperFunctions as H
from speciesprimer import CoreGenes
from speciesprimer import CoreGeneSequences
from speciesprimer import BlastPrep
from speciesprimer import Blast
from speciesprimer import BlastParser
from speciesprimer import PrimerDesign
from speciesprimer import PrimerQualityControl
from speciesprimer import PipelineStatsCollector
from speciesprimer import Summary


class Subgroups(CoreGenes):
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.ffn_dir = os.path.join(self.target_dir, "ffn_files")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "subgroups")
        self.subgroup_dir = self.results_dir
        self.fasta_dir = os.path.join(self.results_dir, "fasta")
        self.singlecopy = os.path.join(
                self.subgroup_dir, "subgroup_genes.csv")
        self.ffn_seqs = os.path.join(self.pangenome_dir, "ffn_sequences.csv")
        self.blast_dir = os.path.join(self.subgroup_dir, "blast")
        self.subgroup_seqs = []
        self.group_dict = {}
        self.subgroup = self.config.subgroup

    def get_subgroup_genes(self, mode):
        G.create_directory(self.subgroup_dir)
        filepath = os.path.join(
                self.pangenome_dir, "gene_presence_absence.csv")
        header = [
            "Gene", "Non-unique Gene name", "Annotation", "No. isolates", 
            "No. sequences", "Avg sequences per isolate", "Genome Fragment", 
            "Order within Fragment", "Accessory Fragment",
            "Accessory Order with Fragment", "QC", 
            "Min group size nuc", "Max group size nuc", "Avg group size nuc"]

        col_selection = header + self.subgroup
        df = pd.read_csv(filepath, dtype=str)      
        df = df.filter(regex=('|'.join(col_selection)))
        df = df.astype({
                "No. isolates": int, "No. sequences": int, 
                "Avg sequences per isolate": float }, copy=False)
        col = df.columns
        accession = col[len(header)::]
        df.dropna(axis=0, inplace=True, subset=accession)
        df = df.loc[df['Avg sequences per isolate'] == 1]
        col = df.columns
        accession = list(col[len(header)::])
        export = ["Gene"] + accession
        df = df.loc[:, export]
#        df["Gene"] = df.Gene.replace({"'":""}, regex=True)
#        df["Gene"] = df.Gene.replace({"/":"-"}, regex=True)
#        df["Gene"] = df.Gene.replace({" ":"_"}, regex=True)
        
        groups = df[df.Gene.str.contains("group_")]
        genenames = df[~df.Gene.str.contains("_")]
        print(groups.shape)
        print(genenames.shape)
        df = pd.concat([genenames, groups], copy=False)
        print(df.shape)
        df.set_index(["Gene"], inplace=True)

        if df.empty:
            print("No subgroup genes were identified")
        else:
            if mode == "normal":
                df.to_csv(self.singlecopy, header=False)
        self.print_gene_stats(len(df.index))
        return len(df.index)

    def print_gene_stats(self, seq_count):
        all_genes = (
            "Continue with " + str(seq_count)
            + " single copy core genes")
        G.logger("> " + all_genes)
        print("\n" + all_genes)
        stats = PipelineStatsCollector(self.target_dir)
        stats.write_stat("single copy core genes: " + str(seq_count))

    def coregene_extract(self, mode="normal"):
        info = "Run: core_gene_extract(" + self.target + ")"
        print(info)
        G.logger(info)
        self.get_subgroup_genes(mode)      
        if mode == "normal":
            locustags = self.get_sequences_from_ffn()
            self.get_fasta(locustags)

class SubgroupCoreGeneSequences(CoreGeneSequences):
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "subgroups")
        self.fasta_dir = os.path.join(self.results_dir, "fasta")
        self.alignments_dir = os.path.join(self.results_dir, "alignments")
        self.consensus_dir = os.path.join(self.results_dir, "consensus")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.conserved_dict = {}    


class SubgroupBlastParser(BlastParser):
    def __init__(self, configuration, results="seqs"):
        self.exception = configuration.exception
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "subgroups")
        self.p3dict = {}
        self.subgroup_dir = self.results_dir
        self.nontargetlist = configuration.nontargetlist
        self.selected = []
        self.mode = "normal"
        self.start = time.time()
        self.blast_dir = os.path.join(self.subgroup_dir, "blast")
        if results == "primer":
            self.mode = "primer"
            self.primer_dir = os.path.join(self.results_dir, "primer")
            self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
            self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
            self.maxgroupsize = 25000


class SubgroupPrimerDesign(PrimerDesign):
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "subgroups")
        self.p3dict = {}
        self.subgroup_dir = self.results_dir
        self.blast_dir = os.path.join(self.subgroup_dir, "blast")
        self.primer_dir = os.path.join(self.subgroup_dir, "primer")


class SubgroupPrimerQualityControl(PrimerQualityControl):
    def __init__(self, configuration, primer3_dict):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.subgroup_dir = os.path.join(self.results_dir, "subgroups")
        self.blast_dir = os.path.join(self.subgroup_dir, "blast")
        self.primer_dir = os.path.join(self.subgroup_dir, "primer")
        self.summ_dir = os.path.join(self.config.path, "Summary", self.target)
        self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
        self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
        self.mfold_dir = os.path.join(self.primer_dir, "mfold")
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.primer3_dict = primer3_dict
        self.call_blastparser = SubgroupBlastParser(self.config, "primer")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.primerlist = []
        self.start = time.time()
        self.mfethreshold = self.config.mfethreshold
        self.referencegenomes = 10
        self.dbinputfiles = []
        self.groupdbfiles = []

    def get_primerinfo(self, selected_seqs, mode):
        G.logger("Run: get_primerinfo(" + self.target + ")")
        val_list = []
        for item in selected_seqs:
            try:
                if len(item) == 2:
                    item = item[0]
                if (item.endswith("_F") or item.endswith("_R")):
                    primer_name = "_".join(item.split("_")[0:-1])
                else:
                    primer_name = item

                shortname = primer_name.split(H.abbrev(self.target) + "_")[1]
                target_id = "_".join(shortname.split("_")[0:-1])
                primerpair = "Primer_pair_" + primer_name.split("_P")[-1]
                template_seq = self.primer3_dict[target_id]["template_seq"]
                annotation = self.primer3_dict[target_id]["annotation"]
                x = self.primer3_dict[target_id][primerpair]
                pp_penalty = round(x["primer_P_penalty"], 2)
                pp_prodsize = x["product_size"]
                pp_prodTM = round(x["product_TM"], 2)
                amp_seq = x["amplicon_seq"]
                lseq = x["primer_L_sequence"]
                rseq = x["primer_R_sequence"]
                lpen = round(x["primer_L_penalty"], 2)
                rpen = round(x["primer_R_penalty"], 2)
                lTM = round(x["primer_L_TM"], 2)
                rTM = round(x["primer_R_TM"], 2)
                if mode == "mfeprimer":
                    info = [
                        primer_name + "_F", lseq, primer_name + "_R", rseq,
                        template_seq]
                    if info not in val_list:
                        val_list.append(info)
                if mode == "mfold":
                    info = [target_id, primerpair, amp_seq, primer_name]
                    if info not in val_list:
                        val_list.append(info)
                if mode == "dimercheck":
                    info = [
                        primer_name, lseq, rseq]
                    if info not in val_list:
                        val_list.append(info)
                if mode == "results":
                    ppc = x['PPC']
                    if self.config.probe:
                        iseq = x["primer_I_sequence"]
                        ipen = round(x["primer_I_penalty"], 2)
                        iTM = round(x["primer_I_TM"], 2)
                    else:
                        iseq = "None"
                        ipen = "None"
                        iTM = "None"

                    info = [
                        primer_name, ppc, pp_penalty, annotation,
                        lseq, lTM, lpen,
                        rseq, rTM, rpen,
                        iseq, iTM, ipen,
                        pp_prodsize, pp_prodTM, amp_seq,
                        template_seq]
                    if info not in val_list:
                        val_list.append(info)

            except Exception:
                G.logger("error in get_primerinfo()" + str(sys.exc_info()))

        return val_list


    def create_assembly_db_files(self):
        strain_dir = os.path.join(self.primer_qc_dir, "straindbs")
        for files in os.listdir(self.fna_dir):
            DIR = "_".join(files.split("_")[0:-1])
            targetdir = os.path.join(strain_dir, DIR)
            G.create_directory(targetdir)
            tofile = os.path.join(targetdir, "nontarget_strain.fas")
            fromfile = os.path.join(self.fna_dir, files)
            if not os.path.isfile(tofile):
                shutil.copy(fromfile, tofile)
            self.groupdbfiles.append(tofile)


    def prepare_MFEprimer_Dbs(self, primerinfos):
        from speciesprimer import errors
        G.logger("Run: prepare_MFEprimer_Dbs(" + self.target + ")")
        G.create_directory(self.primer_qc_dir)
        self.create_template_db_file(primerinfos)
        self.create_assembly_db_files()
        templatefilepath = os.path.join(
                self.primer_qc_dir, "template.sequences")
        dblist = [templatefilepath] + self.groupdbfiles
        for db in self.dbinputfiles:
            dblist.append(db)

        # parallelization try
        pool = multiprocessing.Pool()
        results = [
            pool.apply_async(P.index_database, args=(inputfilepath,))
            for inputfilepath in dblist]
        output = [p.get() for p in results]
        for item in output:
            if item:
                errors.append([self.target, item])

    def MFEprimer_QC(self, primerinfos):
        # option: also allow user provided non-target database created with
        # MFEprimer for primer QC
        G.logger("Run: MFEprimer_QC(" + self.target + ")")
        info_msg = "Start primer quality control(" + self.target + ")"
        print(info_msg)
        G.logger("> " + info_msg)
        os.chdir(self.primer_qc_dir)

        info0 = str(len(primerinfos)) + " primer pair(s) to check"
        print("\n" + info0 + "\n")
        G.logger("> " + info0)

        print("Start MFEprimer with template DB\n")
        G.logger("> Start MFEprimer with template DB")
        template_list = G.run_parallel(
                P.MFEprimer_template, primerinfos,
                [self.primer_qc_dir, self.mfethreshold])
        check_nontarget = self.write_MFEprimer_results(
                template_list, "template")
        msg = " primer pair(s) with good target binding"
        info1 = str(len(check_nontarget)) + msg
        print("\n\n" + info1 + "\n")
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs with good target binding: "
            + str(len(check_nontarget)))
        G.logger("> " + info1)

        nontarget_lists = []
        print("\nStart MFEprimer with nontarget DB\n")
        G.logger("> Start MFEprimer with nontarget DB")
        for index, dbfile in enumerate(self.dbinputfiles):
            info_msg = (
                "nontarget DB " + str(index+1) + "/"
                + str(len(self.dbinputfiles)))
            print(info_msg)
            G.logger(info_msg)
            nontarget_list = G.run_parallel(
                P.MFEprimer_nontarget, check_nontarget,
                [dbfile, self.primer_qc_dir])

            nontarget_lists = list(
                itertools.chain(nontarget_lists, nontarget_list))

        # if the MFEprimer_nontarget.csv has only the table header
        # and no results, then no primer binding was detected
        # and the primers passed the QC
        check_assembly = self.write_MFEprimer_results(
                                                nontarget_lists, "nontarget")
        msg = " primer pair(s) passed non-target PCR check"
        info2 = str(len(check_assembly)) + msg
        print("\n\n" + info2 + "\n")
        G.logger("> " + info2)
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after non-target QC: "
            + str(len(check_assembly)))

        assembly_lists = []
        print("\nStart MFEprimer with assembly DB\n")
        G.logger("> Start MFEprimer with assembly DB")


        for index, dbfile in enumerate(self.groupdbfiles):
            info_msg = (
                "strain assemblies DB " + str(index+1) + "/"
                + str(len(self.groupdbfiles)))
            print(info_msg)
            G.logger(info_msg)
            short = H.abbrev(self.target) + "_"

            assembly_list = G.run_parallel(
                    P.MFEprimer_subgroup, check_assembly,
                    [self.primer_qc_dir, dbfile, self.mfethreshold, short])

            assembly_lists = list(
                    itertools.chain(assembly_lists, assembly_list))

        check_final = self.write_MFEprimer_results(assembly_lists, "assembly")
        msg = " primer pair(s) passed secondary PCR amplicon check\n"
        info3 = str(len(check_final)) + msg
        print("\n\n" + info3 + "\n")
        G.logger("> " + info3)
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after secondary amplicon QC: "
            + str(len(check_final)))
        os.chdir(self.primer_dir)

        # new 07.11.2018 add PPC to results file
        for item in template_list:
            nameF = item[0][0]
            if nameF is not None:
                pname = "_".join(nameF.split("_")[0:-1])
                ppc = item[0][5] + float(self.mfethreshold)
                target_id = "_".join(pname.split("_")[-3:-1])
                shortname = pname.split(H.abbrev(self.target) + "_")[1]
                target_id = "_".join(shortname.split("_")[0:-1])
                primerpair = "Primer_pair_" + pname.split("_P")[-1]
                self.primer3_dict[target_id][primerpair].update({"PPC": ppc})

        primername_list = []
        for primerinfo in check_final:
            primername = "_".join(primerinfo[0].split("_")[0:-1])
            primername_list.append(primername)

        return primername_list

    def prep_mfold(self, mfoldinput, abbr):
        target_id, primerpair, amplicon_seq, primer_name = mfoldinput
        # This removes the Genus species string to shorten the
        # name for mfold (especially for subspecies names). mfold has
        # problems with too long filenames / paths
        short_name = primer_name.split(abbr + "_")[1]
        shortername = "_".join(short_name.split("_")[-2::])
        dir_path = os.path.join(self.mfold_dir, target_id)
        subdir_path = os.path.join(dir_path, primerpair)
        pcr_name = shortername + "_PCR"
        if len(amplicon_seq) >= self.config.minsize:
            G.create_directory(dir_path)
            G.create_directory(subdir_path)
            self.run_mfold(subdir_path, pcr_name, amplicon_seq)

    def read_files(self, filename):
        results = []
        pnum_dir = os.path.split(os.path.dirname(filename))
        pp = pnum_dir[1].split("_")[-1]
        p_dir = os.path.split(pnum_dir[0])
        name = p_dir[1] + "_P" + pp
        primername = H.abbrev(self.target) + "_" + name
        with open(filename, "r", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if re.search("Structure", line):
                    name = "".join(islice(f, 1, 2)).strip()
                    mfoldvalues = self.parse_values(f, line)
                    results.append(
                        self.interpret_values(
                            name, primername, mfoldvalues))
        return results

    def write_results(self, choice):
        G.logger("Run: write_results(" + self.target + ")")
        results = []
        header = [
            "Primer name", "PPC", "Primer penalty", "Annotation",
            "Primer fwd seq", "Primer fwd TM", "Primer fwd penalty",
            "Primer rev seq", "Primer rev TM", "Primer rev penalty",
            "Probe seq)", "Probe TM", "Probe penalty",
            "Amplicon size", "Amplicon TM", "Amplicon sequence",
            "Template sequence"]
        if len(choice) > 0:
            results = self.get_primerinfo(choice, mode="results")
            results.sort(key=lambda x: float(x[1]), reverse=True)
            file_path = os.path.join(
                    self.subgroup_dir,
                    H.abbrev(self.target) + "_primer.csv")
            G.csv_writer(file_path, results, header)

        return results

    def run_primer_qc(self):
        G.logger("Run: run_primer_qc(" + self.target + ")")
        total_results = []
        if self.collect_primer() == 0:
            G.create_directory(self.primerblast_dir)
            blastsum = os.path.join(self.primerblast_dir, "nontargethits.json")
            if not os.path.isfile(blastsum):
                prep = BlastPrep(
                    self.primerblast_dir, self.primerlist,
                    "primer", self.config.blastseqs)
                use_cores, inputseqs = prep.run_blastprep()

                bla = Blast(self.config, self.primerblast_dir, "primer")
                bla.run_blast("primer", use_cores)

            self.call_blastparser.run_blastparser("primer")

            primer_qc_list = self.get_primerinfo(self.primerlist, "mfeprimer")

            for files in os.listdir(self.primer_qc_dir):
                if (
                    files.startswith("BLASTnontarget")
                    and files.endswith(".sequences")
                ):
                    filepath = os.path.join(self.primer_qc_dir, files)
                    self.dbinputfiles.append(filepath)

            self.prepare_MFEprimer_Dbs(primer_qc_list)

            survived_MFEp = self.MFEprimer_QC(primer_qc_list)

            mfoldinput = self.get_primerinfo(survived_MFEp, "mfold")

            self.mfold_analysis(mfoldinput)

            selected_primer, excluded_primer = self.mfold_parser()

            dimercheck = self.dimercheck_primer(
                selected_primer, excluded_primer)

            choice = self.check_primerdimer(dimercheck)

            total_results = self.write_results(choice)

            if total_results != []:
                info = (
                    "Found " + str(len(total_results))
                    + " primer pair(s) for " + self.target)
                print("\n" + info + "\n")
                G.logger("> " + info)
                duration = time.time() - self.start
                G.logger(
                    "> PrimerQC time: "
                    + str(timedelta(seconds=duration)).split(".")[0])
                return total_results

        from speciesprimer import errors
        error_msg = "No compatible primers found"
        duration = time.time() - self.start
        G.logger(
            "> PrimerQC time: "
            + str(timedelta(seconds=duration)).split(".")[0])
        print(error_msg)
        G.logger("> " + error_msg)
        errors.append([self.target, error_msg])

        return total_results


class SubgroupSummary(Summary):
    def __init__(self, configuration, total_results):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "subgroup")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.primer_dir = os.path.join(self.results_dir, "primer")
        self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
        self.mfold_dir = os.path.join(self.primer_dir, "mfold")
        self.summ_dir = os.path.join(
                self.config.path, "Summary", self.target, "strain_specific")
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.aka = H.abbrev(self.target)
        self.g_info_dict = {}
        self.total_results = total_results


def main(config):
    Subgroups(config).run_CoreGenes()
    conserved_seq_dict = SubgroupCoreGeneSequences(config).run_coregeneanalysis()
#    SI.coregene_extract()
#    group_dict = SI.run_groupseqs()
    
#    print(group_dict)
    
#    str_unique = subgroupBlastParser(
#            config).run_blastparser(group_dict)
#    primer_dict = subgroupPrimerDesign(config).run_primerdesign()
#    total_results = subgroupPrimerQualityControl(config, primer_dict).run_primer_qc()
#    subgroupSummary(config, total_results).run_summary(mode="last")
    print(conserved_seq_dict.keys())
    
    if not conserved_seq_dict == 1:
        conserved = SubgroupBlastParser(
                config).run_blastparser(conserved_seq_dict)
        if conserved == 0:
            primer_dict = SubgroupPrimerDesign(config).run_primerdesign()
            total_results = SubgroupPrimerQualityControl(
                config, primer_dict).run_primer_qc()
            SubgroupSummary(config, total_results).run_summary(mode="last")
        else:
            SubgroupSummary(config, total_results).run_summary(mode="last")
    else:
        SubgroupSummary(config, total_results).run_summary(mode="last")


if __name__ == "__main__":
    print(
        "Start this script with speciesprimer.py and the "
        "'--runmode strain' option")
