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
from basicfunctions import GeneralFunctions as G
from basicfunctions import ParallelFunctions as P
from basicfunctions import HelperFunctions as H
from speciesprimer import CoreGenes
from speciesprimer import BlastPrep
from speciesprimer import Blast
from speciesprimer import BlastParser
from speciesprimer import PrimerDesign
from speciesprimer import PrimerQualityControl
from speciesprimer import PipelineStatsCollector
from speciesprimer import Summary



class Singletons(CoreGenes):
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.ffn_dir = os.path.join(self.target_dir, "ffn_files")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.single_dir = os.path.join(self.results_dir, "singletons")
        self.singleton = os.path.join(
                self.single_dir, "singleton_genes.csv")
        self.ffn_seqs = os.path.join(self.pangenome_dir, "ffn_sequences.csv")
        self.blast_dir = os.path.join(self.single_dir, "blast")
        self.singleton_seqs = []
        self.single_dict = {}
        self.strains = self.config.strains

    def get_singleton_genes(self):
        G.create_directory(self.single_dir)
        singleton_count = []
        filepath = os.path.join(
                self.pangenome_dir, "gene_presence_absence.csv")
        newtabledata = []
        with open(filepath, "r") as f:
            reader = csv.reader(f)
            header = next(reader)
            accessions = header[14:]
            genomes = len(accessions)
            rows = []
            if self.strains:
                for s in self.strains:
                    for index, item in enumerate(header):
                        if "_".join(item.split("_")[0:-1]) == s:
                            rows.append(index)
            for row in reader:
                data_row = []
                gene_name = row[0]
                number_isolates = int(row[3])
                number_sequences = int(row[4])
                if self.strains:
                    loci = []
                    for item in rows:
                        loci.append(row[item])
                else:
                    loci = row[14:]
                if number_isolates == 1:
                    singleton_count.append(gene_name)
                    data_row.append(gene_name)
                    for locus in loci:
                        if not locus == "":
                            data_row.append(locus)
                    if len(data_row) > 1:
                        newtabledata.append(data_row)

        if len(newtabledata) > 0:
            G.csv_writer(self.singleton, newtabledata)
        else:
            print("No singleton genes were identified")
            # Add warning
        return len(singleton_count)

    def check_genename(self, gene):
        if "/" in gene:
            gene_name = "-".join(gene.split("/"))
        elif " " in gene:
            gene_name = "-".join(gene.split(" "))
        elif "'" in gene:
            gene_name = str(gene)[0:-1]
            gene = gene + "'"
        else:
            gene_name = gene

        return gene_name

    def get_fasta(self, locustags):
        G.create_directory(self.blast_dir)
        with open(self.singleton, "r") as f:
            reader = csv.reader(f)
            outfile = os.path.join(self.blast_dir, "singleton_sequences.fas")
            with open(outfile, "w") as r:
                for row in reader:
                    gene = self.check_genename(row[0])
                    items = row[1]
                    if "\t" in items:
                        item_s = items.split("\t")
                    else:
                        item_s = [items]

                    for item in item_s:

                        name = "_".join(locustags[item]["name"].split("_")[0:-1])
                        seq = locustags[item]["seq"]
                        ident = '{}_{}'.format(name, gene)
                        record = SeqRecord(
                            Seq(seq),
                            name=item,
                            id=ident,
                            description="")
                        SeqIO.write(record, r, "fasta")
                        self.single_dict.update({ident: seq})
                        self.singleton_seqs.append([ident, seq])

    def coregene_extract(self):
        info = "Run: core_gene_extract(" + self.target + ")"
        print(info)
        G.logger(info)
        self.get_singleton_genes()
        locustags = self.get_sequences_from_ffn()
        self.get_fasta(locustags)

    def run_singleseqs(self):
        name = "singleton"
        blastsum = os.path.join(self.blast_dir, "nontargethits.json")
        if not os.path.isfile(blastsum):
            use_cores, inputseqs = BlastPrep(
                self.blast_dir, self.singleton_seqs, name,
                self.config.blastseqs).run_blastprep()
            Blast(
                self.config, self.blast_dir, "conserved"
            ).run_blast(name, use_cores)
        return self.single_dict


class SingletonBlastParser(BlastParser):
    def __init__(self, configuration, results="seqs"):
        self.exception = configuration.exception
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "singletons")
        self.p3dict = {}
        self.single_dir = self.results_dir
        self.nontargetlist = configuration.nontargetlist
        self.selected = []
        self.mode = "normal"
        self.start = time.time()
        self.blast_dir = os.path.join(self.single_dir, "blast")
        if results == "primer":
            self.mode = "primer"
            self.primer_dir = os.path.join(self.results_dir, "primer")
            self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
            self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
            self.maxgroupsize = 25000


class SingletonPrimerDesign(PrimerDesign):
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "singletons")
        self.p3dict = {}
        self.single_dir = self.results_dir
        self.blast_dir = os.path.join(self.single_dir, "blast")
        self.primer_dir = os.path.join(self.single_dir, "primer")


class SingletonPrimerQualityControl(PrimerQualityControl):
    def __init__(self, configuration, primer3_dict):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.single_dir = os.path.join(self.results_dir, "singletons")
        self.blast_dir = os.path.join(self.single_dir, "blast")
        self.primer_dir = os.path.join(self.single_dir, "primer")
        self.summ_dir = os.path.join(self.config.path, "Summary", self.target)
        self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
        self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
        self.mfold_dir = os.path.join(self.primer_dir, "mfold")
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.primer3_dict = primer3_dict
        self.call_blastparser = SingletonBlastParser(self.config, "primer")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.primerlist = []
        self.start = time.time()
        self.mfethreshold = self.config.mfethreshold
        self.referencegenomes = 10
        self.dbinputfiles = []
        self.singledbfiles = []

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
            self.singledbfiles.append(tofile)


    def prepare_MFEprimer_Dbs(self, primerinfos):
        from speciesprimer import errors
        G.logger("Run: prepare_MFEprimer_Dbs(" + self.target + ")")
        G.create_directory(self.primer_qc_dir)
        self.create_template_db_file(primerinfos)
        self.create_assembly_db_files()
        templatefilepath = os.path.join(
                self.primer_qc_dir, "template.sequences")
        dblist = [templatefilepath] + self.singledbfiles
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


        for index, dbfile in enumerate(self.singledbfiles):
            info_msg = (
                "strain assemblies DB " + str(index+1) + "/"
                + str(len(self.singledbfiles)))
            print(info_msg)
            G.logger(info_msg)
            short = H.abbrev(self.target) + "_"

            assembly_list = G.run_parallel(
                    P.MFEprimer_singleton, check_assembly,
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


class SingletonSummary(Summary):
    def __init__(self, configuration, total_results):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(
                self.pangenome_dir, "results", "singleton")
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
    SI = Singletons(config)
    SI.coregene_extract()
    single_dict = SI.run_singleseqs()
    str_unique = SingletonBlastParser(
            config).run_blastparser(single_dict)
    primer_dict = SingletonPrimerDesign(config).run_primerdesign()
    total_results = SingletonPrimerQualityControl(config, primer_dict).run_primer_qc()
    SingletonSummary(config, total_results).run_summary(mode="last")


if __name__ == "__main__":
    print(
        "Start this script with speciesprimer.py and the "
        "'--runmode strain' option")
