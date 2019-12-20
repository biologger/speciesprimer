#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import time
import itertools
import multiprocessing
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



class Singletons(CoreGenes):
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
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
            for row in reader:
                data_row = []
                gene_name = row[0]
                number_isolates = int(row[3])
                number_sequences = int(row[4])
                loci = row[14:]
                if number_isolates == 1:
                    singleton_count.append(gene_name)
                    data_row.append(gene_name)
                    for locus in loci:
                        if not locus == "":
#                            print(locus)
                            data_row.append(locus)
                    newtabledata.append(data_row)


        G.csv_writer(self.singleton, newtabledata)
        return len(singleton_count)

    def get_fasta(self, locustags):

        def check_genename(gene):
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
        G.create_directory(self.blast_dir)
        with open(self.singleton, "r") as f:
            reader = csv.reader(f)
            outfile = os.path.join(self.blast_dir, "singleton_sequences.fas")
            with open(outfile, "w") as r:
                for row in reader:
                    gene = check_genename(row[0])
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

    def prepare_MFEprimer_Dbs(self, primerinfos):
        from speciesprimer import errors
        G.logger("Run: prepare_MFEprimer_Dbs(" + self.target + ")")
        G.create_directory(self.primer_qc_dir)
        self.create_template_db_file(primerinfos)
        self.create_assembly_db_file()
        assemblyfilepath = os.path.join(
            self.primer_qc_dir,
            H.abbrev(self.target) + ".genomic")
        templatefilepath = os.path.join(
                self.primer_qc_dir, "template.sequences")
        dblist = [assemblyfilepath, templatefilepath]
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

        print("\nStart MFEprimer with assembly DB\n")
        G.logger("> Start MFEprimer with assembly DB")

        dbfile = H.abbrev(self.target) + ".genomic"
        assembly_list = G.run_parallel(
                P.MFEprimer_singleton, check_assembly,
                [self.primer_qc_dir, dbfile, self.mfethreshold])

        check_final = self.write_MFEprimer_results(assembly_list, "assembly")
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
                primerpair = "Primer_pair_" + pname.split("_P")[-1]
                self.primer3_dict[target_id][primerpair].update({"PPC": ppc})

        primername_list = []
        for primerinfo in check_final:
            primername = "_".join(primerinfo[0].split("_")[0:-1])
            primername_list.append(primername)

        return primername_list


class SingletonSummary:
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
        self.summ_dir = os.path.join(self.config.path, "Summary", self.target)
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.aka = H.abbrev(self.target)
        self.g_info_dict = {}
        if total_results is None:
            self.total_results = []
        else:
            self.total_results = total_results


def main(config):
    SI = Singletons(config)
    SI.coregene_extract()
    single_dict = SI.run_singleseqs()
    str_unique = SingletonBlastParser(
            config).run_blastparser(single_dict)
    primer_dict = SingletonPrimerDesign(config).run_primerdesign()
    total_results = SingletonPrimerQualityControl(config, primer_dict).run_primer_qc()
#    SingletonSummary(config, total_results).run_summary(mode="last")


if __name__ == "__main__":
    main()
