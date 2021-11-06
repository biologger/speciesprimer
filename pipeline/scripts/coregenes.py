#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import csv
import shutil
import multiprocessing
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from ipywidgets import widgets
from scripts.configuration import RunConfig
from scripts.configuration import errors
from scripts.configuration import PipelineStatsCollector
from scripts.blastscripts import Blast
from scripts.blastscripts import BlastPrep
from scripts.blastscripts import BlastParser
import basicfunctions
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H

# paths
script_dir = os.path.dirname(os.path.abspath(__file__))
pipe_dir, tail = os.path.split(script_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')
if not pipe_dir in sys.path:
    sys.path.append(pipe_dir)

class PangenomeAnalysis(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def get_numberofgenomes(self):
        inputgenomes = len(os.listdir(self.gff_dir))
        msg = ["genome assemblies for pan-genome analysis:", str(inputgenomes)]
        PipelineStatsCollector(self.config).write_stat(" ".join(msg))
        return inputgenomes

    def run_roary(self):
        num_cpus = str(multiprocessing.cpu_count())
        G.logger("Run: run_roary(" + self.target + ")")
        roary_cmd = [
            "roary", "-f", self.pangenome_dir, "-s", "-p", num_cpus,
            "-cd", "100", "-e", "-n", str(Path(self.gff_dir, "*.gff"))]
        if self.config.skip_tree:
            if "-e" in roary_cmd:
                roary_cmd.remove("-e")
            if "-n" in roary_cmd:
                roary_cmd.remove("-n")
        if self.config.virus:
           roary_cmd.insert(-1, "-t")
           roary_cmd.insert(-1, str(1))
        try:
            G.run_shell(
                " ".join(roary_cmd), printcmd=True, logcmd=True, log=True)
        except (KeyboardInterrupt, SystemExit):
            G.keyexit_rollback("pan-genome analysis", dp=self.pangenome_dir)
            raise

        return roary_cmd

    def run_fasttree(self):
        G.logger("Run: run_fasttree(" + self.target + ")")
        os.chdir(self.pangenome_dir)
        coregenealn = "core_gene_alignment.aln"
        if os.path.isfile(coregenealn):
            tree = H.abbrev(self.target) + "_tree.nwk"
            fasttree_cmd = [
                "FastTreeMP", "-nt", "-gtr", "-nopr", coregenealn, ">", tree]
            try:
                G.run_shell(
                    " ".join(fasttree_cmd), printcmd=True, logcmd=True,
                    log=True)
            except (KeyboardInterrupt, SystemExit):
                G.keyexit_rollback(
                    "fasttree run", dp=self.pangenome_dir, fn=tree)
                raise
        os.chdir(self.target_dir)

    def main(self):
        with self.output:
            G.comm_log("Run: run_pangenome_analysis(" + self.target + ")")
            if os.path.isdir(self.pangenome_dir):
                filepath = os.path.join(
                    self.pangenome_dir, "gene_presence_absence.csv")
                if os.path.isfile(filepath):
                    info = (
                        "Pangenome directory already exists\n"
                        "Continue with existing Pangenome data")
                    G.comm_log("> " + info)
                    self.progress.value = 1.0
                    return 2
                shutil.rmtree(self.pangenome_dir)

            genomes = self.get_numberofgenomes()
            if len(genomes) == 0:
                error_msg = "Error: No genomes survived QC"
                G.comm_log("> " + error_msg)
                errors.append([self.target, error_msg])
                # add some error code to add tipps and tricks messages
                return 1

            if len(genomes) < 2:
                error_msg = "Error: Less than two genomes survived QC"
                G.comm_log("> " + error_msg)
                errors.append([self.target, error_msg])
                # add some error code to add tipps and tricks messages
                return 1

            self.progress.value = 0.1
            self.run_roary()
            self.progress.value = 0.6
            self.run_fasttree()
            self.progress.value = 1.0
            return 0


class CoreGenes(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def get_singlecopy_genes(self):
        filepath = os.path.join(self.pangenome_dir, "gene_presence_absence.csv")
        data = pd.read_csv(filepath)
        genomes = len(data.iloc[:, 14:].columns)
        data = data.query('`No. isolates` == @genomes')
        data = data.query('`No. sequences` == @genomes')
        groupmask = data["Gene"].str.contains("group")
        groupdata = data[groupmask]
        namedata = data[~groupmask]
        uniquenamedata = namedata[~namedata["Gene"].str.contains("_")]
        coregenes = pd.concat([uniquenamedata, groupdata])
        coregenes = coregenes.set_index("Gene")
        coregenes = coregenes.iloc[:, 13:]
        coregenes.to_csv(self.singlecopy, header=False)
        stats = PipelineStatsCollector(self.config)
        stats.write_stat("core genes: " + str(len(data.index)))
        stats.write_stat("single copy core genes: " + str(len(coregenes)))
        G.comm_log("core genes: " + str(len(data.index)))
        G.comm_log("single copy core genes: " + str(len(coregenes)))
        return coregenes

    def get_sequences_from_ffn(self, coregenes):
        if not os.path.isfile(self.ffn_seqs):
            ffn_seqs = []
            accessions = coregenes.columns.to_list()
            total = len(accessions)
            for i, acc in enumerate(accessions):
                ffn_path = Path(self.ffn_dir, acc + ".ffn")
                loci = coregenes[acc].to_list()
                self.progress.value = i+1/total
                with open(ffn_path) as f:
                    records = SeqIO.parse(f, "fasta")
                    for rec in records:
                        if rec.id in loci:
                            ffn_seqs.append([acc, rec.id, str(rec.seq)])

            ffn_df = pd.DataFrame(ffn_seqs)
            ffn_df.to_csv(self.ffn_seqs, index=False, header=False)
        else:
            ffn_df = pd.read_csv(self.ffn_seqs, header=None)
            self.progress.value = 1.0

        return ffn_df

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

    def get_fasta(self, ffn_df):
        ffn_df.columns = ["Accession", "Locus", "Sequence"]
        locustags = ffn_df.set_index(["Locus"])
        with open(self.singlecopy, "r") as f:
            lines = f.readlines()
            total = len(lines)
            reader = csv.reader(lines)
            for i, row in enumerate(reader):
                gene = self.check_genename(row[0])
                outfile = os.path.join(self.fasta_dir, gene + ".fasta")
                self.progress.value = i+1/total
                if not os.path.isdir(outfile):
                    with open(outfile, "w") as r:
                        for item in row[1:]:
                            name = locustags.loc[item, "Accession"]
                            seq = locustags.loc[item, "Sequence"]
                            record = SeqRecord(
                                Seq(seq), name=item,
                                id='{}|{}|{}'.format(name, item, gene),
                                description="")
                            SeqIO.write(record, r, "fasta")

    def main(self):
        with self.output:
            G.comm_log("> Collect results of pan-genome analysis")
            G.create_directory(self.fasta_dir)
            coregenes = self.get_singlecopy_genes()
            if len(coregenes.index) == 0:
                error_msg = "No single copy core genes found"
                G.comm_log(error_msg)
                errors.append([self.target, error_msg])
                return 1
            ffn_df = self.get_sequences_from_ffn(coregenes)
            self.get_fasta(ffn_df)
            return 0


class CoreGeneAlignments(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def run_prank(self):
        G.comm_log("> Start alignment of core gene sequences")
        G.create_directory(self.alignments_dir)
        inputfiles = [f.split(".fasta")[0] for f in os.listdir(self.fasta_dir)]
        outputfiles = [f.split(".best.fas")[0] for f in os.listdir(self.alignments_dir)]
        todo = list(set(inputfiles) - set(outputfiles))
        if len(todo) == 0:
            G.comm_log("> Found " + str(len(outputfiles)) + " alignments")
            self.progress.value = 1.0
            return 0
        cmds = [
            ["prank", "-d=" + os.path.join(self.fasta_dir, f + ".fasta"),
            "-o=" + os.path.join(self.alignments_dir, f)] for f in todo]
        try:
            G.run_parallel(basicfunctions.GeneralFunctions.read_shelloutput, cmds, self.progress)
        except (KeyboardInterrupt, SystemExit):
            G.keyexit_rollback(
                "Prank MSA run", dp=self.alignments_dir)
            raise

    def main(self):
        with self.output:
            G.logger("Run: CoreGeneAlignments(" + self.target + ")")
            self.run_prank()
            return 0

class CoreGeneSequences(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def run_consesus_sequences(self):
        G.comm_log("> Find consensus sequence for aligned core gene sequences")
        G.create_directory(self.consensus_dir)
        inputfiles = [
            f.split(".best.fas")[0] for f in os.listdir(self.alignments_dir)]
        outputfiles = [
            f.split(
                "_consensus.fas")[0] for f in os.listdir(self.consensus_dir)]
        todo = list(set(inputfiles) - set(outputfiles))
        if len(todo) == 0:
            self.progress.value = 1.0
            return 0
        cmds = [
            ["consambig", "-sequence",
             os.path.join(self.alignments_dir,f + ".best.fas"),
            "-outseq",
            os.path.join(self.consensus_dir, f + "_consensus.fas"), "-name",
             f + "_consensus", "-auto"] for f in todo]
        try:
            G.run_parallel(basicfunctions.GeneralFunctions.read_shelloutput, cmds, self.progress)
        except (KeyboardInterrupt, SystemExit):
            G.keyexit_rollback(
                "consensus run", dp=self.consensus_dir)
            raise

    def split_conserved_sequences(self, record):
        split_records = []
        # replace any non-GATC character with N
        sequence = re.sub(
            "[^GATC]", "N", re.sub("[a-z]", "N", str(record.seq)))
        split_seq = re.sub(
            "N{5,50}", "*",
            re.sub("N.{1,20}N", "*", str(sequence)))
        count = 1
        split_list = split_seq.split("*")
        for seq in split_list:
            if len(seq) >= self.config.minsize:
                desc = record.id
                if "group_" in desc:
                    seq_name = (
                        "group_" + desc.split("_")[-2:-1][0]
                        + "_" + str(count))
                else:
                    seq_name = (
                        desc.split("_")[-2:-1][0] + "_" + str(count))
                count += 1
                rec = SeqRecord(Seq(seq), id=seq_name, description="")
                split_records.append(rec)
        return split_records

    def find_conserved_sequences(self):
        conserved_recs = []
        G.comm_log("> Search conserved regions in consensus sequences")
        files = [
            os.path.join(self.consensus_dir, f)
            for f in os.listdir(self.consensus_dir)]
        for filepath in files:
            for record in SeqIO.parse(filepath, "fasta"):
                split_records = self.split_conserved_sequences(record)
                conserved_recs.extend(split_records)

        G.create_directory(self.blast_dir)
        f = os.path.join(self.blast_dir, H.abbrev(self.target) + "_conserved_seqs.fas")
        SeqIO.write(conserved_recs, f, "fasta")
        info = "Number of conserved sequences: " + str(len(conserved_recs))
        PipelineStatsCollector(self.config).write_stat(info)
        G.comm_log("> " + info)
        if len(conserved_recs) == 0:
            error_msg = "Error: no conserved target sequences found"
            G.comm_log(error_msg)
            errors.append([self.target, error_msg])
            return 1
        return 0

    def main(self):
        with self.output:
            G.logger("Run: run_coregeneanalysis(" + self.target + ")")
            self.run_consesus_sequences()
            exitstatus = self.find_conserved_sequences()
            return exitstatus

class ConservedBlast(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def run_blastsearch(self, conserved_seqs):
        use_cores = BlastPrep(
             self.blast_dir, conserved_seqs, "conserved",
             self.config.blastseqs).run_blastprep()
        blast = Blast(
             self.config, self.blast_dir, "conserved"
         )
        self.progress = blast.progress
        self.output = blast.output
        blast.run_blast("conserved", use_cores)

    def main(self):
        with self.output:
            blastsum = os.path.join(self.blast_dir, "BLAST_results_summary.csv")
            if not os.path.isfile(blastsum):
                f = os.path.join(self.blast_dir, H.abbrev(self.target) + "_conserved_seqs.fas")
                conserved_seqs = list(SeqIO.parse(f, "fasta"))
                self.run_blastsearch(conserved_seqs)

            exitstatus = BlastParser(self.config).run_blastparser()
            self.progress.value = 1.0
            return exitstatus
