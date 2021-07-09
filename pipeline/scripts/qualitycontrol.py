#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import multiprocessing
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from speciesprimer import RunConfig
from speciesprimer import Blast
from speciesprimer import BlastParser
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H

class QualityControl(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)

    def move_to_excluded(self, dir_tree, filename=None):
        if self.config.ignore_qc is False:
            if filename is None:
                from_dir = dir_tree
                head, tail = os.path.split(dir_tree)
                to_dir = Path(self.ex_dir, os.path.basename(head), tail)
            else:
                from_dir = Path(dir_tree, filename)
                sub_dir = Path(self.ex_dir, os.path.basename(dir_tree))
                to_dir = Path(sub_dir, filename)
                if not os.path.isdir(sub_dir):
                    G.create_directory(sub_dir)
            G.comm_log(
                "Moved file to excluded_genomes: " + os.path.basename(to_dir))
            shutil.move(from_dir, to_dir)

    def remove_qc_failures(self, accession):
        dirs = [self.fna_dir, self.gff_dir, self.ffn_dir, self.genomic_dir]
        outdir = str(Path(self.annotation_dir, accession))
        self.move_to_excluded(outdir, filename=None)
        for i, d in enumerate(dirs):
            if i == 3:
                accession = H.genomicversion_from_accession(accession)
            for filename in os.listdir(d):
                if filename.startswith(accession):
                    self.move_to_excluded(d, filename)

    def count_contigs(self, filepath):
        contigs = list(SeqIO.parse(filepath, "fasta"))
        if len(contigs) >= self.contiglimit:
            dir_tree, filename = os.path.split(filepath)
            self.move_to_excluded(dir_tree, filename)
            return "max contig"
        return "passed QC"

    def start_prokka(self, inputfile, outputname):
        if self.config.virus:
            genus = 'Genus'
            kingdom = "Viruses"
        else:
            genus = self.target.split("_")[0]
            kingdom = "Bacteria"
        outdir = str(Path(self.annotation_dir, outputname))
        prokka_cmd = [
            "prokka",
            "--kingdom", kingdom,
            "--outdir", outdir,
            "--genus", genus,
            "--locustag", outputname,
            "--prefix", outputname,
            "--cpus", "0",
            str(inputfile)
        ]
        G.comm_log(outputname + " annotation required")
        try:
            G.run_subprocess(prokka_cmd, True, True, False)
        except (KeyboardInterrupt, SystemExit):
            G.keyexit_rollback(
                "annotation", dp=outdir)
            raise

    def annotation(self, filepath, accession):
        fmts = ["fna", "gff", "ffn"]
        dirs = [self.fna_dir, self.gff_dir, self.ffn_dir ]
        testlist = []
        for d in dirs:
            for filename in os.listdir(d):
                if filename.startswith(accession):
                    testlist.append(1)
        if testlist != [1, 1, 1]:
            if os.path.isfile(filepath):
                self.start_prokka(filepath, accession)
                outdir = str(Path(self.annotation_dir, accession))
                for filename in os.listdir(outdir):
                    for i, fmt in enumerate(fmts):
                        if filename.endswith(fmt):
                            fromdir = Path(outdir, filename)
                            todir = Path(dirs[i], filename)
                            shutil.copy(fromdir, todir)
            else:
                return "failed"
        return "annotated"

    def find_qcgene_annotations(self, gff_filepath, qc_gene):
        searchdict = {
            "rRNA": "product=16S ribosomal RNA",
            "tuf": "product=Elongation factor Tu",
            "dnaK": "product=Chaperone protein DnaK",
            "recA": "product=Protein RecA",
            "pheS": "product=Phenylalanine--tRNA ligase alpha subunit"}
        annotations = []
        with open(gff_filepath) as f:
            for line in f:
                if searchdict[qc_gene] in line:
                    gene = line.split("ID=")[1].split(";")[0].split(" ")[0]
                    annotations.append(gene)

        dir_tree, gff_filepath = os.path.split(gff_filepath)
        accession = gff_filepath.split(".")[0]

        if len(annotations) == 0:
            G.comm_log("No sequence(s) for "+ qc_gene +" found in " + str(gff_filepath))
            self.remove_qc_failures(accession)
            record = SeqRecord(Seq(""), id="failed", name="failed", description="failed")
        else:
            ffn_p = Path(self.ffn_dir, accession + ".ffn")
            recs = list(SeqIO.parse(ffn_p, "fasta"))
            qc_rec_data = [[rec.id, str(rec.seq)] for rec in recs if rec.id in annotations]
            max_len = sorted(qc_rec_data, key=lambda x: len(x[1]), reverse=True)[0]
            record = SeqRecord(Seq(max_len[1]), id=max_len[0], name=max_len[0], description=qc_gene)
        return record

    def qc_blast(self, qc_gene):
        qc_dir = os.path.join(self.genomedata_dir, qc_gene + "_QC")
        cores = multiprocessing.cpu_count()
        Blast(self.config, qc_dir, "quality_control").run_blast(qc_gene, cores)
        qc_df = BlastParser(self.config, results="quality_control").bp_parse_results(qc_dir)
        fp = os.path.join(qc_dir, qc_gene + "_QC_report.csv")
        qc_df.to_csv(fp, index=False)
        return qc_df

    def remove_offtarget_genomes(self, qc_report):
        offtarget = qc_report.query('`QC status` == "failed QC"')
        gene_locus = offtarget["Query Seq-id"].to_list()
        for locus in gene_locus:
            accession = "_".join(locus.split("_")[0:-1])
            self.remove_qc_failures(accession)
        G.comm_log("> Quality control removed " + str(len(gene_locus)) + " genomes")
        return offtarget

    def collect_qc_infos(self):
        annotation_df = pd.read_csv(
            Path(self.genomedata_dir, "annotation_report.csv"), index_col=0)
        maxcontigs = annotation_df["Max_contigs"]
        maxcontigs.index = [
            H.genomicversion_from_accession(x) for x in maxcontigs.index.to_list()]
        meta_data = pd.read_csv(
            Path(self.genomedata_dir, "genomes_metadata.csv"), index_col=0)
        qc_meta = meta_data[['AssemblyName', 'Strain', 'AssemblyStatus']]
        reports = [maxcontigs]
        for qc_gene in self.config.qc_gene:
            qc_dir = os.path.join(self.genomedata_dir, qc_gene + "_QC")
            fp = os.path.join(qc_dir, qc_gene + "_QC_report.csv")
            qc_df = pd.read_csv(fp, index_col=0)
            qc_data = qc_df[["QC status", "Species"]]
            qc_data.index = [
            H.genomicversion_from_accession(x) for x in qc_data.index.to_list()]
            qc_data.columns = [qc_gene, qc_gene + " Blast"]
            reports.append(qc_data)
        qc_infos = qc_meta.join(reports, how="outer")
        qc_infos.to_csv(Path(self.target_dir, "QC_report.csv"))

    def quality_control(self):
        qc_files = []
        max_contigs = []
        annotations = []
        accessions = []
        for filename in os.listdir(self.genomic_dir):
            accession = H.accession_from_filename(filename)
            filepath = Path(self.genomic_dir, filename)
            accessions.append(accession)
            qc_files.append(filename)
            max_contigs.append(self.count_contigs(filepath))
            annotations.append(self.annotation(filepath, accession))

        annotation_dict = {"filename": qc_files, "Max_contigs": max_contigs, "Annotated": annotations}

        for qc_gene in self.config.qc_gene:
            qc_sequences = []
            qc_dir = os.path.join(self.genomedata_dir, qc_gene + "_QC")
            G.create_directory(qc_dir)
            for acc in accessions:
                gff_filepath = Path(self.gff_dir, acc + ".gff")
                if os.path.isfile(gff_filepath):
                    gene_locus = self.find_qcgene_annotations(gff_filepath, qc_gene)
                    if gene_locus.id != "failed":
                        qc_sequences.append(gene_locus)

            qcblast_input = Path(qc_dir, qc_gene + ".part-0")
            SeqIO.write(qc_sequences, qcblast_input, "fasta")
            qc_df = self.qc_blast(qc_gene)
            self.remove_offtarget_genomes(qc_df)

        annotation_df = pd.DataFrame(
            data = annotation_dict,
            index=accessions)
        report_path = Path(self.genomedata_dir, "annotation_report.csv")
        annotation_df.to_csv(report_path)
        self.collect_qc_infos()
