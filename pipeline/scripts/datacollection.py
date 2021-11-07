#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import wget
import json
import shutil
import urllib
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio import Entrez
from ipywidgets import widgets
from scripts.configuration import errors
from scripts.configuration import RunConfig
from scripts.configuration  import PipelineStatsCollector
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H

# paths
script_dir = os.path.dirname(os.path.abspath(__file__))
pipe_dir, tail = os.path.split(script_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')

Entrez.tool = "SpeciesPrimer pipeline"

class GenomeDownload(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def get_taxid(self, target):
        Entrez.email = H.get_email_for_Entrez()
        taxid, syn = H.check_input(target, Entrez.email)
        return syn, taxid

    def prepare_dirs(self):
        G.create_directory(self.target_dir)
        G.create_directory(self.config_dir)
        G.create_directory(self.genomic_dir)

    def create_GI_list(self):
        G.logger("Run: create_GI_list(" + self.target + ")")
        removed_gis = []
        with open(os.path.join(dict_path, "no_blast.gi"), "r") as f:
            for line in f:
                if "#" not in line:
                    gi = line.strip()
                    if gi not in removed_gis:
                        removed_gis.append(gi)
        if len(removed_gis) > 0:
            with open(
                os.path.join(self.config_dir, "no_blast.gi"), "w"
            ) as f:
                for gi in removed_gis:
                    f.write(gi + "\n")

    def collect_genomedata(self, taxid, email, maxrecords=2000):
        # 02/12/18 overwrite file because assembly level can change
        Entrez.email = email
        assembly_search = Entrez.esearch(
            db="assembly",
            term="txid" + str(taxid) + "[Orgn]",
            retmax=maxrecords)

        assembly_record = Entrez.read(assembly_search, validate=False)
        uidlist = assembly_record["IdList"]
        assembly_efetch = Entrez.efetch(
            db="assembly",
            id=uidlist,
            rettype="docsum",
            retmode="xml")

        assembly_records = Entrez.read(assembly_efetch, validate=False)
        metadata = self.get_genome_infos(assembly_records)
        return metadata

    def get_genome_infos(self, records):
        genome_data = []
        keys = [
            "AssemblyAccession", 'AssemblyName', "AssemblyStatus",
            'FromType', 'RefSeq_category', 'Taxid', 'Organism',
            'SpeciesTaxid', 'SpeciesName', 'FtpPath_GenBank',
            'FtpPath_RefSeq', 'ExclFromRefSeq']

        for result in records['DocumentSummarySet']['DocumentSummary']:
            data = []
            for k in keys:
                try:
                    value = result[k]
                    if value == []:
                        value = ""
                except KeyError:
                    value = "na"
                data.append(value)
                strain = "unknown"
                infraspecieslist = result['Biosource']['InfraspeciesList']
                for i, item in enumerate(infraspecieslist):
                    if len(item) > 0:
                        strain = infraspecieslist[i]['Sub_value']

            data.insert(2, strain)
            genome_data.append(data)
        keys.insert(2, "Strain")
        df = pd.DataFrame(genome_data, columns=keys)
        metadatafile = (
            os.path.join(self.genomedata_dir, "genomes_metadata.csv"))
        df.to_csv(metadatafile, index=False)

        return df

    def select_assemblies(self, df):
        df = df.replace("", np.nan)
        if self.config.assemblylevel == "offline":
            return pd.DataFrame()
        if self.config.assemblylevel != ["all"]:
            df = df[df['AssemblyStatus'].str.lower().isin(
                [x for x in self.config.assemblylevel])]
        if self.config.genbank:
            df.loc[df["FtpPath_RefSeq"].isna(), "FtpPath_RefSeq"] = df["FtpPath_GenBank"]
        if not self.config.ignore_qc:
            df.loc[:, "Accession"] = df.iloc[:, 0].str.split(".").str[0]
            df.loc[:, "version"] = df.iloc[:, 0].str.split(".").str[1]
            df = df.astype({'version': 'int32'})
            df.sort_values("version", inplace=True, ascending=False)
            df.drop_duplicates("Accession")

        df = df[~df["FtpPath_RefSeq"].isna()]
        return df.sort_index()

    def get_links(self, df):
        if df.empty is False:
            links = df["FtpPath_RefSeq"].to_list()
            urls = [x + "/" + x.split("/")[-1] + "_genomic.fna.gz" for x in links]
            outfile = os.path.join(self.genomedata_dir, "genomic_links.txt")
            with open(outfile, "w") as f:
                for url in urls:
                    f.write(url + "\n")
        else:
            urls = []

        return urls

    def get_ncbi_links(self, taxid, maxrecords=2000):
        G.logger("Run: get_ncbi_links(" + self.target + ")")
        email = H.get_email_for_Entrez()
        metadata = self.collect_genomedata(taxid, email)
        link_df = self.select_assemblies(metadata)
        link_list = self.get_links(link_df)
        statmsg = "genome assemblies from NCBI: " + str(len(link_list))
        if self.config.assemblylevel == ["offline"]:
            statmsg = "genome assemblies from NCBI: 0 (offline/skip download)"
        G.comm_log(statmsg, self.output)
        PipelineStatsCollector(self.config).write_stat(statmsg)
        return statmsg

    def download_is_required(self, filename):
        # check if the genome is already in genomic_fna
        ex_genomic_dir = os.path.join(self.ex_dir, "genomic_fna")
        gdirs = [self.genomic_dir, ex_genomic_dir]
        for gdir in gdirs:
            if os.path.isdir(gdir):
                accession = H.accession_from_filename(filename, version=False)
                for genomic_file in os.listdir(gdir):
                    if genomic_file.startswith(accession):
                        return False

        # check if the genome has been annotated
        accession = H.accession_from_filename(filename, version=True)
        dirs = [self.fna_dir, self.gff_dir, self.ffn_dir]
        exdirs = [os.path.join(self.ex_dir, os.path.split(d)[1]) for d in dirs]
        for testdir in [dirs, exdirs]:
            testlist = []
            for d in testdir:
                if os.path.isdir(d):
                    for annot_file in os.listdir(d):
                        if annot_file.startswith(accession):
                            testlist.append(1)
                if testlist == [1, 1, 1]:
                    return False

        return True

    def download_genomes(self, filename, URL, tries=3):
        for i in range(0, tries):
            G.comm_log("Download..." + filename, newline=True)
            try:
                wget.download(URL)
                break
            except urllib.error.HTTPError:
                G.comm_log("Retry download..." + filename)
                if i == tries-1:
                    error_msg = (
                                "> SpeciesPrimer in unable to "
                                "connect to the NCBI FTP server. "
                                "Please check internet connection "
                                "and NCBI FTP server status")
                    G.comm_log(error_msg)
                    errors.append([self.target, error_msg])

    def ncbi_download(self):
        G.logger("Run: ncbi_download(" + self.target + ")")
        G.create_directory(self.gff_dir)
        G.create_directory(self.ffn_dir)
        G.create_directory(self.fna_dir)
        os.chdir(self.genomic_dir)
        with open(os.path.join(self.genomedata_dir, "genomic_links.txt")) as f:
            URLS = f.readlines()
            total = len(URLS)
            for i, URL in enumerate(URLS):
                self.progress.value = (i+1)/total
                URL = URL.strip()
                filename = Path(URL).parts[-1]
                if self.download_is_required(filename):
                    self.download_genomes(filename, URL, tries=3)
        for files in os.listdir(self.genomic_dir):
            if files.endswith(".gz"):
                G.run_subprocess(["gunzip", files], False, True, False)
        os.chdir(self.target_dir)


    def add_synonym_exceptions(self, syn):
        for item in syn:
            if item not in self.config.exception:
                self.config.exception.append(item)
        conffile = os.path.join(self.config_dir, "config.json")
        with open(conffile) as f:
            for line in f:
                config_dict = json.loads(line)
        config_dict.update({"exception": self.config.exception})
        with open(conffile, "w") as f:
            f.write(json.dumps(config_dict))
        return self.config.exception

    def main(self):
        with self.output:
            G.logger("Run: collect data(" + self.target + ")")
            self.prepare_dirs()
            pan = os.path.join(self.pangenome_dir, "gene_presence_absence.csv")
            if os.path.isfile(pan):
                G.comm_log("> Found pangenome directory, skip data collection")
                self.progress.value = 1.0
                return 0

            if not self.config.offline:
                syn, taxid = self.get_taxid(self.target)
                if syn:
                    self.add_synonym_exceptions(syn)

                self.get_ncbi_links(taxid)
                if not self.config.skip_download:
                    self.ncbi_download()

            G.create_directory(self.gff_dir)
            G.create_directory(self.ffn_dir)
            G.create_directory(self.fna_dir)
            for files in os.listdir(self.genomic_dir):
                if files.endswith(".gz"):
                    filepath = os.path.join(self.genomic_dir, files)
                    G.run_subprocess(
                        ["gunzip", filepath], False, True, False)
            os.chdir(self.target_dir)
            self.create_GI_list()
            return 0

class Annotation(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)
        self.annot_report_path = Path(self.reports_dir, "annotation_report.csv")

    def count_contigs(self, filepath):
        contigs = list(SeqIO.parse(filepath, "fasta"))
        if len(contigs) >= self.contiglimit:
            from scripts.qualitycontrol import QualityControl
            dir_tree, filename = os.path.split(filepath)
            QualityControl(self.config).move_to_excluded(dir_tree, filename)
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
        G.create_directory(self.annotation_dir)
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

    def main(self):
        with self.output:
            G.comm_log("Start Annotation")
            if os.path.isdir(self.pangenome_dir):
                G.comm_log("> Found pangenome directory, skip Annotation ")
                self.progress.value = 1.0
                return 0
            else:
                qc_files = []
                max_contigs = []
                annotations = []
                accessions = []
                genomic_files = os.listdir(self.genomic_dir)
                total = len(genomic_files)
                for i, filename in enumerate(genomic_files):
                    self.progress.value = (i+1)/total
                    accession = H.accession_from_filename(filename)
                    filepath = Path(self.genomic_dir, filename)
                    accessions.append(accession)
                    qc_files.append(filename)
                    max_contigs.append(self.count_contigs(filepath))
                    annotations.append(self.annotation(filepath, accession))

                annotation_dict = {
                    "filename": qc_files,
                    "Max_contigs": max_contigs,
                    "Annotated": annotations}

                annotation_df = pd.DataFrame(
                    data = annotation_dict,
                    index=accessions)

                G.create_directory(self.reports_dir)
                annotation_df.to_csv(self.annot_report_path)

                return 0
