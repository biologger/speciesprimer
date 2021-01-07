#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import time
import logging
import csv
import fnmatch
import signal
import re
import shutil
import multiprocessing
import wget
import json
import urllib
import pandas as pd
import numpy as np
import itertools
from itertools import islice
from datetime import timedelta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio import Entrez
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H
from basicfunctions import ParallelFunctions as P
from basicfunctions import BlastDBError

# paths
pipe_dir = os.path.dirname(os.path.abspath(__file__))
dict_path = os.path.join(pipe_dir, "dictionaries")
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')
errors = []

# general info
systemdirs = [
    "genomic_fna", "config", "ffn_files", "gff_files", "Pangenome",
    "rRNA_QC", "recA_QC", "tuf_QC",
    "dnaK_QC", "pheS_QC"]

Entrez.tool = "SpeciesPrimer pipeline"


class Config:
    def __init__(self, mode="man", config_dict=None):
        if mode == "man":
            self.config_dict = self.generate_config_files()
        else:
            if config_dict:
                self.config_dict = config_dict

    def generate_config_files(self):
        import batchassist
        config_dict = batchassist.Input().initiate()
        return config_dict

    def get_targets(self):
        targets = []
        for key in self.config_dict:
            targets.append(key)
        return targets

    def get_config(self, target):
        minsize = self.config_dict[target]["minsize"]
        maxsize = self.config_dict[target]["maxsize"]
        mpprimer = self.config_dict[target]["mpprimer"]
        exception = self.config_dict[target]["exception"]
        path = self.config_dict[target]["path"]
        intermediate = self.config_dict[target]["intermediate"]
        qc_gene = self.config_dict[target]["qc_gene"]
        mfold = self.config_dict[target]["mfold"]
        skip_download = self.config_dict[target]["skip_download"]
        assemblylevel = self.config_dict[target]["assemblylevel"]
        skip_tree = self.config_dict[target]["skip_tree"]
        nolist = self.config_dict[target]["nolist"]
        offline = self.config_dict[target]["offline"]
        ignore_qc = self.config_dict[target]["ignore_qc"]
        mfethreshold = self.config_dict[target]["mfethreshold"]
        customdb = self.config_dict[target]["customdb"]
        blastseqs = self.config_dict[target]["blastseqs"]
        probe = self.config_dict[target]["probe"]
        virus = self.config_dict[target]["virus"]
        genbank = self.config_dict[target]["genbank"]
        evalue = self.config_dict[target]["evalue"]
        nuc_identity = self.config_dict[target]["nuc_identity"]
        runmode = self.config_dict[target]["runmode"]
        strains = self.config_dict[target]["strains"]

        return (
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, skip_tree, nolist, offline, ignore_qc, mfethreshold,
            customdb, blastseqs, probe, virus, genbank,
            evalue, nuc_identity, runmode, strains)


class CLIconf:
    def __init__(
            self, minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold,
            skip_download, assemblylevel,
            nontargetlist, skip_tree, nolist, offline, ignore_qc, mfethreshold,
            customdb, blastseqs, probe, virus, genbank,
            evalue, nuc_identity, runmode, strains):
        self.minsize = minsize
        self.maxsize = maxsize
        self.mpprimer = mpprimer
        self.exception = exception
        self.target = target
        self.path = path
        self.intermediate = intermediate
        self.qc_gene = qc_gene
        self.mfold = mfold
        self.skip_download = skip_download
        self.assemblylevel = assemblylevel
        self.nontargetlist = nontargetlist
        self.skip_tree = skip_tree
        self.nolist = nolist
        self.offline = offline
        self.ignore_qc = ignore_qc
        self.mfethreshold = mfethreshold
        self.customdb = customdb
        self.blastseqs = blastseqs
        self.probe = probe
        self.virus = virus
        self.genbank = genbank
        self.evalue = evalue
        self.nuc_identity = nuc_identity
        self.runmode = runmode
        self.strains = strains
        self.save_config()

    def save_config(self):
        config_dict = {}
        config_dict.update({"minsize": self.minsize})
        config_dict.update({"maxsize": self.maxsize})
        config_dict.update({"mpprimer": self.mpprimer})
        config_dict.update({"exception": self.exception})
        config_dict.update({"target": self.target})
        config_dict.update({"path": self.path})
        config_dict.update({"intermediate": self.intermediate})
        config_dict.update({"qc_gene": self.qc_gene})
        config_dict.update({"mfold": self.mfold})
        config_dict.update({"skip_download": self.skip_download})
        config_dict.update({"assemblylevel": self.assemblylevel})
        config_dict.update({"skip_tree": self.skip_tree})
        config_dict.update({"nolist": self.nolist})
        config_dict.update({"offline": self.offline})
        config_dict.update({"ignore_qc": self.ignore_qc})
        config_dict.update({"mfethreshold": self.mfethreshold})
        config_dict.update({"customdb": self.customdb})
        config_dict.update({"blastseqs": self.blastseqs})
        config_dict.update({"probe": self.probe})
        config_dict.update({"virus": self.virus})
        config_dict.update({"genbank": self.genbank})
        config_dict.update({"evalue": self.evalue})
        config_dict.update({"nuc_identity": self.nuc_identity})
        config_dict.update({"runmode": self.runmode})
        config_dict.update({"strains": self.strains})

        dir_path = os.path.join(self.path, self.target)
        config_path = os.path.join(self.path, self.target, "config")
        file_path = os.path.join(config_path, "config.json")
        G.create_directory(dir_path)
        G.create_directory(config_path)
        # To do
        # Add warning if a previous file is overwritten
        with open(file_path, "w") as f:
            f.write(json.dumps(config_dict))


class DataCollection():
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.exception = self.config.exception
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.genomic_dir = os.path.join(self.target_dir, "genomic_fna")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.ffn_dir = os.path.join(self.target_dir, "ffn_files")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.contiglimit = 500
        self.ex_dir = os.path.join(
            self.config.path, "excludedassemblies", self.target)

    def get_taxid(self, target):
        Entrez.email = H.get_email_for_Entrez()
        taxid = H.check_input(target, Entrez.email)
        if taxid:
            syn = H.check_species_syn(taxid, Entrez.email, target)
        else:
            syn = None

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
                if "#" in line:
                    pass
                else:
                    gi = line.strip()
                    if gi not in removed_gis:
                        removed_gis.append(gi)
        if len(removed_gis) > 0:
            with open(
                os.path.join(self.config_dir, "no_blast.gi"), "w"
            ) as f:
                for gi in removed_gis:
                    f.write(gi + "\n")

    def get_ncbi_links(self, taxid, maxrecords=2000):

        def collect_genomedata(taxid, email):
            genomedata = []
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

            with open("genomicdata.json", "w") as f:
                f.write(json.dumps(assembly_records))

            with open("genomicdata.json", "r") as f:
                for line in f:
                    gen_dict = json.loads(line)

            for assembly in gen_dict['DocumentSummarySet']['DocumentSummary']:
                accession = assembly["AssemblyAccession"]
                name = assembly["AssemblyName"]
                status = assembly["AssemblyStatus"]
                ftp_path = assembly["FtpPath_RefSeq"]

                if ftp_path == "":
                    if self.config.genbank is True:
                        ftp_path = assembly["FtpPath_GenBank"]
                        if ftp_path == "":
                            ftp_link = "None"
                        else:
                            ftp_link = (
                                ftp_path + "/" + ftp_path.split("/")[-1]
                                + "_genomic.fna.gz"
                            )
                    else:
                        ftp_link = "None"
                else:
                    ftp_link = (
                        ftp_path + "/" + ftp_path.split("/")[-1]
                        + "_genomic.fna.gz"
                    )
                data = [accession, name, status, ftp_link]
                genomedata.append(data)

            return genomedata

        def get_links(linkdata):
            statusdict = {
                "complete": "Complete Genome", "chromosome": "Chromosome",
                "scaffold": "Scaffold", "contig": "Contig",
                "offline": "Offline"}
            links = []
            if self.config.assemblylevel == ["all"]:
                for item in linkdata:
                    if not item[3] == "None":
                        links.append(item[3])
            else:
                for status in self.config.assemblylevel:
                    for item in linkdata:
                        if item[2] == statusdict[status]:
                            if not item[3] == "None":
                                links.append(item[3])

            return links

        def write_links(link_list):
            with open("genomic_links.txt", "w") as f:
                for link in link_list:
                    f.write(link + "\n")

        G.logger("Run: ncbi_genome_links(" + self.target + ")")
        os.chdir(self.config_dir)
        email = H.get_email_for_Entrez()
        genomedata = collect_genomedata(taxid, email)
        time.sleep(1)
        link_list = get_links(genomedata)
        msg = " genome assemblies are available for download"
        statmsg = "genome assemblies from NCBI: " + str(len(link_list))
        if self.config.assemblylevel == ["offline"]:
            msg = "Skip download / Working offline"
            statmsg = "genome assemblies from NCBI: 0 (offline/skip download)"
            info = msg
        else:
            info = str(len(link_list)) + msg
        print(info)
        G.logger("> " + info)
        PipelineStatsCollector(self.target_dir).write_stat(statmsg)

        write_links(link_list)
        os.chdir(self.target_dir)
        return info

    def check_download_files(self, input_line):
        line = input_line.strip()
        zip_file = line.split("/")[-1]
        name_start = (
            zip_file.split(".")[0]
            + "v" + zip_file.split(".")[1].split("_")[0])
        directories = [self.gff_dir, self.ffn_dir]
        gff = False
        ffn = False
        genomic = False
        for directory in directories:
            for files in [f for f in os.listdir(directory)]:
                if files.startswith(name_start):
                    if directory == self.gff_dir:
                        if files.endswith(".gff"):
                            gff = True
                    elif directory == self.ffn_dir:
                        if files.endswith(".ffn"):
                            ffn = True
        for files in [f for f in os.listdir(self.genomic_dir)]:
            if files == zip_file.split(".gz")[0]:
                genomic = True
        if (gff and ffn) is True:
            status = True
        elif genomic is True:
            status = "Extracted"
        else:
            status = False
        return status

    def get_excluded_assemblies(self):
        excluded = []
        # a list of excluded Genomes to keep information
        # after directories are (manually) removed (to save diskspace)
        if os.path.isdir(self.ex_dir):
            filepath = os.path.join(self.ex_dir, "excluded_list.txt")
            if os.path.isfile(filepath):
                with open(filepath, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line not in excluded:
                            excluded.append(line)
        return excluded

    def ncbi_download(self):
        G.logger("Run: ncbi_download(" + self.target + ")")
        G.create_directory(self.gff_dir)
        G.create_directory(self.ffn_dir)
        G.create_directory(self.fna_dir)
        excluded = self.get_excluded_assemblies()
        os.chdir(self.genomic_dir)
        with open(os.path.join(self.config_dir, "genomic_links.txt")) as r:
            for line in r:
                zip_file = line.split("/")[-1].strip()
                target_path = os.path.join(self.genomic_dir, zip_file)
                ftp_path = line.strip()
                file_status = self.check_download_files(line)
                if file_status is True:
                    info = "all required files found for " + zip_file
                    G.logger(info)
                elif file_status == "Extracted":
                    info = "extracted fna file already exists for " + zip_file
                    G.logger(info)
                else:
                    if os.path.isfile(target_path):
                        info = "File already downloaded " + zip_file
                        G.logger(info)
                    else:
                        if [
                                ex for ex in excluded
                                if zip_file.startswith(".".join(ex.split("v")))
                        ]:
                            msg = " already in excludedassemblies"
                            print(zip_file + msg)
                            G.logger(zip_file + msg)
                        else:
                            print_msg = "\n\nDownload..." + zip_file + "\n"
                            info = "Downloaded " + zip_file
                            try:
                                print(print_msg)
                                wget.download(ftp_path)
                                print("\n")
                                G.logger(info)
                            except urllib.error.HTTPError:
                                try:
                                    print("\ntry again\n")
                                    print(print_msg)
                                    wget.download(ftp_path)
                                    print("\n")
                                    G.logger(info)
                                except urllib.error.HTTPError:
                                    try:
                                        print("\ntry a last time\n")
                                        print(print_msg)
                                        wget.download(ftp_path)
                                        print("\n")
                                        G.logger(info)
                                    except urllib.error.HTTPError:
                                        error_msg = (
                                            "SpeciesPrimer in unable to "
                                            "connect to the NCBI FTP server. "
                                            "Please check internet connection "
                                            "and NCBI FTP server status")
                                        print(error_msg)
                                        G.logger("> " + error_msg)
                                        errors.append([self.target, error_msg])
                                        raise

        for files in os.listdir(self.genomic_dir):
            if files.endswith(".gz"):
                G.run_subprocess(["gunzip", files], False, True, False)
        os.chdir(self.target_dir)

    def copy_genome_files(self):
        G.logger("Run: copy_genome_files(" + self.target + ")")
        for root, dirs, files in os.walk(self.target_dir):
            for file_name in files:
                if file_name.endswith(".ffn"):
                    ffn_file = os.path.join(self.ffn_dir, file_name)
                    if file_name not in os.listdir(self.ffn_dir):
                        shutil.copy(os.path.join(root, file_name), ffn_file)
                if file_name.endswith(".gff"):
                    gff_file = os.path.join(self.gff_dir, file_name)
                    if file_name not in os.listdir(self.gff_dir):
                        shutil.copy(os.path.join(root, file_name), gff_file)
                if file_name.endswith(".fna"):
                    if not file_name.endswith("genomic.fna"):
                        if "genomic_fna" not in root.split("/"):
                            fna_file = os.path.join(self.fna_dir, file_name)
                            if file_name not in os.listdir(self.fna_dir):
                                shutil.copy(
                                    os.path.join(root, file_name), fna_file)

    def run_prokka(self):
        annotation_dirs = []

        def excluded_assemblies():
            qc_fail_dirs = []
            # list of failed assemblies
            if os.path.isfile(os.path.join(self.ex_dir, "excluded_list.txt")):
                with open(
                    os.path.join(self.ex_dir, "excluded_list.txt"), "r"
                ) as f:
                    for line in f:
                        line = line.strip()
                        qc_fail_dirs.append(line)
            return qc_fail_dirs

        def get_annotation_dirs():
            dirs = []
            # list of folders of annotated genomes
            for directories in [
                d for d in os.listdir(self.target_dir)
                if os.path.isdir(os.path.join(self.target_dir, d))
            ]:
                name = "_".join(directories.split("_")[:-1])
                dirs.append(name)
            return dirs

        def get_genomic_files():
            genomic_files = []
            # list of genomic files
            if os.path.isdir(self.genomic_dir):
                for file_name in os.listdir(self.genomic_dir):
                    genomic_files.append(file_name)
            return genomic_files

        def get_filenames(directory):
            filenames = []
            for files in [
                f for f in os.listdir(directory)
                if os.path.isfile(os.path.join(directory, f))
            ]:
                name = "_".join(files.split("_")[:-1])
                filenames.append(name)
            return filenames

        def check_incomplete(annot_names, filetype_names, directory):
            incomplete = set(filetype_names) - set(annot_names)
            for name in incomplete:
                for filename in os.listdir(directory):
                    if filename.startswith(name):
                        filepath = os.path.join(directory, filename)
                        os.remove(filepath)
                        info = "Removed " + filepath
                        G.logger(info)

        def start_prokka(filename, fna):
            date = time.strftime("%Y%m%d")
            if self.config.virus:
                genus = 'Genus'
                kingdom = "Viruses"
            else:
                genus = self.target.split("_")[0]
                kingdom = "Bacteria"
            outdir = filename + "_" + date
            prokka_cmd = [
                "prokka",
                "--kingdom", kingdom,
                "--outdir", outdir,
                "--genus", genus,
                "--locustag", filename,
                "--prefix", filename + "_" + date,
                "--cpus", "0",
                "genomic_fna/" + fna
            ]
            info = file_name + " annotation required"
            G.logger(info)
            print("\n" + info)
            try:
                G.run_subprocess(prokka_cmd, True, True, False)
            except (KeyboardInterrupt, SystemExit):
                G.keyexit_rollback(
                    "annotation", dp=os.path.join(self.target_dir, outdir))
                raise

            annotation_dirs.append(filename + "_" + date)

        annotated = []
        excluded = []
        G.logger("Run: run_prokka(" + self.target + ")")

        dirs = get_annotation_dirs()
        qc_fail_dirs = excluded_assemblies()
        genomic_files = get_genomic_files()
        # set intersection to identify the required annotation files
        fna_names = get_filenames(self.fna_dir)
        gff_names = get_filenames(self.gff_dir)
        ffn_names = get_filenames(self.ffn_dir)
        annotatedname = set(fna_names) & set(gff_names) & set(ffn_names)
        for name in annotatedname:
            if name not in dirs:
                dirs.append(name)

        check_incomplete(annotatedname, fna_names, self.fna_dir)
        check_incomplete(annotatedname, gff_names, self.gff_dir)
        check_incomplete(annotatedname, ffn_names, self.ffn_dir)

        # write annotation command
        for fna in genomic_files:
            if fna.endswith("_genomic.fna"):
                file_name = (
                    fna.split(".")[0] + "v" + fna.split(".")[1].split("_")[0])
            else:
                name = fna.split(".fna")[0]
                if "_" in name:
                    file_name = (
                        "_".join("-".join(name.split(".")).split("_")[0:-1]))
                else:
                    file_name = "-".join(name.split("."))

            if file_name in dirs:
                if file_name != '':
                    annotated.append(file_name)
            elif file_name in qc_fail_dirs:
                excluded.append(file_name)
            else:
                start_prokka(file_name, fna)

        if len(annotated) > 0:
            info = "Already annotated: "
            G.logger(info.strip())
            G.logger(annotated)
            print("\n" + info)
            print(annotated)

        if len(excluded) > 0:
            info = "Already in excludedassemblies:"
            G.logger("> " + info)
            G.logger(excluded)
            print("\n" + info)
            print(excluded)

        return annotation_dirs, annotated

    def remove_max_contigs(self):
        maxcontigs = []
        for files in os.listdir(self.genomic_dir):
            if files.endswith(".fna"):
                contigcount = 0
                filepath = os.path.join(self.genomic_dir, files)
                for line in open(filepath).readlines():
                    if ">" in line:
                        contigcount += 1
                if contigcount >= self.contiglimit:
                    if self.config.ignore_qc is False:
                        if files.endswith("_genomic.fna"):
                            file_name = (
                                files.split(".")[0] + "v" +
                                files.split(".")[1].split("_")[0])
                        else:
                            name = files.split(".fna")[0]
                            file_name = (
                                "_".join(
                                    "-".join(name.split(".")).split("_")[0:-1])
                            )
                        maxcontigs.append(file_name)
                        msg = (
                            files + " has more than " + str(self.contiglimit)
                            + " contigs and will be removed before annotation "
                            "to include it in the run use the ignore_qc option"
                        )
                        print(msg)
                        G.logger(msg)
                        excl_path = os.path.join(self.ex_dir, "genomic_fna")
                        if not os.path.isdir(excl_path):
                            G.create_directory(excl_path)
                        shutil.move(filepath, os.path.join(excl_path, files))

        if len(maxcontigs) > 0:
            if os.path.isfile(os.path.join(self.ex_dir, "excluded_list.txt")):
                with open(
                    os.path.join(self.ex_dir, "excluded_list.txt"), "a"
                ) as f:
                    for item in maxcontigs:
                        f.write(item + "\n")
            else:
                with open(
                    os.path.join(self.ex_dir, "excluded_list.txt"), "w"
                ) as f:
                    for item in maxcontigs:
                        f.write(item + "\n")

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

    def collect(self):
        G.logger("Run: collect data(" + self.target + ")")
        self.prepare_dirs()
        pan = os.path.join(self.pangenome_dir, "gene_presence_absence.csv")
        if os.path.isfile(pan):
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

        self.remove_max_contigs()
        self.create_GI_list()
        annotation_dirs, annotated = self.run_prokka()
        self.copy_genome_files()

        if self.config.intermediate is False:
            for directory in annotation_dirs:
                dirpath = os.path.join(self.target_dir, directory)
                if os.path.isdir(dirpath):
                    shutil.rmtree(dirpath)

        return self.config


class QualityControl:
    # dictionary containing the search word for annotated genes in gff files
    # check / update when Prokka version is updated
    searchdict = {
        "rRNA": "product=16S ribosomal RNA",
        "tuf": "product=Translation initiation factor IF-1",
        "dnaK": "product=Chaperone protein DnaK",
        "recA": "product=Protein RecA",
        "pheS": "product=Phenylalanine--tRNA ligase alpha subunit"}

    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.exception = configuration.exception
        self.ex_dir = os.path.join(
            self.config.path, "excludedassemblies", self.target)
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.config_dir = os.path.join(self.target_dir, "config")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.ffn_dir = os.path.join(self.target_dir, "ffn_files")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.qc_gene_search = []
        self.ffn_list = []
        self.contiglimit = 500
        self.no_seq = []
        self.double = []
        self.contig_ex = []
        self.problems = []
        self.passed = []

    def get_excluded_gis(self):
        excluded_gis = []
        gi_file = os.path.join(self.config_dir, "no_blast.gi")
        if os.path.isfile(gi_file):
            if os.stat(gi_file).st_size > 0:
                with open(gi_file, "r") as f:
                    for line in f:
                        gi = line.strip()
                        excluded_gis.append(str(gi))
        return excluded_gis

    def search_qc_gene(self, file_name, qc_gene):
        with open(os.path.join(self.gff_dir, file_name), "r") as f:
            for line in f:
                if self.searchdict[qc_gene] in line:
                    gene = line.split("ID=")[1].split(";")[0].split(" ")[0]
                    if gene not in self.qc_gene_search:
                        self.qc_gene_search.append(gene)

    def count_contigs(self, gff_list, contiglimit):
        exclude = []
        for dirs in os.listdir(self.target_dir):
            if dirs not in systemdirs:
                path = os.path.join(self.target_dir, dirs)
                if os.path.isdir(path):
                    for files in os.listdir(path):
                        if files.endswith(".fna"):
                            contigcount = 0
                            filepath = os.path.join(path, files)
                            file = files.split(".fna")[0]
                            for line in open(filepath).readlines():
                                if ">" in line:
                                    contigcount += 1
                            if contigcount >= contiglimit:
                                exclude.append(file)

        if len(exclude) > 0:
            for item in exclude:
                if item + ".gff" in gff_list:
                    gff_list.remove(item + ".gff")
                data = [item, "", "", "", "", "Max contigs"]
                if data not in self.contig_ex:
                    self.contig_ex.append(data)
            info = (
                "skip " + str(len(self.contig_ex))
                + " Genome(s) with more than " + str(self.contiglimit)
                + " contigs")
            print(info)
            G.logger("> " + info)

        return gff_list

    def identify_duplicates(self, gff_list):
        duplicate = []
        duplicate_test = []
        keep = []
        remove_older_version = []

        def find_potential_duplicates():
            for item in gff_list:
                name = '_'.join(item.split(".gff")[0].split("_")[0:-1])
                if (("GCA" or "GCF") and "v") in name:
                    version = name.split("_")[-1].split("v")[1]
                    common = name.split("v")[0]
                    if int(version) > 1:
                        if common not in duplicate:
                            duplicate.append(common)

        def test_if_duplicate(duplicate):
            for y in duplicate:
                del duplicate_test[:]
                for x in gff_list:
                    x = "_".join(x.split(".")[0].split("_")[0:-1])
                    if str(y) in str(x):
                        if x not in duplicate_test:
                            duplicate_test.append(x)

                if len(duplicate_test) > 0:
                    maxi = max(
                        duplicate_test,
                        key=lambda item: int(item.split("v")[1])
                    )
                    if maxi not in keep:
                        keep.append(maxi)
                    duplicate_test.remove(maxi)
                    for item in duplicate_test:
                        if item not in remove_older_version:
                            remove_older_version.append(item)

        find_potential_duplicates()
        test_if_duplicate(duplicate)

        if len(remove_older_version) > 0:
            for item in remove_older_version:
                for gff_file in gff_list:
                    if item in gff_file:
                        if gff_file in gff_list:
                            gff_list.remove(gff_file)
                        data = [
                            gff_file.split(".gff")[0],
                            "", "", "", "", "Duplicate"]
                        if data not in self.double:
                            self.double.append(data)

            info = (
                "skip " + str(len(self.double)) + " duplicate Genome(s) ")
            print(info)
            G.logger("> " + info)

        return gff_list

    # 12.02.2018 change to generate one QC file
    def check_no_sequence(self, qc_gene, gff):
        ffn_list = []
        sub_gff = []
        sub_gene_search = []
        for file_name in gff:
            name = "_".join(file_name.split("_")[:-1])
            sub_gff.append(name)
        for seq_id in self.qc_gene_search:
            seq_name = "_".join(seq_id.split("_")[:-1])
            sub_gene_search.append(seq_name)
        no_seq_found = set(sub_gff) - set(sub_gene_search)

        if len(no_seq_found) > 0:
            for item in no_seq_found:
                for file_name in gff:
                    if item in file_name:
                        gff.remove(file_name)
                        self.no_seq.append([
                            file_name.split(".gff")[0],
                            "", "", "", "", "QC gene missing"])

            info = (
                "skip " + str(len(self.no_seq)) + " Genome(s) without "
                + qc_gene + " sequence")
            print(info)
            G.logger("> " + info)

        for item in gff:
            ffn = item.split(".gff")[0] + ".ffn"
            ffn_list.append(ffn)

        return ffn_list

    def get_qc_seqs(self, qc_gene):
        G.logger("Run: get_qc_seqs(" + qc_gene + ")")
        G.logger("> Starting QC with " + qc_gene)
        print("Starting QC with " + qc_gene)
        gff = []
        qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
        G.create_directory(qc_dir)
        # find annotation of gene in gff files and store file name
        for files in os.listdir(self.gff_dir):
            if files not in gff:
                gff.append(files)
        info = "found " + str(len(gff)) + " gff files"
        G.logger(info)
        print(info)

        if len(gff) > 0:
            # look for annotations
            if self.contiglimit > 0:
                contig_gff_list = self.count_contigs(gff, self.contiglimit)
                gff_list = self.identify_duplicates(contig_gff_list)
            else:
                gff_list = self.identify_duplicates(gff)

            for item in gff_list:
                self.search_qc_gene(item, qc_gene)

            info = (
                "found " + str(len(self.qc_gene_search)) + " "
                + qc_gene + " annotations in gff files")
            G.logger(info)
            print(info)

            ffn_check = self.check_no_sequence(qc_gene, gff_list)

            # search sequences in ffn files
            for files in os.listdir(self.ffn_dir):
                if files in ffn_check:
                    if files not in self.ffn_list:
                        self.ffn_list.append(files)
            info = (
                    "selected " + str(len(self.ffn_list)) + " "
                    + qc_gene + " sequences from ffn files")
            G.logger(info)
            print(info)

        else:
            error_msg = (
                    "Error: No .gff files found for QualityControl " + qc_gene)
            print(error_msg)
            G.logger("> " + error_msg)
            errors.append([self.target, error_msg])
            return 1
        return 0

    def choose_sequence(self, qc_gene):
        """ find files and choose the longest sequence
        if a sequence was found more than once """
        G.logger("Run: choose_sequence(" + qc_gene + ")")
        qc_seqs = []
        qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
        with open(os.path.join(qc_dir, qc_gene + "_seq"), "w") as o:
            for file_name in self.ffn_list:
                recid = []
                recseq = []
                for seq_id in self.qc_gene_search:
                    match = "_".join(seq_id.split("_")[:-1])
                    if file_name.startswith(match):
                        file_path = os.path.join(self.ffn_dir, file_name)
                        for record in SeqIO.parse(file_path, "fasta"):
                            if seq_id in record.id:
                                if len(str(record.seq)) != 0:
                                    if record.id not in recid:
                                        recid.append(record.id)
                                    if str(record.seq) not in recseq:
                                        recseq.append(str(record.seq))

                # get longest sequence and write to file
                if len(recseq) != 0:
                    q = max(recseq, key=len)
                    z = dict(zip(tuple(recid), tuple(recseq)))
                    a = {v: k for k, v in z.items()}
                    o.write(">" + str(a[q]) + "\n" + str(q) + "\n")
                    qc_seqs.append([str(a[q]), str(q)])

        return qc_seqs

    def qc_blast_parser(self, qc_gene):
        qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
        G.logger("Run: qc_blast_parser(" + qc_gene + ")")

        def delete_blastreport(xmlblastresults):
            os.chdir(qc_dir)
            for file_name in xmlblastresults:
                os.remove(file_name)
                n = file_name.split("_")[-2]
                os.remove(qc_gene + ".part-" + n)
            os.remove(qc_gene + "_seq")

        def parse_blastresults():
            wrote = []
            exceptions = []
            if not self.exception == []:
                for item in self.exception:
                    exception = ' '.join(item.split("_"))
                    exceptions.append(exception)
            gi_list = []
            excluded_gis = self.get_excluded_gis()
            expected = " ".join(self.target.split("_"))
            # collect the blast results
            os.chdir(qc_dir)
            blapa = BlastParser(self.config)
            xmlblastresults = blapa.blastresult_files(qc_dir)
            for filename in xmlblastresults:
                blast_record_list = blapa.parse_BLASTfile(filename)
                for index, blast_record in enumerate(blast_record_list):
                    i = 0
                    try:
                        alignment = blast_record.alignments[i]
                        aln_data = blapa.get_alignmentdata(alignment, [])
                        spec, gi, db_id = aln_data[0], aln_data[1], aln_data[2]
                        query = blast_record.query
                        if str(gi) in excluded_gis:
                            if gi not in gi_list:
                                gi_list.append(gi)

                            while i < len(blast_record.alignments) - 1:
                                i = i+1
                                alignment = blast_record.alignments[i]
                                aln_data = blapa.get_alignmentdata(
                                                                alignment, [])
                                spec, gi, db_id = (
                                        aln_data[0], aln_data[1], aln_data[2])
                                if str(gi) not in excluded_gis:
                                    break
                                else:
                                    if gi not in gi_list:
                                        gi_list.append(gi)

                    except IndexError:
                        query = blast_record.query
                        spec, gi, db_id = "no match", "", ""

                    if expected in spec:
                        if query not in wrote:
                            wrote.append(query)
                            success = [
                                query, gi, db_id, spec,
                                expected, "passed QC"]
                            self.passed.append(success)

                    elif spec in exceptions:
                        if query not in wrote:
                            wrote.append(query)
                            success = [
                                query, gi, db_id, spec,
                                expected, "passed QC"]
                            self.passed.append(success)
                    else:
                        if query not in wrote:
                            wrote.append(query)
                            fail = [
                                query, gi, db_id, spec, expected, "failed QC"]
                            self.problems.append(fail)

            if self.config.intermediate is False:
                delete_blastreport(xmlblastresults)
            if len(gi_list) > 0:
                info = "removed GI's in excluded GI list from results"
                print(info)
                print(gi_list)
                G.logger(info)
                G.logger(gi_list)
            os.chdir(self.target_dir)

        def write_blastresults():
            results = []
            if len(self.passed) > 0:
                for item in self.passed:
                    results.append(item)
            if len(self.problems) > 0:
                for item in self.problems:
                    results.append(item)
            if len(self.no_seq) > 0:
                for item in self.no_seq:
                    results.append(item)
            if len(self.contig_ex) > 0:
                for item in self.contig_ex:
                    results.append(item)
            if len(self.double) > 0:
                for item in self.double:
                    results.append(item)
            # write files
            report = os.path.join(qc_dir, qc_gene + "_QC_report.csv")
            if len(results) > 0:
                header = [
                    "Query", "GI", "DB ID", "Species",
                    "Target species", "QC status"]
                G.csv_writer(report, results, header)
            else:
                error_msg = "No Quality Control results found."
                print(error_msg)
                G.logger("> " + error_msg)
                errors.append([self.target, error_msg])


        parse_blastresults()
        write_blastresults()
        return self.passed

    def remove_qc_failures(self, qc_gene):
        G.logger("Run: remove_qc_failures(" + qc_gene + ")")
        print("Run: remove_qc_failures(" + qc_gene + ")")
        qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
        delete = []
        with open(os.path.join(qc_dir, qc_gene + "_QC_report.csv"), "r") as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if "passed QC" not in row[5]:
                    name = '_'.join(row[0].split("_")[0:-1])
                    if name not in delete:
                        delete.append(name)

        if len(delete) > 0:
            G.create_directory(self.ex_dir)
            self.delete_failed_assemblies(delete)
            info = "Quality control removed " + str(len(delete)) + " Genome(s)"
            G.logger("> " + info)
            print("\n" + info)

        return delete

    def delete_failed_assemblies(self, delete):
        filepath = os.path.join(self.ex_dir, "excluded_list.txt")
        with open(filepath, "a+") as f:
            for item in delete:
                f.write(item + "\n")

        def move_dirs(from_dir, to_dir):
            try:
                shutil.move(from_dir, to_dir)
                info = "Moved: " + from_dir + " to " + to_dir
                G.logger(info)
            except shutil.Error:
                if os.path.isdir(from_dir):
                    shutil.rmtree(from_dir)

        def move_file(from_file, to_file):
            try:
                shutil.move(from_file, to_file)
                info = "Moved: " + from_file + " to " + to_file
                G.logger(info)
            except shutil.Error:
                pass

        def remove_from_list(file_format, file_name):
            if file_format == "ffn":
                if file_name in self.ffn_list:
                    self.ffn_list.remove(file_name)

        def delete_files(item, directory):
            search_str = str(item + "_*")
            for files in os.listdir(directory):
                if re.search(search_str, files):
                    x = '_'.join(files.split("_")[0:-1])
                    if x == item:
                        filepath = os.path.join(directory, files)
                        file_format = files.split(".")[1]
                        from_file = filepath
                        ff_dir = os.path.join(
                                self.ex_dir, file_format + "_files")
                        G.create_directory(ff_dir)
                        to_file = os.path.join(ff_dir, files)
                        move_file(from_file, to_file)
                        info = "remove: " + files
                        G.logger(info)
                        remove_from_list(file_format, files)
                        try:
                            os.remove(filepath)
                        except FileNotFoundError:
                            pass

        def remove_directories():
            for item in delete:
                search_str = str(item + "_*")
                for root, dirs, files in os.walk(self.target_dir):
                    for directory in dirs:
                        if directory not in systemdirs:
                            if re.search(search_str, directory):
                                dir_path = os.path.join(root, directory)
                                file_names = os.listdir(dir_path)
                                if len(file_names) == 0:
                                    info = "Remove empty directory " + dir_path
                                    G.logger(info)
                                    os.rmdir(dir_path)
                                else:
                                    x = '_'.join(
                                        file_names[0].split("_")[0:-1])
                                    if x == item:
                                        from_dir = os.path.join(
                                            root, directory)
                                        to_dir = os.path.join(
                                            self.ex_dir, directory)
                                        move_dirs(from_dir, to_dir)

        def remove_files():
            for item in delete:
                search_str = str(item + "_*")
                if "GCF" and "v" in item or "GCA" and "v" in item:
                    accession = ".".join(item.split("v"))
                    search_str = str(accession + "_*")
                # genomic_fna
                genome_dir = os.path.join(self.target_dir, "genomic_fna")
                if os.path.isdir(genome_dir):
                    for files in os.listdir(genome_dir):
                        if re.search(search_str, files):
                            G.create_directory(
                                os.path.join(self.ex_dir, "genomic_fna"))
                            from_file = os.path.join(genome_dir, files)
                            to_file = os.path.join(
                                self.ex_dir, "genomic_fna", files)
                            move_file(from_file, to_file)

                    delete_files(item, self.ffn_dir)
                    delete_files(item, self.gff_dir)
                    delete_files(item, self.fna_dir)

                # clean the remaining directories
                for root, dirs, files in os.walk(self.target_dir):
                    for file_name in files:
                        if re.search(search_str, file_name):
                            x = '_'.join(file_name.split("_")[0:-1])
                            if x == item:
                                info = "remove " + str(file_name)
                                G.logger(info)
                                os.remove(os.path.join(root, file_name))
        remove_directories()
        remove_files()

    def check_passed_list(self, passed_list, qc_gene):
        if self.config.ignore_qc:
            info = (
                "--ignore_qc option: Check quality control"
                " files in the Summary directory")
            G.logger("> " + info)
            print("\n" + info)
            return 0

        if len(passed_list) == 0:
            error_msg = "Error: No genomes survived QC"
            print(error_msg)
            G.logger("> " + error_msg)
            errors.append([self.target, error_msg])
            self.remove_qc_failures(qc_gene)
            return 1

        if len(passed_list) < 2:
            error_msg = "Error: Less than two genomes survived QC"
            print(error_msg)
            G.logger("> " + error_msg)
            errors.append([self.target, error_msg])
            self.remove_qc_failures(qc_gene)
            return 1

        self.remove_qc_failures(qc_gene)
        return 0


    def quality_control(self, qc_gene):
        pan = os.path.join(self.pangenome_dir, "gene_presence_absence.csv")
        qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
        if os.path.isfile(pan):
            info = "Found Pangenome directory, skip QC " + qc_gene
            G.logger("> " + info)
            print(info)
            return 0

        print("\nRun: quality_control(" + qc_gene + ")")
        G.logger("Run: quality_control(" + qc_gene + ")")

        if self.get_qc_seqs(qc_gene) == 0:
            qc_seqs = self.choose_sequence(qc_gene)
            if len(qc_seqs) > 0:
                use_cores, inputseqs = BlastPrep(
                        qc_dir, qc_seqs, qc_gene,
                        self.config.blastseqs).run_blastprep()
                Blast(
                    self.config, qc_dir, "quality_control"
                ).run_blast(qc_gene, use_cores)

            passed_list = self.qc_blast_parser(qc_gene)
            exitcode = self.check_passed_list(passed_list, qc_gene)
        else:
            exitcode = 1
        return exitcode


class PangenomeAnalysis:
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")

    def get_numberofgenomes(self):
        inputgenomes = os.listdir(self.gff_dir)
        number = len(inputgenomes)
        stats = PipelineStatsCollector(self.target_dir)
        msg = ["genome assemblies for pan-genome analysis:", str(number)]
        stats.write_stat(" ".join(msg))
        return number

    def run_roary(self):
        num_cpus = str(multiprocessing.cpu_count())
        G.logger("Run: run_roary(" + self.target + ")")
        os.chdir(self.target_dir)
        roary_cmd = [
            "roary", "-f", "./Pangenome", "-s", "-p", num_cpus,
            "-cd", "100", "-e", "-n", "./gff_files/*.gff"]
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
                "fasttree", "-nt", "-gtr", "-nopr", coregenealn, ">", tree]
            try:
                G.run_shell(
                    " ".join(fasttree_cmd), printcmd=True, logcmd=True,
                    log=True)
            except (KeyboardInterrupt, SystemExit):
                G.keyexit_rollback(
                    "fasttree run", dp=self.pangenome_dir, fn=tree)
                raise
        os.chdir(self.target_dir)

    def run_pangenome_analysis(self):
        G.logger("> Start pan-genome analysis")
        print("Start pan-genome analysis")
        G.logger("Run: run_pangenome_analysis(" + self.target + ")")
        exitstat = 0
        if os.path.isdir(self.pangenome_dir):
            filepath = os.path.join(
                self.pangenome_dir, "gene_presence_absence.csv")
            if os.path.isfile(filepath):
                info = (
                    "Pangenome directory already exists\n"
                    "Continue with existing Pangenome data")
                print(info)
                G.logger("> " + info)
                exitstat = 2
                return exitstat

            shutil.rmtree(self.pangenome_dir)

        self.get_numberofgenomes()
        self.run_roary()
        self.run_fasttree()
        return exitstat


class CoreGenes:
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.ffn_dir = os.path.join(self.target_dir, "ffn_files")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.fasta_dir = os.path.join(self.results_dir, "fasta")
        self.all_core_path = os.path.join(self.pangenome_dir, "allcoregenes")
        self.multi_path = os.path.join(self.pangenome_dir, "multiannotated")
        self.singlecopy = os.path.join(
                self.pangenome_dir, "singlecopy_genes.csv")
        self.ffn_seqs = os.path.join(self.pangenome_dir, "ffn_sequences.csv")

    def get_singlecopy_genes(self, mode):
        total_count = []
        single_count = []
        all_core = []
        multi_annotated = []
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
                average_seq_per_isolate = float(row[5])
                loci = row[14:]
                if number_isolates == genomes:
                    total_count.append(gene_name)
                    if number_sequences == genomes:
                        if average_seq_per_isolate == 1:
                            single_count.append(gene_name)
                            data_row.append(gene_name)
                            for locus in loci:
                                data_row.append(locus)
                            if "group" in gene_name:
                                all_core.append(gene_name)
                                newtabledata.append(data_row)
                            if "group" not in gene_name:
                                if len(gene_name.split("_")) > 1:
                                    multi_annotated.append(gene_name)
                                else:
                                    all_core.append(gene_name)
                                    newtabledata.append(data_row)

        if mode == "normal":
            G.csv_writer(self.singlecopy, newtabledata)

        self.print_gene_stats(all_core, total_count)
        coregenesummary = [
            len(all_core), len(total_count),
            len(multi_annotated), len(newtabledata)]
        return coregenesummary

    def print_gene_stats(self, all_core, total_count):
        all_genes = (
            "Continue with " + str(len(all_core))
            + " single copy core genes")
        print("\n# of core genes: " + str(len(total_count)))
        G.logger("> " + all_genes)
        print("\n" + all_genes)
        stats = PipelineStatsCollector(self.target_dir)
        stats.write_stat("core genes: " + str(len(total_count)))
        stats.write_stat("single copy core genes: " + str(len(all_core)))

    def get_sequences_from_ffn(self):
        locustags = {}
        if os.path.isfile(self.ffn_seqs):
            msg = "Found ffn_sequences.csv file, read from file"
            print("\n" + msg + "\n")
            G.logger(msg)
            with open(self.ffn_seqs, "r") as f:
                r = csv.reader(f)
                for row in r:
                    locustags.update({row[1]: {"name": row[0], "seq": row[2]}})
        else:
            ffn_data = []
            for files in os.listdir(self.ffn_dir):
                if files.endswith(".ffn"):
                    filepath = os.path.join(self.ffn_dir, files)
                    with open(filepath) as f:
                        records = SeqIO.parse(f, "fasta")
                        for record in records:
                            name = files.split(".ffn")[0]
                            recid = record.id
                            locus = recid.split(" ")[0]
                            seq = str(record.seq)
                            locustags.update(
                                    {locus: {"name": name, "seq": seq}})
                            ffn_data.append([name, locus, seq])

            G.csv_writer(self.ffn_seqs, ffn_data)

        return locustags

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
        with open(self.singlecopy, "r") as f:
            reader = csv.reader(f)
            for row in reader:
                gene = self.check_genename(row[0])
                outfile = os.path.join(self.fasta_dir, gene + ".fasta")
                with open(outfile, "w") as r:
                    for item in row[1:]:
                        name = locustags[item]["name"]
                        seq = locustags[item]["seq"]
                        record = SeqRecord(
                            Seq(seq),
                            name=item,
                            id='{}|{}|{}'.format(name, item, gene),
                            description="")
                        SeqIO.write(record, r, "fasta")

    def coregene_extract(self, mode="normal"):
        info = "Run: core_gene_extract(" + self.target + ")"
        print(info)
        G.logger(info)
        self.get_singlecopy_genes(mode)
        if mode == "normal":
            locustags = self.get_sequences_from_ffn()
            self.get_fasta(locustags)

    def remove_intermediatefiles(self):
        filelist = [self.singlecopy, self.ffn_seqs]
        for item in filelist:
            filepath = os.path.join(self.pangenome_dir, item)
            os.remove(filepath)

    def run_CoreGenes(self):
        count = 0
        print("\nCollect results of pan-genome analysis\n")
        G.logger("Run: run_CoreGenes(" + self.target + ")")
        G.logger("> Collect results of pan-genome analysis")
        os.chdir(self.pangenome_dir)
        coregenes = os.path.join(self.results_dir, "fasta", "coregenes.txt")
        if not os.path.isdir(os.path.join(self.results_dir, "fasta")):
            G.create_directory(self.results_dir)
            G.create_directory(os.path.join(self.results_dir, "fasta"))
            self.coregene_extract()
            if self.config.intermediate is False:
                self.remove_intermediatefiles()
        elif os.path.isfile(coregenes):
            print("found: " + coregenes)
            with open(coregenes, "r") as f:
                for line in f:
                    count = count + 1
            if count <= 1:
                shutil.rmtree(os.path.join(self.results_dir, "fasta"))
                G.create_directory(self.results_dir)
                G.create_directory(os.path.join(self.results_dir, "fasta"))
                self.coregene_extract()
                if self.config.intermediate is False:
                    self.remove_intermediatefiles()
            else:
                print("Single copy fasta files: " + str(count))
                self.coregene_extract(mode="statistics")
        else:
            if len(os.listdir(os.path.join(self.results_dir, "fasta"))) <= 1:
                shutil.rmtree(os.path.join(self.results_dir, "fasta"))
                G.create_directory(self.results_dir)
                G.create_directory(os.path.join(self.results_dir, "fasta"))
                self.coregene_extract()
                if self.config.intermediate is False:
                    self.remove_intermediatefiles()
            else:
                self.coregene_extract(mode="statistics")


class CoreGeneSequences:
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.fasta_dir = os.path.join(self.results_dir, "fasta")
        self.alignments_dir = os.path.join(self.results_dir, "alignments")
        self.consensus_dir = os.path.join(self.results_dir, "consensus")
        self.blast_dir = os.path.join(self.results_dir, "blast")

    def seq_alignments(self):

        def get_run_commands(run_file):
            cmds = []
            for files in os.listdir(self.fasta_dir):
                if files.endswith(".fasta"):
                    file_name = files.split(".")[0]
                    fasta_files.append(files)
                    result_path = os.path.join(
                        self.alignments_dir, file_name + ".best.fas")
                    if not os.path.isfile(result_path):
                        fasta_path = os.path.join("fasta", files)
                        aligned_path = os.path.join("alignments", file_name)
                        cmds.append(
                            "prank -d=" + fasta_path + " -o=" + aligned_path)
            if len(cmds) > 0:
                with open(run_file, "w") as w:
                    for cmd in cmds:
                        w.write(cmd + "\n")

        fasta_files = []
        print("Start alignment of core gene sequences")
        G.logger("> " + "Start alignment of core gene sequences")
        G.logger("Run: seq_alignments(" + self.target + ")")
        G.create_directory(self.alignments_dir)
        run_file = os.path.join(self.results_dir, "run_prank")
        coregenes = os.path.join(self.results_dir, "fasta", "coregenes.txt")
        if os.path.isfile(run_file):
            os.remove(run_file)
        if os.path.isfile(coregenes):
            info = "Skip " + " ".join(["parallel", "-a", run_file])
            print(info)
            G.logger(info)
            return

        os.chdir(self.results_dir)
        get_run_commands(run_file)
        if os.path.isfile(run_file):
            try:
                G.run_subprocess(
                    ["parallel", "-a", run_file], True, True, False)
            except (KeyboardInterrupt, SystemExit):
                G.keyexit_rollback(
                    "Prank MSA run", dp=self.alignments_dir, fp=run_file)
                raise

        with open(coregenes, "w") as f:
            for fastafile in fasta_files:
                f.write(fastafile + "\n")
                if self.config.intermediate is False:
                    filepath = os.path.join(self.fasta_dir, fastafile)
                    os.remove(filepath)

        if self.config.intermediate is False:
            if os.path.isfile(run_file):
                os.remove(run_file)
        os.chdir(self.target_dir)

    def get_consensus_input(self, cons_summary):
        input_files = []
        output_files = []
        aln_summary = os.path.join(
                                self.alignments_dir, "alignments_summary.txt")
        if not os.path.isfile(aln_summary):
            with open(aln_summary, "w") as f:
                for files in os.listdir(self.alignments_dir):
                    if files.endswith(".best.fas"):
                        filename = files.split(".")[0]
                        input_files.append(filename)
                        f.write(files + "\n")

        if os.path.isfile(cons_summary):
            target = self.target + "_"
            with open(cons_summary, "r") as f:
                for line in f:
                    if ">" in line:
                        li = line.strip()
                        filename = li.split(target)[1].split("_consensus")[0]
                        output_files.append(filename)
        else:
            for files in os.listdir(self.consensus_dir):
                if files.endswith("_consens.fasta"):
                    filename = files.split("_consens.fasta")[0]
                    output_files.append(filename)

        task =  set(input_files) - set(output_files)

        return list(task), input_files

    def write_consensus_data(self, cons_summary):
        recordlist = []
        if os.path.isfile(cons_summary):
            records = list(SeqIO.parse(cons_summary, "fasta"))
            recordlist.append(records)
        for files in os.listdir(self.consensus_dir):
            if files.endswith("_consens.fasta"):
                filepath = os.path.join(self.consensus_dir, files)
                records = list(SeqIO.parse(filepath, "fasta"))
                recordlist.append(records)
                if self.config.intermediate is False:
                    os.remove(filepath)

        with open(cons_summary, "w") as out:
            for records in recordlist:
                SeqIO.write(records, out, "fasta")

    def write_consensus_commands(self, run_file, task):
        with open(run_file, "w") as w:
            for file_name in task:
                seqname = self.target + "_" + file_name + "_consensus"
                result_path = os.path.join(
                    "consensus", file_name + "_consens.fasta")
                aligned_path = os.path.join(
                    "alignments", file_name + ".best.fas")
                w.write(
                    "consambig -sequence " + aligned_path
                    + " -outseq " + result_path
                    + " -name " + seqname + " -auto\n")

    def seq_consensus(self):

        info = "Find consensus sequence for aligned core gene sequences"
        print(info)
        G.logger("> " + info)
        G.logger("Run: seq_consensus(" + self.target + ")")
        cons_summary = os.path.join(
                self.consensus_dir, "consensus_summary.txt")
        G.create_directory(self.consensus_dir)
        run_file = os.path.join(self.results_dir, "run_consensus")
        if os.path.isfile(run_file):
            os.remove(run_file)
        os.chdir(self.results_dir)
        task, input_files = self.get_consensus_input(cons_summary)
        if len(task) == 0:
            info = "Skip " + " ".join(["parallel", "-a", run_file])
            print(info)
            G.logger(info)
            return

        self.write_consensus_commands(run_file, task)

        try:
            G.run_subprocess(
                ["parallel", "-a", run_file], True, True, False)
        except (KeyboardInterrupt, SystemExit):
            G.keyexit_rollback(
                "consensus run", dp=self.consensus_dir, fp=run_file)
            raise

        self.write_consensus_data(cons_summary)

        if self.config.intermediate is False:
            for filename in input_files:
                filepath = os.path.join(
                    self.alignments_dir, filename + ".best.fas")
                if os.path.isfile(filepath):
                    os.remove(filepath)

        os.chdir(self.target_dir)

    def conserved_seqs(self):
        info = "Search conserved regions in consensus sequences"
        print(info)
        filepaths = []
        G.logger("> " + info)
        G.logger("Run: conserved(" + self.target + ")")
        G.create_directory(self.blast_dir)
        conserv_seqs = []
        cons_summary = os.path.join(
            self.consensus_dir, "consensus_summary.txt")
        result_path = os.path.join(
            self.blast_dir,
            H.abbrev(self.target) + "_conserved")

        if os.path.isfile(cons_summary):
            filepaths.append(cons_summary)
        else:
            for files in os.listdir(self.consensus_dir):
                if files.endswith("_consens.fasta"):
                    file_path = os.path.join(self.consensus_dir, files)
                    filepaths.append(file_path)

        for filepath in filepaths:
            for record in SeqIO.parse(filepath, "fasta"):
                # replace any non-GATC character with N
                sequence = re.sub(
                    "[^GATC]", "N", re.sub("[a-z]", "N", str(record.seq)))
                split_seq = re.sub(
                    "N{5,50}", "*",
                    re.sub("N.{1,20}N", "*", str(sequence)))

                count = 1
                split_list = split_seq.split("*")
                for item in split_list:
                    if len(item) >= self.config.minsize:
                        desc = record.id
                        if "group_" in desc:
                            seq_name = (
                                "group_" + desc.split("_")[-2:-1][0]
                                + "_" + str(count))
                        else:
                            seq_name = (
                                desc.split("_")[-2:-1][0] + "_" + str(count))
                        seq = item
                        count += 1
                        conserv_seqs.append([seq_name, seq])

        if len(conserv_seqs) == 0:
            error_msg = "Error: no conserved target sequences found"
            info = "Number of conserved sequences: 0"
            print(error_msg)
            errors.append([self.target, error_msg])
            G.logger("> " + error_msg)
            G.logger("> " + info)
            PipelineStatsCollector(self.target_dir).write_stat(info)
            return 1

        with open(result_path, "w") as f:
            for seq in conserv_seqs:
                f.write(">" + seq[0] + "\n")
                f.write(seq[1] + "\n")
        info = "Number of conserved sequences: " + str(len(conserv_seqs))
        PipelineStatsCollector(self.target_dir).write_stat(info)
        print(info)
        G.logger("> " + info)
        return conserv_seqs


    def run_coregeneanalysis(self):
        G.logger("Run: run_coregeneanalysis(" + self.target + ")")
        self.seq_alignments()
        self.seq_consensus()
        conserved_seqs = self.conserved_seqs()
        if conserved_seqs == 1:
            return 1
        name = "conserved"
        blastsum = os.path.join(self.blast_dir, "BLAST_results_summary.csv")
        if not os.path.isfile(blastsum):
            use_cores, inputseqs = BlastPrep(
                self.blast_dir, conserved_seqs, "conserved",
                self.config.blastseqs).run_blastprep()
            Blast(
                self.config, self.blast_dir, "conserved"
            ).run_blast(name, use_cores)


class BlastPrep():
    def __init__(self, directory, input_list, name, maxpartsize):
        self.list_dict = {}
        self.input_list = input_list
        self.maxpartsize = maxpartsize
        self.filename = name
        self.directory = directory

    def create_listdict(self):
        groups = len(self.input_list) // self.maxpartsize
        if len(self.input_list) % self.maxpartsize > 0:
            groups = groups + 1
        for i in range(0, groups):
            if i not in self.list_dict.keys():
                self.list_dict.update({i: []})

    def get_equalgroups(self):
        self.input_list.sort(key=lambda x: int(len(x[1])), reverse=True)
        list_start = 0
        list_end = len(self.input_list) - 1
        removed_key = []
        key = 0
        i = list_start
        while key in self.list_dict.keys():
            if key in removed_key:
                key = key + 1
            else:
                item = self.input_list[i]
                if len(self.list_dict[key]) < self.maxpartsize:
                    self.list_dict[key].append(item)
                    key = key + 1
                    if key not in self.list_dict.keys():
                        key = 0
                    if i == list_end:
                        break
                    else:
                        i = i + 1
                else:
                    removed_key.append(key)

    def write_blastinput(self):
        inputsequences = []
        for key in self.list_dict.keys():
            if len(self.list_dict[key]) > 0:
                file_name = os.path.join(
                    self.directory, self.filename + ".part-"+str(key))
                with open(file_name, "w") as f:
                    for item in self.list_dict[key]:
                        f.write(">" + item[0] + "\n")
                        inputsequences.append(item[0])
                        f.write(item[1] + "\n")
        return inputsequences

    def run_blastprep(self):
        G.logger("Run: run_blastprep - Preparing files for BLAST")
        print("\nPreparing files for BLAST")
        self.create_listdict()
        self.get_equalgroups()
        cores = multiprocessing.cpu_count()
        inputseqs = self.write_blastinput()
        return cores, inputseqs


class Blast:
    def __init__(self, configuration, directory, mode):
        print("Start BLAST")
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.config.path, self.target, "config")
        self.directory = directory
        self.mode = mode

    def get_blast_cmd(self, blastfile, filename, cores):
        fmt_file = os.path.join(dict_path, "blastfmt6.csv")
        fmts = pd.read_csv(fmt_file, header=None).dropna()
        blast_fmt = " ".join(["6"] + list(fmts[0]))

        if self.mode == "quality_control":
            blast_cmd = [
                "blastn", "-task", "megablast", "-num_threads",
                str(cores), "-query", blastfile, "-max_target_seqs", "5",
                "-max_hsps", "1", "-out", filename, "-outfmt", blast_fmt]

        if self.mode == "conserved":
            blast_cmd = [
                "blastn", "-task", "dc-megablast", "-num_threads",
                str(cores), "-query", blastfile, "-max_target_seqs",
                "2000", "-evalue", "500", "-out", filename, "-outfmt", blast_fmt]

        if self.mode == "primer":
            blast_cmd = [
                "blastn", "-task", "blastn-short", "-num_threads",
                str(cores), "-query", blastfile,
                "-evalue", "500", "-out", filename, "-outfmt", blast_fmt]

        blast_cmd.append("-db")
        if self.config.customdb:
            blast_cmd.append(self.config.customdb)
        else:
            blast_cmd.append("nt")

        return blast_cmd

    def search_blastfiles(self, directory):
        blast_files = []
        for files in os.listdir(directory):
            if fnmatch.fnmatch(files, "*.part*"):
                blast_files.append(files)
        return blast_files

    def run_blast(self, name, use_cores):
        G.logger("Run: run_blast - Start BLAST")
        blastfiles = self.search_blastfiles(self.directory)
        if len(blastfiles) > 0:
            blastfiles.sort(key=lambda x: int(x.split("part-")[1]))
            start = time.time()
            os.chdir(self.directory)
            for blastfile in blastfiles:
                blast_cmd = False
                part = str(blastfile).split("-")[1]
                filename = name + "_" + part + "_results.csv"
                results_path = os.path.join(self.directory, filename)
                if self.mode == "quality_control":
                    blast_cmd = self.get_blast_cmd(
                        blastfile, filename, use_cores)
                elif not os.path.isfile(results_path):
                    blast_cmd = self.get_blast_cmd(
                        blastfile, filename, use_cores)
                else:
                    if os.stat(results_path).st_size == 0:
                        blast_cmd = self.get_blast_cmd(
                            blastfile, filename, use_cores)
                    else:
                        info = "Skip Blast step for " + blastfile
                        print("\n" + info)
                        G.logger("> " + info)

                if blast_cmd:
                    try:
                        G.run_subprocess(blast_cmd)
                    except (KeyboardInterrupt, SystemExit):
                        G.keyexit_rollback(
                            "BLAST search", dp=self.directory, fn=filename)
                        raise

            duration = time.time() - start
            G.logger(
                "> Blast duration: "
                + str(timedelta(seconds=duration)).split(".")[0])
            os.chdir(self.target_dir)

class BlastParser:
    def __init__(self, configuration, results="seqs"):
        self.exception = configuration.exception
        self.config = configuration
        self.target = configuration.target
        self.evalue = self.config.evalue
        self.nuc_identity = self.config.nuc_identity
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.nontargetlist = configuration.nontargetlist
        self.selected = []
        self.mode = "normal"
        self.start = time.time()
        if results == "primer":
            self.mode = "primer"
            self.primer_dir = os.path.join(self.results_dir, "primer")
            self.blast_dir = os.path.join(self.primer_dir, "primerblast")
            #self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
            self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
            self.maxgroupsize = 25000

    def blastresult_files(self):
        blastresults = []
        print("Run: blastresults_files(" + self.target + ")")
        for files in [f for f in os.listdir(self.blast_dir) if f.endswith("results.csv")]:
            file_path = os.path.join(self.blast_dir, files)
            if file_path not in blastresults:
                blastresults.append(file_path)
        blastresults.sort()
        return blastresults

    def check_blastdb_errors(self, blastdf, filename):    
        if len(blastdf.index) == 0:
            error_msg = " ".join([
                "A problem with the BLAST results file",
                filename, "was detected.",
                "Please check if the file was removed and start the run again"])
            
        elif len(blastdf[blastdf["Subject Seq-id"].str.contains("gnl|BL_ORD_ID|", regex=False)]) > 0:
            error_msg = (
                "Problem with custom DB, Please use the '-parse_seqids'"
                " option for the makeblastdb command")

        elif len(blastdf[blastdf["Subject Title"].str.contains("No definition line", regex=False)]) > 0:
            error_msg = (
                "Error: No definition line in Subject Title"
                "\nData is missing in the custom BLAST DB. At least "
                "a unique sequence identifier and the species name "
                "is required for each entry.\nExpected format: "
                ">seqid species name optional description")           
        else:
            return

        print("\n" + error_msg + "\n")
        logging.error("> " + error_msg, exc_info=True)
        errors.append([self.target, error_msg])
        os.remove(filename)
        print("removed " + filename)
        raise BlastDBError(error_msg)

    def check_seq_ends(self, rejected):
        #w_mode = 'a'
        #if os.path.isfile("BLAST_results_summary.csv") is True and filenum == 0:
            #w_mode = 'w'    
        # Filter non-aligned endings
        rejected.loc[:, 'overhang'] = (
                            rejected.loc[:,'Query sequence length'] - rejected.loc[:, 'End of alignment in query'])
        partials_max = rejected.sort_values('overhang', ascending=True).drop_duplicates(['Query Seq-id'])
        keep_max = partials_max[partials_max['overhang'] >= self.config.minsize]
        partials_min = rejected.sort_values('Start of alignment in query', ascending=True).drop_duplicates(['Query Seq-id'])
        keep_min = partials_min[partials_min["Start of alignment in query"] >= self.config.minsize]
        keep_min = keep_min.assign(Start=1)

        # Write sequence range data
        mindata = keep_min[['Query Seq-id', "Start", "Start of alignment in query"]]
        maxdata = keep_max[['Query Seq-id', 'End of alignment in query', 'Query sequence length']]
        mindata.columns = ["ID", "Start", "Stop"]
        maxdata.columns = ["ID", "Start", "Stop"]
            
        return pd.concat([mindata, maxdata], sort=False)
   
    def get_excluded_gis(self):
        excluded_gis = []
        gi_file = os.path.join(self.config_dir, "no_blast.gi")
        if os.path.isfile(gi_file):
            if os.stat(gi_file).st_size > 0:
                with open(gi_file, "r") as f:
                    for line in f:
                        gi = line.strip()
                        excluded_gis.append(str(gi))
        return excluded_gis
        
    def parse_blastrecords(self, blastdf, excluded_gis):
        exceptions = [H.subspecies_handler(self.target, mode="space")]
        if self.exception != []:
            for item in self.exception:
                exception = ' '.join(item.split("_"))
                if exception not in exceptions:
                    exceptions.append(exception)
                    
        target_filter = blastdf["Subject Title"].str.contains("|".join(exceptions))
        target_hits = blastdf[target_filter]
        offtarget_hits = blastdf[~target_filter]

        if self.config.nolist is True:
            # if no target list is defined this is the result
            accept = target_hits.copy()
            reject = offtarget_hits.copy()                
        else:
            offtarget_filter = offtarget_hits["Subject Title"].str.contains("|".join(self.nontargetlist))
            listed = offtarget_hits[offtarget_filter]                                        
            unlisted = offtarget_hits[~offtarget_filter]
            accept = pd.concat([target_hits, unlisted], sort=False)
            reject = listed.copy()

        # remove excluded sequences from the results
        reject = reject[~reject["Subject GI"].isin(excluded_gis)]
        reject = reject[~reject["Subject accession"].isin(excluded_gis)]
        
        if self.mode == "normal":
            # Filter according to configuration
            reject = reject[reject['Percentage of identical matches'] >= self.config.nuc_identity]
            reject = reject[reject['Expect value'] <= self.config.evalue]        
            # check if sequence ends without off-target matches exist
            partial = self.check_seq_ends(reject)
            speciesdata = blastdf[['Query Seq-id', 'Subject accession', "Subject Title"]]            
            reject = reject[['Query Seq-id', "Subject accession"]]
        else:
            partial = pd.DataFrame()
            speciesdata = pd.DataFrame()
            reject = reject[[
                            "Query Seq-id", "Subject accession", "Start of alignment in subject", 
                            "End of alignment in subject", "Subject sequence length", "Subject Title"]]

        accept = accept[['Query Seq-id', "Subject accession"]]
        return accept, reject, partial, speciesdata

    def blastresults_summary(self, accept, reject, partial):       
        # write summary
        ac_sum = accept.groupby(['Query Seq-id'])["Subject accession"].count()
        re_sum = reject.groupby(['Query Seq-id'])["Subject accession"].count()
        summary = pd.concat([ac_sum, re_sum], axis=1, sort=False).replace(np.nan, 0)
        summary.columns = ["Target hits", "Off-target hits"]
        summary.index.name = "Query"
        sum_file = os.path.join(self.blast_dir, "BLAST_results_summary.csv")
        summary.to_csv(sum_file)
        specific_seqs = summary[summary["Off-target hits"] == 0].index.to_list()
        #unspecific_seqs = summary[summary["Off-target hits"] != 0].index.to_list()
        if len(partial.index) != 0:
            part_file = os.path.join(self.blast_dir, "partialseqs.csv")
            partial.to_csv(part_file, index=False, header=False)
            
        return specific_seqs

    def write_mostcommonhits(self, df):
        to_file = os.path.join(self.blast_dir, "mostcommonhits.csv")
        total = len(df.index)
        queries = len(set(df["Query Seq-id"]))
        df = self.get_species_names_from_title(df)
        mostcommon = pd.DataFrame(df.drop_duplicates(["Query Seq-id", "Species"])["Species"].value_counts())
        mostcommon.index.name ="Species"
        mostcommon.columns = ["BLAST hits [count]"]
        mostcommon["BLAST hits [% of queries]"] = mostcommon["BLAST hits [count]"].apply(                                                                    lambda x: round(100/queries*x, 1))
        mostcommon.sort_values("BLAST hits [% of queries]", ascending=False, inplace=True)
        f_head = str("Total BLAST hits,Number of queries\n" + str(total) + "," + str(queries) + "\n")
        with open(to_file, "w") as f:
            f.write(f_head)
        mostcommon.to_csv(to_file, mode='a')

        
    def get_species_names_from_title(self, df):
        if self.config.virus is True:
            df.loc[:, "Species"] = df.loc[:, "Subject Title"].str.split(",").str[0]
        else:    
            subsp_filter = df["Subject Title"].str.contains("|".join(["subsp.", "pv."]))
            df.loc[subsp_filter, "Species"] = df.loc[
                                                    subsp_filter, "Subject Title"
                                                        ].str.split(" ").str[0:4].apply(
                                                                    lambda x: ' '.join(x))
            df.loc[~subsp_filter, "Species"] = df.loc[
                                                    ~subsp_filter, "Subject Title"
                                                        ].str.split(" ").str[0:2].apply(
                                                                    lambda x: ' '.join(x))
        return df

    def bp_parse_results(self):
        target_dfs = []
        offtarget_dfs = []
        partial_dfs = []
        speciesdata_dfs = []
        blastresults = self.blastresult_files()
        excluded_gis = self.get_excluded_gis()
        print("Excluded GI(s):", excluded_gis)
        fmt_file = os.path.join(dict_path, "blastfmt6.csv")
        header = list(pd.read_csv(fmt_file, header=None)[1].dropna())
        tmp_file = os.path.join(self.blast_dir, "tmp_mostcommon.csv")
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)
        for i, filename in enumerate(blastresults):
            print(
                "\nopen BLAST result file " + str(i+1)
                + "/" + str(len(blastresults)))
            try:
                blastdf = pd.read_csv(filename, sep="\t", header=None)
                blastdf.columns = header
            except pd.errors.EmptyDataError:
                blastdf = pd.DataFrame()
            
            self.check_blastdb_errors(blastdf, filename)

            accept, reject, partial, speciesdata = self.parse_blastrecords(blastdf, excluded_gis)
            target_dfs.append(accept)
            offtarget_dfs.append(reject)
            partial_dfs.append(partial)
            speciesdata_dfs.append(speciesdata)
        
        target = pd.concat(target_dfs)
        offtarget = pd.concat(offtarget_dfs)
        partial = pd.concat(partial_dfs)
        speciesdata = pd.concat(speciesdata_dfs)
        
        specific_seqs = self.blastresults_summary(target, offtarget, partial)
       
        if self.mode == "primer":
            return offtarget
        
        self.write_mostcommonhits(speciesdata)
        return specific_seqs
        

    def changed_primer3_input(self, file_path, controlfile_path):

        def find_difference():
            new = []
            old = []
            with open(file_path) as n:
                for line in n:
                    if "SEQUENCE_ID=" in line:
                        if line.strip() not in new:
                            new.append(line.strip())
                    if "PRIMER_PICK_INTERNAL_OLIGO=" in line:
                        if line.strip() not in new:
                            new.append(line.strip())

            with open(controlfile_path) as o:
                for line in o:
                    if "SEQUENCE_ID=" in line:
                        if line.strip() not in old:
                            old.append(line.strip())
                    if "PRIMER_PICK_INTERNAL_OLIGO=" in line:
                        if line.strip() not in old:
                            old.append(line.strip())

            diff = list(set(new) ^ set(old))
            return diff

        if os.path.isfile(controlfile_path):
            diff = find_difference()
            if len(diff) > 0:
                info1 = (
                    "Due to changed settings primer design "
                    "and quality control will start from scratch")
                info2 = "Differences in primer3 input:"
                for info in [info1, info2, diff]:
                    G.logger(info)
                    print(info)
                primer_dir = os.path.join(self.results_dir, "primer")
                if os.path.isdir(primer_dir):
                    G.logger("Delete primer directory")
                    print("Delete primer directory")
                    shutil.rmtree(primer_dir)

                shutil.copy(file_path, controlfile_path)
        else:
            shutil.copy(file_path, controlfile_path)
    
    def write_primer3_input(self, selected_seqs):        
        #G.logger("Run: write_primer3_input(" + self.target + ")")
        conserved_seqs = os.path.join(self.blast_dir, H.abbrev(self.target) + "_conserved")
        conserved_seq_dict = SeqIO.to_dict(SeqIO.parse(conserved_seqs, "fasta"))
        file_path = os.path.join(self.results_dir, "primer3_input")
        controlfile_path = os.path.join(self.results_dir, ".primer3_input")
        if self.config.probe is True:
            probe = "\nPRIMER_PICK_INTERNAL_OLIGO=1"
        else:
            probe = ""

        with open(file_path, "w") as f:
            for item in selected_seqs:
                f.write(
                        "SEQUENCE_ID=" + item + "\nSEQUENCE_TEMPLATE="
                        + str(conserved_seq_dict[item].seq)
                        + "\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize) + "-"
                        + str(self.config.maxsize) + probe + "\n=\n")
        
            partial_file = os.path.join(os.getcwd(), "partialseqs.csv")
            if os.path.isfile(partial_file):
                parts = pd.read_csv(partial_file, header=None)
                seq_id = parts[0].to_list()
                start = parts[1].to_list()
                end = parts[2].to_list()
                for i, idx in enumerate(seq_id):
                    f.write(
                        "SEQUENCE_ID=" + idx + "\nSEQUENCE_TEMPLATE="
                        + str(conserved_seq_dict[idx].seq)[start[i]:end[i]]
                        + "\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize) + "-"
                        + str(self.config.maxsize) + probe + "\n=\n")
        self.changed_primer3_input(file_path, controlfile_path)

    def find_primerbinding_offtarget_seqs(self, df):
        df.loc[:, "Primer pair"] = df.loc[:, "Query Seq-id"].str.split("_").str[0:-1].apply(
                                                                        lambda x: '_'.join(x))
        df.sort_values(['Start of alignment in subject'], inplace=True)

        fwd_df = df[df["Query Seq-id"].str.endswith("_F")]
        rev_df =  df[df["Query Seq-id"].str.endswith("_R")]
        int_df = pd.merge(
                    fwd_df, rev_df, how ='inner', 
                    on =['Subject accession', 'Primer pair'], suffixes=("_F", "_R"))

        f = int_df[[
            'Subject accession', 'Start of alignment in subject_F', 
            'End of alignment in subject_F', 'Subject sequence length_F']] 
        r = int_df[[
            'Subject accession', 'Start of alignment in subject_R', 
            'End of alignment in subject_R', 'Subject sequence length_R']]
        std_cols = [
            'Subject accession', 'Start of alignment in subject', 
            'End of alignment in subject', 'Subject sequence length']

        f.columns, r.columns = std_cols, std_cols
        common = pd.concat([f, r], sort=False)
        common.reset_index(drop=True, inplace=True)
        return common
        
    def get_primerBLAST_DBIDS(self, offtarget):
        print("\nGet sequence accessions of BLAST hits\n")
        G.logger("> Get sequence accessions of BLAST hits")
        G.create_directory(self.primer_qc_dir)
        overhang=2000
        output_path = os.path.join(self.primer_qc_dir, "primerBLAST_DBIDS.csv")
        if os.path.isfile(output_path):
            return 0
        # data manipulation
        strandfilter = offtarget['Start of alignment in subject'] > offtarget['End of alignment in subject']
        offtarget.loc[strandfilter, "Start overhang"] = offtarget.loc[strandfilter, 'End of alignment in subject'] - overhang
        offtarget.loc[~strandfilter, "Start overhang"] = offtarget.loc[~strandfilter,'Start of alignment in subject'] - overhang
        offtarget.loc[strandfilter, "End overhang"] = offtarget.loc[strandfilter, 'End of alignment in subject'] + overhang
        offtarget.loc[~strandfilter, "End overhang"] = offtarget.loc[~strandfilter,'Start of alignment in subject'] + overhang

        overfilter = offtarget["End overhang"] > offtarget['Subject sequence length']
        offtarget.loc[overfilter, "End overhang"] = offtarget.loc[overfilter, 'Subject sequence length']
        lowfilter = offtarget["Start overhang"] < 1
        offtarget.loc[lowfilter, "Start overhang"] = 1
        # datatype to int
        offtarget["Start overhang"] = offtarget["Start overhang"].astype('Int64')
        offtarget["End overhang"] = offtarget["End overhang"].astype('Int64')

        max_range = offtarget["End overhang"].max()
        stepsize = self.config.maxsize + overhang*2 + 1

        collection = []
        for i in range(1, max_range,  stepsize):            
            j = i + overhang*2 + self.config.maxsize
            sub = offtarget[offtarget["Start overhang"].between(i, j, inclusive=True)]
            mini = sub.groupby(["Subject accession"])["Start overhang"].min()
            maxi = sub.groupby(["Subject accession"])["End overhang"].max()
            submax = pd.concat([mini, maxi], axis=1)
            submax.columns = ["Start", "Stop"]
            
            submax.sort_values(["Start", "Stop"], inplace=True, ascending=False)
            submax.drop_duplicates(inplace=True)             
            
            collection.append(submax)

        results = pd.concat(collection)
        
        if len(results.index) == 0:
            msg = (
                "Error did not find any sequences for non-target DB. "
                + "Please check the species list and/or BLAST database")
            print(msg)
            G.logger("> " + msg)
            errors.append([self.target, msg])
            return 1

        results.index.name = "accession"
        results.to_csv(output_path, header=None)
        return 0

    def write_nontarget_sequences(self, nonreddata):
        maxsize = 25000
        if self.config.customdb is not None:
            db = self.config.customdb
        else:
            db = "nt"
        seqcount = len(nonreddata)
        info = "Found " + str(seqcount) + " sequences for the non-target DB"
        G.logger(info)
        print("\n" + info + "\n")
        parts = seqcount//maxsize + 1
        
        for part in range(0, parts):
            print("Working on part " + str(part + 1) + "/" + str(parts))
            end = (part+1)*maxsize
            if end > seqcount:
                end = seqcount
            filename = "BLASTnontarget" + str(part) + ".sequences"
            filepath = os.path.join(self.primer_qc_dir, filename)
            if not os.path.isfile(filepath):
                info = "Start writing " + filename
                print(info)
                G.logger(info)
                print("Start DB extraction")
                G.logger("Start DB extraction")
                data = nonreddata[part*maxsize:end]
                fasta_seqs = G.run_parallel(
                        P.get_seq_fromDB, data, db)
                try:
                    with open(filepath, "w") as r:
                        for fastainfo in fasta_seqs:
                            name = fastainfo[0].split(":")
                            fastainfo[0] = (
                                name[0] + "_" + "_".join(name[1].split("-")))
                            fastadata = "\n".join(fastainfo) + "\n"
                            r.write(fastadata)
                except (KeyboardInterrupt, SystemExit):
                    G.keyexit_rollback("DB extraction", fp=filepath)
                    raise

            else:
                info2 = "Skip writing " + filename
                print(info2)
                G.logger(info2)
                return
            info3 = "Finished writing " + filename
            print(info3)
            G.logger(info3)
        
    def run_blastparser(self):
        if self.mode == "primer":
            print("Start primer blast parser")
            offtarget = self.bp_parse_results()
            db_seqs = self.find_primerbinding_offtarget_seqs(offtarget)
            self.get_primerBLAST_DBIDS(db_seqs)
            
            dbids = os.path.join(self.primer_qc_dir, "primerBLAST_DBIDS.csv")
            nonreddata = pd.read_csv(dbids, header=None).values.tolist()
            self.write_nontarget_sequences(nonreddata)

            duration = time.time() - self.start
            G.logger(
                "> Primer blast parser time: "
                + str(timedelta(seconds=duration)).split(".")[0])
            print(timedelta(seconds=duration))
        else:
            specific_ids = self.bp_parse_results()
            self.write_primer3_input(specific_ids)
            duration = time.time() - self.start
            info = ("species specific conserved sequences: "
                    + str(len(specific_ids)))
            G.logger(
                "> Blast parser time: "
                + str(timedelta(seconds=duration)).split(".")[0])
            print(timedelta(seconds=duration))
            G.logger("> " + info)
            PipelineStatsCollector(self.target_dir).write_stat(info)

            if len(specific_ids) == 0:
                msg = "No conserved sequences without non-target match found"
                print(msg)
                G.logger("> " + msg)
                errors.append([self.target, msg])
                return 1

            return 0


class BlastParser_dev:
    def __init__(self, configuration, results="seqs"):
        self.exception = configuration.exception
        self.config = configuration
        self.target = configuration.target
        self.evalue = self.config.evalue
        self.nuc_identity = self.config.nuc_identity
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.nontargetlist = configuration.nontargetlist
        self.selected = []
        self.mode = "normal"
        self.start = time.time()
        if results == "primer":
            self.mode = "primer"
            self.primer_dir = os.path.join(self.results_dir, "primer")
            self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
            self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
            self.maxgroupsize = 25000

    def blastresult_files(self, blast_dir):
        blastresults = []
        for files in [f for f in os.listdir(blast_dir) if f.endswith("results.csv")]:
            file_path = os.path.join(blast_dir, files)
            if file_path not in blastresults:
                blastresults.append(file_path)
        return blastresults

    def check_blastdb_errors(self, blastdf, filename):
        if len(blastdf.index) == 0:
            error_msg = " ".join([
                "A problem with the BLAST results file",
                filename, "was detected.",
                "Please check if the file was removed and start the run again"])
            os.remove(filename)
            print("removed " + filename)

        elif len(blastdf[blastdf["Subject Seq-id"].str.contains("gnl|BL_ORD_ID|", regex=False)]) > 0:
            error_msg = (
                "Problem with custom DB, Please use the '-parse_seqids'"
                " option for the makeblastdb command")

        elif len(blastdf[blastdf["Subject Title"].str.contains("No definition line", regex=False)]) > 0:
            error_msg = (
                "Error: No definition line in Subject Title"
                "\nData is missing in the custom BLAST DB. At least "
                "a unique sequence identifier and the species name "
                "is required for each entry.\nExpected format: "
                ">seqid species name optional description")

        else:
            return

        print("\n" + error_msg + "\n")
        logging.error("> " + error_msg, exc_info=True)
        errors.append([self.target, error_msg])
        raise BlastDBError(error_msg)

    def check_seq_ends(self, rejected, filenum):
        w_mode = 'a'
        if os.path.isfile("partialseqs.csv") is True and filenum == 0:
            w_mode = 'w'

        # Filter non-aligned endings
        rejected.loc[:, 'overhang'] = (
                            rejected.loc[:,'Query sequence length'] - rejected.loc[:, 'End of alignment in query'])
        partials_max = rejected.sort_values('overhang', ascending=True).drop_duplicates(['Query Seq-id'])
        keep_max = partials_max[partials_max['overhang'] >= self.config.minsize]
        partials_min = rejected.sort_values('Start of alignment in query', ascending=True).drop_duplicates(['Query Seq-id'])
        keep_min = partials_min[partials_min["Start of alignment in query"] >= self.config.minsize]
        keep_min = keep_min.assign(Start=1)

        # Write sequence range data
        mindata = keep_min[['Query Seq-id', "Start", "Start of alignment in query"]]
        maxdata = keep_max[['Query Seq-id', 'End of alignment in query', 'Query sequence length']]
        mindata.columns = ["ID", "Start", "Stop"]
        maxdata.columns = ["ID", "Start", "Stop"]
        write_partial = pd.concat([mindata, maxdata])
        if len(write_partial.index) != 0:
            write_partial.to_csv("partialseqs.csv", mode=w_mode, index=False, header=False)

    def get_excluded_gis(self):
        excluded_gis = []
        gi_file = os.path.join(self.config_dir, "no_blast.gi")
        if os.path.isfile(gi_file):
            if os.stat(gi_file).st_size > 0:
                with open(gi_file, "r") as f:
                    for line in f:
                        gi = line.strip()
                        excluded_gis.append(str(gi))
        return excluded_gis

    def parse_blastrecords(self, blastdf, filenum, excluded_gis):
        exceptions = [H.subspecies_handler(self.target, mode="space")]
        if self.exception != []:
            for item in self.exception:
                exception = ' '.join(item.split("_"))
                if exception not in exceptions:
                    exceptions.append(exception)

        target_filter = blastdf["Subject Title"].str.contains("|".join(exceptions))
        target_hits = blastdf[target_filter]
        offtarget_hits = blastdf[~target_filter]

        if self.config.nolist is True:
            # if no target list is defined this is the result
            accept = target_hits.copy()
            reject = offtarget_hits.copy()
        else:
            offtarget_filter = offtarget_hits["Subject Title"].str.contains("|".join(self.nontargetlist))
            listed = offtarget_hits[offtarget_filter]
            unlisted = offtarget_hits[~offtarget_filter]
            accept = pd.concat([target_hits, unlisted])
            reject = listed.copy()

        # remove excluded sequences from the results

        reject = reject[~reject["Subject GI"].isin(excluded_gis)]

        # write summary
        ac_sum = accept.groupby(['Query Seq-id'])["Subject Seq-id"].count()
        re_sum = reject.groupby(['Query Seq-id'])["Subject Seq-id"].count()
        summary = pd.concat([ac_sum, re_sum], axis=1).replace(np.nan, 0)
        summary.columns = ["Target hits", "Off-target hits"]
        summary.to_csv("BLAST_results_summary_" + str(filenum) + ".csv")
        specific_seqs = summary[summary["Off-target hits"] == 0].index.to_list()
        unspecific_seqs = summary[summary["Off-target hits"] != 0].index.to_list()

        if self.mode == "normal":
            # check if sequence ends without off-target matches exist
            self.check_seq_ends(reject, filenum)

        return specific_seqs, unspecific_seqs


    def bp_parse_results(self, blast_dir):
        specific_ids = []
        unspecific_ids = []
        blastresults = self.blastresult_files(blast_dir)
        excluded_gis = self.get_excluded_gis()
        print("Excluded GI(s):", excluded_gis)
        fmt_file = os.path.join(dict_path, "blastfmt6.csv")
        header = list(pd.read_csv(fmt_file, header=None)[1].dropna())
        for i, filename in enumerate(blastresults):
            print(
                "\nopen BLAST result file " + str(i+1)
                + "/" + str(len(blastresults)))

            blastdf = pd.read_csv(filename, sep="\t", header=None)
            blastdf.columns = header
            self.check_blastdb_errors(blastdf, filename)
            specific_seqs, unspecific_seqs = self.parse_blastrecords(blastdf, i, excluded_gis)
            specific_ids.extend(specific_seqs)
            unspecific_ids.extend(unspecific_seqs)
        return specific_ids, unspecific_ids

    def changed_primer3_input(self, file_path, controlfile_path):

        def find_difference():
            new = []
            old = []
            with open(file_path) as n:
                for line in n:
                    if "SEQUENCE_ID=" in line:
                        if line.strip() not in new:
                            new.append(line.strip())
                    if "PRIMER_PICK_INTERNAL_OLIGO=" in line:
                        if line.strip() not in new:
                            new.append(line.strip())

            with open(controlfile_path) as o:
                for line in o:
                    if "SEQUENCE_ID=" in line:
                        if line.strip() not in old:
                            old.append(line.strip())
                    if "PRIMER_PICK_INTERNAL_OLIGO=" in line:
                        if line.strip() not in old:
                            old.append(line.strip())

            diff = list(set(new) ^ set(old))
            return diff

        if os.path.isfile(controlfile_path):
            diff = find_difference()
            if len(diff) > 0:
                info1 = (
                    "Due to changed settings primer design "
                    "and quality control will start from scratch")
                info2 = "Differences in primer3 input:"
                for info in [info1, info2, diff]:
                    G.logger(info)
                    print(info)
                primer_dir = os.path.join(self.results_dir, "primer")
                if os.path.isdir(primer_dir):
                    G.logger("Delete primer directory")
                    print("Delete primer directory")
                    shutil.rmtree(primer_dir)

                shutil.copy(file_path, controlfile_path)
        else:
            shutil.copy(file_path, controlfile_path)

    def write_primer3_input(self, selected_seqs):
        G.logger("Run: write_primer3_input(" + self.target + ")")
        conserved_seqs = os.path.join(self.blast_dir, H.abbrev(self.target) + "_conserved")
        conserved_seq_dict = SeqIO.to_dict(SeqIO.parse(conserved_seqs, "fasta"))
        file_path = os.path.join(self.results_dir, "primer3_input")
        controlfile_path = os.path.join(self.results_dir, ".primer3_input")
        if self.config.probe is True:
            probe = "\nPRIMER_PICK_INTERNAL_OLIGO=1"
        else:
            probe = ""

        with open(file_path, "w") as f:
            for item in selected_seqs:
                f.write(
                        "SEQUENCE_ID=" + item + "\nSEQUENCE_TEMPLATE="
                        + str(conserved_seq_dict[item].seq)
                        + "\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize) + "-"
                        + str(self.config.maxsize) + probe + "\n=\n")

            partial_file = os.path.join(os.getcwd(), "partialseqs.csv")
            if os.path.isfile(partial_file):
                parts = pd.read_csv(partial_file, header=None)
                seq_id = parts[0].to_list()
                start = parts[1].to_list()
                end = parts[2].to_list()
                for i, idx in enumerate(seq_id):
                    f.write(
                        "SEQUENCE_ID=" + idx + "\nSEQUENCE_TEMPLATE="
                        + str(conserved_seq_dict[idx].seq)[start[i]:end[i]]
                        + "\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize) + "-"
                        + str(self.config.maxsize) + probe + "\n=\n")

        self.changed_primer3_input(file_path, controlfile_path)

class BlastParser_old:
    def __init__(self, configuration, results="seqs"):
        self.exception = configuration.exception
        self.config = configuration
        self.target = configuration.target
        self.evalue = self.config.evalue
        self.nuc_identity = self.config.nuc_identity
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.nontargetlist = configuration.nontargetlist
        self.selected = []
        self.mode = "normal"
        self.start = time.time()
        #self.taxid = self.config.taxid
        self.taxid = None
        if results == "primer":
            self.mode = "primer"
            self.primer_dir = os.path.join(self.results_dir, "primer")
            self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
            self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
            self.maxgroupsize = 25000

    def blastresult_files(self, blast_dir):
        blastresults = []
        G.logger("Run: blastresults_files(" + self.target + ")")
        for files in [f for f in os.listdir(blast_dir) if f.endswith(".csv")]:
            file_path = os.path.join(blast_dir, files)
            if file_path not in blastresults:
                blastresults.append(file_path)
        return blastresults

    def get_alignmentdata_old(self, alignment, exceptions):
        if "gnl|BL_ORD_ID|" in alignment.hit_id:
            error_msg = (
                "Problem with custom DB, Please use the '-parse_seqids'"
                " option for the makeblastdb command")
            print("\n" + error_msg + "\n")
            G.logger("> " + error_msg)
            errors.append([self.target, error_msg])
            raise BlastDBError(error_msg)

        if alignment.hit_def == "No definition line":
            error_msg = (
                "Error: No definition line in " + alignment.title +
                "\nData is missing in the custom BLAST DB. At least "
                "a unique sequence identifier and the species name "
                "is required for each entry.\nExpected format: "
                ">seqid species name optional description")
            print("\n" + error_msg + "\n")
            G.logger("> " + error_msg)
            errors.append([self.target, error_msg])
            raise BlastDBError(error_msg)


        if "gi|" in alignment.hit_id:
            gi = alignment.hit_id.split("gi|")[1].split("|")[0]
        else:
            gi = alignment.accession

        db_id = alignment.accession
        lname = alignment.hit_def
        name = lname.split(" ")
        if self.config.virus:
            identity = "unknown"
            for ex in exceptions:
                if ex in lname:
                    test = str(" ".join(name[0:2]))
                    if " ".join(ex.split(" ")[0:2]) == test:
                        identity = " ".join(self.target.split("_"))
            if identity == "unknown":
                if "," in lname:
                    identity = lname.split(",")[0]
                else:
                    identity = lname

        elif not "PREDICTED" in lname:
            if len(name) >= 3:
                if "subsp" in str(" ".join(name)):
                    identity = str(" ".join(name[0:4]))
                elif "pv" in str(" ".join(name)):
                    identity = str(" ".join(name[0:4]))
                else:
                    identity = str(" ".join(name[0:2]))
            else:
                identity = str(" ".join(name[0:2]))

        else:
            identity = None

        for hsp in alignment.hsps:
            score = hsp.score
            e_value = hsp.expect
            query = hsp.query
            match = hsp.match
            subject = hsp.sbjct
            subject_start = hsp.sbjct_start
            align_length = hsp.align_length
            nuc_ident = hsp.identities

            if (identity and gi and db_id and score and e_value) is not None:
                return (
                    identity, gi, db_id, score, e_value, query, match,
                    subject, subject_start, align_length, nuc_ident)

    def get_seq_ends(self, blast_record, alignment, query_start, query_end):
        for hsp in alignment.hsps:
            # check start of sequence for not aligned parts
            qstart = hsp.query_start
            aln_len = hsp.align_length
            query_start.append(qstart)
            # check end of sequence for not aligned parts
            end_len = (blast_record.query_letters - (qstart + aln_len - 1))
            query_end.append(end_len)

    def check_seq_ends_old(self, blast_record, query_start, query_end):
        seq_ends = []
        if len(query_start) > 0:
            if min(query_start) >= self.config.minsize:
                seq_range = "[1:" + str(min(query_start)) + "]"
                seq_ends.append([blast_record.query, seq_range])
        if len(query_end) > 0:
            if min(query_end) >= self.config.minsize:
                qletters = str(blast_record.query_letters)
                seq_range = "[" + str(min(query_end)) + ":" + qletters + "]"
                seqlen = int(qletters) - min(query_end)
                if seqlen >= self.config.minsize:
                    seq_ends.append([blast_record.query, seq_range])
        if len(seq_ends) > 0:
            filename = os.path.join(self.blast_dir, "partialseqs.txt")
            with open(filename, "w") as f:
                for end in seq_ends:
                    gi = str(end[0])
                    s_range = str(end[1])
                    f.write(gi + " " + s_range + "\n")

    def write_nontargethits(
            self, dir_path, align_dict, result_format="json",
            name="nontargethits"):
        G.logger("Run: write_nontargethits(" + self.target + ")")
        file_path = os.path.join(dir_path, name + ".json")
        if result_format == "json":
            with open(file_path, "w") as q:
                q.write(json.dumps(align_dict))
        else:
            pass
            # csv

    def changed_primer3_input(self, file_path, controlfile_path):

        def find_difference():
            new = []
            old = []
            with open(file_path) as n:
                for line in n:
                    if "SEQUENCE_ID=" in line:
                        if line.strip() not in new:
                            new.append(line.strip())
                    if "PRIMER_PICK_INTERNAL_OLIGO=" in line:
                        if line.strip() not in new:
                            new.append(line.strip())

            with open(controlfile_path) as o:
                for line in o:
                    if "SEQUENCE_ID=" in line:
                        if line.strip() not in old:
                            old.append(line.strip())
                    if "PRIMER_PICK_INTERNAL_OLIGO=" in line:
                        if line.strip() not in old:
                            old.append(line.strip())

            diff = list(set(new) ^ set(old))
            return diff

        if os.path.isfile(controlfile_path):
            diff = find_difference()
            if len(diff) > 0:
                info1 = (
                    "Due to changed settings primer design "
                    "and quality control will start from scratch")
                info2 = "Differences in primer3 input:"
                for info in [info1, info2, diff]:
                    G.logger(info)
                    print(info)
                primer_dir = os.path.join(self.results_dir, "primer")
                if os.path.isdir(primer_dir):
                    G.logger("Delete primer directory")
                    print("Delete primer directory")
                    shutil.rmtree(primer_dir)

                shutil.copy(file_path, controlfile_path)
        else:
            shutil.copy(file_path, controlfile_path)

    def write_primer3_input(self, selected_seqs, conserved_seq_dict):
        G.logger("Run: write_primer3_input(" + self.target + ")")
        file_path = os.path.join(self.results_dir, "primer3_input")
        controlfile_path = os.path.join(self.results_dir, ".primer3_input")
        if self.config.probe is True:
            probe = "\nPRIMER_PICK_INTERNAL_OLIGO=1"
        else:
            probe = ""
        with open(file_path, "w") as f:
            for item in selected_seqs:
                if "complete" in item[1]:
                    f.write(
                        "SEQUENCE_ID=" + item[0] + "\nSEQUENCE_TEMPLATE="
                        + str(conserved_seq_dict[item[0]])
                        + "\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize) + "-"
                        + str(self.config.maxsize) + probe + "\n=\n")
                else:
                    x = str(conserved_seq_dict[item[0]])
                    i = item[1].strip("[]").split(":")
                    subseq = x[int(i[0])-1:int(i[1])-1]
                    f.write(
                        "SEQUENCE_ID="+item[0] + "\nSEQUENCE_TEMPLATE="
                        + str(subseq)
                        + "\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize) + "-"
                        + str(self.config.maxsize) + probe + "\n=\n")

        self.changed_primer3_input(file_path, controlfile_path)


    def get_excluded_gis(self):
        excluded_gis = []
        gi_file = os.path.join(self.config_dir, "no_blast.gi")
        if os.path.isfile(gi_file):
            if os.stat(gi_file).st_size > 0:
                with open(gi_file, "r") as f:
                    for line in f:
                        gi = line.strip()
                        excluded_gis.append(str(gi))
        return excluded_gis



    def parse_blastrecords_old(self, blast_record):

        def add_align_dict_data(aln_data, alignment):
            [
                identity, gi, db_id, score, e_value, query, match,
                subject, subject_start, align_length, nuc_ident] = aln_data
            perc_coverage = round(100/query_length * align_length, 0)
            perc_ident = round(100/align_length * nuc_ident, 0)
            ids = {identity: {
                "gi": gi, "db_id": db_id, "score": score,
                "e_value": e_value, "query": query,
                "match": match, "subject": subject,
                "subject_start": subject_start,
                "perc_coverage": perc_coverage,
                "perc_ident": perc_ident,
                "query_length": query_length}}
            if ids not in hits:
                hits.append(ids)
            if self.mode == "normal":
                self.get_seq_ends(
                    blast_record, alignment,
                    query_start, query_end)

        align_dict = {}
        hits = []
        query_start = []
        query_end = []
        query_length = blast_record.query_length
        exceptions = []
        if not self.exception == []:
            for item in self.exception:
                exception = ' '.join(item.split("_"))
                if exception not in exceptions:
                    exceptions.append(exception)

        if len(blast_record.alignments) == 0:
            align_dict.update({blast_record.query: {}})

        for alignment in blast_record.alignments:
            align_dict.update({blast_record.query: {}})

            aln_data = self.get_alignmentdata(alignment, exceptions)

            if aln_data:
                if self.config.nolist:
                    targetspecies = " ".join(str(self.target).split("_"))
                    if "subsp" in self.target:
                        targetspecies = (
                            "subsp.".join(targetspecies.split("subsp")))
                    if "pv" in self.target:
                        targetspecies = (
                            "pv.".join(targetspecies.split("pv")))

                    if not (
                        str(aln_data[0]) == str(targetspecies) or
                        str(aln_data[0]) in exceptions
                    ):
                        add_align_dict_data(aln_data, alignment)
                else:
                    if str(aln_data[0]) not in exceptions:
                        for species in self.nontargetlist:
                            if str(aln_data[0]) == str(species):
                                add_align_dict_data(aln_data, alignment)

                align_dict.update({blast_record.query: hits})
        if self.mode == "normal":
            self.check_seq_ends(blast_record, query_start, query_end)
        return align_dict

    def get_selected_sequences(self, align_dict):
        seq_ends = []
        selected_seqs = []
        excluded_seqs = []
        G.logger("Run: get_selected_sequences( " + self.target + ")")
        filename = os.path.join(self.blast_dir, "partialseqs.txt")
        if os.path.isfile(filename):
            with open(filename, "r") as f:
                for line in f:
                    line = line.strip()
                    gi = line.split(" ")[0]
                    s_range = line.split(" ")[1]
                    seq_ends.append([gi, s_range])
        for key in align_dict:
            if len(align_dict[key]) == 0:
                selected_seqs.append([key, "complete"])
            else:
                excluded_seqs.append(key)
            for item in seq_ends:
                if not (
                        [item[0], "complete"] in selected_seqs or
                        item in selected_seqs):
                    selected_seqs.append(item)
                    if item[0] in excluded_seqs:
                        excluded_seqs.remove(item[0])
        info1 = "selected sequences: " + str(len(selected_seqs))
        info2 = "excluded sequences: " + str(len(excluded_seqs))
        G.logger("> " + info1)
        G.logger("> " + info2)
        print("\n" + info1)
        print("\n" + info2)
        return selected_seqs

    def create_posdict(self, nonred_dict):
        from collections import defaultdict
        posdict = defaultdict(list)
        for key in nonred_dict.keys():
            if not len(nonred_dict[key]) == 0:
                for species in nonred_dict[key]:
                    poskey = nonred_dict[key][species]['main_id']
                    pos = int(
                        nonred_dict[key][species]["subject_start"])
                    if pos not in posdict[poskey]:
                        posdict[poskey].append(pos)
        return posdict

    def sort_nontarget_sequences(self, nonred_dict):
        nonreddata = []
        filepath = os.path.join(self.primer_qc_dir, "primerBLAST_DBIDS.csv")
        overhang = 2000
        if os.path.isfile(filepath):
            with open(filepath, "r") as f:
                reader = csv.reader(f)
                next(reader, None)
                for row in reader:
                    nonreddata.append(row)
            return nonreddata

        posdict = self.create_posdict(nonred_dict)
        for key in posdict.keys():
            posdict[key].sort()
            inrange = []
            for index, item in enumerate(posdict[key]):
                if index == 0:
                    stop = item + overhang
                    if item > overhang:
                        start = item - overhang
                    else:
                        start = 1
                    inrange.append([key, start, stop])
                else:
                    if item < stop:
                        stop = item + overhang
                        inrange.append([key, start, stop])
                        if index == len(posdict[key]) - 1:
                            nonreddata.append(inrange[-1])
                    else:
                        nonreddata.append(inrange[-1])
                        inrange = []
                        if item > overhang:
                            start = item - overhang
                        else:
                            start = 1
                        stop = item + overhang
                        inrange.append([key, start, stop])

        if len(nonreddata) > 0:
            header = ["Accession", "Start pos", "Stop pos"]
            G.csv_writer(filepath, nonreddata, header)

        return nonreddata

    def write_nontarget_sequences(self, nonreddata):
        maxsize = 25000
        if self.config.customdb is not None:
            db = self.config.customdb
        else:
            # the version 5 databases will no longer have
            # "_v5" as part of the database names
            # https://ncbiinsights.ncbi.nlm.nih.gov/2020/01/28/blast-db-ftp/
            db = "nt"
        seqcount = len(nonreddata)
        info = "Found " + str(seqcount) + " sequences for the non-target DB"
        G.logger(info)
        print("\n" + info + "\n")
        parts = seqcount//maxsize + 1

        for part in range(0, parts):
            print("Working on part " + str(part + 1) + "/" + str(parts))
            end = (part+1)*maxsize
            if end > seqcount:
                end = seqcount
            filename = "BLASTnontarget" + str(part) + ".sequences"
            filepath = os.path.join(self.primer_qc_dir, filename)
            if not os.path.isfile(filepath):
                info = "Start writing " + filename
                print(info)
                G.logger(info)
                print("Start DB extraction")
                G.logger("Start DB extraction")
                data = nonreddata[part*maxsize:end]
                fasta_seqs = G.run_parallel(
                        P.get_seq_fromDB, data, db)
                try:
                    with open(filepath, "w") as r:
                        for fastainfo in fasta_seqs:
                            name = fastainfo[0].split(":")
                            fastainfo[0] = (
                                name[0] + "_" + "_".join(name[1].split("-")))
                            fastadata = "\n".join(fastainfo) + "\n"
                            r.write(fastadata)
                except (KeyboardInterrupt, SystemExit):
                    G.keyexit_rollback("DB extraction", fp=filepath)
                    raise

            else:
                info2 = "Skip writing " + filename
                print(info2)
                G.logger(info2)
                return
            info3 = "Finished writing " + filename
            print(info3)
            G.logger(info3)

    def get_primerBLAST_DBIDS(self, nonred_dict):
        print("\nGet sequence accessions of BLAST hits\n")
        G.logger("> Get sequence accessions of BLAST hits")
        G.create_directory(self.primer_qc_dir)
        nonreddata = self.sort_nontarget_sequences(nonred_dict)
        if len(nonreddata) == 0:
            msg = (
                "Error did not find any sequences for non-target DB. "
                + "Please check the species list and/or BLAST database")
            print(msg)
            G.logger("> " + msg)
            errors.append([self.target, msg])
            return 1

        self.write_nontarget_sequences(nonreddata)
        return 0

    def remove_redundanthits(self, align_dict):
        nonred_dict = {}
        for key in align_dict.keys():
            n = 0
            conv_seq = []
            summarydict = {}
            for hit in align_dict[key]:
                species = list(hit.keys())[0]
                (
                    gi, db_id, coverage, identity, query, sub_start,
                    evalue, score, subject, match, query_length
                    ) = (
                        hit[species]['gi'], hit[species]['db_id'],
                        hit[species]['perc_coverage'],
                        hit[species]['perc_ident'],
                        hit[species]['query'], hit[species]['subject_start'],
                        hit[species]['e_value'], hit[species]['score'],
                        hit[species]['subject'], hit[species]['match'],
                        hit[species]['query_length'])
                if [species, subject] in conv_seq:
                    summarydict[species]["gi_ids"].append(gi)
                else:
                    conv_seq.append([species, subject])
                    if species in summarydict.keys():
                        mainspecies = species + "_" + str(n)
                        n = n + 1
                    else:
                        mainspecies = species

                    summ_info = {
                        mainspecies: {
                            "title": species, "query": query, "score": score,
                            "main_id": gi, "identity": identity,
                            "coverage": coverage, "sequence": subject,
                            "query_length": query_length, "evalue": evalue,
                            'subject_start': sub_start, "gi_ids": []}}
                    summarydict.update(summ_info)
            nonred_dict.update({key: summarydict})
        return nonred_dict

    def bp_read_nontarget_hits(self, file_path, excluded_gis):
        align_dict = {}
        info = "Read nontargethits"
        print(info)
        G.logger("> " + info)
        with open(file_path, 'r') as f:
            for line in f:
                align_dict = json.loads(line)
        if len(excluded_gis) > 0:
            ex_gis = []
            deletekey = []
            for key in align_dict.keys():
                for species in align_dict[key].keys():
                    gi = align_dict[key][species]['main_id']
                    if str(gi) in excluded_gis:
                        if gi not in ex_gis:
                            ex_gis.append(gi)
                        deletekey.append([key, species])
            for item in deletekey:
                del align_dict[item[0]][item[1]]

            if len(ex_gis) > 0:
                info = "removed GI's in excluded GI list from results"
                G.logger(info, ex_gis)
                print("\n" + info, ex_gis)

        return align_dict

    def bp_parse_results(self, blast_dir):
        blastresults = self.blastresult_files(blast_dir)
        nr = 1
        for filename in blastresults:
            print(
                "\nopen BLAST result file " + str(nr)
                + "/" + str(len(blastresults)))

            blastdf = pd.read_csv(filename, sep="\t", header=None)
            blastdf.columns = self.header
            check_blastdb_errors(blastdf)

    def check_blastdb_errors(blastdf):
        error_msg = None
        if len(blastdf.index) == 0:
            error_msg = " ".join([
                "A problem with the BLAST results file",
                filename, "was detected.",
                "Please check if the file was removed and start the run again"])
            os.remove(filename)
            print("removed " + filename)

        if blastdf["Subject Seq-id"].str.contains("gnl|BL_ORD_ID|"):
            error_msg = (
                "Problem with custom DB, Please use the '-parse_seqids'"
                " option for the makeblastdb command")

        if blastdf["Subject Title"].str.contains("No definition line"):
            error_msg = (
                "Error: No definition line in " + alignment.title +
                "\nData is missing in the custom BLAST DB. At least "
                "a unique sequence identifier and the species name "
                "is required for each entry.\nExpected format: "
                ">seqid species name optional description")

        if error_msg:
            print("\n" + error_msg + "\n")
            logging.error("> " + error_msg, exc_info=True)
            errors.append([self.target, error_msg])
            raise BlastDBError(error_msg)

    def check_seq_ends(rejected):
        rejected["Query overhang"] = rejected['End of alignment in query'] - rejected['Alignment length']
        rejected[rejected['Start of alignment in query'] >= self.config.minsize]

        partials_min = rejected.query(
                'Start of alignment in query >= "@self.config.minsize"')
        partials_max = rejected.query('Query overhang >= "@self.config.minsize"')

        partials_min["range"] = "[1:" + partials_min['Start of alignment in query'] + "]"
        partials_max["range"] = (
                "[" + partials_max['End of alignment in query']
                + ":" + partials_max['Query sequence length'])

        partial_seqs = pd.concat([partials_min, partials_max])
        filename = os.path.join(self.blast_dir, "partialseqs.txt")
        if len(partial_seqs.index) != 0:
            partial_seqs[["Query Seq-id", "range"]].to_csv(filename, sep=" ")

    def parse_blastrecords(self, blastdf):
        align_dict = {key: {} for key in set(blastdf["Query Seq-id"])}
        exceptions = [self.target]
        if self.exception != []:
            for item in self.exception:
                exception = ' '.join(item.split("_"))
                if exception not in exceptions:
                    exceptions.append(exception)

        if self.taxid is None:
            target_hits = blastdf[
                            blastdf["Subject Title"].str.contains("|".join(exceptions))]
            offtarget_hits = blastdf[
                        ~blastdf["Subject Title"].str.contains("|".join(exceptions))]
            accept = target_hits
            reject = offtarget_hits

            if self.config.nolist is False:
                listed_offtarget = offtarget_hits[
                                        offtarget_hits["Subject Title"].str.contains(
                                                            "|".join(self.nontargetlist))]
                unlisted = offtarget_hits[
                                        ~offtarget_hits["Subject Title"].str.contains(
                                                            "|".join(self.nontargetlist))]
                accept = pd.concat([target_hits, unlisted])
                reject = listed_offtarget

            excluded_gis = self.get_excluded_gis()
            reject = reject[~reject["sgi"].str.contains("|".join(excluded_gis))]

        else:
            accept = blastdf.query('Subject Taxonomy ID(s) == "@self.taxid"')
            reject = blastdf.query('Subject Taxonomy ID(s) != "@self.taxid"')

            if len(self.target.split(" ")) == 2:
                sub_sp = reject[reject['Subject Title'].str.contains("|".join(["subsp.", "pv."]))]
                sub_check = sub_sp[sub_sp['Subject Title'].str.contains(" ".join(self.target[0:2]))]
                accept = pd.concat([accept, sub_check])
                reject = pd.concat([accept, reject]).drop_duplicates(keep=False)

        if self.mode == "normal":
            # check if sequence ends without off-target matches exist
            self.check_seq_ends(reject)

    def bp_parse_xml_files(self, blast_dir):
        align_dict = {}
        xmlblastresults = self.blastresult_files(blast_dir)
        nr = 1
        for filename in xmlblastresults:
            print(
                "\nopen BLAST result file " + str(nr)
                + "/" + str(len(xmlblastresults)))
            blastrecords = self.parse_BLASTfile(filename)
            print(
                "read BLAST results file " + str(nr)
                + "/" + str(len(xmlblastresults)))
            nr = nr + 1
            total = len(blastrecords)
            rec = 1
            for record in blastrecords:
                print(
                    '\r read record ' + str(rec) + "/" + str(total), end=''
                )
                result = self.parse_blastrecords(record)
                align_dict.update(result)
                rec = rec + 1

        if self.config.intermediate is False:
            for file_name in xmlblastresults:
                os.remove(file_name)
                name = os.path.basename(file_name)
                n = name.split("_")
                filename = n[0] + ".part-" + n[1]
                filepath = os.path.join(blast_dir, filename)
                os.remove(filepath)

        return align_dict

    def blast_parser(self, blast_dir):
        G.logger("Run: blast_parser")
        file_path = os.path.join(blast_dir, "nontargethits.json")
        excluded_gis = self.get_excluded_gis()
        if os.path.isfile(file_path):
            align_dict = self.bp_read_nontarget_hits(file_path, excluded_gis)
        else:
            align_dict = self.bp_parse_xml_files(blast_dir)

            if len(excluded_gis) > 0:
                ex_gis = []
                for key in align_dict.copy().keys():
                    for index, item in enumerate(
                        align_dict.copy()[key]
                    ):
                        for species in item.copy().keys():
                            gi = item[species]['gi']
                            if str(gi) in excluded_gis:
                                if gi not in ex_gis:
                                    ex_gis.append(gi)
                                del align_dict[key][index]

                if len(ex_gis) > 0:
                    info = "removed GI's in excluded GI list from results"
                    info2 = ex_gis
                    G.logger(info)
                    G.logger(info2)
                    print("\n" + info)
                    print(info2)

        return align_dict

    def commonhit_counter(self, nonred):
        import collections
        uniq_count = []
        hitlist = []
        c = collections.Counter()
        output = os.path.join(self.blast_dir, "mostcommonhits.csv")
        keycount = 0
        for key in nonred.keys():
            keycount = keycount + 1
            for species in nonred[key]:
                try:
                    duplicate = int(species.split("_")[-1])
                    hitspec = species.split("_")[0:-1][0]
                except ValueError:
                    hitspec = species
                main = nonred[key][species]['main_id']
                hitlist.append([main, hitspec])
                for gi in nonred[key][species]['gi_ids']:
                    hitlist.append([gi, hitspec])

        for item in hitlist:
            c.update([item[0] + " " + item[1]])
        c_sorted = sorted(c.most_common())
        for key, val in c_sorted:
            uniq_count.append([key, val])

        uniq_count.sort(key=lambda x: int(x[1]), reverse=True)
        results = []
        total_hits = sum(x[1] for x in uniq_count)
        for uniq in uniq_count:
            perc = round(100 / keycount * int(uniq[1]), 2)
            if float(perc) > 1:
                gi = uniq[0].split(" ")[0]
                spec = " ".join(uniq[0].split(" ")[1::])
                results.append([gi, spec, perc, uniq[1]])

        header = [
            ["Total BLAST hits", "Number of queries"], [total_hits, keycount],
            ["GI", "Species", "BLAST hits [%]", "BLAST hits [count]"]]
        G.csv_writer(output, results, header)


    def filter_nonreddict(self, nonred_dict):
        filtdict = {}
        for key, val in nonred_dict.items():
            filtdict.update({key: {}})
            for item in val:
                nident = val[item]["identity"]
                evalue = val[item]["evalue"]
                if (
                    int(nident) >= self.config.nuc_identity and
                    float(evalue) <= self.config.evalue
                ):
                    filtdict[key] = val

        return filtdict


    def run_blastparser(self, conserved_seq_dict):
        if self.mode == "primer":
            filepath = os.path.join(self.primerblast_dir, "nontargethits.json")
            G.logger("Run: run_blastparser(" + self.target + "), primer")
            align_dict = self.blast_parser(self.primerblast_dir)
            if not os.path.isfile(filepath):
                nonred_dict = self.remove_redundanthits(align_dict)
                self.write_nontargethits(
                        self.primerblast_dir, nonred_dict,
                        "json")
            else:
                nonred_dict = align_dict

            self.get_primerBLAST_DBIDS(nonred_dict)

            duration = time.time() - self.start
            G.logger(
                "> Primer blast parser time: "
                + str(timedelta(seconds=duration)).split(".")[0])
        else:
            filepath = os.path.join(self.blast_dir, "nontargethits.json")
            G.logger(
                "Run: run_blastparser(" + self.target
                + "), conserved sequences")
            align_dict = self.blast_parser(self.blast_dir)
            if not os.path.isfile(filepath):
                nonred_dict = self.remove_redundanthits(align_dict)
                self.commonhit_counter(nonred_dict)
                self.write_nontargethits(
                        self.blast_dir, nonred_dict,
                        "json")
            else:
                nonred_dict = align_dict

            filtered_dict = self.filter_nonreddict(nonred_dict)

            selected_seqs = self.get_selected_sequences(filtered_dict)
            self.write_primer3_input(selected_seqs, conserved_seq_dict)
            duration = time.time() - self.start
            info = ("species specific conserved sequences: "
                    + str(len(selected_seqs)))
            G.logger(
                "> Blast parser time: "
                + str(timedelta(seconds=duration)).split(".")[0])
            G.logger("> " + info)
            PipelineStatsCollector(self.target_dir).write_stat(info)

            if len(selected_seqs) == 0:
                msg = "No conserved sequences without non-target match found"
                print(msg)
                G.logger("> " + msg)
                errors.append([self.target, msg])
                return 1

            return 0


class PrimerDesign():
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.primer_dir = os.path.join(self.results_dir, "primer")
        self.p3dict = {}

    def run_primer3(self):
        G.logger("Run: run_primer3(" + self.target + ")")
        input_file = os.path.join(self.results_dir, "primer3_input")
        output_file = os.path.join(self.primer_dir, "primer3_output")
        settings_file = os.path.join(dict_path, "p3parameters")
        if not os.path.isfile(output_file):
            primer3cmd = [
                "primer3_core", "-p3_settings_file=" + settings_file,
                "-echo_settings_file", "-output=" + output_file, input_file]
            try:
                G.run_subprocess(primer3cmd, True, True, True)
            except (KeyboardInterrupt, SystemExit):
                G.keyexit_rollback("primer3 run", fp=output_file)
                raise
        else:
            info = "Skip primerdesign with primer3"
            G.logger("> " + info)
            print(info)

    def parse_Primer3_output(self, p3_output):
        settings_file = os.path.join(pipe_dir, "p3parameters")

        def parseSeqId(key, value):
            if key.startswith("SEQUENCE_ID"):
                if "group" in value:
                    if value.startswith("group_"):
                        value = "g" + value.split("group_")[1]
                    else:
                        spval = value.split("group_")
                        value = spval[0] + "g" + spval[1]

                seq_id = value
                self.p3dict.update({seq_id: {"Primer_pairs": None}})
                p3list.append(seq_id)

        def parseTemplate(key, value):
            if key.startswith("SEQUENCE_TEMPLATE"):
                self.p3dict[p3list[-1]].update(
                    {"template_seq": value})

        def countPrimer(key, value):
            if key.startswith("PRIMER_PAIR_NUM_RETURNED"):
                self.p3dict[p3list[-1]]["Primer_pairs"] = int(value)
                for i in range(0, int(value)):
                    self.p3dict[p3list[-1]].update({"Primer_pair_"+str(i): {}})

        def parsePrimer(key, value, pos):
            pp = "Primer_pair_"
            if key.startswith("PRIMER_" + pos):
                i = key.split("_")[2]
                if key.endswith("_PENALTY"):
                    primer_rpen = value
                    self.p3dict[p3list[-1]][pp + str(i)].update(
                        {"primer_" + pos[0] + "_penalty": float(primer_rpen)})
                if key.endswith("_SEQUENCE"):
                    primer_rseq = value
                    self.p3dict[p3list[-1]][pp + str(i)].update(
                        {"primer_" + pos[0] + "_sequence": primer_rseq})
                if key.endswith("_TM"):
                    right_TM = value
                    self.p3dict[p3list[-1]][pp + str(i)].update(
                        {"primer_" + pos[0] + "_TM": float(right_TM)})

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

        def read_primeroutput(p3_output):
            settings = 1
            primerdatasets = []
            primerdata = []
            primererror = []
            with open(p3_output) as f:
                for line in f:
                    line = line.strip()
                    if "PRIMER_ERROR=Cannot open " in line:
                        msg = (
                            "primer3 cannot open the settingsfile: "
                            + settings_file)
                        G.logger(">" + msg)
                        errors.append([self.target, msg])
                        raise Exception
                    if "P3_SETTINGS_FILE_END=" in line:
                        settings = 0
                    if settings == 0:
                        if not (
                            line == "=" or "P3_SETTINGS_FILE_END=" in line
                        ):
                            if "PRIMER_ERROR=" in line:
                                primerdata.append(line)
                                primererror.append(primerdata)
                                primerdata = []
                            else:
                                primerdata.append(line)
                    if line == "=":
                        primerdatasets.append(primerdata)
                        primerdata = []

            if primererror != []:
                errfile = os.path.join(self.primer_dir, "primer3_errors.csv")
                msg = "Detected errors during primer3 run, check:\n" + errfile
                G.logger(">" + msg)
                print("\n" + msg + "\n")
                errors.append([self.target, msg])
                G.csv_writer(errfile, primererror)
            return primerdatasets

        p3list = []
        info = "Run: parse_Primer3_output(" + self.target + ")"
        print(info)
        G.logger(info)
        primerdatasets = read_primeroutput(p3_output)
        for primerdata in primerdatasets:
            for item in primerdata:
                key = item.split("=")[0]
                value = item.split("=")[1]
                parseSeqId(key, value)
                parseTemplate(key, value)
                countPrimer(key, value)
                parsePrimer(key, value, "RIGHT")
                parsePrimer(key, value, "LEFT")
                parsePrimer(key, value, "INTERNAL")
                parsePrimerPair(key, value)

    def get_amplicon_seq(self):
        def PCR(left, rc_right, temp):
            pcr_product = (
                temp[temp.index(left):temp.index(rc_right)] + rc_right)
            return pcr_product

        for key in self.p3dict.keys():
            if self.p3dict[key]["Primer_pairs"] > 0:
                template = self.p3dict[key]["template_seq"]
                for pp in self.p3dict[key].keys():
                    if "Primer_pair_" in pp:
                        lprimer = self.p3dict[key][pp]["primer_L_sequence"]
                        rprimer = self.p3dict[key][pp]["primer_R_sequence"]
                        rev_compl = str(Seq(rprimer).reverse_complement())
                        pcr_product = PCR(lprimer, rev_compl, template)
                        self.p3dict[key][pp].update(
                            {"amplicon_seq": pcr_product})

    def get_annotation_info(self):
        annot_dict = {}
        gpa_file = os.path.join(
                self.pangenome_dir, "gene_presence_absence.csv")
        with open(gpa_file, "r") as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                gene = row[0]
                annotation = row[2]
                if not gene in annot_dict.keys():
                    if "group" in gene:
                        if gene.startswith("group_"):
                            gene = "g" + gene.split("group_")[1]
                    annot_dict.update({gene: annotation})

        return annot_dict

    def update_annotation_info(self):
        annot_dict = self.get_annotation_info()
        for key in self.p3dict.keys():
            genename = "_".join(key.split("_")[0:-1])
            try:
                annotation = annot_dict[genename]
            except KeyError:
                for i in range(1, len(key.split("_"))):
                    geneid = "_".join(key.split("_")[i::])
                    if geneid in annot_dict.keys():
                        annotation = annot_dict[geneid]
            self.p3dict[key]["annotation"] = annotation

    def write_primer3_data(self):
        file_path = os.path.join(self.primer_dir, "primer3_summary.json")
        with open(file_path, "w") as f:
            f.write(json.dumps(self.p3dict))

    def run_primerdesign(self):
        info = "Start primer design"
        print(info)
        G.logger("> " + info)
        G.logger("Run: run_primerdesign(" + self.target + ")")
        G.create_directory(self.primer_dir)
        self.run_primer3()
        p3_output = os.path.join(self.primer_dir, "primer3_output")
        self.parse_Primer3_output(p3_output)
        self.get_amplicon_seq()
        self.update_annotation_info()
        self.write_primer3_data()
        return self.p3dict


class PrimerQualityControl:
    def __init__(self, configuration, primer3_dict):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.primer_dir = os.path.join(self.results_dir, "primer")
        self.summ_dir = os.path.join(self.config.path, "Summary", self.target)
        self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
        self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
        self.mfold_dir = os.path.join(self.primer_dir, "mfold")
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.primer3_dict = primer3_dict
        self.call_blastparser = BlastParser(self.config, "primer")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.primerlist = []
        self.start = time.time()
        self.mfethreshold = self.config.mfethreshold
        self.referencegenomes = 10
        self.dbinputfiles = []

    def collect_primer(self):
        G.logger("Run: collect_primer(" + self.target + ")")
        primer_sorted = []
        for key in self.primer3_dict:
            if not self.primer3_dict[key]["Primer_pairs"] == 0:
                for item in self.primer3_dict[key]:
                    if "Primer_pair_" in item:
                        primer_sorted.append([
                            key, item,
                            self.primer3_dict[key][item]["primer_P_penalty"]])

        primer_sorted.sort(key=lambda x: float(x[2]))
        for count, item in enumerate(primer_sorted):
            self.get_blast_input(item)

        if len(self.primerlist) == 0:
            info = "potential Primer pair(s): 0"
            error_msg = "No primers found for conserved sequences"
            PipelineStatsCollector(self.target_dir).write_stat(info)
            print("\n" + error_msg)
            G.logger("> " + error_msg)
            G.logger("> " + info)
            errors.append([self.target, error_msg])
            return 1

        info = (
            "Number of potential primer pair(s) found "
            + str(len(self.primerlist)//2))
        print("\n" + info + "\n")
        G.logger("> " + info)
        PipelineStatsCollector(self.target_dir).write_stat(
            "potential primer pair(s): "
            + str(len(self.primerlist)//2))
        return 0

    def get_blast_input(self, item):
        # start new primername definition here
        # start 15.12.2017
        # get primername without direction split("_")
        p_fwd_name = (
            H.abbrev(self.target) + "_" + item[0]
            + "_P" + item[1].split("_")[-1]) + "_F"
        p_rev_name = (
            "_".join(p_fwd_name.split("_")[0:-1]) + "_R")
        p_fwd_seq = self.primer3_dict[item[0]][item[1]]['primer_L_sequence']
        p_rev_seq = self.primer3_dict[item[0]][item[1]]['primer_R_sequence']
        self.primerlist.append([p_fwd_name, p_fwd_seq])
        self.primerlist.append([p_rev_name, p_rev_seq])

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
                target_id = "_".join(primer_name.split("_")[-3:-1])
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
                G.logger(
                    "error in get_primerinfo()"
                    + str(sys.exc_info()))

        return val_list

    def create_template_db_file(self, primerinfos):
        wrote = []
        file_path = os.path.join(self.primer_qc_dir, "template.sequences")
        with open(file_path, "w") as f:
            for nameF, seqF, nameR, seqR, templ_seq in primerinfos:
                if templ_seq not in wrote:
                    primername = "_".join(nameF.split("_")[0:-2])
                    f.write(">" + primername + "\n" + templ_seq + "\n")
                    wrote.append(templ_seq)

    def get_QC_data(self):
        qc_data = []
        if os.path.isdir(self.summ_dir):
            for files in os.listdir(self.summ_dir):
                if files.endswith("_qc_sequences.csv"):
                    file_path = os.path.join(self.summ_dir, files)
                    with open(file_path) as f:
                        reader = csv.reader(f)
                        next(reader, None)
                        for row in reader:
                            qc_data.append(row)
        return qc_data

    def find_QC_assemblies(self):
        qc_data = self.get_QC_data()
        remove = []
        qc_acc = []
        assembly_dict = {}
        for item in qc_data:
            accession, assembly_stat, rRNA, tuf, recA, dnaK, pheS = (
                item[0], item[2], item[4], item[6],
                item[8], item[10], item[12])
            qc_list = rRNA, tuf, recA, dnaK, pheS
            assembly_dict.update({accession: assembly_stat})

            if self.config.ignore_qc is True:
                for qc_gene in qc_list:
                    if accession not in qc_acc:
                        qc_acc.append(accession)
            else:
                for qc_gene in qc_list:
                    if "passed QC" == qc_gene:
                        if accession not in qc_acc:
                            qc_acc.append(accession)
                    elif "" == qc_gene:
                        if accession not in qc_acc:
                            qc_acc.append(accession)
                    else:
                        if accession not in remove:
                            remove.append(accession)
        check = list(set(qc_acc) - set(remove))
        check.sort()

        return assembly_dict, check

    def create_assembly_db_file(self):
        # add option to choose a folder for Reference genomes?
        def assembly_selection(stat):
            for item in check:
                if len(ref_assembly) < self.referencegenomes:
                    if stat == "":
                        if item not in ref_assembly:
                            ref_assembly.append(item)
                    else:
                        if assembly_dict[item] == stat:
                            if item not in ref_assembly:
                                ref_assembly.append(item)
                else:
                    break

        ref_assembly = []

        if self.config.virus is True:
            assemblies = os.listdir(self.fna_dir)
            for assem in assemblies:
                if assem.endswith(".fna"):
                    ref_assembly.append("_".join(assem.split("_")[0:-1]))

        else:
            assembly_dict, check = self.find_QC_assemblies()
            assembly_stats = ["Complete Genome", "Chromosome", "Scaffold", ""]
            for stat in assembly_stats:
                assembly_selection(stat)

        ref_assembly.sort()
        target_fasta = []
        for files in os.listdir(self.fna_dir):
            for item in ref_assembly:
                acc = "v".join(item.split("."))
                if acc in files:
                    if files.endswith(".fna"):
                        with open(os.path.join(self.fna_dir, files)) as f:
                            records = SeqIO.parse(f, "fasta")
                            target_fasta.append(list(records))
        db_name = H.abbrev(self.target) + ".genomic"
        file_path = os.path.join(self.primer_qc_dir, db_name)
        with open(file_path, "w") as fas:
            for item in target_fasta:
                SeqIO.write(item, fas, "fasta")

    def prepare_MFEprimer_Dbs(self, primerinfos):
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

    def write_MFEprimer_results(self, input_list, name):
        outputlist = []
        val_list = []
        for item in input_list:
            # item[0] is the primerinfo
            if len(item[0]) > 1:
                if not item[0] in outputlist:
                    outputlist.append(item[0])
            # item[1] are the results of MFEprimer
            if len(item[1]) > 1:
                for values in item[1]:
                    val = values.split("\t")
                    val_list.append(val)

        G.csv_writer("MFEprimer_" + name + ".csv", val_list)

        return outputlist

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

        dbfile = os.path.join(
                self.primer_qc_dir, H.abbrev(self.target) + ".genomic")
        assembly_list = G.run_parallel(
                P.MFEprimer_assembly, check_assembly,
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

    def mfold_analysis(self, mfoldinputlist):
        info = "Start mfold analysis of PCR products"
        print(info)
        G.logger("> " + info)
        G.logger("Run: mfold_analysis(" + self.target + ")")
        print("Run: mfold_analysis(" + self.target + ")")

        if os.path.isdir(self.mfold_dir):
            shutil.rmtree(self.mfold_dir)
        G.create_directory(self.mfold_dir)
        os.chdir(self.mfold_dir)

        for mfoldinput in mfoldinputlist:
            abbr = H.abbrev(self.target)
            self.prep_mfold(mfoldinput, abbr)
        os.chdir(self.target_dir)

    def prep_mfold(self, mfoldinput, abbr):
        target_id, primerpair, amplicon_seq, primer_name = mfoldinput
        # This removes the Genus species string to shorten the
        # name for mfold (especially for subspecies names). mfold has
        # problems with too long filenames / paths
        short_name = primer_name.split(abbr + "_")[1]
        dir_path = os.path.join(self.mfold_dir, target_id)
        subdir_path = os.path.join(dir_path, primerpair)
        pcr_name = short_name + "_PCR"
        if len(amplicon_seq) >= self.config.minsize:
            G.create_directory(dir_path)
            G.create_directory(subdir_path)
            self.run_mfold(subdir_path, pcr_name, amplicon_seq)

    def run_mfold(self, subdir_path, seq_name, description):
        file_path = os.path.join(subdir_path, seq_name)
        with open(file_path, "w") as f:
            f.write("> " + seq_name + "\n" + description)
        if not os.path.isfile(file_path + ".det"):
            os.chdir(subdir_path)
            mfold_cmd = [
                "mfold", "SEQ=" + seq_name, "NA=DNA", "T=60", "MG_CONC=0.003"]
            G.run_subprocess(mfold_cmd, False, True, False)
            os.chdir(self.mfold_dir)

    def mfold_parser(self):

        def multiple_structures():
            selected = test[0][0]
            selected_primer.append(selected)
            for index, structure_nr in enumerate(mfold):
                selected_name, passed, excluded, failed = structure_nr
                filename, structure, dG, dH, dS, Tm = passed
                if index == 0:
                    first = selected
                else:
                    first = ""
                pos_results.append(
                    [first, structure, dG, dH, dS, Tm])

        def multiple_failed_structures():
            if mfold[0][2] is None:
                excluded = mfold[0][0]
            else:
                excluded = mfold[0][2]
            excluded_primer.append(excluded)
            for index, structure_nr in enumerate(mfold):
                selected, passed, excluded_name, failed = structure_nr
                if passed is None:
                    filename, structure, dG, dH, dS, Tm = failed
                    if index == 0:
                        first = excluded
                    else:
                        first = ""
                    neg_results.append(
                        [first, structure, dG, dH, dS, Tm])

                if failed is None:
                    filename, structure, dG, dH, dS, Tm = passed
                    if index == 0:
                        first = excluded
                    else:
                        first = ""
                    neg_results.append(
                        [first, structure, dG, dH, dS, Tm])

        selected_primer = []
        excluded_primer = []
        pos_results = []
        neg_results = []
        info = "Run: mfold_parser(" + self.target + ")\n"
        print(info)
        G.logger(info)
        file_list = self.find_mfold_results()
        for mfoldfiles in file_list:
            mfold = self.read_files(mfoldfiles)

            if len(mfold) == 1:
                selected, passed, excluded, failed = mfold[0]
                if not (selected is None or passed is None):
                    selected_primer.append(selected)
                    filename, structure, dG, dH, dS, Tm = passed
                    pos_results.append([selected, structure, dG, dH, dS, Tm])
                else:
                    excluded_primer.append(excluded)
                    filename, structure, dG, dH, dS, Tm = failed
                    neg_results.append([excluded, structure, dG, dH, dS, Tm])

            elif len(mfold) > 1:
                test = []
                for structure_nr in mfold:
                    selected, passed, excluded, failed = structure_nr
                    if not (selected is None or passed is None):
                        test.append([selected, passed, excluded, failed])

                if len(test) == len(mfold):
                    multiple_structures()
                else:
                    multiple_failed_structures()

        passfile = os.path.join(self.mfold_dir, "mfold_passed.csv")
        failfile = os.path.join(self.mfold_dir, "mfold_failed.csv")
        header = ["primer", "structure", "dG", "dH", "dS", "Tm"]
        G.csv_writer(passfile, pos_results, header)
        G.csv_writer(failfile, neg_results, header)

        ex = str(len(excluded_primer))
        pas = str(len(selected_primer))
        ex_info = ex + " primer pair(s) excluded by mfold"
        pass_info = pas + " primer pair(s) to continue"
        print("\n\n" + ex_info + "\n")
        print(pass_info)
        G.logger("> " + ex_info)
        G.logger("> " + pass_info)
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after mfold: "
            + str(len(selected_primer)))

        if self.config.intermediate is False:
            for dirs in os.listdir(self.mfold_dir):
                dir_path = os.path.join(self.mfold_dir, dirs)
                if os.path.isdir(dir_path):
                    shutil.rmtree(dir_path)

        return selected_primer, excluded_primer

    def find_mfold_results(self):
        file_list = []
        for root, dirs, files in os.walk(self.mfold_dir):
            for item in files:
                if item.endswith(".det"):
                    if os.path.join(root, item) not in file_list:
                        file_list.append(os.path.join(root, item))

        return file_list

    def parse_values(self, file, line):
        structure = line
        v = "".join(islice(file, 0, 1)).strip()
        dG = v.split("=")[1].split("d")[0].strip()
        dH = v.split("=")[2].split("d")[0].strip()
        dS = v.split("=")[3].split("T")[0].strip()
        Tm = v.split("=")[4].strip(" ").split(" ", 1)[0]
        try:
            y = float(Tm)
        except ValueError:
            Tm = "999"

        return structure, dG, dH, dS, Tm

    def interpret_values(self, name, primername, mfoldvalues):
        structure, dG, dH, dS, Tm = mfoldvalues
        mfold_output = [name, structure, dG, dH, dS, Tm]

        if float(dG) <= float(self.config.mfold):
            return [None, None, primername, mfold_output]

        return [primername, mfold_output, None, None]

    def get_primername(self, name):
        # adds the genus species info again to the
        # primername after mfold
        primer_name = (
            H.abbrev(self.target) + "_"
            + "_".join(name.split("_")[0:-1]))
        return primer_name

    def read_files(self, filename):
        results = []
        with open(filename, "r", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if re.search("Structure", line):
                    name = "".join(islice(f, 1, 2)).strip()
                    primername = self.get_primername(name)
                    mfoldvalues = self.parse_values(f, line)
                    results.append(
                        self.interpret_values(
                            name, primername, mfoldvalues))
        return results

    def dimercheck_primer(self, selected_seqs, excluded_primer):
        G.logger("Run: dimercheck_primer(" + self.target + ")")
        dimercheck = []
        for primer in selected_seqs:
            primername = "_".join(primer[0].split("_")[0:-1])
            if primername not in excluded_primer:
                if primer not in dimercheck:
                    dimercheck.append(primer)

        return dimercheck

    def check_primerdimer(self, dimercheck):
        info = "Start analysis of primer dimers"
        print(info)
        G.logger("> " + info)
        G.logger("Run: check_primerdimer(" + self.target + ")")

        def write_dimercheck_input(item, primer_name, lseq, rseq):
            fwd_name = primer_name + "_F"
            rev_name = primer_name + "_R"
            fwd_file = os.path.join(self.dimercheck_dir, fwd_name)
            rev_file = os.path.join(self.dimercheck_dir, rev_name)
            p_file = os.path.join(self.dimercheck_dir, primer_name)
            with open(fwd_file, "w") as f:
                f.write(
                    ">"+fwd_name+"\n"+lseq+"\n>"+fwd_name + "_d\n"+lseq+"\n")
            with open(rev_file, "w") as f:
                f.write(
                    ">"+rev_name+"\n"+rseq+"\n>"+rev_name + "_d\n"+rseq+"\n")
            with open(p_file, "w") as f:
                f.write(">"+fwd_name+"\n"+lseq+"\n>"+rev_name+"\n"+rseq+"\n")

        def get_dimercheck_output(item, primer_name, lseq, rseq, choice):
            test = []
            summary_data = []
            fwd_name = primer_name + "_F"
            rev_name = primer_name + "_R"
            fwd_file = os.path.join(
                self.dimercheck_dir, fwd_name + "_dimer_out")
            rev_file = os.path.join(
                self.dimercheck_dir, rev_name + "_dimer_out")
            p_file = os.path.join(
                self.dimercheck_dir, primer_name + "_dimer_out")
            names = [fwd_name, rev_name, primer_name]
            files = [fwd_file, rev_file, p_file]
            summary_data.append(primer_name)
            for file_name in files:
                with open(file_name, "r") as f:
                    for line in f:
                        line = line.strip().split("\t")
                        val = float(line[2].strip())
                        summary_data.append(val)
                        if file_name == p_file:
                            if val >= float(self.config.mpprimer):
                                name = file_name.split("/")[-1].split(
                                    "_dimer_out")[0]
                                test.append(name)
                        else:
                            if val >= float(self.config.mpprimer) - 1:
                                name = file_name.split("/")[-1].split(
                                    "_dimer_out")[0]
                                test.append(name)

                if self.config.intermediate is False:
                    outfile = file_name
                    name = os.path.basename(outfile).split("_dimer_out")[0]
                    infile = os.path.join(self.dimercheck_dir, name)
                    os.remove(outfile)
                    os.remove(infile)

            if test == names:
                if primer_name not in choice:
                    choice.append(primer_name)

            return summary_data

        def run_primerdimer_check():
            choice = []
            if os.path.isdir(self.dimercheck_dir):
                shutil.rmtree(self.dimercheck_dir)
            G.create_directory(self.dimercheck_dir)
            primer_data = self.get_primerinfo(
                dimercheck, mode="dimercheck")
            for item in primer_data:
                write_dimercheck_input(item, item[0], item[1], item[2])

            for files in os.listdir(self.dimercheck_dir):
                if not fnmatch.fnmatch(files, "*_dimer_out"):
                    input_file = os.path.join(self.dimercheck_dir, files)
                    output_file = os.path.join(
                        self.dimercheck_dir, files + "_dimer_out")
                    dimer_cmd = [
                        "MPprimer_dimer_check.pl", "-f", input_file, "-d", "3",
                        ">", output_file]
                    G.run_shell(
                        " ".join(dimer_cmd), printcmd=False, logcmd=False,
                        log=False)

            dimer_summary = []
            for item in primer_data:
                summary_data = get_dimercheck_output(
                    item, item[0], item[1], item[2], choice)
                dimer_summary.append(summary_data)
            filepath = os.path.join(
                self.dimercheck_dir, "dimercheck_summary.csv")
            if len(dimer_summary) > 0:
                header = ["primer", "fwd-fwd dG", "rev-rev dG", "fwd-rev dG"]
                G.csv_writer(filepath, dimer_summary, header)
            return choice

        choice = run_primerdimer_check()
        info = str(len(choice)) + " primer left after primer-dimer check"
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after primer QC: "
            + str(len(choice)))
        print("\n" + info)
        G.logger("> " + info)
        return choice

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
                    self.results_dir,
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
                duration = time.time() - self.start
                G.logger(
                    "> PrimerQC time: "
                    + str(timedelta(seconds=duration)).split(".")[0])
                return total_results

        error_msg = "No compatible primers found"
        duration = time.time() - self.start
        G.logger(
            "> PrimerQC time: "
            + str(timedelta(seconds=duration)).split(".")[0])
        print(error_msg)
        G.logger("> " + error_msg)
        errors.append([self.target, error_msg])
        return total_results


class Summary:
    def __init__(self, configuration, total_results):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.primer_dir = os.path.join(self.results_dir, "primer")
        self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
        self.mfold_dir = os.path.join(self.primer_dir, "mfold")
        self.summ_dir = os.path.join(self.config.path, "Summary", self.target)
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.aka = H.abbrev(self.target)
        self.g_info_dict = {}
        self.total_results = total_results

    def collect_qc_infos(self, qc_gene):

        def genome_info_init(row):
            accessiondata = {
                "name": "", "assemblystatus": "",
                "strain": "",
                "rRNA": {
                    "status": "", "hit": "",
                    "GI": "", "DB_id": ""},
                "tuf": {
                    "status": "", "hit": "",
                    "GI": "", "DB_id": ""},
                "recA": {
                    "status": "", "hit": "",
                    "GI": "", "DB_id": ""},
                "dnaK": {
                    "status": "", "hit": "",
                    "GI": "", "DB_id": ""},
                "pheS": {
                    "status": "", "hit": "",
                    "GI": "", "DB_id": ""}}

            accession = "_".join(row[0].split("_")[0:-1])
            if "GCF" in accession or "GCA" in accession:
                accession = ".".join(accession.split("v"))

            if accession not in self.g_info_dict.keys():
                add_accession = {
                    accession: accessiondata}
                self.g_info_dict.update(add_accession)

            self.g_info_dict[accession][qc_gene]["GI"] = row[1]
            self.g_info_dict[accession][qc_gene]["DB_id"] = row[2]
            self.g_info_dict[accession][qc_gene]["hit"] = row[3]
            self.g_info_dict[accession][qc_gene]["status"] = row[5]

        if not qc_gene == "None":
            G.logger(
                "Run: collect_qc_infos(" + self.target + " " + qc_gene + ")")
            qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
            qc_report = os.path.join(qc_dir, qc_gene + "_QC_report.csv")
            try:
                with open(qc_report, "r") as r:
                    reader = csv.reader(r)
                    next(reader, None)
                    for row in reader:
                        genome_info_init(row)
            except FileNotFoundError:
                pass
        else:
            qc_gene = "rRNA"
            filelist = os.listdir(self.fna_dir)
            for filename in filelist:
                row = [filename, "", "", "", "", ""]
                genome_info_init(row)

    def get_genome_infos(self):
        G.logger("Run: get_genome_infos(" + self.target + ")")
        genome_data = []
        genomedata = os.path.join(self.config_dir, "genomicdata.json")
        if os.path.isfile(genomedata):
            with open(genomedata) as q:
                for line in q:
                    records = json.loads(line)

            for result in records['DocumentSummarySet']['DocumentSummary']:
                accession = result["AssemblyAccession"]
                name = result["AssemblyName"]
                status = result["AssemblyStatus"]
                strain = "unknown"
                infraspecieslist = result['Biosource']['InfraspeciesList']
                for i, item in enumerate(infraspecieslist):
                    if len(item) > 0:
                        strain = infraspecieslist[i]['Sub_value']

                data = [accession, name, status, strain]
                genome_data.append(data)

            for item in genome_data:
                accession = item[0]
                if accession in self.g_info_dict.keys():
                    ncbi_info = {
                        "name": item[1], "strain": item[3],
                        "assemblystatus": item[2]}
                    self.g_info_dict[item[0]].update(ncbi_info)

    def write_genome_info(self):
        G.logger("Run: write_genome_info(" + self.target + ")")
        file_name = self.aka + "_qc_sequences.csv"
        filepath = os.path.join(self.summ_dir, file_name)
        header = [
                "Assembly accession", "Assembly name", "Assembly status",
                "Strain", "rRNA", "rRNA Blast", "tuf", "tuf Blast", "recA",
                "recA Blast", "dnaK", "dnaK Blast", "pheS", "pheS Blast"]
        genomeinfo = []
        for key in self.g_info_dict:
            k = self.g_info_dict[key]
            infos = [
                key, k["name"], k['assemblystatus'], k["strain"],
                k["rRNA"]["status"], k["rRNA"]["hit"],
                k["tuf"]["status"], k["tuf"]["hit"],
                k["recA"]["status"], k["recA"]["hit"],
                k["dnaK"]["status"], k["dnaK"]["hit"],
                k["pheS"]["status"], k["pheS"]["hit"]]
            genomeinfo.append(infos)
        G.csv_writer(filepath, genomeinfo, header)

        file_name = self.aka + "_qc_sequences_details.csv"
        filepath = os.path.join(self.summ_dir, file_name)
        header = [
                "Assembly accession", "Assembly name", "Assembly status",
                "Strain", "rRNA", "rRNA Blast", "Hit GI", "Hit DB_id",
                "tuf", "tuf Blast", "Hit GI", "Hit DB_id",
                "recA", "recA Blast", "Hit GI", "Hit DB_id",
                "dnaK", "dnaK Blast", "Hit GI", "Hit DB_id",
                "pheS", "pheS Blast", "Hit GI", "Hit DB_id"]
        detailinfo = []
        for key in self.g_info_dict:
            k = self.g_info_dict[key]
            infos = [
                key, k["name"], k['assemblystatus'], k["strain"],
                k["rRNA"]["status"], k["rRNA"]["hit"],
                k["rRNA"]["GI"], k["rRNA"]["DB_id"],
                k["tuf"]["status"], k["tuf"]["hit"],
                k["tuf"]["GI"], k["tuf"]["DB_id"],
                k["recA"]["status"], k["recA"]["hit"],
                k["recA"]["GI"], k["recA"]["DB_id"],
                k["dnaK"]["status"], k["dnaK"]["hit"],
                k["dnaK"]["GI"], k["dnaK"]["DB_id"],
                k["pheS"]["status"], k["pheS"]["hit"],
                k["pheS"]["GI"], k["pheS"]["DB_id"]]
            detailinfo.append(infos)
        G.csv_writer(filepath, detailinfo, header)

    def write_results(self, total_results):
        today = time.strftime("%Y_%m_%d", time.localtime())
        if not total_results == []:
            G.logger("Run: write_results(" + self.target + ")")
            wrote = []
            outputlist = []
            header = [
                "Primer name", "PPC", "Primer penalty", "Annotation",
                "Primer fwd seq", "Primer fwd TM", "Primer fwd penalty",
                "Primer rev seq", "Primer rev TM", "Primer rev penalty",
                "Probe seq", "Probe TM", "Probe penalty",
                "Amplicon size", "Amplicon TM", "Amplicon sequence",
                "Template sequence"]
            filepath = os.path.join(self.summ_dir, self.aka + "_primer.csv")
            if os.path.isfile(filepath):
                filepath = os.path.join(
                    self.summ_dir, self.aka + "_primer" + today + ".csv")
            for result in total_results:
                if result not in wrote:
                    wrote.append(result)
                    outputlist.append(result)
            G.csv_writer(filepath, outputlist, header)

    def copy_pangenomeinfos(self):
        for filename in os.listdir(self.pangenome_dir):
            filepath = os.path.join(self.pangenome_dir, filename)
            if filename.endswith("core_gene_alignment.aln"):
                shutil.copy(filepath, self.summ_dir)
            if filename.endswith("_tree.nwk"):
                shutil.copy(filepath, self.summ_dir)
            if filename.endswith("Rplots.pdf"):
                shutil.copy(filepath, self.summ_dir)

    def copy_config(self):
        today = time.strftime("%Y_%m_%d", time.localtime())
        for filename in os.listdir(self.config_dir):
            if filename.startswith("config.json"):
                filepath = os.path.join(self.config_dir, filename)
                if os.path.isfile(os.path.join(self.config_dir, filename)):
                    targetpath = os.path.join(
                        self.summ_dir, filename.split(".json")[0] + "_"
                        + today + ".json")
                    shutil.copy(filepath, targetpath)
                else:
                    shutil.copy(filepath, self.summ_dir)

    def copy_pipelinestats(self):
        for files in os.listdir(self.target_dir):
            if files.startswith("pipeline_stats_") and files.endswith(".txt"):
                filename = files
                filepath = os.path.join(self.target_dir, filename)
                abbr = H.abbrev(self.target)
                targetpath = os.path.join(
                    self.summ_dir, abbr + "_" + filename)
                shutil.copy(filepath, targetpath)

    def copy_mostcommon_hits(self):
        filepath = os.path.join(self.blast_dir, "mostcommonhits.csv")
        targetpath = os.path.join(self.summ_dir, "mostcommonhits.csv")
        if os.path.isfile(filepath):
            try:
                shutil.copy(filepath, targetpath)
            except OSError:
                pass

    def last_summary_nolist(self):
        specieslist = []
        blast_path = os.path.join(
                self.blast_dir, "nontargethits.json")
        primerblast_path = os.path.join(
                self.primerblast_dir, "nontargethits.json")
        filelist = [blast_path, primerblast_path]
        for file_path in filelist:
            if os.path.isfile(file_path):
                with open(file_path, 'r') as f:
                    for line in f:
                        blast_dict = json.loads(line)

                for key in blast_dict.keys():
                    for speciesname in blast_dict[key]:
                        if speciesname not in specieslist:
                            specieslist.append(speciesname)

        listpath = os.path.join(
            self.summ_dir, "potential_specieslist.txt")
        with open(listpath, "w") as sp_l:
            for speciesname in specieslist:
                genus = speciesname.split(" ")[0]
                if not (
                    "." in genus or "-" in speciesname or
                    len(re.findall(r'[A-Z]', genus)) == len(genus)
                ):
                    if len(re.findall(r'[0-9]', speciesname)) == 0:
                        sp_l.write(speciesname + "\n")


    def run_summary(self, mode="normal"):
        G.logger("Run: run_summary(" + self.target + ")")
        G.create_directory(self.summ_dir)
        self.write_results(self.total_results)
        if not self.config.qc_gene is None:
            for qc_gene in self.config.qc_gene:
                self.collect_qc_infos(qc_gene)
        else:
            self.collect_qc_infos("None")
        self.copy_mostcommon_hits()
        self.get_genome_infos()
        self.write_genome_info()
        if mode == "last":
            if self.config.nolist:
                self.last_summary_nolist()
            # copy coregenealignment, trees to summary_dir
            self.copy_pangenomeinfos()
            self.copy_config()
            PipelineStatsCollector(
                self.target_dir).write_stat("End: " + str(time.ctime()))
            self.copy_pipelinestats()
            msg = [
                "SpeciesPrimer run finished for", self.target,
                "\n", "End:", time.ctime(), "\n", "See results in",
                self.summ_dir]
            print(" ".join(msg))
            G.logger(" ".join(msg))


class PipelineStatsCollector():
    def __init__(self, target_dir):
        self.target_dir = target_dir
        today = time.strftime("%Y_%m_%d", time.localtime())
        self.statfile = os.path.join(
            target_dir, "pipeline_stats_" + today + ".txt")

    def write_stat(self, info):
        if os.path.isfile(self.statfile):
            with open(self.statfile, "a") as f:
                f.write(info + "\n")
        else:
            with open(self.statfile, "w") as f:
                f.write(info + "\n")


def commandline():
    # command line interface
    parser = argparse.ArgumentParser(
        prog="SpeciesPrimer",
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="This Pipeline is designed to find prokaryotic species "
        "specific primer pairs for PCR",
        allow_abbrev=False)
    # define targets
    parser.add_argument(
        "-t", "--target", nargs="*", type=str, help="Bacterial species in "
        "format: 'Genus_species' e.g. 'Lactobacillus_casei'"
        " use spaces to separate different species, Virus species in format "
        "Species_genus e.g. Wuhan_coronavirus")
    # path
    parser.add_argument(
        "-p", "--path", type=str, help="Absolute path of working directory, "
        "Default = current working directory", default=os.getcwd())
    # Download
    parser.add_argument(
        "--skip_download", action="store_true",
        help="Skip download of genomes")
    # Quality Control choose a qc gene
    parser.add_argument(
       "--qc_gene", nargs="*", type=str,
       choices=["rRNA", "tuf", "recA", "dnaK", "pheS"], default=["rRNA"],
       help="Use Gene to Blast during initial quality control step "
       "(default = '16S rRNA')")
    # Define an exception
    parser.add_argument(
       "-x", "--exception", nargs="*", type=str, default=[],
       help="Define a bacterial species for which assay target sequence "
       "similarity is tolerated; format: 'Genus_species'")
    # skip_tree option
    parser.add_argument(
       "--skip_tree", action="store_true",
       help="Faster Pangenome analysis, no core gene alignment and tree for "
       "troubleshooting is generated.")
    # primer3 amplicon size
    parser.add_argument(
       "--minsize", "-min", type=int, default=70,
       help="Minimal Amplicon size (Default = 70)")
    parser.add_argument(
       "--maxsize", "-max", type=int, default=200,
       help="Maximal Amplicon size (Default = 200)")
    # mfold-threshold
    parser.add_argument(
        "--mfold", type=float, default=-3.0, help="Delta G threshold for "
        " secondary structures in PCR products at 60 degree Celsius"
        "calculated by mfold (default = -3.0)")
    # mpprimer-threshold
    parser.add_argument(
        "--mpprimer", type=float, default=-3.5, help="Delta G threshold for "
        "3'-end primer-primer binding calculated by MPprimer (default = -3.5)")
    parser.add_argument(
        "--mfethreshold", type=int, default=90, help="Threshold for "
        " MFEprimer PPC for nontarget sequences (default = 90)")
    parser.add_argument(
        "--assemblylevel", "-l", nargs="*", type=str, default=["all"],
        choices=[
            "complete", "chromosome", "scaffold", "contig", "all", "offline"],
        help="Limit downloads of Genomes to assembly status")
    # intermediate
    parser.add_argument(
       "--intermediate", action="store_true", default=False,
       help="Do not delete the intermediate files")
    parser.add_argument(
        "--blastseqs", type=int, choices=[100, 500, 1000, 2000, 5000],
        help="Set the number of sequences per BLAST search. "
        "Decrease the number of sequences if BLAST slows down due to low "
        "memory, default=1000", default=1000)
    parser.add_argument(
        "--probe", action="store_true", help="Primer3 designs also an "
        "internal oligo [Experimental!]")
    parser.add_argument(
        "--nolist", action="store_true", help="Species list is not used"
        " and only sequences without blast hits are used for primer design "
        "[Experimental, not recommended for nt DB!]")
    parser.add_argument(
        "--offline", action="store_true", help="Work offline no data from"
        " NCBI is collected, use your own Genomic sequences")
    parser.add_argument(
        "--ignore_qc", action="store_true", help="Genomes which do not"
        " pass quality control are included in the analysis")
    parser.add_argument(
        "--customdb", type=str, default=None,
        help="Absolute filepath of a custom database for blastn")
    parser.add_argument(
        "-e", "--email", type=str, default=None,
        help="A valid email address to make use of NCBI's E-utilities"
        " default=None (user input is required later for online functions)")
    parser.add_argument(
        "--configfile", type=str, default=None,
        help="Provide the path to a JSON format inputfile, "
        "with keys and a list of "
        "settings or a path to a custom settings file. "
        "Key and example file name: "
        '["genus_abbrev", "genus_abbrev.csv"], '
        '["species_list","species_list.txt"], '
        '["p3settings", "p3parameters"], '
        '["excludedgis", "no_blast.gi"]'
        "The current settings files will be overwritten")
    parser.add_argument(
        "--runmode", "-m", type=str, default=["species"],
        choices=["species", "strain"], help="Singleton is a new feature "
        "under development")
    parser.add_argument(
        "--strains", nargs="*", type=str, help="Start of filename of annotated "
        "fna file, GCF_XYZXYZXYZv1, will only search for singletons for this "
        "genome", default = [])

    parser.add_argument(
        "-g", "--genbank", action="store_true",
        help="Download genome assemblies from Genbank"
            )
    parser.add_argument(
        "--evalue", type=float, default=500.0,
        help="E-value threshold for BLAST search, "
        "all results with a lower value pass")
    parser.add_argument(
        "--nuc_identity", type=int, default=0,
        help="Nucleotide identity % threshold for BLAST search, "
        "all results with a lower value pass")
    parser.add_argument(
        "-v", "--virus", action="store_true",
        help="Design primers for viruses")
    # Version
    parser.add_argument(
        "-V", "--version", action="version", version="%(prog)s 2.2.0-dev")

    return parser


def citation():
    citation = """
    Please cite SpeciesPrimer if you use any of the results it produces:
    Dreier M, Berthoud H, Shani N, Wechsler D, Junier P. 2020.
    SpeciesPrimer: a bioinformatics pipeline dedicated to the design
    of qPCR primers for the quantification of bacterial species.
    PeerJ 8:e8544 https://doi.org/10.7717/peerj.8544
    """
    print(citation)
    return citation


def auto_run():
    tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')
    with open(tmp_db_path, 'r') as f:
        for line in f:
            tmp_db = json.loads(line)
    if tmp_db["new_run"]['modus'] == 'continue':
        mode = 'continue'
        data = tmp_db["new_run"]['path'], tmp_db["new_run"]['targets']
    else:
        mode = "new"
        data = tmp_db["new_run"]["targets"]

    import batchassist
    config_dict = batchassist.Input().gui_runner(mode, data)
    targets = []
    for key in config_dict:
        targets.append(key)
    conf_from_file = Config(mode="auto", config_dict=config_dict)
    use_configfile = True
    return targets, conf_from_file, use_configfile


def get_configuration_from_file(target, conf_from_file):
    (
        minsize, maxsize, mpprimer, exception, target, path,
        intermediate, qc_gene, mfold, skip_download,
        assemblylevel, skip_tree, nolist,
        offline, ignore_qc, mfethreshold, customdb,
        blastseqs, probe, virus, genbank,
        evalue, nuc_identity, runmode, strains
    ) = conf_from_file.get_config(target)
    if nolist:
        nontargetlist = []
    else:
        nontargetlist = H.create_non_target_list(target)

    config = CLIconf(
        minsize, maxsize, mpprimer, exception, target, path,
        intermediate, qc_gene, mfold, skip_download,
        assemblylevel, nontargetlist, skip_tree, nolist,
        offline, ignore_qc, mfethreshold, customdb,
        blastseqs, probe, virus, genbank,
        evalue, nuc_identity, runmode, strains)

    return config


def run_pipeline_for_target(target, config):
    print("\nStart searching primer for " + target)
    G.logger("> Start searching primer for " + target)
    target_dir = os.path.join(config.path, target)
    PipelineStatsCollector(target_dir).write_stat(
        target + " pipeline statistics:")
    PipelineStatsCollector(target_dir).write_stat(
        "Start: " + str(time.ctime()))
    newconfig = DataCollection(config).collect()
    if newconfig != 0:
        config = newconfig
    if config.virus is True:
        config.qc_gene = []
    qc_count = []
    for qc_gene in config.qc_gene:
        qc = QualityControl(config).quality_control(qc_gene)
        qc_count.append(qc)
    if not sum(qc_count) == 0:
        # writes QC summary in summary directory, run is finished
        total_results = []
        Summary(config, total_results).run_summary()
    else:
        # writes QC summary in summary directory, run continues
        try:
            total_results = []
            Summary(config, total_results).run_summary()
        except FileNotFoundError:
            pass

        PangenomeAnalysis(config).run_pangenome_analysis()

        if "strain" in config.runmode:
            import strainprimer
            msg = "Start searching for strain specific primers"
            print(msg)
            G.logger(msg)
            strainprimer.main(config)

        if "species" in config.runmode:
            msg = "Start searching for species specific primers"
            print(msg)
            G.logger(msg)
            CoreGenes(config).run_CoreGenes()
            conserved_seq_dict = CoreGeneSequences(
                    config).run_coregeneanalysis()
            if not conserved_seq_dict == 1:
                conserved = BlastParser(
                        config).run_blastparser(conserved_seq_dict)
                if conserved == 0:
                    primer_dict = PrimerDesign(config).run_primerdesign()
                    total_results = PrimerQualityControl(
                        config, primer_dict).run_primer_qc()
                    Summary(config, total_results).run_summary(mode="last")
                else:
                    Summary(config, total_results).run_summary(mode="last")
            else:
                Summary(config, total_results).run_summary(mode="last")


def get_configuration_from_args(target, args):
    if args.nolist:
        nontargetlist = []
    else:
        nontargetlist = H.create_non_target_list(target)

    config = CLIconf(
        args.minsize, args.maxsize, args.mpprimer, args.exception,
        target, args.path, args.intermediate,
        args.qc_gene, args.mfold, args.skip_download,
        args.assemblylevel, nontargetlist,
        args.skip_tree, args.nolist, args.offline,
        args.ignore_qc, args.mfethreshold, args.customdb,
        args.blastseqs, args.probe, args.virus, args.genbank,
        args.evalue, args.nuc_identity, args.runmode, args.strains)

    if args.configfile:
        exitstat = H.advanced_pipe_config(args.configfile)
        if exitstat:
            sys.exit()

    return config

def exitatsigterm(signalNumber, frame):
    raise SystemExit('GUI stop')


def main(mode=None):
    today = time.strftime("%Y_%m_%d", time.localtime())
    use_configfile = False

    if mode == "auto":
        signal.signal(signal.SIGTERM, exitatsigterm)
        os.chdir(os.path.join("/", "primerdesign"))
        logfile = os.path.join(os.getcwd(), "speciesprimer_" + today + ".log")
        logging.basicConfig(
            filename=logfile, level=logging.DEBUG, format="%(message)s")
        targets, conf_from_file, use_configfile = auto_run()

    else:
        parser = commandline()
        # required for testing (ignore pytest arguments)
        args, unknown = parser.parse_known_args()
        logfile = os.path.join(os.getcwd(), "speciesprimer_" + today + ".log")
        logging.basicConfig(
            filename=logfile, level=logging.DEBUG, format="%(message)s")
        if args.target is None:
            conf_from_file = Config()
            targets = conf_from_file.get_targets()
            use_configfile = True
        else:
            targets = args.target

        if not os.path.isabs(args.path):
            args.path = os.path.join(os.getcwd(), args.path)

        if args.configfile:
            if not os.path.isabs(args.configfile):
                args.configfile = os.path.join(os.getcwd(), args.configfile)

        if args.email:
            H.get_email_for_Entrez(args.email)

        if unknown:
            msg = "\t!!! Warning !!! the following arguments are not known:"
            print("\n" + msg)
            print("\n".join(unknown))
            G.logger(msg)
            G.logger(unknown)

    G.logger(citation())

    for target in targets:
        target = target.capitalize()
        if use_configfile:
            config = get_configuration_from_file(target, conf_from_file)
        else:
            config = get_configuration_from_args(target, args)

        today = time.strftime("%Y_%m_%d", time.localtime())
        G.logger("> Start log: " + target + " " + today)
        H.BLASTDB_check(config)
        G.logger(config.__dict__)

        try:
            run_pipeline_for_target(target, config)
        except Exception as exc:
            msg = [
                "fatal error while working on", target,
                "check logfile", logfile]
            target_dir = os.path.join(config.path, config.target)
            PipelineStatsCollector(target_dir).write_stat(
                "Error: " + str(time.ctime()))
            print(" ".join(msg))
            print(exc)
            import traceback
            traceback.print_exc()
            G.logger(" ".join(msg))
            errors.append([target, " ".join(msg)])
            logging.error(
                "fatal error while working on " + target, exc_info=True)

        except (KeyboardInterrupt, SystemExit):
            logging.error(
                "SpeciesPrimer was stopped while working on " + target,
                exc_info=True)
            raise

    if len(errors) > 0:
        print("Error report: ")
        G.logger("> Error report: ")
        for index, error in enumerate(errors):
            error_nr = "Error " + str(index + 1) + ":"
            print("for target " + error[0])
            print(error_nr)
            print(error[1])
            G.logger("> for target " + error[0])
            G.logger("> " + str(error_nr))
            G.logger("> " + str(error[1]))


if __name__ == "__main__":
    main()
