#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import time
import logging
import csv
import fnmatch
import re
import shutil
import multiprocessing
from multiprocessing import Process
import wget
import json
import tempfile
import sqlite3
import urllib
import itertools
from itertools import islice
from datetime import timedelta
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastnCommandline
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H

# paths
pipe_bin = os.path.abspath(__file__)
pipe_dir = pipe_bin.split("bin")[0]
dict_path = os.path.join(pipe_dir, "dictionaries")
errors = []

# general info
systemdirs = [
    "genomic_fna", "config", "ffn_files", "gff_files", "Pangenome",
    "rRNA_control", "recA_control", "tuf_control",
    "dnaK_control", "pheS_control"]


Entrez.tool = "SpeciesPrimer pipeline"


class Config:

    def __init__(self):
        self.config_dict = self.generate_config_files()

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
        verbosity = self.config_dict[target]["verbosity"]
        qc_gene = self.config_dict[target]["qc_gene"]
        mfold = self.config_dict[target]["mfold"]
        skip_download = self.config_dict[target]["skip_download"]
        assemblylevel = self.config_dict[target]["assemblylevel"]
        skip_tree = self.config_dict[target]["skip_tree"]
        nolist = self.config_dict[target]["nolist"]
        offline = self.config_dict[target]["offline"]
        ignore_qc = self.config_dict[target]["ignore_qc"]
        mfethreshold = self.config_dict[target]["mfethreshold"]
        remoteblast = self.config_dict[target]["remoteblast"]
        blastseqs = self.config_dict[target]["blastseqs"]
        probe = self.config_dict[target]["probe"]

        return (
            minsize, maxsize, mpprimer, exception, target, path,
            verbosity, qc_gene, mfold, skip_download,
            assemblylevel, skip_tree, nolist, offline, ignore_qc, mfethreshold,
            remoteblast, blastseqs, probe)


class CLIconf:
    def __init__(
            self, minsize, maxsize, mpprimer, exception, target, path,
            verbosity, qc_gene, mfold,
            skip_download, assemblylevel,
            nontargetlist, skip_tree, nolist, offline, ignore_qc, mfethreshold,
            remoteblast, blastseqs, probe):
        self.minsize = minsize
        self.maxsize = maxsize
        self.mpprimer = mpprimer
        self.exception = exception
        self.target = target
        self.path = path
        self.verbosity = verbosity
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
        self.remoteblast = remoteblast
        self.blastseqs = blastseqs
        self.probe = probe
        self.save_config()

    def save_config(self):
        config_dict = {}
        config_dict.update({"minsize": self.minsize})
        config_dict.update({"maxsize": self.maxsize})
        config_dict.update({"mpprimer": self.mpprimer})
        config_dict.update({"exception": self.exception})
        config_dict.update({"target": self.target})
        config_dict.update({"path": self.path})
        config_dict.update({"verbosity": self.verbosity})
        config_dict.update({"qc_gene": self.qc_gene})
        config_dict.update({"mfold": self.mfold})
        config_dict.update({"skip_download": self.skip_download})
        config_dict.update({"assemblylevel": self.assemblylevel})
        config_dict.update({"skip_tree": self.skip_tree})
        config_dict.update({"nolist": self.nolist})
        config_dict.update({"offline": self.offline})
        config_dict.update({"ignore_qc": self.ignore_qc})
        config_dict.update({"mfethreshold": self.mfethreshold})
        config_dict.update({"remoteblast": self.remoteblast})
        config_dict.update({"blastseqs": self.blastseqs})
        config_dict.update({"probe": self.probe})

        dir_path = os.path.join(self.path, self.target)
        config_path = os.path.join(self.path, self.target, "config")
        file_path = os.path.join(config_path, "config.json")
        G.create_directory(dir_path)
        G.create_directory(config_path)
        with open(file_path, "w") as f:
            f.write(json.dumps(config_dict))


class DataCollection():

    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.genomic_dir = os.path.join(self.target_dir, "genomic_fna")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.ffn_dir = os.path.join(self.target_dir, "ffn_files")
        self.fna_dir = os.path.join(self.target_dir, "fna_files")
        self.ex_dir = os.path.join(self.config.path, "excludedassemblies", self.target)

    def get_email_for_Entrez(self):
        mailfile = os.path.join(pipe_dir, "user_email.txt")
        if os.path.isfile(mailfile):
            with open(mailfile) as f:
                for line in f:
                    if "@" and "." in line:
                        email = line.strip()
        else:
            email = input(
                "To make use of NCBI's E-utilities, NCBI requires you to specify "
                "your email address. \n"
                "Your email address will be stored locally in this directory: \n"
                + mailfile + "\n")

            if "@" and "." in email:
                with open(mailfile, "w") as f:
                    f.write(email.strip())

        return email

    def get_taxid(self, target):
        Entrez.email = self.get_email_for_Entrez()

        taxids = H.taxiddict(dict_path)
        taxid = H.check_input(target, taxids, Entrez.email)
        return taxid, Entrez.email

    def prepare_dirs(self):
        G.create_directory(self.target_dir)
        G.create_directory(self.config_dir)
        G.create_directory(self.genomic_dir)

    def create_GI_list(self):
        G.logger("Run: create_GI_list(" + self.target + ")")
        removed_gis = []

        def create_bin(input_file, output_file):
            cmd_GI_bin = [
                'blastdb_aliastool', '-gi_file_in', input_file,
                '-gi_file_out', output_file]
            if os.path.isfile(output_file):
                info = "GI bin already exists"
                G.logger(info)
                if self.config.verbosity == 2:
                    print(info)
            else:
                if self.config.verbosity == 2:
                    G.run_subprocess(cmd_GI_bin, True, True, False, True)
                else:
                    G.run_subprocess(cmd_GI_bin, True, True, True, False)
                info = "Created GI bin"
                G.logger(info)
                if self.config.verbosity == 2:
                    print(info)

        def removed_geneidentifier():
            for files in os.listdir(os.path.join(pipe_dir, "NO_Blast")):
                if files.endswith(".gi"):
                    with open(os.path.join(pipe_dir, "NO_Blast", files), "r") as f:
                        for line in f:
                            if "#" in line:
                                pass
                            else:
                                gi = line.strip()
                                if gi not in removed_gis:
                                    removed_gis.append(gi)
            if len(removed_gis) > 0:
                with open(os.path.join(self.config_dir, "NO_Blast.gi"), "w") as f:
                    for gi in removed_gis:
                        f.write(gi+"\n")

        if os.path.isdir(os.path.join(pipe_dir, "NO_Blast")):
            removed_geneidentifier()
            gi_file = os.path.join(self.config_dir, "NO_Blast.gi")
            gi_bin = os.path.join(self.config_dir,"NO_Blast.bin")
            if os.path.isfile(gi_file):
                if os.stat(gi_file).st_size > 0:
                    create_bin(gi_file, gi_bin)
        else:
            pass

    def get_ncbi_links(self, taxid, email):

        def collect_genomedata(taxid, email):
            genomedata = []
            # 02/12/18 overwrite file because assembly level can change
            Entrez.email = email
            assembly_search = Entrez.esearch(
                db="assembly",
                term="txid" + str(taxid) + "[Orgn]",
                retmax=2000)
            assembly_record = Entrez.read(assembly_search)
            uidlist = assembly_record["IdList"]
            assembly_efetch = Entrez.efetch(
                db="assembly",
                id=uidlist,
                rettype="docsum",
                retmode="xml")
            assembly_records = Entrez.read(assembly_efetch)

            with open("genomicdata.json", "w") as f:
                f.write(json.dumps(assembly_records))

            with open("genomicdata.json", "r") as f:
                for line in f:
                    genomic_dict = json.loads(line)

            for assembly in genomic_dict['DocumentSummarySet']['DocumentSummary']:
                accession = assembly["AssemblyAccession"]
                name = assembly["AssemblyName"]
                status = assembly["AssemblyStatus"]
                ftp_path = assembly["FtpPath_RefSeq"]
                if ftp_path is "":
                    ftp_link = "None"
                else:
                    ftp_link = (
                        ftp_path + "/" + ftp_path.split("/")[-1] + "_genomic.fna.gz")
                data = [accession, name, status, ftp_link]
                genomedata.append(data)

            return genomedata

        def get_links(linkdata):
            statusdict = {
                "complete": "Complete Genome", "chromosome": "Chromosome",
                "scaffold": "Scaffold", "contig": "Contig"}
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
        genomedata = collect_genomedata(taxid, email)
        link_list = get_links(genomedata)
        info = str(len(link_list)) + " genome assemblies are available for download"
        print(info)
        PipelineStatsCollector(self.target_dir).write_stat(
            "genome assemblies from NCBI: " + str(len(link_list)) )
        G.logger(info)
        write_links(link_list)
        os.chdir(self.target_dir)

    def ncbi_download(self):
        G.logger("Run: ncbi_download(" + self.target + ")")
        G.create_directory(self.gff_dir)
        G.create_directory(self.ffn_dir)
        G.create_directory(self.fna_dir)
        excluded = []
        # a list of excluded Genomes to keep information
        # after directories are (manually) removed (to save diskspace)
        if os.path.isdir(self.ex_dir):
            try:
                filepath = os.path.join(self.ex_dir, "excluded_list.txt")
                with open(filepath, "r") as f:
                    for line in f:
                        line = line.strip()
                        if not line in excluded:
                            excluded.append(line)
            except FileNotFoundError:
                pass

        def check_files(input_line):
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
            if (gff and ffn) == True:
                return True
            elif genomic == True:
                return "Extracted"
            else:
                return False

        os.chdir(self.genomic_dir)
        with open(os.path.join(self.config_dir, "genomic_links.txt")) as r:
            for line in r:
                zip_file = line.split("/")[-1].strip()
                target_path = os.path.join(self.genomic_dir, zip_file)
                ftp_path = line.strip()
                file_status = check_files(line)
                if file_status == True:
                    info = "all required files found for " + zip_file
                    G.logger(info)
                    if self.config.verbosity == 2:
                        print(info)
                elif file_status == "Extracted":
                    info = "extracted fna file already exists for " + zip_file
                    G.logger(info)
                    if self.config.verbosity == 2:
                        print(info)
                else:
                    if os.path.isfile(target_path):
                        info = "File already downloaded " + zip_file
                        G.logger(info)
                        if self.config.verbosity == 2:
                            print(info)
                    else:
                        if [
                            ex for ex in excluded
                            if zip_file.startswith(".".join(ex.split("v")))]:
                                print(zip_file + " already in excludedassemblies")
                        else:
                            try:
                                print("\n\nDownload..." + zip_file + "\n")
                                wget.download(ftp_path)
                                print("\n")
                                info = "Downloaded " + zip_file
                                G.logger(info)
                                if self.config.verbosity == 2:
                                    print(info)
                            except urllib.error.HTTPError:
                                try:
                                    print("\ntry again\n")
                                    print("\n\nDownload..." + zip_file + "\n")
                                    wget.download(ftp_path)
                                    print("\n")
                                    info = "Downloaded " + zip_file
                                    G.logger(info)
                                    if self.config.verbosity == 2:
                                        print(info)
                                except urllib.error.HTTPError:
                                    try:
                                        print("\ntry a last time\n")
                                        print("\n\nDownload..." + zip_file + "\n")
                                        wget.download(ftp_path)
                                        print("\n")
                                        info = "Downloaded " + zip_file
                                        G.logger(info)
                                        if self.config.verbosity == 2:
                                            print(info)
                                    except urllib.error.HTTPError:
                                        raise


        for files in os.listdir(self.genomic_dir):
            if files.endswith(".gz"):
                G.run_subprocess(["gunzip", files], False, True, False, False)
        os.chdir(self.target_dir)

    def copy_genome_files(self):
        G.logger("Run: copy_genome_files(" + self.target + ")")
        for root, dirs, files in os.walk(self.target_dir):
            for file_name in files:
                if file_name.endswith(".ffn"):
                    ffn_file = os.path.join(self.ffn_dir, file_name)
                    if not file_name in os.listdir(self.ffn_dir):
                        shutil.copy(os.path.join(root, file_name), ffn_file)
                if file_name.endswith(".gff"):
                    gff_file = os.path.join(self.gff_dir, file_name)
                    if not file_name in os.listdir(self.gff_dir):
                        shutil.copy(os.path.join(root, file_name), gff_file)

                if file_name.endswith(".fna"):
                    if not file_name.endswith("genomic.fna"):
                        if not "genomic_fna" in root.split("/"):
                            fna_file = os.path.join(self.fna_dir, file_name)
                            if not file_name in os.listdir(self.fna_dir):
                                shutil.copy(os.path.join(root, file_name), fna_file)


    def run_prokka(self):
        annotated = []
        excluded = []
        G.logger("Run: run_prokka(" + self.target + ")")
        date = time.strftime("%Y%m%d")
        genus = self.target.split("_")[0]
        fna_files = []
        dirs = []
        qc_fail_dir = []
        # list of genomic files
        if os.path.isdir(self.genomic_dir):
            for file_name in os.listdir(self.genomic_dir):
                if file_name.endswith(".fna"):
                    fna_files.append(file_name)

        # list of folders of annotated genomes
        for directories in [
            d for d in os.listdir(self.target_dir)
            if os.path.isdir(os.path.join(self.target_dir, d))]:
            name = "_".join(directories.split("_")[:-1])
            dirs.append(name)

        if os.path.isdir(self.fna_dir):
            for files in [
                f for f in os.listdir(self.fna_dir)
                if os.path.isfile(os.path.join(self.fna_dir, f))]:
                name = "_".join(files.split("_")[:-1])
                if not name in dirs:
                    dirs.append(name)
        if os.path.isdir(self.fna_dir):
            for files in [
                f for f in os.listdir(self.gff_dir)
                if os.path.isfile(os.path.join(self.gff_dir, f))]:
                name = "_".join(files.split("_")[:-1])
                if not name in dirs:
                    dirs.append(name)
        if os.path.isdir(self.fna_dir):
            for files in [
                f for f in os.listdir(self.ffn_dir)
                if os.path.isfile(os.path.join(self.ffn_dir, f))]:
                name = "_".join(files.split("_")[:-1])
                if not name in dirs:
                    dirs.append(name)

        # list of failed assemblies
        if os.path.isfile(os.path.join(self.ex_dir, "excluded_list.txt")):
            with open(os.path.join(self.ex_dir, "excluded_list.txt"), "r") as f:
                for line in f:
                    line = line.strip()
                    qc_fail_dir.append(line)

        # write annotation command
        for fna in fna_files:
            if fna.endswith("_genomic.fna"):
                file_name = (
                    fna.split(".")[0] + "v" + fna.split(".")[1].split("_")[0])
            else:
                name = fna.split(".fna")[0]
                file_name = (
                    "_".join("-".join(name.split(".")).split("_")[0:-1]))

            if file_name in dirs:
                if not file_name == '':
                    annotated.append(file_name)

            elif file_name in qc_fail_dir:

                excluded.append(file_name)
            else:
                prokka_cmd = [
                    "prokka", "--kingdom", "Bacteria", "--outdir",
                    file_name + "_" + date, "--genus", genus, "--locustag",
                    file_name, "--prefix", file_name + "_" + date, "genomic_fna/"
                    + fna]
                info = "\n" + file_name + " annotation required"
                G.logger(info)
                if not self.config.verbosity == 0:
                    print(info)
                if self.config.verbosity == 2:
                    G.run_subprocess(prokka_cmd, True, True, False, True)
                else:
                    G.run_subprocess(prokka_cmd, True, True, False, False)

        if len(annotated) > 0:
            info = "\nAlready annotated: "
            G.logger(info.strip())
            G.logger(annotated)
            if not self.config.verbosity == 0:
                print(info)
                print(annotated)

        if len(excluded) > 0:
            info = "\nAlready in excludedassemblies:"
            G.logger(info)
            G.logger(excluded)
            if not self.config.verbosity == 0:
                print(info)
                print(excluded)

    def collect(self):
        G.logger("Run: collect data(" + self.target + ")")
        self.prepare_dirs()
        pan = os.path.join(self.pangenome_dir, "gene_presence_absence.csv")
        if os.path.isfile(pan):
            return 0

        if not self.config.offline:
            if not self.config.skip_download:
                taxid, email = self.get_taxid(self.target)
                self.get_ncbi_links(taxid, email)
                self.ncbi_download()

            else:
                G.create_directory(self.gff_dir)
                G.create_directory(self.ffn_dir)
                G.create_directory(self.fna_dir)
                os.chdir(self.target_dir)
        else:
            G.create_directory(self.gff_dir)
            G.create_directory(self.ffn_dir)
            G.create_directory(self.fna_dir)
            os.chdir(self.target_dir)

        self.create_GI_list()
        self.run_prokka()
        self.copy_genome_files()


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
        self.ex_dir = os.path.join(self.config.path, "excludedassemblies", self.target)
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

    # make a subclass of quality control?
    def get_qc_seqs(self, qc_gene):
        G.logger("Run: get_qc_seqs(" + qc_gene + ")")
        G.logger("Starting QC with " + qc_gene)
        print("Starting QC with " + qc_gene)
        gff = []
        qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")

        def search_qc_gene(file_name):
            with open(os.path.join(self.gff_dir, file_name), "r") as f:
                for line in f:
                    if self.searchdict[qc_gene] in line:
                        gene = line.split("ID=")[1].split(";")[0].split(" ")[0]
                        if gene not in self.qc_gene_search:
                            self.qc_gene_search.append(gene)

        def count_contigs(gff_list, contiglimit):
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
                G.logger(info)

            return gff_list

        def identify_duplicates(gff_list):
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
                        maxi = max(duplicate_test, key=lambda item: int(item.split("v")[1]))

                        if not maxi in keep:
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
                            data = [gff_file.split(".gff")[0], "","","","", "Duplicate"]
                            if data not in self.double:
                                self.double.append(data)

                info = ("skip " + str(len(self.double)) + " duplicate Genome(s) ")
                print(info)
                G.logger(info)

            return gff_list

        # 12.02.2018 change to generate one QC file
        def check_no_sequence(qc_gene, gff):
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
                G.logger(info)

            for item in gff:
                ffn = item.split(".gff")[0] + ".ffn"
                ffn_list.append(ffn)

            return ffn_list

        def run_get_qc_seqs():
            G.create_directory(qc_dir)
            # find annotation of gene in gff files and store file name
            for files in os.listdir(self.gff_dir):
                if files not in gff:
                    gff.append(files)
            info = "found " + str(len(gff)) + " gff files"
            G.logger(info)
            if not self.config.verbosity == 0:
                print(info)

            if len(gff) > 0:
                # look for annotations
                if self.contiglimit > 0:
                    contig_gff_list = count_contigs(gff, self.contiglimit)
                    gff_list = identify_duplicates(contig_gff_list)
                else:
                    gff_list = identify_duplicates(gff)

                for item in gff_list:
                    search_qc_gene(item)

                info = (
                    "found " + str(len(self.qc_gene_search)) + " "
                    + qc_gene + " annotations in gff files")
                G.logger(info)
                if not self.config.verbosity == 0:
                    print(info)

                ffn_check = check_no_sequence(qc_gene, gff_list)

             # search sequences in ffn files
                for files in os.listdir(self.ffn_dir):
                    if files in ffn_check:
                        if files not in self.ffn_list:
                            self.ffn_list.append(files)
                info = (
                        "selected " + str(len(self.ffn_list)) + " "
                        + qc_gene + " sequences from ffn files")
                G.logger(info)
                if not self.config.verbosity == 0:
                    print(info)

            else:
                error_msg = "Error: No .gff files found"
                print(error_msg)
                G.logger(error_msg)
                errors.append([self.target, error_msg])
                return 1
            return 0

        status = run_get_qc_seqs()
        return status

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
                        for record in SeqIO.parse(
                            os.path.join(self.ffn_dir, file_name), "fasta"):
                            if seq_id in record.id:
                                if not len(str(record.seq)) == 0:
                                    if record.id not in recid:
                                        recid.append(record.id)
                                    if str(record.seq) not in recseq:
                                        recseq.append(str(record.seq))

                # get longest sequence and write to file
                if not len(recseq) == 0:
                    q = max(recseq, key=len)
                    z = dict(zip(tuple(recid), tuple(recseq)))
                    a = {v: k for k, v in z.items()}
                    o.write(">" + str(a[q]) + "\n" + str(q) + "\n")
                    qc_seqs.append([">" + str(a[q]) + "\n", str(q) + "\n"])

        return qc_seqs


    def qc_blast_parser(self, qc_gene):
        """ reads the blast results """
        qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
        if self.config.verbosity == 2:
            print("Run: qc_blast_parser(" + qc_gene + ")")
        G.logger("Run: qc_blast_parser(" + qc_gene + ")")
        xmlblastresults = []
        wrote = []
        problems = []
        passed = []

        def find_blastresults():
            if os.path.isdir(qc_dir):
                for filename in os.listdir(qc_dir):
                    if filename.endswith("_results.xml"):
                        if not os.stat(os.path.join(qc_dir, filename)).st_size == 0:
                            xmlblastresults.append(filename)
            if len(xmlblastresults) > 0:
                return True

        def get_blastresults_info(blast_record):
            short = str(blast_record.alignments[0]).split("|")[4].strip(" ").split(" ")
            gi = str(blast_record.alignments[0]).split("|")[1].strip(" ")
            db_id = str(blast_record.alignments[0]).split("|")[3].strip(" ")
            if "subsp." in str(blast_record.alignments[0]).split("|")[4].strip(" "):
                spec = " ".join(short[0:2]) + " " + short[2].split(".")[0] + " " + short[3]
            else:
                spec = short[0]+" "+short[1]

            return spec, gi, db_id

        def parse_blastresults():
            expected = " ".join(self.target.split("_"))
            # collect the blast results
            os.chdir(qc_dir)
            for file_name in xmlblastresults:
                result_handle = open(file_name)
                blast_records = NCBIXML.parse(result_handle)
                blast_records = list(blast_records)
                for blast_record in blast_records:
                    spec, gi, db_id = get_blastresults_info(blast_record)
                    if expected in spec:
                        if blast_record.query not in wrote:
                            wrote.append(blast_record.query)
                            passed.append([blast_record.query, gi, db_id, spec, expected, "passed QC"])
                    else:
                        if blast_record.query not in wrote:
                            wrote.append(blast_record.query)
                            problems.append([blast_record.query, gi, db_id, spec, expected, "failed QC"])

            os.chdir(self.target_dir)

        def write_blastresults():
            # write files
            report = os.path.join(qc_dir, qc_gene + "_QC_report.csv")
            with open(report, "w") as f:
                writer = csv.writer(f)
                writer.writerow(["Query", "GI", "DB ID", "Species", "Target species", "QC status"])
                for item in passed:
                    writer.writerow(item)
                if len(problems) > 0:
                    for item in problems:
                        writer.writerow(item)
                if len(self.no_seq) > 0:
                    for item in self.no_seq:
                        writer.writerow(item)
                if len(self.contig_ex) > 0:
                    for item in self.contig_ex:
                        writer.writerow(item)
                if len(self.double) > 0:
                    for item in self.double:
                        writer.writerow(item)


        def run_blast_parser():
            if find_blastresults():
                parse_blastresults()
                write_blastresults()
                return passed

        passed = run_blast_parser()
        return passed

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

        # remove_qc_failures()
        if len(delete) > 0:
            G.create_directory(self.ex_dir)
            self.delete_failed_assemblies(delete)
            info = "Quality control removed " + str(len(delete)) + " Genome(s)"
            G.logger(info)
            if self.config.verbosity > 0:
                print("\n" + info)


    def delete_failed_assemblies(self, delete):
        filepath = os.path.join(self.ex_dir, "excluded_list.txt")
        with open(filepath, "a+") as f:
            for item in delete:
                f.write(item + "\n")

        def move_dirs(from_dir, to_dir):
            try:
                shutil.move(from_dir, to_dir)
                info = "Moved: "+ from_dir + " to " + to_dir
                G.logger(info)
                if self.config.verbosity == 2:
                    print(info)
            except shutil.Error:
                if os.path.isdir(from_dir):
                    shutil.rmtree(from_dir)
                else:
                    raise

        def move_file(from_file, to_file):
            try:
                shutil.move(from_file, to_file)
                info = "Moved: "+ from_file + " to " + to_file
                G.logger(info)
                if self.config.verbosity == 2:
                    print(info)
            except shutil.Error:
                if os.path.isfile(from_file):
                    os.remove(from_file)
                else:
                    raise

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
                        info = "remove: " + files
                        if self.config.verbosity == 2:
                            print(info)
                        file_format = files.split(".")[1]
                        remove_from_list(file_format, files)
                        os.remove(os.path.join(directory, files))

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
                                    info = "Remove empty directory "+ dir_path
                                    G.logger(info)
                                    if self.config.verbosity == 2:
                                        print(info)
                                    os.rmdir(dir_path)
                                else:
                                    x = '_'.join(file_names[0].split("_")[0:-1])
                                    if x == item:
                                        from_dir = os.path.join(root, directory)
                                        to_dir = os.path.join(self.ex_dir)
                                        move_dirs(from_dir, to_dir)

        def remove_files():
            for item in delete:
                search_str = str(item + "_*")
                if (("GCF" or "GCA") and "v") in item:
                    accession = ".".join(item.split("v"))
                    genomicsearch_str = str(accession + "_*")
                # genomic_fna
                genome_dir = os.path.join(self.target_dir, "genomic_fna")
                if os.path.isdir(genome_dir):
                    for files in os.listdir(genome_dir):
                        if re.search(genomicsearch_str, files):
                            G.create_directory(os.path.join(self.ex_dir, "genomic_fna"))
                            from_file = os.path.join(genome_dir, files)
                            to_file = os.path.join(self.ex_dir, "genomic_fna", files)
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
                                info = "remove" + str(file_name)
                                G.logger(info)
                                if self.config.verbosity == 2:
                                    print(info)
                                os.remove(os.path.join(root, file_name))
        remove_directories()
        remove_files()

    def quality_control(self, qc_gene):
        pan = os.path.join(self.pangenome_dir, "gene_presence_absence.csv")
        if os.path.isfile(pan):
            info = "Found Pangenome directory, skip QC " + qc_gene
            G.logger(info)
            if self.config.verbosity == 2:
                print(info)
            return 0
        else:
            print("\nRun: quality_control(" + qc_gene + ")")
            G.logger("Run: quality_control(" + qc_gene + ")")
            qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
            if self.get_qc_seqs(qc_gene) == 0:
                qc_seqs = self.choose_sequence(qc_gene)
                use_cores, inputseqs = BlastPrep(qc_dir, qc_seqs, qc_gene, self.config.blastseqs).run_blastprep()
                Blast(self.config, qc_dir, "quality_control").run_blast(qc_gene, use_cores)
                passed_list = self.qc_blast_parser(qc_gene)

                if self.config.ignore_qc:
                    info = (
                        "--ignore_qc option: Check quality control"
                        " directories")
                    G.logger(info)
                    if not self.config.verbosity == 0:
                        print("\n" + info)
                    return 0

                elif passed_list:
                    if len(passed_list) < 2:
                        error_msg = "Error: Less than two genomes survived QC"
                        print(error_msg)
                        G.logger(error_msg)
                        errors.append([self.target, error_msg])
                        return 1
                    else:
                        self.remove_qc_failures(qc_gene)
                        return 0
                else:
                    error_msg = "Error: No genomes survived QC"
                    print(error_msg)
                    G.logger(error_msg)
                    errors.append([self.target, error_msg])
                    return 1
            else:
                return 1


class PangenomeAnalysis:
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")


    def run_roary(self):
        G.logger("Run: run_roary(" + self.target + ")")
        os.chdir(self.target_dir)
        if self.config.skip_tree:
            roary_cmd = (
                "roary -f ./Pangenome -s -p " + str(multiprocessing.cpu_count())
                + " -cd 100 ./gff_files/*.gff")
        else:
            roary_cmd = (
                "roary -f ./Pangenome -e -n -s -p " + str(multiprocessing.cpu_count())
                + " -cd 100 ./gff_files/*.gff")

        G.run_shell(roary_cmd, printcmd=True, logcmd=True, log=True, printoption=False)

    def run_fasttree(self):
        G.logger("Run: run_fasttree(" + self.target + ")")
        os.chdir(self.pangenome_dir)
        coregenealn = "core_gene_alignment.aln"
        if os.path.isfile(coregenealn):
            tree = H.abbrev(self.target, dict_path) + "_tree.newick"
            treecmd = "fasttree -nt -gtr -nopr " + coregenealn + " > " + tree
            G.run_shell(treecmd, printcmd=True, logcmd=True, log=True, printoption=False)
        os.chdir(self.target_dir)

    def run_pangenome_analysis(self):
        G.logger("Run: run_pangenome_analysis(" + self.target + ")")
        if not os.path.isdir(self.pangenome_dir):
            if self.config.skip_tree:
                self.run_roary()
            else:
                self.run_roary()
                self.run_fasttree()
        else:
            if os.path.isfile(os.path.join(self.pangenome_dir, "gene_presence_absence.csv")):
                info = (
                    "Pangenome directory already exists\n"
                    "Continue with existing Pangenome data")
                if not self.config.verbosity == 0:
                    print(info)
                G.logger(info)
            else:
                shutil.rmtree(self.pangenome_dir)
                if self.config.skip_tree:
                    self.run_roary()
                else:
                    self.run_roary()
                    self.run_fasttree()

class CoreGenes:
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.gff_dir = os.path.join(self.target_dir, "gff_files")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.all_core_path = os.path.join(self.pangenome_dir, "allcoregenes")
        self.multi_path = os.path.join(self.pangenome_dir, "multiannotated")
        self.iterator = 0

    def copy_DBGenerator(self):
        src = os.path.join(pipe_dir, "ext-scripts", "DBGenerator.py")
        dst = os.path.join(self.pangenome_dir, "DBGenerator.py")
        shutil.copy(src, dst)

    def run_DBGenerator(self):
        G.logger("Run: run_DBGenerator(" + self.target + ")")
        if not os.path.isfile(self.target + ".db"):
            genome_locus = os.path.join(self.pangenome_dir, "genomas_locus.csv")
            locus_seq = os.path.join(self.pangenome_dir, "locus_sequence.csv")
            pangen_locus = os.path.join(self.pangenome_dir, "pangenoma_locus.csv")
            pangen = os.path.join(self.pangenome_dir, "pangenoma.csv")
            file_paths = [genome_locus, locus_seq, pangen_locus, pangen]
            for path in file_paths:
                if os.path.isfile(path):
                    os.remove(path)
            DBcmd = self.pangenome_dir + "/DBGenerator.py ../ffn_files"
            G.run_shell(DBcmd, printcmd=True, logcmd=True, log=True, printoption=True)


    def create_sqldb(self):
        G.logger("Run: create_sqldb(" + self.target + ")")
        if not os.path.isfile(self.target + ".db"):
            connection = sqlite3.connect(self.target + ".db")
            cursor = connection.cursor()
            sql_tables = [
                """create table genomas_locus (cod text, locus text);""",
                """create table pangenoma (gene text, non_unique_gene_name text,
                annotation text, no_isolates integer, no_sequences integer,
                avg_sequences_per_isolate integer, genome_fragment integer,
                order_within_fragment integer, accessory_fragment integer,
                accessory_order_with_fragment integer,
                qc text, min_group_size_nuc integer,
                max_group_size_nuc integer, avg_group_size_nuc integer);""",
                """create table pangenoma_locus (gene text, locus text);""",
                """create table locus_sequence (locus text, sequence text);""",
            ]
            for cmd in sql_tables:
                cursor.execute(cmd)
            db_data = []
            with open("genomas_locus.csv","r") as f:
                fieldnames = ["cod text", "locus text"]
                csv_file = csv.DictReader(f, delimiter="|", fieldnames=fieldnames)
                db_data = [(i["cod text"], i["locus text"]) for i in csv_file]
            cursor.executemany("INSERT INTO genomas_locus VALUES (?, ?);", db_data)

            db_data = []
            with open("pangenoma_locus.csv","r") as f:
                fieldnames = ["gene text", "locus text"]
                csv_file = csv.DictReader(f, delimiter="|", fieldnames=fieldnames)
                db_data = [(i["gene text"], i["locus text"]) for i in csv_file]
            cursor.executemany(
                "INSERT INTO pangenoma_locus VALUES (?, ?);", db_data)

            db_data = []
            with open("locus_sequence.csv","r") as f:
                fieldnames = ["locus text", "sequence text"]
                csv_file = csv.DictReader(f, delimiter="|", fieldnames=fieldnames)
                db_data = [(i["locus text"], i["sequence text"]) for i in csv_file]
            cursor.executemany(
                "INSERT INTO locus_sequence VALUES (?, ?);", db_data)

            db_data = []
            with open("pangenoma.csv","r") as f:
                fieldnames = [
                    "gene text", "non_unique_gene_name text",
                    "annotation text", "no_isolates integer",
                    "no_sequences integer",
                    "avg_sequences_per_isolate integer", "genome_fragment integer",
                    "order_within_fragment integer", "accessory_fragment integer",
                    "accessory_order_with_fragment integer", "qc text",
                    "min_group_size_nuc integer", "max_group_size_nuc integer",
                    "avg_group_size_nuc integer"]
                csv_file = csv.DictReader(f, delimiter="|", fieldnames=fieldnames)
                db_data = [(
                    i["gene text"], i["non_unique_gene_name text"],
                    i["annotation text"], i["no_isolates integer"],
                    i["no_sequences integer"],
                    i["avg_sequences_per_isolate integer"],
                    i["genome_fragment integer"],
                    i["order_within_fragment integer"],
                    i["accessory_fragment integer"],
                    i["accessory_order_with_fragment integer"], i["qc text"],
                    i["min_group_size_nuc integer"],
                    i["max_group_size_nuc integer"],
                    i["avg_group_size_nuc integer"]
                ) for i in csv_file]
            cursor.executemany(
                "INSERT INTO pangenoma VALUES (?, ?, ?, ?, ?, ?, ?,"
                    "?, ?, ?, ?, ?, ?, ?);", db_data)

            sql_index = [
                """create index genomas_locus_index on genomas_locus
                (cod, locus);""",
                """create index pangenoma_index on pangenoma(gene,
                non_unique_gene_name, annotation, no_isolates, no_sequences,
                avg_sequences_per_isolate, genome_fragment, order_within_fragment,
                accessory_fragment, accessory_order_with_fragment, qc,
                min_group_size_nuc, max_group_size_nuc, avg_group_size_nuc);""",
                """create index pangenoma_locus_index on
                pangenoma_locus(gene, locus);""",
                """create index locus_sequence_index on
                locus_sequence(locus, sequence);"""
            ]

            for cmd in sql_index:
                cursor.execute(cmd)
            connection.commit()
            connection.close()
        else:
            info ="DB exists"
            if not self.config.verbosity == 0:
                print(info)
            G.logger(info)

    def extract_sql_data(self, gene, cursor):
        if "/" in gene:
            gene_name = "-".join(gene.split("/"))
        elif " " in gene:
            gene_name = "-".join(gene.split(" "))
        else:
            gene_name = gene
        results_file = os.path.join(self.results_dir, "fasta", gene_name + ".fasta")
        with open(results_file, "w") as f:
            cursor.execute(
                "select '>' || cod ||"
                " '|' || locus_sequence.locus || '|' || pangenoma.gene || "
                "x'0a' || sequence from locus_sequence\n"
                "inner join pangenoma_locus on locus_sequence.locus = "
                "pangenoma_locus.locus inner join pangenoma on "
                "pangenoma_locus.gene = pangenoma.gene inner join "
                "genomas_locus on locus_sequence.locus = "
                "genomas_locus.locus  where pangenoma.gene = '" + gene + "';")

            for x in cursor.fetchall():
                f.write(x[0] + "\n")

    def coregene_extract(self, mode="normal"):
        info = "Run: core_gene_extract(" + self.target + ")"
        print(info)
        G.logger(info)
        total_count = []
        single_count = []
        all_core = []
        multi_annotated = []
        gff = []
        def nr_genomes():
            # count number of genomes used for pangenome analysis
            for files in [f for f in os.listdir(self.gff_dir) if f.endswith(".gff")]:
                if files not in gff:
                    gff.append(files)
            genomes = int(len(gff) - int(self.iterator))
            return genomes

        def find_genes():
            # open gene_presence_absence.csv to find single copy genes with gene names
            genomes = nr_genomes()
            PipelineStatsCollector(self.target_dir).write_stat(
                "genome assemblies for pan-genome analysis: " + str(genomes))
            filepath = os.path.join(self.pangenome_dir, "gene_presence_absence.csv")
            with open(filepath, "r") as csvfile:
                reader = csv.reader(csvfile, delimiter=",", quotechar='"')
                next(reader, None)
                for row in reader:
                    gene_name = row[0]
                    number_isolates = int(row[3])
                    number_sequences = int(row[4])
                    average_seq_per_isolate = float(row[5])

                    if number_isolates == genomes:
                        total_count.append(gene_name)
                        if number_sequences == genomes:
                            if average_seq_per_isolate == 1:
                                single_count.append(gene_name)
                                if "group" in gene_name:
                                    all_core.append(gene_name)
                                if "group" not in gene_name:
                                    if len(gene_name.split("_")) > 1:
                                        multi_annotated.append(gene_name)
                                    else:
                                        all_core.append(gene_name)

        def print_gene_stats():
            all_genes = (
                "\nContinue with " + str(len(all_core))
                + " single copy core genes")
            print("\n# of core genes: " + str(len(total_count)))
            G.logger(all_genes)
            print(all_genes)
            stats = PipelineStatsCollector(self.target_dir)
            stats.write_stat("core genes: " + str(len(total_count)))
            stats.write_stat("single copy core genes: " + str(len(all_core)))

        def write_files():
            all_core_path = os.path.join(self.results_dir, "allcoregenes")
            multi_path = os.path.join(self.results_dir, "multiannotated")
            if not len(all_core) == 0:
                with open(all_core_path, "w") as f:
                    for gene in all_core:
                        f.write(gene+"\n")
            if not len(multi_annotated) == 0:
                multi_annotated.sort()
                with open(multi_path, "w") as f:
                    for gene in multi_annotated:
                        f.write(gene+"\n")

        def get_fasta_fromDB():
            if not len(all_core) == 0:
                db_path = os.path.join(self.pangenome_dir, self.target + ".db")
                if os.path.isfile(db_path):
                    connection = sqlite3.connect(db_path)
                    cursor = connection.cursor()
                    for gene in all_core:
                        self.extract_sql_data(gene, cursor)
            else:
                error_msg = "Error: No core genes found"
                print(error_msg)
                errors.append([self.target, error_msg])
                G.logger(error_msg)

        def run_coregene_extract(mode):
            if mode == "normal":
                find_genes()
                print_gene_stats()
                write_files()
                get_fasta_fromDB()
            else:
                find_genes()
                print_gene_stats()


        run_coregene_extract(mode)

    def run_CoreGenes(self):
        print("\nCollect Results of Pangenome analysis\n")
        G.logger("Run: run_CoreGenes(" + self.target + ")")
        G.logger("Collect Results of Pangenome analysis")
        os.chdir(self.pangenome_dir)
        if not os.path.isdir(os.path.join(self.results_dir, "fasta")):
            self.copy_DBGenerator()
            self.run_DBGenerator()
            self.create_sqldb()
            G.create_directory(self.results_dir)
            G.create_directory(os.path.join(self.results_dir,"fasta"))
            self.coregene_extract()
        else:
            if len(os.listdir(os.path.join(self.results_dir, "fasta"))) <= 1:
                shutil.rmtree(os.path.join(self.results_dir, "fasta"))
                self.copy_DBGenerator()
                self.run_DBGenerator()
                self.create_sqldb()
                G.create_directory(self.results_dir)
                G.create_directory(os.path.join(self.results_dir,"fasta"))
                self.coregene_extract()
            else:
                self.coregene_extract("statistics")


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
        self.conserved_dict = {}

    def seq_alignments(self):
        G.logger("Run: seq_alignments(" + self.target + ")")
        G.create_directory(self.alignments_dir)
        run_file = os.path.join(self.results_dir, "run_prank")
        if os.path.isfile(run_file):
            os.remove(run_file)
        os.chdir(self.results_dir)
        with open(run_file, "w") as w:
            for files in os.listdir(self.fasta_dir):
                if files.endswith(".fasta"):
                    file_name = files.split(".")[0]
                    result_path = (
                        os.path.join(self.alignments_dir, file_name + ".best.fas"))
                    if not os.path.isfile(result_path):
                        fasta_path = os.path.join("fasta", files)
                        aligned_path = os.path.join("alignments", file_name)
                        w.write("prank -d=" + fasta_path + " -o=" + aligned_path + "\n")

        if os.stat(run_file).st_size == 0:
            info = "Skip " + " ".join(["parallel", "-a", run_file])
            print("\n" + info)
            G.logger(info)
        else:
            G.run_subprocess(["parallel", "-a", run_file], True, True, False, False)

        os.chdir(self.target_dir)

    def seq_consensus(self):
        G.logger("Run: seq_consensus(" + self.target + ")")
        G.create_directory(self.consensus_dir)
        run_file = os.path.join(self.results_dir, "run_consensus")
        if os.path.isfile(run_file):
            os.remove(run_file)

        os.chdir(self.results_dir)

        with open(run_file, "w") as w:
            for files in os.listdir(self.alignments_dir):
                if files.endswith(".best.fas"):
                    file_name = files.split(".")[0]
                    result_path = (
                        os.path.join("consensus", file_name + "_consens.fasta"))
                    if not os.path.isfile(result_path):
                        aligned_path = os.path.join("alignments", files)
                        w.write(
                            "consambig -sequence " + aligned_path
                            + " -outseq " + result_path +" -name "
                            + self.target + "_" + file_name + "_consensus -auto\n")

        if os.stat(run_file).st_size == 0:
            info = "Skip " + " ".join(["parallel", "-a", run_file])
            print("\n" + info)
            G.logger(info)
        else:
            G.run_subprocess(["parallel", "-a", run_file], True, True, True, False)
        os.chdir(self.target_dir)

    def conserved_seqs(self):
        G.logger("Run: conserved_seqs(" + self.target + ")")
        G.create_directory(self.blast_dir)
        conserv_seqs = []
        result_path = os.path.join(
            self.blast_dir, H.abbrev(self.target, dict_path) + "_conserved_seqs")
        for files in os.listdir(self.consensus_dir):
            if files.endswith("_consens.fasta"):
                file_path = os.path.join(self.consensus_dir, files)
                for record in SeqIO.parse(file_path, "fasta"):
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
                            seq_name = files.split("_consens")[0] + "_" + str(count)
                            seq = item
                            count += 1
                            conserv_seqs.append(["> " + seq_name + "\n", seq + "\n"])
                            self.conserved_dict.update({seq_name: seq})

        if len(conserv_seqs) > 0:
            with open(result_path, "w") as f:
                for seq in conserv_seqs:
                    f.write(seq[0])
                    f.write(seq[1])
            PipelineStatsCollector(self.target_dir).write_stat(
                "Number of conserved sequences: " + str(len(conserv_seqs)))
            return conserv_seqs

        else:
            error_msg = "Error: no conserved target sequences found"
            print(error_msg)
            errors.append([self.target, error_msg])
            G.logger(error_msg)
            PipelineStatsCollector(self.target_dir).write_stat("Number of conserved sequences: 0")
            return 1

    def run_coregeneanalysis(self):
        G.logger("Run: run_coregeneanalysis(" + self.target + ")")
        self.seq_alignments()
        self.seq_consensus()
        conserved_seqs = self.conserved_seqs()
        if conserved_seqs == 1:
            return 1
        name = "conserved"
        use_cores, inputseqs = BlastPrep(self.blast_dir, conserved_seqs, "conserved_seqs", self.config.blastseqs).run_blastprep()
        Blast(self.config, self.blast_dir, "conserved_seqs").run_blast(name, use_cores)
        return self.conserved_dict

class BlastPrep():
    def __init__(self, directory, input_list, name, maxpartsize):
        self.list_dict = {}
        self.input_list = input_list
        self.maxpartsize = maxpartsize
        self.filename = name
        self.directory = directory

    def create_listdict(self):
        groups = len(self.input_list)//self.maxpartsize
        if len(self.input_list)%self.maxpartsize > 0:
           groups = groups + 1
        for i in range(0, groups):
            if not i in self.list_dict.keys():
                self.list_dict.update({i:[]})

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
                    if not key in self.list_dict.keys():
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
                file_name = os.path.join(self.directory, self.filename + ".part-"+str(key))
                with open(file_name, "w") as f:
                    for item in self.list_dict[key]:
                        f.write(item[0])
                        inputsequences.append(item[0].split(">")[1].strip())
                        f.write(item[1])
        return inputsequences

    def run_blastprep(self):
        G.logger("Run: run_blastprep")
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

    def run_blast(self, name, use_cores):
        G.logger("Run: run_blast")

        def get_blast_cmd(filename, cores):
            gi_bin = os.path.join(self.config_dir, "NO_Blast.bin")

            if self.mode == "quality_control":
                if self.config.remoteblast:
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="megablast", remote=True,
                        max_target_seqs=1, max_hsps=1, query=blastfile, db="nt",
                        out=filename, outfmt=5)
                elif os.path.isfile(gi_bin):
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="megablast", num_threads=cores,
                        negative_gilist=gi_bin, max_target_seqs=1, max_hsps=1,
                        query=blastfile, db="nt", out=filename, outfmt=5)
                else:
                        blast_cmd = NcbiblastnCommandline(
                            cmd="blastn", task="megablast", num_threads=cores,
                            max_target_seqs=1, max_hsps=1, query=blastfile, db="nt",
                            out=filename, outfmt=5)
                return blast_cmd

            if self.mode == "conserved_seqs":
                if self.config.remoteblast:
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="dc-megablast", query=blastfile,
                        db="nt", max_target_seqs=2000, evalue=500,
                        out=filename, outfmt=5, remote=True)
                elif os.path.isfile(gi_bin):
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="dc-megablast", query=blastfile,
                        db="nt", max_target_seqs=2000, evalue=500,
                        out=filename, outfmt=5, negative_gilist=gi_bin,
                        num_threads=cores)
                else:
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="dc-megablast", query=blastfile,
                        db="nt", max_target_seqs=2000, evalue=500,
                        out=filename, outfmt=5,
                        num_threads=cores)
                return blast_cmd

            if self.mode == "primer":
                if self.config.remoteblast:
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="blastn-short", query=blastfile,
                        db="nt", evalue=500, out=filename, outfmt=5,
                        remote=True)
                elif os.path.isfile(gi_bin):
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="blastn-short", query=blastfile,
                        db="nt", evalue=500, out=filename, outfmt=5,
                        negative_gilist=gi_bin, num_threads=cores)
                else:
                    blast_cmd = NcbiblastnCommandline(
                        cmd="blastn", task="blastn-short", query=blastfile,
                        db="nt", evalue=500, out=filename, outfmt=5,
                        num_threads=cores)
                return blast_cmd

        def search_blastfiles(directory):
            blast_files = []
            for files in os.listdir(directory):
                if fnmatch.fnmatch(files, "*.part*"):
                    blast_files.append(files)
            return blast_files

        blastfiles = search_blastfiles(self.directory)
        if len(blastfiles) > 0:
            blastfiles.sort(key=lambda x: int(x.split("part-")[1]))
            start = time.time()
            os.chdir(self.directory)
            for blastfile in blastfiles:
                blast_cmd = False
                filename = name + "_" +str(blastfile).split("-")[1] + "_results.xml"
                results_path = os.path.join(self.directory, filename)
                if not os.path.isfile(results_path):
                    blast_cmd = get_blast_cmd(filename, use_cores)
                else:
                    if os.stat(results_path).st_size == 0:
                        blast_cmd = get_blast_cmd(filename, use_cores)
                    else:
                        info = "Skip Blast step for " + blastfile
                        print("\n" + info)
                        G.logger(info)

                if blast_cmd:
                    info = "Run: " + str(blast_cmd)
                    print(info)
                    G.logger(info)
                    stdout, stderr = blast_cmd()

            duration = time.time() - start
            G.logger("Blast duration:")
            G.logger(str(timedelta(seconds=duration)))
            os.chdir(self.target_dir)


class BlastParser:
    def __init__(self, configuration, results="seqs"):
        self.exception = configuration.exception
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.nontargetlist = configuration.nontargetlist
        self.start_match = []
        self.end_match = []
        self.selected = []
        self.start_match = []
        self.end_match = []
        self.mode = "normal"
        self.start = time.time()
        if results == "primer":
            self.mode = "primer"
            self.primer_dir = os.path.join(self.results_dir, "primer")
            self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
            self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
            self.maxgroupsize = 25000
            self.primer_align_dict = {}

    def blastresult_files(self, blast_dir):
        xmlblastresults = []
        G.logger("Run: blastresults_files(" + self.target + ")")
        for files in [f for f in os.listdir(blast_dir) if f.endswith(".xml")]:
            file_path = os.path.join(blast_dir, files)
            if file_path not in xmlblastresults:
                xmlblastresults.append(file_path)
        return xmlblastresults

    def get_alignmentdata(self, alignment):
        for hsp in alignment.hsps:
            score = hsp.score
            e_value = hsp.expect
            query = hsp.query
            match = hsp.match
            subject = hsp.sbjct
            subject_start = hsp.sbjct_start
            align_length = hsp.align_length
            nuc_ident = hsp.identities

            identity = False
            gi = alignment.title.split("|")[1].strip()
            db_id = alignment.title.split("|")[3].strip()
            name_long = (alignment.title.split("|")[4].split(",")[0].strip(" "))
            if re.search("PREDICTED", name_long):
                pass
            else:
                name = name_long.split(" ")
                if len(name) >= 3:
                    if "subsp" in str(" ".join(name)):
                        identity = str(" ".join(name[0:4]))
                    else:
                        identity = str(" ".join(name[0:2]))
                else:
                    identity = str(" ".join(name[0:2]))
            if (identity and gi and db_id and score and e_value) is not None:
                return identity, gi, db_id, score, e_value, query, match, subject, subject_start, align_length, nuc_ident

    def get_seq_ends(self, blast_record, alignment, query_start, query_end):
        for hsp in alignment.hsps:
            # check start of sequence for not aligned parts
            hsp.query_start
            query_start.append(hsp.query_start)
            # check end of sequence for not aligned parts
            end_len = (blast_record.query_letters - (hsp.query_start + hsp.align_length -1))
            query_end.append(end_len)

    def check_seq_ends(self, blast_record, query_start, query_end):
        if len(query_start) > 0:
            if min(query_start) >= self.config.minsize:
                seq_range = "[1:" + str(min(query_start)) + "]"
                self.start_match.append([blast_record.query, seq_range])
        if len(query_end) > 0:
            if min(query_end) >= self.config.minsize:
                seq_range = "[" + str(min(query_end)) + ":" + str(blast_record.query_letters) + "]"
                self.end_match.append([blast_record.query, seq_range])

    def write_blastsummary(self, dir_path, align_dict, result_format="json"):
        G.logger("Run: write_blastsummary(" + self.target + ")")
        file_path = os.path.join(dir_path, "blastsummary.json")
        if result_format == "json":
            with open(file_path, "w") as q:
                q.write(json.dumps(align_dict))
        else:
            pass
            # csv

    def write_primer3_input(self, selected_seqs, conserved_seq_dict):
        G.logger("Run: write_primer3_input(" + self.target + ")")
        file_path = os.path.join(self.results_dir, "primer3_input")
        controlfile_path = os.path.join(self.results_dir, ".primer3_input")
        if self.config.probe == True:
            probe = "\nPRIMER_PICK_INTERNAL_OLIGO=1"
        else:
            probe = ""
        with open(file_path, "w") as f:
            for item in selected_seqs:
                if "complete" in item[1]:
                    f.write(
                        "SEQUENCE_ID="+item[0]+"\nSEQUENCE_TEMPLATE="
                        + str(conserved_seq_dict[item[0]])
                        +"\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize)+ "-"
                        + str(self.config.maxsize) + probe + "\n=\n")
                else:
                    x = str(conserved_seq_dict[item[0]])
                    i = item[1].strip("[]").split(":")
                    subseq = x[int(i[0])-1:int(i[1])-1]
                    f.write(
                        "SEQUENCE_ID="+item[0]+"\nSEQUENCE_TEMPLATE="
                        + str(subseq)
                        + "\nPRIMER_PRODUCT_SIZE_RANGE="
                        + str(self.config.minsize)+ "-"
                        + str(self.config.maxsize) + probe + "\n=\n")

        if os.path.isfile(controlfile_path):
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
            if len(diff) > 0:
                info1 = (
                    "Due to changed settings primer design "
                    "and quality control will start from scratch")
                info2 = "Differences in primer3 input:"
                G.logger(info1)
                G.logger(info2)
                G.logger(diff)
                print(info1)
                print(info2)
                print(diff)
                primer_dir = os.path.join(self.results_dir, "primer")
                if os.path.isdir(primer_dir):
                    G.logger("Delete primer directory")
                    print("Delete primer directory")
                    shutil.rmtree(primer_dir)

                shutil.copy(file_path, controlfile_path)

        else:
            shutil.copy(file_path, controlfile_path)

    def parse_BLASTfile(self, filename):
        record_list = []
        result_handle = open(filename)
        blast_records = NCBIXML.parse(result_handle)
        record_list = list(blast_records)
        return record_list

    def get_excluded_gis(self):
        excluded_gis = []
        gi_file = os.path.join(self.config_dir, "NO_Blast.gi")
        if os.path.isfile(gi_file):
            if os.stat(gi_file).st_size > 0:
                with open(gi_file, "r") as f:
                    for line in f:
                        gi = line.strip()
                        excluded_gis.append(str(gi))
        return excluded_gis

    def parse_blastrecords(self, blast_record):
        align_dict = {}
        hits = []
        if self.exception:
            exception = ' '.join(self.exception.split("_"))
        else:
            exception = None
        query_start = []
        query_end = []
        query_length = blast_record.query_length

        for alignment in blast_record.alignments:
            align_dict.update({blast_record.query: {}})
            (
                identity, gi, db_id, score, e_value, query, match,
                subject, subject_start, align_length, nuc_ident

            ) = self.get_alignmentdata(alignment)

            perc_coverage = round(100/query_length * align_length,0)
            perc_ident = round(100/align_length * nuc_ident, 0)

            if self.config.nolist:
                if "subsp" in self.target:
                    targetspecies = "subsp.".join(" ".join(str(self.target).split("_")).split("subsp"))
                    if not (
                        str(identity) == str(targetspecies) or
                        str(identity) == str(exception)):
                        ids = {identity: {
                            "gi": gi, "db_id": db_id, "score": score,
                            "e_value": e_value, "query": query,
                            "match": match, "subject": subject,
                            "subject_start": subject_start,
                            "perc_coverage": perc_coverage,
                            "perc_ident": perc_ident}}
                        if ids not in hits:
                            hits.append(ids)
                        if self.mode == "normal":
                            self.get_seq_ends(blast_record, alignment, query_start, query_end)
                else:
                    targetspecies = " ".join(str(self.target).split("_"))
                    if not (
                        str(identity) == str(targetspecies) or
                        str(identity) == str(exception)):
                        ids = {identity: {
                            "gi": gi, "db_id": db_id, "score": score,
                            "e_value": e_value, "query": query,
                            "match": match, "subject": subject,
                            "subject_start": subject_start,
                            "perc_coverage": perc_coverage,
                            "perc_ident": perc_ident}}
                        if ids not in hits:
                            hits.append(ids)
                        if self.mode == "normal":
                            self.get_seq_ends(blast_record, alignment, query_start, query_end)
            else:
                if not str(identity) == str(exception):
                    for species in self.nontargetlist:
                        if str(identity) == str(species):
                            ids = {identity: {
                                "gi": gi, "db_id": db_id, "score": score,
                                "e_value": e_value, "query": query,
                                "match": match, "subject": subject,
                                "subject_start": subject_start,
                                "perc_coverage": perc_coverage,
                                "perc_ident": perc_ident}}
                            if ids not in hits:
                                hits.append(ids)
                            if self.mode == "normal":
                                self.get_seq_ends(blast_record, alignment, query_start, query_end)

            align_dict.update({blast_record.query: hits})
        if self.mode == "normal":
            self.check_seq_ends(blast_record, query_start, query_end)
        return align_dict

    def get_selected_sequences(self, align_dict):
        selected_seqs = []
        excluded_seqs = []
        G.logger("Run: get_selected_sequences( " + self.target + ")")
        for key in align_dict:
            if len(align_dict[key]) == 0:
                selected_seqs.append([key, "complete"])
            else:
                excluded_seqs.append(key)
            for item in self.start_match:
                if not (
                        [item[0], "complete"] in selected_seqs or
                        item in selected_seqs):
                    selected_seqs.append(item)
                    if item[0] in excluded_seqs:
                        excluded_seqs.remove(item[0])
            for item in self.end_match:
                if not (
                        [item[0], "complete"] in selected_seqs or
                        item in selected_seqs):
                    selected_seqs.append(item)
                    if item[0] in excluded_seqs:
                        excluded_seqs.remove(item[0])
        info = (
                "\nselected sequences: "+ str(len(selected_seqs)) +
                "\nexcluded sequences: "+ str(len(excluded_seqs)))
        G.logger(info)
        print(info)
        return selected_seqs

    def get_seq_range(self, db_id, overhang=2000):
        accession = db_id[0]
        seq_start = int(db_id[1])
        blastdbcmd = "blastdbcmd -db nt -entry " +  accession + " -outfmt %l"
        seqlen_cmd = G.read_shelloutput(blastdbcmd, printcmd=False, logcmd=False, printoption=False)
        seq_len = int(seqlen_cmd[0])
        if seq_start > overhang:
            start = seq_start - overhang
        else:
            start = 1
        if seq_start + overhang < seq_len - 1:
            stop = seq_start + overhang
        else:
            stop = seq_len - 1
        seq_cmd = (
            " ".join(["blastdbcmd", "-db", "nt", "-entry", accession,
            "-range", str(start) + "-" + str(stop), "-outfmt", "%f"]))
        fasta = G.read_shelloutput(seq_cmd, printcmd=False, logcmd=False, printoption=False)

        return fasta

    def write_nontarget_sequences(self, files):
        statuscount = 0
        part = files.split(".txt")[0][-1:]
        with open(os.path.join(self.primer_qc_dir, files), "r") as f:
            filename = "BLASTnontarget" + str(part) + ".sequences"
            print("Start " + filename)
            G.logger("Start " + filename)
            if not os.path.isfile(os.path.join(self.primer_qc_dir,filename)):
                with open(os.path.join(self.primer_qc_dir,filename), "w") as r:
                    for line in f:
                        statuscount = statuscount + 1
                        line = line.strip("\n")
                        db_id = line.split(" ")[0]
                        sbjct_start = line.split(" ")[1]
                        fasta_seq = self.get_seq_range([db_id, sbjct_start], 2000)
                        fasta_info = "\n".join(fasta_seq) + "\n"
                        r.write(fasta_info)
        print("Finished " + filename)
        G.logger("Finished " + filename)
        return

    def write_DBIDS(self, prefix, part, suffix, info):
        filename = prefix + str(part) + suffix
        filepath = os.path.join(self.primer_qc_dir, filename)
        with open(filepath, "a+") as f:
            f.write(info)

    def create_primerBLAST_DBIDS(self):
        print("\nGet sequence accessions of BLAST hits\n")
        G.create_directory(self.primer_qc_dir)
        DBID_files = []
        idcount = 0
        part = 0
        for files in os.listdir(self.primer_qc_dir):
            if files.startswith("primerBLAST_DBIDS"):
                DBID_files.append(files)
                with open(os.path.join(self.primer_qc_dir, files), "r") as f:
                    for line in f:
                        idcount = idcount + 1
                part = part + 1

        if len(DBID_files) == 0:
            prefix = "primerBLAST_DBIDS"
            suffix = ".txt"
            for key in self.primer_align_dict.keys():
                if not len(self.primer_align_dict[key]) == 0:
                    for i in range(0, len(self.primer_align_dict[key])):
                        for species in self.primer_align_dict[key][i]:
                            db_id = self.primer_align_dict[key][i][species]['db_id']
                            sbjct_start = self.primer_align_dict[key][i][species]["subject_start"]
                            idcount = idcount + 1
                            part = idcount//self.maxgroupsize
                            info = str(db_id) + " " + str(sbjct_start) + "\n"
                            self.write_DBIDS(prefix, part, suffix, info)
                            filename = prefix + str(part) + suffix
                            if not filename in DBID_files:
                                DBID_files.append(filename)

        if idcount == 0:
            print("Error did not find any sequences for non-target DB")
            G.logger("Error did not find any sequences for non-target DB")
            return

        else:
            print("\n" + str(idcount) + " sequences for non-target DB")
            print("Start extraction of sequences from BLAST DB")
            G.logger("found " + str(idcount) + " sequences for non-target DB")
            pool = multiprocessing.Pool(processes=4)
            results = [pool.apply_async(self.write_nontarget_sequences, args=(files,)) for files in DBID_files]
            output = [p.get() for p in results]


    def conserved_blast_parser(self, conserved_seq_dict):
        align_dict = {}
        excluded_gis = self.get_excluded_gis()
        file_path = os.path.join(self.blast_dir, "blastsummary.json")
        if os.path.isfile(file_path):
            print("Read blastsummary")
            with open(file_path, 'r') as f:
                for line in f:
                    align_dict = json.loads(line)

            for key in align_dict.copy().keys():
                for index, item in enumerate(align_dict.copy()[key]):
                    for species in item.copy().keys():
                        gi = item[species]['gi']
#                        evalue = item[species]['e_value']
#                        cover = item[species]['perc_coverage']
#                        ident = item[species]['perc_ident']
                        if str(gi) in excluded_gis:
                            del align_dict[key][index]
#                       # option for a e value or coverage cut off
#                        elif (cover < float(70) and ident < float(80)):
#                            del align_dict[key][index]
#                        elif evalue > 10:
#                            del align_dict[key][index]

            selected_seqs = self.get_selected_sequences(align_dict)
            self.write_primer3_input(selected_seqs, conserved_seq_dict)
        else:
            xmlblastresults = self.blastresult_files(self.blast_dir)
            nr = 1
            for filename in xmlblastresults:
                print("\nopen BLAST result file " + str(nr) + "/"  + str(len(xmlblastresults)))
                blastrecords = self.parse_BLASTfile(filename)
                print("read BLAST results file " + str(nr) + "/"  + str(len(xmlblastresults)))
                nr = nr + 1
                total = len(blastrecords)
                rec = 1
                for record in blastrecords:
                    print('\r read record ' + str(rec) + "/" + str(total), end='')
                    result = self.parse_blastrecords(record)
                    align_dict.update(result)
                    rec = rec + 1

            for key in align_dict.copy().keys():
                for index, item in enumerate(align_dict.copy()[key]):
                    for species in item.copy().keys():
                        gi = item[species]['gi']
#                        evalue = item[species]['e_value']
#                        cover = item[species]['perc_coverage']
#                        ident = item[species]['perc_ident']
                        if str(gi) in excluded_gis:
                            del align_dict[key][index]
#        option for a e value or coverage cut off
#                        elif (cover < float(70) and ident < float(80)):
#                            del align_dict[key][index]
#                        elif evalue > 10:
#                            del align_dict[key][index]

            self.write_blastsummary(self.blast_dir, align_dict, "json")
            selected_seqs = self.get_selected_sequences(align_dict)
            self.write_primer3_input(selected_seqs, conserved_seq_dict)



        duration = time.time() - self.start
        G.logger("Blast parser time:")
        G.logger(str(timedelta(seconds=duration)))
        PipelineStatsCollector(self.target_dir).write_stat(
            "species specific conserved sequences: "
            + str(len(selected_seqs)))

        if len(selected_seqs) == 0:
            error_msg = "No conserved sequences without non-target match found"
            print(error_msg)
            G.logger(error_msg)
            errors.append([self.target, error_msg])
            return 1
        else:
            return 0

    def primer_blast_parser(self):
        summaryfile = os.path.join(self.primerblast_dir, "blastsummary.json")
        if not os.path.isfile(summaryfile):
            xmlblastresults = self.blastresult_files(self.primerblast_dir)
            nr = 1
            for filename in xmlblastresults:
                print("\nopen BLAST result file " + str(nr) + "/"  + str(len(xmlblastresults)))
                blastrecords = self.parse_BLASTfile(filename)
                print("read BLAST results file " + str(nr) + "/"  + str(len(xmlblastresults)))
                nr = nr + 1
                total = len(blastrecords)
                rec = 1
                for record in blastrecords:
                    print('\r read record ' + str(rec) + "/" + str(total), end='')
                    result = self.parse_blastrecords(record)
                    self.primer_align_dict.update(result)
                    rec = rec + 1

            self.write_blastsummary(self.primerblast_dir, self.primer_align_dict, "json")
        else:
            print("\nRead primerblast summary")
            with open(summaryfile, 'r') as f:
                for line in f:
                    self.primer_align_dict = json.loads(line)

        self.create_primerBLAST_DBIDS()

        duration = time.time() - self.start
        G.logger("Primer blast parser time:")
        G.logger(str(timedelta(seconds=duration)))

    def run_blastparser(self, conserved_seq_dict):
        G.logger("Run: run_blastparser(" + self.target + ")")
        if self.mode == "primer":
            self.primer_blast_parser()
        else:
            exitstatus = self.conserved_blast_parser(conserved_seq_dict)
            return exitstatus

class PrimerDesign():
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.primer_dir = os.path.join(self.results_dir, "primer")
        self.primer3_dict = {}

    def run_primer3(self):
        G.logger("Run: run_primer3(" + self.target + ")")
        input_file = os.path.join(self.results_dir, "primer3_input")
        output_file = os.path.join(self.primer_dir, "primer3_output")
        settings_file = os.path.join(pipe_dir, "p3parameters")
        if not os.path.isfile(output_file):
            primer3cmd = [
                "primer3_core", "-p3_settings_file=" + settings_file,
                "-echo_settings_file", "-output=" + output_file, input_file]
            G.run_subprocess(primer3cmd, True, True, False, True)
        else:
            info = "Skip primerdesign with primer3"
            G.logger(info)
            print(info)

    def parse_Primer3_output(self, path_to_file):
        G.logger("Run: parse_Primer3_output(" + self.target + ")")
        def parseSeqId(key, value):
            if key.startswith("SEQUENCE_ID"):
                if "group" in value:
                    value = "g" + value.split("group_")[1]
                seq_id = value
                self.primer3_dict.update({seq_id: {"Primer_pairs": None}})
                primer3_list.append(seq_id)

        def parseTemplate(key, value):
            if key.startswith("SEQUENCE_TEMPLATE"):
                self.primer3_dict[primer3_list[-1]].update({"Template_seq": value})

        def countPrimer(key, value):
            if key.startswith("PRIMER_PAIR_NUM_RETURNED"):
                self.primer3_dict[primer3_list[-1]]["Primer_pairs"] = int(value)
                for count in range(0, int(value)):
                    self.primer3_dict[primer3_list[-1]].update({"Primer_pair_"+str(count): {}})

        def parseRightPrimer(key, value):
            if key.startswith("PRIMER_RIGHT"):
                count = key.split("_")[2]
                if key.endswith("_PENALTY"):
                    primer_rpen = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_R_penalty": float(primer_rpen)})
                if key.endswith("_SEQUENCE"):
                    primer_rseq = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_R_sequence": primer_rseq})
                if key.endswith("_TM"):
                    right_TM = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_R_TM": float(right_TM)})

        def parseLeftPrimer(key, value):
            if key.startswith("PRIMER_LEFT"):
                count = key.split("_")[2]
                if key.endswith("_PENALTY"):
                    primer_lpen = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_L_penalty": float(primer_lpen)})
                if key.endswith("_SEQUENCE"):
                    primer_lseq = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_L_sequence": primer_lseq})
                if key.endswith("_TM"):
                    left_TM = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_L_TM": float(left_TM)})

        def parseInternalProbe(key, value):
            if key.startswith("PRIMER_INTERNAL"):
                count = key.split("_")[2]
                if key.endswith("_PENALTY"):
                    primer_ipen = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_I_penalty": float(primer_ipen)})
                if key.endswith("_SEQUENCE"):
                    primer_iseq = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_I_sequence": primer_iseq})
                if key.endswith("_TM"):
                    int_TM = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_I_TM": float(int_TM)})

        def parsePrimerPair(key, value):
            if key.startswith("PRIMER_PAIR"):
                count = key.split("_")[2]
                if (key.endswith('_PENALTY') and "WT" not in key):
                    primer_ppen = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"primer_P_penalty": float(primer_ppen)})
                if (key.endswith('_PRODUCT_SIZE') and "WT" not in key):
                    prod_size = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"product_size": int(prod_size)})
                if (key.endswith("_PRODUCT_TM") and "WT" not in key):
                    prod_TM = value
                    self.primer3_dict[primer3_list[-1]]["Primer_pair_"+str(count)].update({"product_TM": float(prod_TM)})

        primer3_list = []
        info = "Run: parse_Primer3_output(" + self.target + ")"
        print(info)
        with open(path_to_file, "r") as p:
            for line in p:
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

    def write_primer3_data(self):
        file_path = os.path.join(self.primer_dir, "primer3_summary.json")
        with open(file_path, "w") as f:
            f.write(json.dumps(self.primer3_dict))

    def run_primerdesign(self):
        G.logger("Run: run_primerdesign(" + self.target + ")")
        G.create_directory(self.primer_dir)
        self.run_primer3()
        p3_output = os.path.join(self.primer_dir, "primer3_output")
        self.parse_Primer3_output(p3_output)
        self.write_primer3_data()
        return self.primer3_dict


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
            error_msg = "\nNo primers found for conserved sequences"
            PipelineStatsCollector(self.target_dir).write_stat(
                "potential Primersystems: 0")
            print(error_msg)
            G.logger(error_msg)
            errors.append([self.target, error_msg])
            return 1

        else:
            info = "\nNumber of potential Primersystems found " + str(len(self.primerlist)//2)
            print(info)
            G.logger(info)
            PipelineStatsCollector(self.target_dir).write_stat(
                "potential primersystems: "
                + str(len(self.primerlist)//2))
            return 0

    def get_blast_input(self, item):
        # start new primername definition here
        # start 15.12.2017
        # get primername without direction split("_")
        p_fwd_name = (
            H.abbrev(self.target, dict_path) + "_" + item[0]
            + "_P" + item[1].split("_")[-1]) + "_F"
        p_rev_name = (
            H.abbrev(self.target, dict_path) + "_" + item[0]
            + "_P" + item[1].split("_")[-1]) + "_R"
        p_fwd_seq = self.primer3_dict[item[0]][item[1]]['primer_L_sequence']
        p_rev_seq = self.primer3_dict[item[0]][item[1]]['primer_R_sequence']
        self.primerlist.append([">"+p_fwd_name + "\n", p_fwd_seq + "\n"])
        self.primerlist.append([">"+p_rev_name + "\n", p_rev_seq + "\n"])

    def collect_selected_primerinfo(self, selected_seqs, mode):
        G.logger("Run: collect_selected_primerinfo(" + self.target + ")")
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
                primerpair = "Primer_pair_"+ primer_name.split("_P")[-1]
                template_seq = self.primer3_dict[target_id]["Template_seq"]
                x = self.primer3_dict[target_id][primerpair]
                pp_penalty = round(x["primer_P_penalty"], 2)
                pp_prodsize = x["product_size"]
                pp_prodTM = round(x["product_TM"], 2)

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
                    if not info in val_list:
                        val_list.append(info)
                if mode == "mfold":
                    info = [
                        target_id, primerpair, lseq, rseq, template_seq,
                        primer_name]
                    if not info in val_list:
                        val_list.append(info)
                if mode == "dimercheck":
                    info = [
                        primer_name, lseq, rseq]
                    if not info in val_list:
                        val_list.append(info)
                if mode == "results":
                    amp_seq = x["amplicon_seq"]
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
                        primer_name, ppc, pp_penalty, target_id,
                        lseq, lTM, lpen,
                        rseq, rTM, rpen,
                        iseq, iTM, ipen,
                        pp_prodsize, pp_prodTM, amp_seq,
                        template_seq]
                    if not info in val_list:
                        val_list.append(info)

            except Exception:
                G.logger("error in collect_selected_primerinfo()" + str(sys.exc_info()))

        return val_list

    def make_nontargetDB(self, inputfiles):
        db_name = inputfiles
        db_path = os.path.join(self.primer_qc_dir, inputfiles + ".sqlite3.db")
        if not os.path.isfile(db_path):
            G.logger("Start index non-target DB " + inputfiles)
            print("\nStart index non-target DB " + inputfiles)
            self.index_Database(db_name)
            print("Done indexing non-target DB " + inputfiles)
            G.logger("Done indexing non-target DB " + inputfiles)

    def index_Database(self, db_name):
        start = time.time()
        os.chdir(self.primer_qc_dir)
        cmd = "IndexDb.sh " + db_name + " 9"
        G.run_shell(cmd, printcmd=True, logcmd=True, log=False, printoption="")
        os.chdir(self.primer_dir)
        end = time.time() - start
        G.logger("Run: index_Database(" + db_name + ") time: " + str(timedelta(seconds=end)))

    def prepare_MFEprimer_Dbs(self, primerinfos):
        G.logger("Run: prepare_MFEprimer_Dbs(" + self.target + ")")
        G.create_directory(self.primer_qc_dir)

        def create_template_db(primerinfos):
            wrote = []
            file_path = os.path.join(self.primer_qc_dir, "template.sequences")
            with open(file_path, "w") as f:
                for nameF, seqF, nameR, seqR, templ_seq in primerinfos:
                    if not templ_seq in wrote:
                        primername = "_".join(nameF.split("_")[0:-2])
                        f.write(">" + primername + "\n" + templ_seq + "\n")
                        wrote.append(templ_seq)

        def get_QC_data():
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

        def create_assembly_db(db_name):
            # add option to choose a folder for Reference genomes?
            qc_data = get_QC_data()

            remove = []
            qc_acc = []
            assembly_dict = {}
            ref_assembly = []
            for item in qc_data:
                accession = item[0]
                assembly_stat = item[2]
                rRNA = item[4]
                tuf = item[6]
                recA = item[8]
                dnaK = item[10]
                pheS = item[12]
                qc_list = rRNA, tuf, recA, dnaK, pheS
                assembly_dict.update({accession: assembly_stat})

                for qc_gene in qc_list:
                    if "passed QC" == qc_gene:
                        if not accession in qc_acc:
                            qc_acc.append(accession)
                    elif "" == qc_gene:
                        if not accession in qc_acc:
                            qc_acc.append(accession)
                    else:
                        if not accession in remove:
                            remove.append(accession)

            if self.config.ignore_qc == True:
                check = set(qc_acc)
            else:
                check = set(qc_acc) - set(remove)
            for item in check:
                if len(ref_assembly) < self.referencegenomes:
                    if assembly_dict[item] == "Complete Genome":
                        if not item in ref_assembly:
                            ref_assembly.append(item)

            if len(ref_assembly) < self.referencegenomes:
                for item in check:
                    if len(ref_assembly) < self.referencegenomes:
                        if assembly_dict[item] == "Chromosome":
                            if not item in ref_assembly:
                                ref_assembly.append(item)

            if len(ref_assembly) < self.referencegenomes:
                for item in check:
                    if len(ref_assembly) < self.referencegenomes:
                        if assembly_dict[item] == "Scaffold":
                            if not item in ref_assembly:
                                ref_assembly.append(item)

            if len(ref_assembly) < self.referencegenomes:
                for item in check:
                    if len(ref_assembly) < self.referencegenomes:
                        if not item in ref_assembly:
                            ref_assembly.append(item)

            target_fasta = []
            for files in os.listdir(self.fna_dir):
                for item in ref_assembly:
                    acc = "v".join(item.split("."))
                    if acc in files:
                        if files.endswith(".fna"):
                            with open(os.path.join(self.fna_dir, files)) as f:
                                for line in f:
                                    target_fasta.append(line)

            file_path = os.path.join(self.primer_qc_dir, db_name)
            with open(file_path, "w") as fas:
                for item in target_fasta:
                    fas.write(item)


        def make_templateDB():
            db_name = "template.sequences"
            db_path = os.path.join(self.primer_qc_dir, db_name + ".sqlite3.db")
            if not os.path.isfile(db_path):
                G.logger("create template DB")
                print("\ncreate template DB")
                create_template_db(primerinfos)
                G.logger("index template DB")
                print("index template DB")
                self.index_Database(db_name)
                print("Done indexing template DB")

        def make_assemblyDB():
            db_name = H.abbrev(self.target, dict_path) + ".genomic"
            db_path = os.path.join(self.primer_qc_dir, db_name + ".sqlite3.db")
            if not os.path.isfile(db_path):
                G.logger("create target genome assembly DB")
                print("\ncreate target genome assembly DB")
                create_assembly_db(db_name)
                G.logger("index target genome assembly DB")
                print("index target genome assembly DB")
                self.index_Database(db_name)
                print("Done indexing genome assembly DB")


        process_make_templateDB = Process(target=make_templateDB)
        process_make_assemblyDB = Process(target=make_assemblyDB)
        process_make_templateDB.start()
        process_make_assemblyDB.start()
        process_make_templateDB.join()
        process_make_assemblyDB.join()

        # parallelization try
        pool = multiprocessing.Pool(processes=4)
        results = [pool.apply_async(self.make_nontargetDB, args=(inputfiles,)) for inputfiles in self.dbinputfiles]
        output = [p.get() for p in results]


    def MFEprimer_template(self, primerinfo):
        [nameF, seqF, nameR, seqR, templ_seq] = primerinfo
        with tempfile.NamedTemporaryFile(mode='w+', dir=self.primer_qc_dir, prefix="primer", suffix=".fa", delete=False) as primefile:
            primefile.write(">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")

        db = "template.sequences"
        cmd = (
            "MFEprimer.py -i " + primefile.name + " -d " + db
            + " -k 9 --tab --ppc 10")

        result = G.read_shelloutput(cmd, printcmd=False, logcmd=False, printoption=False)
        os.unlink(primefile.name)
        if len(result) == 2:
            val = result[1].split("\t")
            pp_F = "_".join(val[1].split("_")[0:-1])
            pp_R = "_".join(val[2].split("_")[0:-1])
            p_F = "_".join(val[1].split("_")[0:-2])
            primername = val[3]
            ppc = float(val[4])
            if (
                pp_F == pp_R and
                p_F == primername and
                ppc >= float(self.mfethreshold)):
                ppc_val = ppc - float(self.mfethreshold)
                return [[nameF, seqF, nameR, seqR, templ_seq, ppc_val], result]
            else:
                return [[None], result]
        else:
            return [[None], result]

    def MFEprimer_nontarget(self, primerinfo, inputfiles):
        nameF, seqF, nameR, seqR, templ_seq, ppc_val = primerinfo
        with tempfile.NamedTemporaryFile(mode='w+', dir=self.primer_qc_dir, prefix="primer", suffix=".fa", delete=False) as primefile:
            primefile.write(">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")
        db = inputfiles
        cmd = (
            "MFEprimer.py -i " + primefile.name + " -d " + db
            + " -k 9 --tab --ppc 10")
        result = G.read_shelloutput(cmd, printcmd=False, logcmd=False, printoption=False)
        os.unlink(primefile.name)
        if not len(result) == 1:
            for index, item in enumerate(result):
                if index > 0:
                    val = item.split("\t")
                    result_ppc = float(val[4])
                    if result_ppc > ppc_val:
                        return [[None], result]

        return [primerinfo, result]


    def MFEprimer_assembly(self, primerinfo):
        target_product = []
        nameF, seqF, nameR, seqR, templ_seq, ppc_val = primerinfo
        with tempfile.NamedTemporaryFile(mode='w+', dir=self.primer_qc_dir, prefix="primer", suffix=".fa", delete=False) as primefile:
            primefile.write(">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")
        db = H.abbrev(self.target, dict_path) + ".genomic"
        cmd = (
            "MFEprimer.py -i " + primefile.name + " -d " + db
            + " -k 9 --tab --ppc 10")
        result = G.read_shelloutput(cmd, printcmd=False, logcmd=False, printoption=False)
        os.unlink(primefile.name)
        for index, item in enumerate(result):
            if index > 0:
                val = item.split("\t")
                result_ppc = float(val[4])
                product_len = int(val[5])
                targetID = val[3]
                if result_ppc == ppc_val + self.mfethreshold:  #new commandline
                    target_product.append(targetID)
                elif result_ppc > ppc_val:
                    return [[None], result]
        counts = Counter(target_product)
        for item in counts.keys():
            if counts[item] == 1:
                return [primerinfo, result]
            else:
                return [[None], result]



    def MFEprimer_QC(self, primerinfos):
        # option: also allow user provided non-target database created with
        # MFEprimer for primer QC
        def select_and_write_MFEprimer_results(input_list, name):
            outputlist = []
            with open("MFEprimer_" + name + ".csv", "w") as f:
                writer = csv.writer(f)
                for item in input_list:
                    # item[0] is the primerinfo
                    if len(item[0]) > 1:
                        if not item[0] in outputlist:
                            outputlist.append(item[0])
                    # item[1] are the results of MFEprimer
                    for values in item[1]:
                        val = values.split("\t")
                        writer.writerow(val)

            return outputlist

        G.logger("Run: MFEprimer_QC(" + self.target + ")")
        print("Start primer quality control(" + self.target + ")")
        os.chdir(self.primer_qc_dir)

        info0 = str(len(primerinfos)) + " primer pair(s) to check"
        print("\n" + info0 + "\n")
        G.logger(info0)

        print("Start MFEprimer with template DB\n")
        G.logger("Start MFEprimer with template DB")
        template_list = G.run_parallel(self.MFEprimer_template, primerinfos)
        check_nontarget = select_and_write_MFEprimer_results(template_list, "template")

        info1 = str(len(check_nontarget)) + " primer pair(s) with good target binding"
        print("\n\n" + info1 + "\n")
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs with good target binding: "
            + str(len(check_nontarget)))

        G.logger(info1)

        nontarget_lists = []
        print("\nStart MFEprimer with nontarget DB\n")
        G.logger("Start MFEprimer with nontarget DB")
        for index, inputfiles in enumerate(self.dbinputfiles):
            print("nontarget DB " + str(index+1) + "/" + str(len(self.dbinputfiles)))
            nontarget_list = G.run_parallel(self.MFEprimer_nontarget, check_nontarget, inputfiles)
            nontarget_lists = list(itertools.chain(nontarget_lists, nontarget_list))

        # if the MFEprimer_nontarget.csv has only the table header and no results,
        # then no primer binding was detected and the primers passed the QC
        check_assembly = select_and_write_MFEprimer_results(nontarget_lists, "nontarget")

        info2 = str(len(check_assembly)) + " primer pair(s) passed non-target PCR check"
        print("\n\n" +  info2 + "\n")
        G.logger(info2)
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after non-target QC: "
            + str(len(check_assembly)))

        print("\nStart MFEprimer with assembly DB\n")
        G.logger("Start MFEprimer with assembly DB")

        assembly_list =  G.run_parallel(self.MFEprimer_assembly, check_assembly)
        check_final = select_and_write_MFEprimer_results(assembly_list, "assembly")

        info3 = str(len(check_final)) + " primer pair(s) passed secondary PCR amplicon check\n"
        print("\n\n" +  info3 + "\n")
        G.logger(info3)
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after secondary amplicon QC: "
            + str(len(check_final)))
        os.chdir(self.primer_dir)

        ### new 07.11.2018 add PPC to results file
        for item in template_list:
            nameF = item[0][0]
            if not nameF == None:
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
        G.logger("Run: mfold_analysis(" + self.target + ")")
        if not self.config.verbosity == 0:
            print("Run: mfold_analysis(" + self.target + ")")

        if os.path.isdir(self.mfold_dir):
            shutil.rmtree(self.mfold_dir)
        G.create_directory(self.mfold_dir)
        os.chdir(self.mfold_dir)

        for mfoldinput in mfoldinputlist:
            target_id, primerpair, pcr_product = self.prep_mfold(mfoldinput)
            self.primer3_dict[target_id][primerpair].update({"amplicon_seq": pcr_product})
        os.chdir(self.target_dir)

    def prep_mfold(self, mfoldinput):
        target_id, primerpair, lseq, rseq, template_seq, primer_name = mfoldinput
        # This removes the Genus species string to shorten the
        # name for mfold (especially for subspecies names). mfold has
        # problems with too long filenames / paths
        short_name = primer_name.split(H.abbrev(self.target, dict_path) + "_")[1]
        dir_name = target_id
        subdir = primerpair
        primer_fwd = lseq
        primer_rev_convert = Seq(rseq)
        rev_compl = str(primer_rev_convert.reverse_complement())
        dir_path = os.path.join(self.mfold_dir, dir_name)
        subdir_path = os.path.join(dir_path, subdir)
        pcr_name = short_name + "_PCR"
        pcr_prod = (
            template_seq[template_seq.index(primer_fwd):template_seq.index(rev_compl)])
        pcr_product = pcr_prod+rev_compl
        if len(pcr_product) >= self.config.minsize:
            G.create_directory(dir_path)
            G.create_directory(subdir_path)
            self.run_mfold(subdir_path, pcr_name, pcr_product)

        return target_id, primerpair, pcr_product

    def run_mfold(self, subdir_path, seq_name, description):
        file_path = os.path.join(subdir_path, seq_name)
        with open(file_path, "w") as f:
            f.write("> " + seq_name + "\n" + description)
        if not os.path.isfile(file_path + ".det"):
            os.chdir(subdir_path)
            mfold_cmd = ["mfold", "SEQ=" + seq_name, "NA=DNA", "T=60", "MG_CONC=0.003"]
            G.run_subprocess(mfold_cmd, False, True, False, False)
            os.chdir(self.mfold_dir)

    def mfold_parser(self):
        selected_primer = []
        excluded_primer = []

        info = "Run: mfold_parser(" + self.target + ")\n"
        print(info)
        G.logger(info)
        file_list = self.find_mfold_results()
        pas =  open(os.path.join(self.mfold_dir, "mfold_passed.csv"), "w")
        pos_results = csv.writer(pas)
        pos_results.writerow(["primer", "structure", "dG", "dH", "dS", "Tm"])
        fail = open(os.path.join(self.mfold_dir, "mfold_failed.csv"), "w")
        neg_results = csv.writer(fail)
        neg_results.writerow(["primer", "structure", "dG", "dH", "dS", "Tm"])
        for mfoldfiles in file_list:
            mfold = self.read_files(mfoldfiles)
#        mfold_data = G.run_parallel(self.read_files, file_list, args=False, verbosity="none")

            if len(mfold) == 1:
                selected, passed, excluded, failed = mfold[0]
                if not (selected == None or passed == None):
                    selected_primer.append(selected)
                    filename, structure, dG, dH, dS, Tm = passed
                    pos_results.writerow([selected, structure, dG, dH, dS, Tm])
                else:
                    excluded_primer.append(excluded)
                    filename, structure, dG, dH, dS, Tm = failed
                    neg_results.writerow([excluded, structure, dG, dH, dS, Tm])

            elif len(mfold) > 1:
                test = []
                for structure_nr in mfold:
                    selected, passed, excluded, failed = structure_nr
                    if not (selected == None or passed == None):
                        test.append([selected, passed, excluded, failed])

                if len(test) == len(mfold):
                    selected = test[0][0]
                    selected_primer.append(selected)

                    for index, structure_nr in enumerate(mfold):
                        selected_name, passed, excluded, failed = structure_nr
                        filename, structure, dG, dH, dS, Tm = passed
                        if index == 0:
                            pos_results.writerow([selected, structure, dG, dH, dS, Tm])
                        else:
                            pos_results.writerow(["", structure, dG, dH, dS, Tm])
                else:
                    if mfold[0][2] == None:
                        excluded = mfold[0][0]
                    else:
                        excluded = mfold[0][2]
                    excluded_primer.append(excluded)
                    for index, structure_nr in enumerate(mfold):
                        selected, passed, excluded_name, failed = structure_nr
                        if passed == None:
                            filename, structure, dG, dH, dS, Tm = failed
                            if index == 0:
                                neg_results.writerow([excluded, structure, dG, dH, dS, Tm])
                            else:
                                neg_results.writerow(["", structure, dG, dH, dS, Tm])
                        if failed == None:
                            filename, structure, dG, dH, dS, Tm = passed
                            if index == 0:
                                neg_results.writerow([excluded, structure, dG, dH, dS, Tm])
                            else:
                                neg_results.writerow(["", structure, dG, dH, dS, Tm])

            else:
                pass

        pas.close()
        fail.close()

        ex_info = str(len(excluded_primer)) + " primer pair(s) excluded by mfold"
        pass_info = str(len(selected_primer)) + " primer pair(s) to continue"

        print("\n\n" + ex_info + "\n")
        print(pass_info)
        G.logger(ex_info)
        G.logger(pass_info)
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after mfold: "
            + str(len(selected_primer)))

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
            y
        except ValueError:
            Tm = "1111"

        return structure, dG, dH, dS, Tm

    def interpret_values(self, name, primername, mfoldvalues):
        structure, dG, dH, dS, Tm = mfoldvalues
        mfold_output = [name, structure, dG, dH, dS, Tm]

        if float(dG) <= float(self.config.mfold):
            return [None, None, primername, mfold_output]
        else:
            return [primername, mfold_output, None, None]

    def get_primername(self, name):
        # adds the genus species info again to the
        # primername after mfold
        primer_name = (
            H.abbrev(self.target, dict_path) + "_"
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
                        results.append(self.interpret_values(name, primername, mfoldvalues))

        return results

    def get_primer_for_dimercheck(self, selected_seqs, excluded_primer):
        G.logger("Run: get_primer_for_dimercheck(" + self.target + ")")
        dimercheck = []
        for primer in selected_seqs:
            primername = "_".join(primer[0].split("_")[0:-1])
            if not primername in excluded_primer:
                if not primer in dimercheck:
                    dimercheck.append(primer)

        return dimercheck

    def check_primerdimer(self, dimercheck):
        G.logger("Run: check_primerdimer(" + self.target + ")")
        def get_primer_name(item, primer_direction):
            name = (
                H.abbrev(self.target, dict_path) + "_" + primer_direction + "_"
                + "_".join(item[0].split("_")[-4:]))
            return name

        def write_dimercheck_input(item, primer_name, lseq, rseq):
            fwd_name = primer_name + "_F"
            rev_name = primer_name + "_R"
            fwd_file = os.path.join(self.dimercheck_dir, fwd_name)
            rev_file = os.path.join(self.dimercheck_dir, rev_name)
            p_file = os.path.join(self.dimercheck_dir, primer_name)
            with open(fwd_file, "w") as f:
                f.write(">"+fwd_name+"\n"+lseq+"\n>"+fwd_name + "_d\n"+lseq+"\n")
            with open(rev_file, "w") as f:
                f.write(">"+rev_name+"\n"+rseq+"\n>"+rev_name + "_d\n"+rseq+"\n")
            with open(p_file, "w") as f:
                f.write(">"+fwd_name+"\n"+lseq+"\n>"+rev_name+"\n"+rseq+"\n")

        def get_dimercheck_output(item, primer_name, lseq, rseq, choice):
            test = []
            summary_data = []
            fwd_name = primer_name + "_F"
            rev_name = primer_name + "_R"
            fwd_file = os.path.join(self.dimercheck_dir, fwd_name + "_dimer_out")
            rev_file = os.path.join(self.dimercheck_dir, rev_name + "_dimer_out")
            p_file = os.path.join(self.dimercheck_dir, primer_name + "_dimer_out")
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
                                name = file_name.split("/")[-1].split("_dimer_out")[0]
                                test.append(name)
                        else:
                            if val >= float(self.config.mpprimer) - 1:
                                name = file_name.split("/")[-1].split("_dimer_out")[0]
                                test.append(name)


            if test == names:
                if primer_name not in choice:
                    choice.append(primer_name)

            return summary_data


        def run_primerdimer_check():
            choice = []
            if os.path.isdir(self.dimercheck_dir):
                shutil.rmtree(self.dimercheck_dir)
            G.create_directory(self.dimercheck_dir)
            primer_data = self.collect_selected_primerinfo(dimercheck, mode="dimercheck")
            for item in primer_data:
                write_dimercheck_input(item, item[0], item[1], item[2])

            for files in os.listdir(self.dimercheck_dir):
                if not fnmatch.fnmatch(files, "*_dimer_out"):
                    input_file = os.path.join(self.dimercheck_dir, files)
                    output_file = os.path.join(self.dimercheck_dir, files + "_dimer_out")
                    dimer_cmd = (
                        "MPprimer_dimer_check.pl -f " + input_file + " -d 3 > "
                        + output_file)

                    G.run_shell(dimer_cmd, printcmd=False, logcmd=False, log=False, printoption=False)


            dimer_summary = []
            for item in primer_data:
                summary_data = get_dimercheck_output(item, item[0], item[1], item[2], choice)
                dimer_summary.append(summary_data)

            with open(os.path.join(self.dimercheck_dir, "dimercheck_summary.csv"), "w") as f:
                writer = csv.writer(f)
                writer.writerow(["primer", "fwd-fwd dG", "rev-rev dG", "fwd-rev dG"])
                for dimer in dimer_summary:
                    writer.writerow(dimer)
            return choice

        choice = run_primerdimer_check()
        info = str(len(choice)) + " primer left after primer-dimer check"
        PipelineStatsCollector(self.target_dir).write_stat(
            "primer pairs left after primer QC: "
            + str(len(choice)))
        print("\n" + info)
        G.logger(info)
        return choice

    def write_results(self, choice):
        G.logger("Run: write_results(" + self.target + ")")
        header = [
            "Primer name", "PPC", "Primer penalty", "Gene",
            "Primer fwd seq", "Primer fwd TM", "Primer fwd penalty",
            "Primer rev seq", "Primer rev TM", "Primer rev penalty",
            "Probe seq)", "Probe TM", "Probe penalty",
            "Amplicon size", "Amplicon TM", "Amplicon sequence",
            "Template sequence"]
        if len(choice) > 0:
            results = self.collect_selected_primerinfo(choice, mode="results")
            results.sort(key=lambda x: float(x[1]), reverse=True)
            file_path = os.path.join(
                    self.results_dir,
                    H.abbrev(self.target, dict_path) + "_primer.csv")

            with open(file_path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                for result in results:
                    writer.writerow(result)
            return results
        else:
            results = []
            return results


    def run_primer_qc(self):
        G.logger("Run: run_primer_qc(" + self.target + ")")
        total_results = []
        if self.collect_primer() == 0:
            G.create_directory(self.primerblast_dir)

            prep = BlastPrep(self.primerblast_dir, self.primerlist, "primer", self.config.blastseqs)
            use_cores, inputseqs = prep.run_blastprep()

            bla = Blast(self.config, self.primerblast_dir, "primer")
            bla.run_blast("primer", use_cores)

            self.call_blastparser.run_blastparser("primer")

            primer_qc_list = self.collect_selected_primerinfo(inputseqs, "mfeprimer")

            for files in os.listdir(self.primer_qc_dir):
                if files.startswith("BLASTnontarget") and files.endswith(".sequences"):
                    self.dbinputfiles.append(files)

            self.prepare_MFEprimer_Dbs(primer_qc_list)

            survived_MFEp = self.MFEprimer_QC(primer_qc_list)

            mfoldinput = self.collect_selected_primerinfo(survived_MFEp, "mfold")

            self.mfold_analysis(mfoldinput)

            selected_primer, excluded_primer = self.mfold_parser()

            dimercheck = self.get_primer_for_dimercheck(selected_primer, excluded_primer)

            choice = self.check_primerdimer(dimercheck)

            total_results = self.write_results(choice)
            if not total_results == []:
                info = (
                    "\nFound " + str(len(total_results))
                    + " primer system(s) for " + self.target + "\n")
                print(info)
                G.logger(info)
                duration = time.time() - self.start
                G.logger("PrimerQC time:")
                G.logger(str(timedelta(seconds=duration)))
                return total_results
            else:
                error_msg = "No compatible primers found"
                print(error_msg)
                G.logger(error_msg)
                errors.append([self.target, error_msg])
                duration = time.time() - self.start
                G.logger("PrimerQC time:")
                G.logger(str(timedelta(seconds=duration)))
                return total_results
        else:
            error_msg = "No compatible primers found"
            duration = time.time() - self.start
            G.logger("PrimerQC time:")
            G.logger(str(timedelta(seconds=duration)))
            print(error_msg)
            G.logger(error_msg)
            errors.append([self.target, error_msg])
            return total_results

class Summary:
    def __init__(self, configuration, total_results):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.pangenome_dir = os.path.join(self.target_dir, "Pangenome")
        self.results_dir = os.path.join(self.pangenome_dir, "results")
        self.blast_dir = os.path.join(self.results_dir, "blast")
        self.primer_dir = os.path.join(self.results_dir, "primer")
        self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
        self.mfold_dir = os.path.join(self.primer_dir, "mfold")
        self.summ_dir = os.path.join(self.config.path, "Summary", self.target)
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.aka = H.abbrev(self.target, dict_path)
        self.g_info_dict = {}
        if total_results is None:
            self.total_results = []
        else:
            self.total_results = total_results


    def collect_qc_infos(self, qc_gene):
        if not qc_gene == "None":
            G.logger("Run: collect_qc_infos(" + self.target +  " " +qc_gene + ")")
            qc_dir = os.path.join(self.target_dir, qc_gene + "_QC")
            with open(os.path.join(qc_dir, qc_gene + "_QC_report.csv"), "r") as r:
                reader = csv.reader(r)
                next(reader, None)
                for row in reader:
                    if ("GCF" or "GCA") in row[0]:
                        accession = ".".join("_".join(row[0].split("_")[0:-1]).split("v"))
                    else:
                        accession = "_".join(row[0].split("_")[0:-1])

                    if not accession in self.g_info_dict.keys():
                        add_accession = {
                            accession: {
                                "name": "","assemblystatus":"", "strain": "",
                                "rRNA": {"status": "", "hit": "", "GI": "", "DB_id": ""},
                                "tuf": {"status": "", "hit": "", "GI": "", "DB_id": ""},
                                "recA": {"status": "", "hit": "", "GI": "", "DB_id": ""},
                                "dnaK": {"status": "", "hit": "", "GI": "", "DB_id": ""},
                                "pheS": {"status": "", "hit": "", "GI": "", "DB_id": ""}}}

                        self.g_info_dict.update(add_accession)

                    self.g_info_dict[accession][qc_gene]["GI"] = row[1]
                    self.g_info_dict[accession][qc_gene]["DB_id"] = row[2]
                    self.g_info_dict[accession][qc_gene]["hit"] = row[3]
                    self.g_info_dict[accession][qc_gene]["status"] = row[5]


    def get_genome_infos(self):
        G.logger("Run: get_genome_infos(" + self.target + ")")
        genome_data = []
        genomedata = os.path.join(self.config_dir, "genomicdata.json")
        if os.path.isfile(genomedata):
            with open(genomedata) as q:
                for line in q:
                    assembly_records = json.loads(line)

            for result in assembly_records['DocumentSummarySet']['DocumentSummary']:
                accession = result["AssemblyAccession"]
                name = result["AssemblyName"]
                status = result["AssemblyStatus"]
                strain = "unknown"
                for i, item in enumerate(result['Biosource']['InfraspeciesList']):
                    if len(item) > 0:
                        strain = result['Biosource']['InfraspeciesList'][i]['Sub_value']

                data = [accession, name, status, strain]
                genome_data.append(data)

            for item in genome_data:
                accession = item[0]
                if accession in self.g_info_dict.keys():
                    ncbi_info = {"name": item[1], "strain": item[3], "assemblystatus":item[2]}
                    self.g_info_dict[item[0]].update(ncbi_info)

    def write_genome_info(self):
        G.logger("Run: write_genome_info(" + self.target + ")")
        """write qc infos to csv file"""
        file_name = self.aka + "_qc_sequences.csv"
        with open(os.path.join(self.summ_dir, file_name), "w") as f:
            writer = csv.writer(f)
            writer.writerow([
                "Assembly accession", "Assembly name",
                "Assembly status", "Strain",
                "rRNA", "rRNA Blast",
                "tuf", "tuf Blast",
                "recA", "recA Blast",
                "dnaK", "dnaK Blast",
                "pheS", "pheS Blast"])

            for key in self.g_info_dict:
                k = self.g_info_dict[key]
                infos = [
                    key, k["name"], k['assemblystatus'], k["strain"],
                    k["rRNA"]["status"], k["rRNA"]["hit"],
                    k["tuf"]["status"], k["tuf"]["hit"],
                    k["recA"]["status"], k["recA"]["hit"],
                    k["dnaK"]["status"], k["dnaK"]["hit"],
                    k["pheS"]["status"], k["pheS"]["hit"]]
                writer.writerow(infos)

        file_name = self.aka + "_qc_sequences_details.csv"
        with open(os.path.join(self.summ_dir, file_name), "w") as f:
            writer = csv.writer(f)
            writer.writerow([
                "Assembly accession", "Assembly name",
                "Assembly status", "Strain",
                "rRNA", "rRNA Blast", "Hit GI", "Hit DB_id",
                "tuf", "tuf Blast", "Hit GI", "Hit DB_id",
                "recA", "recA Blast", "Hit GI", "Hit DB_id",
                "dnaK", "dnaK Blast", "Hit GI", "Hit DB_id",
                "pheS", "pheS Blast", "Hit GI", "Hit DB_id"])

            for key in self.g_info_dict:
                k = self.g_info_dict[key]
                infos = [
                    key, k["name"], k['assemblystatus'], k["strain"],
                    k["rRNA"]["status"], k["rRNA"]["hit"], k["rRNA"]["GI"], k["rRNA"]["DB_id"],
                    k["tuf"]["status"], k["tuf"]["hit"], k["tuf"]["GI"], k["tuf"]["DB_id"],
                    k["recA"]["status"], k["recA"]["hit"], k["recA"]["GI"], k["recA"]["DB_id"],
                    k["dnaK"]["status"], k["dnaK"]["hit"], k["dnaK"]["GI"], k["dnaK"]["DB_id"],
                    k["pheS"]["status"], k["pheS"]["hit"], k["pheS"]["GI"], k["pheS"]["DB_id"]]
                writer.writerow(infos)


    def write_results(self, total_results):
        today = time.strftime("%Y_%m_%d", time.localtime())
        if not total_results == []:
            G.logger("Run: write_results(" + self.target + ")")
            wrote = []
            header = [
                "Primer name", "PPC", "Primer penalty", "Gene",
                "Primer fwd seq", "Primer fwd TM", "Primer fwd penalty",
                "Primer rev seq", "Primer rev TM", "Primer rev penalty",
                "Probe seq)", "Probe TM", "Probe penalty",
                "Amplicon size", "Amplicon TM", "Amplicon sequence",
                "Template sequence"]
            path = os.path.join(self.summ_dir, self.aka +"_primer.csv")
            if os.path.isfile(path):
                path = os.path.join(self.summ_dir, self.aka +"_primer" + today + ".csv")
            with open(path, "w") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                for result in total_results:
                    if not result in wrote:
                        wrote.append(result)
                        writer.writerow(result)

    def copy_pangenomeinfos(self):
        for filename in os.listdir(self.pangenome_dir):
            filepath = os.path.join(self.pangenome_dir, filename)
            if filename.endswith("core_gene_alignment.aln"):
                shutil.copy(filepath, self.summ_dir)
            if filename.endswith("_tree.newick"):
                shutil.copy(filepath, self.summ_dir)
            if filename.endswith("Rplots.pdf"):
                shutil.copy(filepath, self.summ_dir)
            if filename.endswith("summary_statistics.txt"):
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
                targetpath = os.path.join(self.summ_dir, H.abbrev(self.target, dict_path) + "_" + filename)
                shutil.copy(filepath, targetpath)


    def run_summary(self, mode="normal"):
        G.logger("Run: run_summary(" + self.target + ")")
        G.create_directory(self.summ_dir)
        self.write_results(self.total_results)
        for qc_gene in self.config.qc_gene:
            self.collect_qc_infos(qc_gene)
        self.get_genome_infos()
        self.write_genome_info()
        if mode == "last":
            if self.config.nolist:
                specieslist = []
                try:
                    file_path = os.path.join(
                        self.blast_dir, "blastsummary.json")
                    with open(file_path, 'r') as f:
                        for line in f:
                            blast_dict = json.loads(line)

                    for key in blast_dict.keys():
                        for species in blast_dict[key]:
                            for speciesname in species.keys():
                                if not speciesname in specieslist:
                                    specieslist.append(speciesname)
                except FileNotFoundError:
                    pass

                try:
                    file_path = os.path.join(
                        self.primerblast_dir, "blastsummary.json")
                    with open(file_path, 'r') as f:
                        for line in f:
                            primerblast_dict = json.loads(line)

                    for key in primerblast_dict.keys():
                        for species in primerblast_dict[key]:
                            for speciesname in species.keys():
                                if not speciesname in specieslist:
                                    specieslist.append(speciesname)
                except FileNotFoundError:
                    pass

                listpath = os.path.join(self.summ_dir, "potential_specieslist.txt")
                with open(listpath, "w") as l:
                    for speciesname in specieslist:
                        if not (
                                "." in speciesname.split(" ")[0] or
                                "-" in speciesname or
                                len(re.findall(r'[A-Z]', speciesname.split(" ")[0])) == len(speciesname.split(" ")[0])):
                                    if len(re.findall(r'[0-9]', speciesname)) == 0:
                                        l.write(speciesname + "\n")

            # copy coregenealignment, trees to summary_dir
            self.copy_pangenomeinfos()
            self.copy_config()
            PipelineStatsCollector(self.target_dir).write_stat("End: " + str(time.ctime()))
            self.copy_pipelinestats()

class PipelineStatsCollector():
    def __init__(self, target_dir):
        self.target_dir = target_dir
        today = time.strftime("%Y_%m_%d", time.localtime())
        self.statfile = os.path.join(target_dir, "pipeline_stats_" + today + ".txt")

    def check_file(self, target_dir):
        if os.path.isfile(self.statfile):
            return True
        else:
            return False

    def write_stat(self, info):
        if self.check_file:
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
        "-t", "--target", nargs="*", type=str, help="Bacterial species in format: "
        "'Genus_species' e.g. 'Lactobacillus_casei'"
        " use space to separate different species")
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
       "-x", "--exception", type=str, default=None,
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
        choices=["complete", "chromosome", "scaffold", "contig", "all", "offline"],
        help="Limit downloads of Genomes to assembly status")
    # verbosity
    parser.add_argument(
       "-v", "--verbosity", type=int, choices=[0, 1, 2], default=1,
       help="increase output verbosity (quiet=0; default=1)")
    parser.add_argument(
        "--blastseqs", type=int, choices=[100, 500, 1000, 2000, 5000],
        help="Set the number of sequences per BLAST search. "
        "Decrease the number of sequences if BLAST slows down due to low "
        "memory or if you use the BLAST remote option, default=1000",
        default=1000)
    parser.add_argument(
        "--probe", action="store_true", help="Use primer3 to design an internal probe"
        "[Experimental, not recommended!]")
    parser.add_argument(
        "--nolist", action="store_true", help="Species list is not used"
        " and only sequences without blast hits are used for primer design "
        "[Experimental, not recommended!]")
    parser.add_argument(
        "--offline", action="store_true", help="Work offline no data from"
        " NCBI is collected, use your own Genomic sequences")
    parser.add_argument(
        "--ignore_qc", action="store_true", help="Genomes which do not"
        " pass quality control are included in the analysis")
    parser.add_argument(
        "--remote", action="store_true", help="Use remote option of BLAST+"
        "(Can be very slow, not recommended)")
    # Version
    parser.add_argument(
        "-V", "--version", action="version", version="%(prog)s 1.0")
    args = parser.parse_args()
    return args

def main():
    use_configfile = False
    args = commandline()
    today = time.strftime("%Y_%m_%d", time.localtime())
    logging.basicConfig(
        filename=os.path.join(os.getcwd(), "speciesprimer_" + today + ".log"),
        level=logging.DEBUG)

    if args.target is None:
        conf_from_file = Config()
        targets = conf_from_file.get_targets()
        use_configfile = True

    else:
        targets = args.target

    for target in targets:
        if use_configfile:
            (minsize, maxsize, mpprimer, exception, target, path,
            verbosity, qc_gene, mfold, skip_download,
            assemblylevel, skip_tree, nolist,
            offline, ignore_qc, mfethreshold, remoteblast,
            blastseqs, probe) = conf_from_file.get_config(target)
            if nolist:
                nontargetlist = []
                config = CLIconf(
                    minsize, maxsize, mpprimer, exception, target, path,
                    verbosity, qc_gene, mfold, skip_download,
                    assemblylevel, nontargetlist, skip_tree,
                    nolist, offline, ignore_qc, mfethreshold, remoteblast,
                    blastseqs, probe)
            else:
                nontargetlist = H.create_non_target_list(target)
                config = CLIconf(
                    minsize, maxsize, mpprimer, exception, target, path,
                    verbosity, qc_gene, mfold, skip_download,
                    assemblylevel, nontargetlist, skip_tree,
                    nolist, offline, ignore_qc, mfethreshold, remoteblast,
                    blastseqs, probe)
        else:
            if args.nolist:
                nontargetlist = []
            else:
                nontargetlist = H.create_non_target_list(target)

            config = CLIconf(
                args.minsize, args.maxsize, args.mpprimer, args.exception,
                target, args.path, args.verbosity,
                args.qc_gene, args.mfold, args.skip_download,
                args.assemblylevel, nontargetlist,
                args.skip_tree, args.nolist, args.offline,
                args.ignore_qc, args.mfethreshold, args.remoteblast,
                args.blastseqs, args.probe)

        today = time.strftime("%Y_%m_%d", time.localtime())
        G.logger("Start log: " + target + " " + today)
        G.logger(config.__dict__)

        try:

            print("\nStart searching primer for " + target)
            target_dir = os.path.join(config.path, target)
            PipelineStatsCollector(target_dir).write_stat(
                target + " pipeline statistics:")
            PipelineStatsCollector(target_dir).write_stat("Start: " + str(time.ctime()))
            DataCollection(config).collect()
            qc_count = []
            for qc_gene in config.qc_gene:
                qc = QualityControl(config).quality_control(qc_gene)
                qc_count.append(qc)
            if not sum(qc_count) == 0:
                total_results = []
                Summary(config, total_results).run_summary()
            else:
                # writes QC summary in summary directory
                try:
                    total_results = []
                    Summary(config, total_results).run_summary()
                except FileNotFoundError:
                    pass
                # end
                PangenomeAnalysis(config).run_pangenome_analysis()
                CoreGenes(config).run_CoreGenes()
                conserved_seq_dict = CoreGeneSequences(config).run_coregeneanalysis()
                if not conserved_seq_dict == 1:
                    conserved = BlastParser(config).run_blastparser(conserved_seq_dict)
                    if conserved == 0:
                        primer_dict = PrimerDesign(config).run_primerdesign()
                        total_results = PrimerQualityControl(config, primer_dict).run_primer_qc()
                        Summary(config, total_results).run_summary(mode="last")

                    else:
                        Summary(config, total_results).run_summary(mode="last")
                else:
                    Summary(config, total_results).run_summary(mode="last")

        except Exception as exc:
            error_msg = "fatal error while working on " + target + " check logfile"
            PipelineStatsCollector(target_dir).write_stat("Error: " + str(time.ctime()))
            print(error_msg)
            print(exc)
            G.logger(error_msg)
            errors.append([target, error_msg])
            logging.error("fatal error while working on " + target, exc_info=True)

        except KeyboardInterrupt:
            logging.error("KeyboardInterrupt while working on " + target, exc_info=True)
            raise

    if len(errors) > 0:
        error_log = os.path.join(config.path, "Errorreport_" + today + ".txt")
        with open(error_log, "w") as f:
            print("Error report: ")
            for index, error in enumerate(errors):
                error_nr = "Error " + str(index + 1) + ":"
                print("for target " + error[0])
                print(error_nr)
                print(error[1])
                G.logger("for target " + error[0])
                G.logger(error_nr)
                G.logger(error[1])
                f.write("for target " + error[0] + "\n")
                f.write(error_nr + "\n")
                f.write(error[1] + "\n")

if __name__ == "__main__":
    main()
