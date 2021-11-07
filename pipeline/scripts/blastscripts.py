#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time
import shutil
import logging
import fnmatch
import multiprocessing
import pandas as pd
import numpy as np
from datetime import timedelta
from ipywidgets import widgets
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# paths
script_dir = os.path.dirname(os.path.abspath(__file__))
pipe_dir, tail = os.path.split(script_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')
if not pipe_dir in sys.path:
    sys.path.append(pipe_dir)

from scripts.configuration import errors
from scripts.configuration import RunConfig
from scripts.configuration import PipelineStatsCollector
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H
from basicfunctions import ParallelFunctions as P
from basicfunctions import BlastDBError


class BlastPrep():
    def __init__(self, directory, input_records, name, maxpartsize):
        self.list_dict = {}
        self.input_records = input_records
        self.maxpartsize = maxpartsize
        self.filename = name
        self.directory = directory

    def create_listdict(self):
        groups = len(self.input_records) // self.maxpartsize
        if len(self.input_records) % self.maxpartsize > 0:
            groups = groups + 1
        for i in range(0, groups):
            if i not in self.list_dict.keys():
                self.list_dict.update({i: []})

    def get_equalgroups(self):
        self.input_records.sort(key=lambda x: int(len(x.seq)), reverse=True)
        list_start = 0
        list_end = len(self.input_records) - 1
        removed_key = []
        key = 0
        i = list_start
        while key in self.list_dict.keys():
            if key in removed_key:
                key = key + 1
            else:
                item = self.input_records[i]
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
        for key in self.list_dict.keys():
            if len(self.list_dict[key]) > 0:
                file_name = os.path.join(
                    self.directory, self.filename + ".part-"+str(key))
                SeqIO.write(self.list_dict[key], file_name, "fasta")


    def run_blastprep(self):
        G.comm_log("Run: run_blastprep - Preparing files for BLAST")
        self.create_listdict()
        self.get_equalgroups()
        cores = multiprocessing.cpu_count()
        self.write_blastinput()
        return cores


class Blast(RunConfig):
    def __init__(self, configuration, directory, mode):
        RunConfig.__init__(self, configuration)
        self.directory = directory
        self.mode = mode
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

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
        total = len(blastfiles)
        if total > 0:
            blastfiles.sort(key=lambda x: int(x.split("part-")[1]))
            start = time.time()
            os.chdir(self.directory)
            for i, blastfile in enumerate(blastfiles):
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
                        G.comm_log("> " + info)
                        self.progress.value = (i+1)/total

                if blast_cmd:
                    try:
                        G.run_subprocess(blast_cmd)
                        self.progress.value = (i+1)/total
                    except (KeyboardInterrupt, SystemExit):
                        G.keyexit_rollback(
                            "BLAST search", dp=self.directory, fn=filename)
                        raise

            duration = time.time() - start
            G.comm_log(
                "> Blast duration: "
                + str(timedelta(seconds=duration)).split(".")[0])
            os.chdir(self.target_dir)


class BlastParser(RunConfig):
    def __init__(self, configuration, results="conserved"):
        RunConfig.__init__(self, configuration)
        self.exception = configuration.exception
        self.evalue = self.config.evalue
        self.nuc_identity = self.config.nuc_identity
        self.nontargetlist = configuration.nontargetlist
        self.selected = []
        self.mode = results
        self.start = time.time()
        self.maxgroupsize = 25000

    def blastresult_files(self, blast_dir):
        blastresults = []
        for files in [f for f in os.listdir(blast_dir) if f.endswith("results.csv")]:
            file_path = os.path.join(blast_dir, files)
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

        G.comm_log("\n" + error_msg + "\n")
        logging.error("> " + error_msg, exc_info=True)
        errors.append([self.target, error_msg])
        os.remove(filename)
        G.comm_log("removed " + filename)
        raise BlastDBError(error_msg)

    def check_seq_ends(self, rejected):
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
        target_sp = " ".join([self.target.split("_")[0], H.subspecies_handler(self.target, mode="space")])
        exceptions = [target_sp]
        if self.exception != []:
            for item in self.exception:
                exception = ' '.join(item.split("_"))
                if exception not in exceptions:
                    exceptions.append(exception)

        if self.mode == "quality_control":
            # remove excluded sequences from the results
            blastdf = blastdf[~blastdf["Subject GI"].isin(excluded_gis)]
            blastdf = blastdf[~blastdf["Subject accession"].isin(excluded_gis)]
            blastdf = blastdf.sort_values(
                ["Query Seq-id", "Bit score"], ascending=False)
            blastdf = blastdf.drop_duplicates(["Query Seq-id"])
            blastdf = self.get_species_names_from_title(blastdf)
            mask = blastdf["Species"].str.contains("|".join(exceptions))
            blastdf.loc[mask, "QC status"] = "passed QC"
            blastdf.loc[~mask, "QC status"] = "failed QC"
            blastdf["Target species"] = target_sp
            QC_results = blastdf[[
                    "Query Seq-id", "Subject GI", "Subject accession",
                    "Species", "Target species", "QC status"]]
            na = pd.DataFrame()
            return QC_results, na, na, na

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

        if self.mode == "conserved":
            # Filter according to configuration
            reject = reject[reject['Percentage of identical matches'] >= self.config.nuc_identity]
            reject = reject[reject['Expect value'] <= self.config.evalue]
            # check if sequence ends without off-target matches exist
            partial = self.check_seq_ends(reject)
            speciesdata = blastdf[['Query Seq-id', 'Subject accession', "Subject Title"]]
            reject = reject[['Query Seq-id', "Subject accession"]]
        else:
            partial, speciesdata = pd.DataFrame(), pd.DataFrame()
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
        mostcommon["BLAST hits [% of queries]"] = mostcommon["BLAST hits [count]"].apply(lambda x: round(100/queries*x, 1))
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

    def bp_parse_results(self, blast_dir):
        target_dfs = []
        offtarget_dfs = []
        partial_dfs = []
        speciesdata_dfs = []
        blastresults = self.blastresult_files(blast_dir)
        excluded_gis = self.get_excluded_gis()
        G.comm_log("Excluded GI(s):", excluded_gis)
        fmt_file = os.path.join(dict_path, "blastfmt6.csv")
        header = list(pd.read_csv(fmt_file, header=None)[1].dropna())
        for i, filename in enumerate(blastresults):
            G.comm_log(
                "\nopen BLAST result file " + str(i+1)
                + "/" + str(len(blastresults)))
            try:
                blastdf = pd.read_csv(filename, sep="\t", header=None)
                blastdf.columns = header
                blastdf = blastdf.astype(
                    {"Subject GI": str, "Subject accession": str})

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

        if self.mode == "quality_control":
            if len(accept.index) == 0:
                error_msg = "No Quality Control results found."
                G.comm_log("> " + error_msg)
                errors.append([self.target, error_msg])

            return accept

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
                    G.comm_log(info)
                if os.path.isdir(self.primer_dir):
                    G.comm_log("Delete primer directory")
                    shutil.rmtree(self.primer_dir)

                shutil.copy(file_path, controlfile_path)
        else:
            shutil.copy(file_path, controlfile_path)

    def write_primer3_input(self, selected_seqs):
        G.create_directory(self.primer_dir)
        conserved_seqs = os.path.join(
            self.blast_dir, H.abbrev(self.target) + "_conserved_seqs.fas")
        conserved_seq_dict = SeqIO.to_dict(SeqIO.parse(conserved_seqs, "fasta"))
        file_path = os.path.join(self.primer_dir, "primer3_input")
        controlfile_path = os.path.join(self.primer_dir, ".primer3_input")
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

            partial_file = os.path.join(self.primer_dir, "partialseqs.csv")
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


    def run_blastparser(self):
        specific_ids = self.bp_parse_results(self.blast_dir)
        self.write_primer3_input(specific_ids)
        duration = time.time() - self.start
        info = ("Species specific conserved sequences: "
                + str(len(specific_ids)))
        G.comm_log(
            "> Blast parser time: "
            + str(timedelta(seconds=duration)).split(".")[0])
        G.comm_log("> " + info)
        PipelineStatsCollector(self.config).write_stat(info)

        if len(specific_ids) == 0:
            msg = "> No conserved sequences without non-target match found"
            G.comm_log(msg)
            errors.append([self.target, msg])
            return 1

        return 0

class PrimerBlastParser(RunConfig):
    def __init__(self, configuration, results="conserved"):
        RunConfig.__init__(self, configuration)
        self.exception = configuration.exception
        self.evalue = self.config.evalue
        self.nuc_identity = self.config.nuc_identity
        self.nontargetlist = configuration.nontargetlist
        self.selected = []
        self.mode = "primer"
        self.start = time.time()
        self.maxgroupsize = 25000
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def find_primerbinding_offtarget_seqs(self, df):
        df.loc[:, "Primer pair"] = df.loc[:, "Query Seq-id"].str.split("_").str[0:-1].apply(
                lambda x: '_'.join(x))
        df.sort_values(['Start of alignment in subject'], inplace=True)

        fwd_df = df[df["Query Seq-id"].str.endswith("_F")]
        rev_df =  df[df["Query Seq-id"].str.endswith("_R")]
        int_df = pd.merge(
                    fwd_df, rev_df, how ='inner',
                    on =['Subject accession', 'Primer pair'],
                    suffixes=("_F", "_R"))

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
        G.comm_log("> Get sequence accessions of BLAST hits")
        G.create_directory(self.primer_qc_dir)
        overhang=2000
        output_path = os.path.join(self.primer_qc_dir, "primerBLAST_DBIDS.csv")
        if os.path.isfile(output_path):
            return 0
        # data manipulation
        strandfilter = offtarget[
            'Start of alignment in subject'] > offtarget['End of alignment in subject']
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

        # data binning
        max_range = offtarget["End overhang"].max()
        stepsize = self.config.maxsize + overhang*2 + 1
        collection = []
        for i in range(1, max_range,  stepsize):
            j = i + overhang*2 + self.config.maxsize
            sub = offtarget[offtarget["Start overhang"].between(i, j, inclusive='both')]
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
            G.comm_log("> " + msg)
            errors.append([self.target, msg])
            return 1

        results.index.name = "accession"
        results.to_csv(output_path, header=None)
        return 0

    def write_nontarget_sequences(self):
        # faster but requires more RAM
        db = self.config.customdb
        if db is None:
            db = "nt"

        dbids = os.path.join(self.primer_qc_dir, "primerBLAST_DBIDS.csv")
        df = pd.read_csv(dbids, header=None)
        df.columns = ["Accession", "Start", "Stop"]
        df.sort_values(["Accession"], inplace=True)

        seqcount = len(df.index)
        G.comm_log("Found " + str(seqcount) + " sequences for the non-target DB", newline=True)
        parts = len(df.index)//self.maxgroupsize + 1
        chunks = np.array_split(df, parts)

        for part, chunk in enumerate(chunks):
            start = chunk.groupby(["Accession"])["Start"].min()
            stop = chunk.groupby(["Accession"])["Stop"].max()
            one_extraction = pd.concat([start, stop], axis=1).reset_index().values.tolist()

            keys = list(set(chunk["Accession"]))
            range_dict = {}
            for k in keys:
                ranges = (chunk[chunk["Accession"] == k][["Start", "Stop"]].values - start[k]).tolist()
                range_dict.update({k: ranges})

            filename = "BLASTnontarget" + str(part) + ".sequences"
            filepath = os.path.join(self.primer_qc_dir, filename)
            if not os.path.isfile(filepath):
                G.comm_log("Start writing " + filename)
                G.comm_log("Start DB extraction")
                fasta_seqs = G.run_parallel(
                        P.get_seq_fromDB, one_extraction, self.progress, args=db)
                try:
                    self.write_sequences(fasta_seqs, range_dict, filepath)
                except (KeyboardInterrupt, SystemExit):
                    G.keyexit_rollback("DB extraction", fp=filepath)
                    raise
                G.comm_log("Finished writing " + filename)
            else:
                G.comm_log("Skip writing " + filename)

    def write_sequences(self, fasta_seqs, range_dict, filename):
        recs = []
        for item in fasta_seqs:
            fulldesc = item[0]
            desc = " ".join(fulldesc.split(" ")[1:])
            acc = fulldesc.split(".")[0][1::]
            acc_desc = fulldesc.split(":")[0][1::]
            seqrange = fulldesc.split(":")[1].split(" ")[0].split("-")
            fullseq = "".join(item[1::])
            for start, stop in range_dict[acc]:
                desc_range = str(int(seqrange[0]) + start) + "_" + str(int(seqrange[0]) + stop)
                acc_id = acc_desc + "_" + desc_range
                seq = fullseq[start:stop]
                rec = SeqRecord(Seq(seq), id=acc_id, description=desc)
                recs.append(rec)

        SeqIO.write(recs, filename, "fasta")

    def main(self):
        with self.output:
            G.comm_log("Start primer blast parser")
            bp = BlastParser(self.config)
            bp.mode = "primer"
            offtarget = bp.bp_parse_results(self.primerblast_dir)
            db_seqs = self.find_primerbinding_offtarget_seqs(offtarget)
            exitstatus = self.get_primerBLAST_DBIDS(db_seqs)
            self.write_nontarget_sequences()
            if exitstatus != 1:
                self.progress.value = 1.0

            duration = time.time() - self.start
            G.comm_log(
                "> Primer blast parser time: "
                + str(timedelta(seconds=duration)).split(".")[0])
            return exitstatus