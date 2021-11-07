#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import shutil
import multiprocessing
import pandas as pd
from datetime import timedelta
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H
from basicfunctions import ParallelFunctions as P
from scripts.configuration import errors
from scripts.configuration import RunConfig
from scripts.configuration  import PipelineStatsCollector
from scripts.blastscripts import BlastPrep
from scripts.blastscripts import Blast
from scripts.blastscripts import BlastParser
from ipywidgets import widgets, VBox, HBox
from IPython.display import display
from IPython.display import HTML

# paths
script_dir = os.path.dirname(os.path.abspath(__file__))
pipe_dir, tail = os.path.split(script_dir)
dict_path = os.path.join(pipe_dir, "dictionaries")
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')

class PrimerDesign(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.p3dict = {}
        self.p3settings = os.path.join(dict_path, "p3parameters")
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def run_primer3(self):
        G.logger("Run: run_primer3(" + self.target + ")")
        input_file = os.path.join(self.primer_dir, "primer3_input")
        output_file = os.path.join(self.primer_dir, "primer3_output")
        if not os.path.isfile(output_file):
            primer3cmd = [
                "primer3_core", "-p3_settings_file=" + self.p3settings,
                "-output=" + output_file, input_file]
            try:
                G.run_subprocess(primer3cmd, True, True, True)
            except (KeyboardInterrupt, SystemExit):
                G.keyexit_rollback("primer3 run", fp=output_file)
                raise
        else:
            info = "Skip primerdesign with primer3"
            G.logger("> " + info)
            print(info)

    def read_primeroutput(self, p3_output):
        primerdatasets = []
        primerdata = []
        primererror = []
        with open(p3_output) as f:
            for line in f:
                line = line.strip()
                if "PRIMER_ERROR=Cannot open " in line:
                    msg = (
                        "primer3 cannot open the settingsfile: "
                        + self.p3settings)
                    G.logger(">" + msg)
                    errors.append([self.target, msg])
                    raise Exception

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
            G.comm_log(">" + msg)
            errors.append([self.target, msg])
            G.csv_writer(errfile, primererror)
        return primerdatasets

    def get_annotation_info(self):
        pangenome_data = os.path.join(
            self.pangenome_dir, "gene_presence_absence.csv")
        pandata = pd.read_csv(pangenome_data, index_col=0)
        annot_dict = pandata["Annotation"].to_dict()
        return annot_dict

    def get_amplicon_sequence(self, fwd_p, rev_p, template):
        rev_compl = str(Seq(rev_p).reverse_complement())
        pcr_product = (
                template[
                    template.index(str(fwd_p)):template.index(rev_compl)
                        ] + rev_compl
                    )
        return pcr_product

    def collect_primerdata(self, rawprimerdata, annot_dict, abbrev):
        results = []
        p_dict = {p.split("=")[0]: p.split("=")[1] for p in rawprimerdata}
        seq_keys = [
            'SEQUENCE_ID', 'PRIMER_PAIR_NUM_RETURNED', 'SEQUENCE_TEMPLATE']
        seq_info = [p_dict[k] for k in seq_keys]
        pairs = int(seq_info[1])
        if pairs != 0:
            for i in range(0, pairs):
                pname = abbrev + "_" + seq_info[0] + "_P" + str(i)
                primerpair_keys = [
                    'PRIMER_PAIR_' + str(i) + '_PENALTY',
                    'PRIMER_PAIR_' + str(i) + '_PRODUCT_SIZE',
                    'PRIMER_PAIR_' + str(i) + "_PRODUCT_TM"]
                positions = ["LEFT", "RIGHT"]
                if self.config.probe:
                    positions.append("INTERNAL")
                for pos in positions:
                    primer_keys = [
                        "PRIMER_" + pos + "_" + str(i) + '_SEQUENCE',
                        "PRIMER_" + pos + "_" + str(i) + '_TM',
                        "PRIMER_" + pos + "_" + str(i) + '_PENALTY']
                    primerpair_keys.extend(primer_keys)
                primer_data = [p_dict[k] for k in primerpair_keys]
                annotation = annot_dict["_".join(seq_info[0].split("_")[0:-1])]
                ampliconseq = self.get_amplicon_sequence(
                    primer_data[3], primer_data[6], seq_info[2])
                if self.config.probe is False:
                    primer_data.extend(["", "", ""])
                primer_data.extend(
                    [ampliconseq, seq_info[2], annotation, "", pname])
                results.append(primer_data)
            return results

    def save_primer3_summary(self, primerinfos):
        raw_header = [
            "Primer penalty", "Amplicon size", "Amplicon TM",
            "Primer fwd seq", "Primer fwd TM", "Primer fwd penalty",
            "Primer rev seq", "Primer rev TM", "Primer rev penalty",
            "Probe seq", "Probe TM", "Probe penalty", "Amplicon sequence",
            "Template sequence", "Annotation", "PPC", "Primer name"]
        header = [
            "Primer name", "PPC", "Primer penalty", "Annotation",
            "Primer fwd seq", "Primer fwd TM", "Primer fwd penalty",
            "Primer rev seq", "Primer rev TM", "Primer rev penalty",
            "Probe seq", "Probe TM", "Probe penalty",
            "Amplicon size", "Amplicon TM", "Amplicon sequence",
            "Template sequence"]
        df = pd.DataFrame(primerinfos, columns=raw_header)
        p3_summary = Path(self.primer_dir, "primer3_summary.csv")
        df[header].to_csv(p3_summary, index=False)

    def main(self):
        with self.output:
            abbrev = H.abbrev(self.config.target)
            G.comm_log("> Start primer design")
            self.run_primer3()
            self.progress.value = 0.8
            p3_output = os.path.join(self.primer_dir, "primer3_output")
            annot_dict = self.get_annotation_info()
            p3data = self.read_primeroutput(p3_output)
            primerinfos = []
            for primerdata in p3data:
                p_info = self.collect_primerdata(primerdata, annot_dict, abbrev)
                if p_info:
                    primerinfos.extend(p_info)

            if len(primerinfos) != 0:
                self.save_primer3_summary(primerinfos)
                self.progress.value = 1.0
                return 0
            else:
                error_msg = "No primers found"
                G.comm_log(error_msg)
                errors.append([self.target, error_msg])
                return 1


class PrimerQualityControl(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.p3_summary = Path(self.primer_dir, "primer3_summary.csv")
        self.primerlist = []

    def get_primerinfo(self, primer_dir, selected_seqs=[], mode="primer"):
        p_df = pd.read_csv(self.p3_summary, index_col=0)
        if mode == "primer":
            keys = [
                'Primer fwd seq', 'Primer rev seq']
        if mode == "mfold":
            keys = ['Amplicon sequence']
        if mode == "results":
            keys = p_df.columns
            ppc_file = os.path.join(
                self.primer_dir,
                "MFEprimer_template_results.csv")
            ppc_df = pd.read_csv(ppc_file)
            ppc_df.loc[:, "primer"] = ppc_df.loc[:, "FpID"].str.split(
                                            "_").str[0:-1].str.join("_")
            ppc_df = ppc_df.set_index("primer")
            ppc_dict = ppc_df["PPC"].to_dict()
            p_df.loc[:, "PPC"] = p_df.index.map(ppc_dict)
            return p_df.loc[selected_seqs, keys].reset_index()

        if selected_seqs == []:
            data = p_df.loc[:, keys].reset_index().to_numpy()
        else:
            data = p_df.loc[selected_seqs, keys].reset_index().to_numpy()
        return data



class PrimerBlast(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)
        self.p3_summary = Path(self.primer_dir, "primer3_summary.csv")
        self.start = time.time()

    def collect_primer(self):
        G.create_directory(self.primerblast_dir)
        p_df = pd.read_csv(self.p3_summary)
        primerpairs = len(p_df.index)
        if primerpairs == 0:
            error_msg = "Potential primer pair(s): 0"
            PipelineStatsCollector(self.config).write_stat(error_msg)
            G.comm_log("> " + error_msg)
            errors.append([self.target, error_msg])
            return 1

        info = "Potential primer pair(s): " + str(primerpairs)
        G.comm_log("> " + info)
        PipelineStatsCollector(self.config).write_stat(info)

        p_df = p_df.sort_values('Primer penalty')
        p_names = p_df["Primer name"].to_list()
        fwdnames = [p + "_F" for p in p_names]
        revnames = [p + "_R" for p in p_names]
        fwdseq = p_df['Primer fwd seq'].to_list()
        revseq = p_df['Primer rev seq'].to_list()
        recs = []
        for i, pn in enumerate(fwdnames):
            f = SeqRecord(Seq(fwdseq[i]), id=pn, description="")
            r = SeqRecord(Seq(revseq[i]), id=revnames[i], description="")
            recs.extend([f, r])
        primerfasta = Path(self.primerblast_dir, "primer_seqs.fas")
        SeqIO.write(recs, primerfasta, "fasta")
        return recs

    def primer_blast(self, primerlist):
        blastsum = os.path.join(
            self.primerblast_dir, "BLAST_results_summary.csv")
        if not os.path.isfile(blastsum):
            prep = BlastPrep(
                self.primerblast_dir, primerlist,
                "primer", self.config.blastseqs)
            use_cores = prep.run_blastprep()
            pblast = Blast(self.config, self.primerblast_dir, "primer")
            self.progress.value = pblast.progress.value
            self.output = pblast.output
            pblast.run_blast("primer", use_cores)
            self.progress.value = 1.0

    def main(self):
        with self.output:
            primer_records = self.collect_primer()
            if primer_records != 1:
                self.primer_blast(primer_records)
                return 0
            return 1


class MFEprimerQC(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.p3_summary = Path(self.primer_dir, "primer3_summary.csv")
        self.outputlayout={
            'border': '1px solid black', 'width': 'auto',
            'height': '350px','overflow': 'auto'}
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)
        self.referencegenomes = 10
        self.mfethreshold = self.config.mfethreshold
        self.offtargetdb = []
        self.primerlist = []

    def create_template_db_file(self):
        df = pd.read_csv(self.p3_summary)[["Primer name", "Template sequence"]]
        primerdata = df.drop_duplicates("Template sequence").to_numpy()
        recs = [
            SeqRecord(
                Seq(p[-1]), id="_".join(p[0].split("_")[0:-1]), description="")
                for p in primerdata]
        file_path = os.path.join(self.primer_qc_dir, "template.sequences")
        SeqIO.write(recs, file_path, "fasta")

    def create_assembly_db_file(self):
        qc_report = Path(self.reports_dir, "inputfiles_report.csv")
        qcdf = pd.read_csv(qc_report, index_col=0)
        if self.config.ignore_qc is False:
            qcdf = qcdf[~qcdf.iloc[:, -2].isna()]
        stat_dict = {
            "Complete Genome": 0, "Chromosome": 1,
            "Scaffold": 2, "Contig": 3, "": 4}
        qcdf.loc[:, "sorter"] = qcdf["AssemblyStatus"].map(stat_dict)
        qcdf = qcdf.sort_values("sorter")
        selection = qcdf.head(self.referencegenomes).index.to_list()
        files = [
            Path(self.genomic_dir, s)
            for s in os.listdir(self.genomic_dir)
            if H.accession_from_filename(s, version=False) in selection]
        recs = []
        for f in files:
            record = SeqIO.parse(f, "fasta")
            recs.extend(record)
        file_path = os.path.join(self.primer_qc_dir, "genomic.sequences")
        SeqIO.write(recs, file_path, "fasta")

    def prepare_MFEprimer_Dbs(self, dbsubprogress):
        self.create_template_db_file()
        dbsubprogress.children[0].value = 1.0
        self.create_assembly_db_file()
        dbsubprogress.children[1].value = 1.0
        self.offtargetdb = ["template.sequences"]
        for filename in os.listdir(self.primer_qc_dir):
            if (
                filename.startswith("BLASTnontarget")
                and filename.endswith(".sequences")
            ):
                self.offtargetdb.append(filename)
        self.offtargetdb.append("genomic.sequences")
        inputdbs = [
            str(Path(self.primer_qc_dir, n)) for n in self.offtargetdb]
        # parallelization try
        pool = multiprocessing.Pool()
        results = [
            pool.apply_async(
                P.index_database, args=(inputfilepath,))
                    for inputfilepath in inputdbs]
        output = [p.get() for p in results]
        for item in output:
            if item:
                errors.append([self.target, item])

        dbsubprogress.children[2].value = 1.0

    def run_MFEprimer(self, db_file, primerinfos, subprogress):
        db_name = db_file.split(".")[0]
        G.comm_log("> Start MFEprimer with " + db_name + " DB")
        mfe_results = G.run_parallel(
                P.MFEprimer_run, primerinfos, subprogress,
                args=[db_file, self.primer_qc_dir])
        header = mfe_results[0][0].split("\t")
        rawdata = [r[1:] for r in mfe_results if len(r) > 1]
        data = [i.split("\t") for item in rawdata for i in item]
        df = pd.DataFrame(data, columns=header)
        fp = os.path.join(
            self.primer_dir,
            "MFEprimer_" + db_name + "_results.csv")
        df.to_csv(fp, index=False)
        return df

    def gui_subprogress(self):
        dbsubprogress = VBox([
            widgets.FloatProgress(value=0, min=0.0, max=1.0),
            widgets.FloatProgress(value=0, min=0.0, max=1.0),
            widgets.FloatProgress(value=0, min=0.0, max=1.0)])
        dbsublabel = VBox([
            widgets.Label("Templates"),
            widgets.Label("Assemblies"),
            widgets.Label("Indexing")])
        subprogress = VBox([
            widgets.FloatProgress(value=0, min=0.0, max=1.0),
            widgets.FloatProgress(value=0, min=0.0, max=1.0),
            widgets.FloatProgress(value=0, min=0.0, max=1.0)])
        sublabel = VBox([
            widgets.Label("Templates"),
            widgets.Label("Off-targets"),
            widgets.Label("Assemblies (Secondary amplicons)")])
        remain_label = VBox([
            widgets.Label(""),
            widgets.Label(""),
            widgets.Label("")])
        display(HTML('<p> Database creation </p>'))
        display(HBox([dbsubprogress, dbsublabel]))
        display(HTML('<p> MFEprimer runs </p>'))
        display(HBox([subprogress, sublabel, remain_label]))
        return dbsubprogress, subprogress, remain_label

    def interpret_MFEprimer_results(self, result, dbname, primerinfos):
        db_name = dbname.split(".")[0]
        ppc_range = 100 - self.mfethreshold
        if result.empty:
            G.comm_log("No data found with MFEprimer DB " + dbname)

        result.loc[:, "primer"] = result.loc[
                            :, "FpID"].str.split("_").str[0:-1].str.join("_")
        result.loc[:, "PPC"] = result.loc[:, "PPC"].astype(float)

        if db_name == "template":
            count = result.loc[:, "primer"].value_counts() == 1
            primernames = [i for i in count.index if count[i]]
            targets = result[result.loc[:, "primer"].isin(primernames)]
            targets = targets.query('PPC >= @self.mfethreshold')
            primerlist = targets["primer"].to_list()

        if "nontarget" in db_name:
            offtarget = result.query('PPC >= @ppc_range')
            offtarget_primer = offtarget["primer"].to_list()
            primerlist = [
                p[0] for p in primerinfos if not p[0] in offtarget_primer]

        if db_name == "genomic":
            strongbinding = result.query('PPC >= @ppc_range')
            count = strongbinding.loc[
                        :, "primer"].value_counts() == self.referencegenomes
            primerlist = [i for i in count.index if count[i]]

        G.comm_log(
            "Continue QC with " + str(len(primerlist)) + " primer pairs", True)
        return primerlist

    def main(self):
        with self.output:
            info_msg = "Start primer quality control(" + self.target + ")"
            G.comm_log("> " + info_msg)
            os.chdir(self.primer_qc_dir)
            dbsubprogress, subprogress, remain_label = self.gui_subprogress()
            self.prepare_MFEprimer_Dbs(dbsubprogress)
            dbsubcount = [s.value for s in dbsubprogress.children]
            self.progress.value = sum(dbsubcount)/6

            for i, dbname in enumerate(self.offtargetdb):
                partial_prog = subprogress.children[i]
                primerinfos = H.get_primerinfo(
                    self.primer_dir, self.primerlist, mode="primer")
                result = self.run_MFEprimer(dbname, primerinfos, partial_prog)
                self.primerlist = self.interpret_MFEprimer_results(
                                    result, dbname, primerinfos)
                remain_label.children[i].value = (
                    "\t# Primer pairs passed: " + str(len(self.primerlist)))
                subcount = [s.value for s in subprogress.children]
                self.progress.value = (sum(dbsubcount) + sum(subcount))/6

            if len(self.primerlist) == 0:
                return 1
            primerdf = pd.DataFrame(self.primerlist)
            fp = Path(self.primer_qc_dir, "mfe_primerlist.csv")
            primerdf.to_csv(fp, index=False, header=False)
            return 0


class MFoldQC(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.p3_summary = Path(self.primer_dir, "primer3_summary.csv")
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def mfold_analysis(self, mfoldinputlist):
        results = []
        info = "Start mfold analysis of PCR products"
        G.comm_log("> " + info)
        if os.path.isdir(self.mfold_dir):
            shutil.rmtree(self.mfold_dir)
        G.create_directory(self.mfold_dir)
        tot = len(mfoldinputlist)
        for i, mfoldinput in enumerate(mfoldinputlist):
            print('\rprogress ' + str(i+1) + "/" + str(tot), end='')
            primername = mfoldinput[0]
            genename = "_".join(primername.split("_")[-3:-1])
            primerpair = mfoldinput[0].split("_")[-1]
            fullname = genename + "_" + primerpair
            directory = Path(self.mfold_dir, genename, primerpair)
            G.create_directory(directory)
            filepath = Path(directory, fullname)
            rec = SeqRecord(
                Seq(mfoldinput[1]), id=fullname, description="")
            SeqIO.write(rec, filepath, "fasta")
            os.chdir(directory)
            mfold_cmd = [
                "mfold", "SEQ=" + fullname, "NA=DNA", "T=60", "MG_CONC=0.003"]
            G.run_subprocess(mfold_cmd, False, True, False)
            result = self.parse_mfold_results(primername, fullname + ".det")
            results.extend(result)
            os.chdir(self.mfold_dir)
            self.progress.value = (i+1)/tot
        print()
        primerlist = self.interpret_mfold_results(results)
        G.comm_log(
            "Continue QC with " + str(len(primerlist)) + " primer pairs", True)
        return primerlist

    def parse_mfold_results(self, primername, filename):
        results = []
        count = 0
        with open(filename, "r", errors="ignore") as f:
            for i, line in enumerate(f):
                line = line.strip()
                if "Structure" in line:
                    get_value = i + 3
                    count += 1
                if count > 0:
                    if i == get_value:
                        params = ["dG", "dH", "dS", "Tm", '℃']
                        splitlst = [
                            v.split(" ") for i, v in enumerate(
                                line.split(" = ")) if i != 0]
                        values = [
                            val for sub in splitlst for val in sub if (
                                val != '' and val not in params)]

                        results.append([primername, count] + values)
        return results

    def interpret_mfold_results(self, results):
        header=["Primer name", "Structure", "dG", "dH", "dS", "Tm"]
        df = pd.DataFrame(results, columns=header)
        df.loc[:, ["dG", "dH", "dS", "Tm"]]  = df.loc[
                                    :, ["dG", "dH", "dS", "Tm"]].astype(float)
        mfoldsummary = Path(self.primer_dir, "mfold_results.csv")
        df.to_csv(mfoldsummary, index=False)
        failed = df.query('dG <= @self.config.mfold')["Primer name"].unique()
        passed = df.query('dG > @self.config.mfold')["Primer name"].unique()
        primerlist = list(set(passed) - set(failed))
        return primerlist

    def main(self):
        with self.output:
            fp = Path(self.primer_qc_dir, "mfe_primerlist.csv")
            df = pd.read_csv(fp, header=None)
            passed_primerlist = df.iloc[:, 0].to_list()
            mfoldinput = H.get_primerinfo(
                self.primer_dir, passed_primerlist, mode="mfold")
            passed_primerlist = self.mfold_analysis(mfoldinput)
            if len(passed_primerlist) == 0:
                return 1
            primerdf = pd.DataFrame(passed_primerlist)
            fp = Path(self.primer_qc_dir, "mfold_primerlist.csv")
            primerdf.to_csv(fp, index=False, header=False)
            return 0

class MPprimerDimerQC(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        self.progress = widgets.FloatProgress(value=0, min=0.0, max=1.0)
        self.output = widgets.Output(layout=self.outputlayout)

    def run_MPprimer_dimer_check(self, primer):
        p_name = primer[0]
        fwd = SeqRecord(
                Seq(primer[1]),
                id=p_name + "_F",
                description="")
        fwd2 = SeqRecord(
                Seq(primer[1]),
                id=p_name + "_F2",
                description="")
        rev = SeqRecord(
                Seq(primer[2]),
                id=p_name + "_R",
                description="")
        rev2 = SeqRecord(
                Seq(primer[2]),
                id=p_name + "_R2",
                description="")

        names = [p_name, p_name + "_F", p_name + "_R"]
        seqs = [[fwd, rev], [fwd, fwd2], [rev, rev2]]
        results = []
        for i, fn in enumerate(names):
            fp = Path(self.dimercheck_dir, fn)
            SeqIO.write(seqs[i], fp, "fasta")
            output = str(fp) + "_dimer"
            dimer_cmd = ["MPprimer_dimer_check.pl", "-f", str(fp), "-d", "3"]
            output = G.read_shelloutput(dimer_cmd)
            data = [i.split("\t") for i in output]
            results.extend(data)

        return results

    def check_primerdimers(self, dimer_input):
        G.comm_log("> Start analysis of primer dimers")
        if os.path.isdir(self.dimercheck_dir):
            shutil.rmtree(self.dimercheck_dir)
        G.create_directory(self.dimercheck_dir)

        results = []
        tot = len(dimer_input)
        for i, primerdata in enumerate(dimer_input):
            print('\rprogress ' + str(i+1) + "/" + str(tot), end='')
            output = self.run_MPprimer_dimer_check(primerdata)
            results.extend(output)
            self.progress.value = (i+1)/tot
        print()
        primerlist = self.interpret_primerdimer_results(results)
        info = "Primer pairs left after primer QC: " + str(len(primerlist))
        G.comm_log("> " + info, True)
        PipelineStatsCollector(self.config).write_stat(info)
        return primerlist

    def interpret_primerdimer_results(self, results):
        df = pd.DataFrame(results, columns=["Primer1", "Primer2", "dG"])
        df.loc[:, "dG"]  = df.loc[:, "dG"].astype(float)
        dimersummary = Path(self.primer_dir, "dimer_results.csv")
        df.to_csv(dimersummary, index=False)
        df.loc[:, "Primer name"] = df.loc[
            :, "Primer1"].str.split("_").str[0:-1].str.join("_")
        failed = df.query('dG <= @self.config.mpprimer')["Primer name"].unique()
        passed = df.query('dG > @self.config.mpprimer')["Primer name"].unique()
        primerlist = list(set(passed) - set(failed))
        return primerlist

    def main(self):
        with self.output:
            fp = Path(self.primer_qc_dir, "mfold_primerlist.csv")
            df = pd.read_csv(fp, header=None)
            passed_primerlist = df.iloc[:, 0].to_list()
            dimer_input = H.get_primerinfo(
                self.primer_dir, passed_primerlist, mode="primer")
            passed_primerlist = self.check_primerdimers(dimer_input)
            results_df = H.get_primerinfo(
                self.primer_dir, passed_primerlist, mode="results")

            if results_df.empty:
                error_msg = "No compatible primers found"
                G.comm_log("> " + error_msg)
                errors.append([self.target, error_msg])
                return 1

            results_df = results_df.sort_values("PPC", ascending=False)
            fp = Path(self.reports_dir, "primer_report.csv")
            results_df.to_csv(fp, index=False)
            return 0


