#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import shutil
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H
from scripts.configuration import RunConfig
from scripts.configuration  import PipelineStatsCollector

class Summary(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)

    def copy_primerresults(self):
        aka = H.abbrev(self.target)
        today = time.strftime("%Y_%m_%d", time.localtime())
        filepath = os.path.join(self.summ_dir, aka + "_primer.csv")
        if os.path.isfile(filepath):
            filepath = os.path.join(
                self.summ_dir, aka + "_primer" + today + ".csv")
        fp = os.path.join(self.target_dir, "primer_report.csv")
        try:
            shutil.copy(fp, filepath)
        except OSError:
            pass

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

    def run_summary(self):
        G.logger("Run: run_summary(" + self.target + ")")
        G.create_directory(self.summ_dir)
        self.copy_primerresults()
        self.copy_mostcommon_hits()
        self.copy_pangenomeinfos()
        self.copy_config()
        PipelineStatsCollector(
            self.config).write_stat("End: " + str(time.ctime()))
        self.copy_pipelinestats()
        msg = [
            "SpeciesPrimer run finished for", self.target,
            "\n", "End:", time.ctime(), "\n", "See results in",
            self.summ_dir]
        G.comm_log(" ".join(msg))