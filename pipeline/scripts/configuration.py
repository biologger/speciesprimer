#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import time
from basicfunctions import GeneralFunctions as G

errors = []

class Config:
    def __init__(self, mode="man", config_dict={}):
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
    ## attention!!!! position of nontargetlist has changed
    def __init__(
            self, minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold,
            skip_download, assemblylevel,
            skip_tree, nolist, offline, ignore_qc, mfethreshold,
            customdb, blastseqs, probe, virus, genbank,
            evalue, nuc_identity, runmode, strains,
            nontargetlist):
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
        self.gui = False
        #self.save_config()
    def set_gui(self):
        self.gui = True

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


class RunConfig():
    def __init__(self, configuration):
        self.config = configuration
        self.target = configuration.target
        self.target_dir = os.path.join(self.config.path, self.target)
        self.config_dir = os.path.join(self.target_dir, "config")
        self.reports_dir = os.path.join(self.target_dir, "reports")
        self.genomedata_dir = os.path.join(self.target_dir, "genomedata")
        self.genomic_dir = os.path.join(self.genomedata_dir, "genomic_fna")
        self.ex_dir = os.path.join(self.genomedata_dir, "excluded_genomes")
        self.annotation_dir = os.path.join(self.genomedata_dir, "annotations")
        self.fna_dir = os.path.join(self.genomedata_dir, "fna_files")
        self.ffn_dir = os.path.join(self.genomedata_dir, "ffn_files")
        self.gff_dir = os.path.join(self.genomedata_dir, "gff_files")
        self.pangenome_dir = os.path.join(self.target_dir, "pangenome")
        self.coregene_dir = os.path.join(self.target_dir, "coregenes")
        self.fasta_dir = os.path.join(self.coregene_dir, "fasta")
        self.singlecopy = os.path.join(self.coregene_dir, "singlecopy_genes.csv")
        self.ffn_seqs = os.path.join(self.coregene_dir, "ffn_sequences.csv")
        self.alignments_dir = os.path.join(self.coregene_dir, "alignments")
        self.consensus_dir = os.path.join(self.coregene_dir, "consensus")
        self.blast_dir = os.path.join(self.coregene_dir, "blast")
        self.primer_dir = os.path.join(self.target_dir, "primer")
        self.primerblast_dir = os.path.join(self.primer_dir, "primerblast")
        self.primer_qc_dir = os.path.join(self.primer_dir, "primer_QC")
        self.mfold_dir = os.path.join(self.primer_dir, "mfold")
        self.dimercheck_dir = os.path.join(self.primer_dir, "dimercheck")
        self.summ_dir = os.path.join(self.config.path, "Summary", self.target)
        self.contiglimit = 500
        self.outputlayout={
            'border': '1px solid black', 'width': 'auto',
            'height': '200px','overflow': 'auto'}


class PipelineStatsCollector(RunConfig):
    def __init__(self, configuration):
        RunConfig.__init__(self, configuration)
        today = time.strftime("%Y_%m_%d", time.localtime())
        self.statfile = os.path.join(
            self.target_dir, "pipeline_stats_" + today + ".txt")

    def write_stat(self, info):
        if os.path.isfile(self.statfile):
            with open(self.statfile, "a") as f:
                f.write(info + "\n")
        else:
            with open(self.statfile, "w") as f:
                f.write(info + "\n")

