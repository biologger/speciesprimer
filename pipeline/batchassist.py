#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import sys
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H

# paths
pipe_dir = os.path.dirname(os.path.abspath(__file__))

class Input:
    def __init__(self):
        self.target_list = []
        self.config_dict = {}
        self.config_paths = []

    def gui_runner(self, mode, data):
        if mode == "continue":
            path, targets = data

            config_dict = Output().run_gui_output(targets, path)
            return config_dict

        elif mode == "new":
            self.config_dict = data
            for target in self.config_dict.keys():
                self.target_list.append(target)
                self.write_config_file(target)
            return self.config_dict

    def initiate(self):
        mode = input(
            "Create new config files or start pipeline with previously "
            "generated files?\ntype (n)ew or (s)tart:\n")
        if (mode.lower() == "n" or mode.lower() == "new"):
            config_dict = self.main()
            return config_dict
        elif (mode.lower() == "s" or mode.lower() == "start"):
            config_dict = Output().run_output()
            return config_dict
        else:
            print("no valid input, type: (n)ew or (s)tart")
            sys.exit()

    def targetinput(self):
        targets = input(
            "Please specify target species (comma separated) "
            "or type help for format examples: \n")
        return targets

    def helpmessage(self, targets):
        print(
            "Accepted input formats are e.g. Lactobacillus_helveticus, "
            "Lactobacillus delbrueckii, Lactococcus lactis subsp lactis, "
            "Lactococcus_lactis_subsp_cremoris")
        targets = input(
            "Please specify target species (comma separated) "
            "or type help for format examples: \n")
        return targets

    def parse_targets(self, targets):
        if isinstance(targets, list):
            items = targets
        else:
            items = list(targets.split(","))
        for i, item in enumerate(items):
            x = item.strip()
            x = x.capitalize()
            self.target_list.insert(i, x)

        for i, item in enumerate(self.target_list):
            target = item.strip()
            if " " in target:
                self.target_list.remove(target)
                new_target = "_".join(target.split(" "))
                self.target_list.insert(i, new_target)
        for target in self.target_list:
            self.config_dict.update({target: {"target": target}})
        return self.target_list

    def get_path(self, target, index, listlen):
        if "path" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            inpath = input(
                "Please specify a path to use as the working directory "
                "or hit return to use the current working directory:\n")
            if inpath:
                path = inpath
            else:
                path = os.getcwd()
            if index == 0:
                if not self.value_for_all("path", path, listlen):
                    self.config_dict[target].update({"path": path})
            else:
                self.config_dict[target].update({"path": path})

            print("path", path)

    def get_skip_tree(self, target, index, listlen):
        if "skip_tree" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            skip_tree = input(
                "Skip the core gene alignment and tree for visualization and "
                "troubleshooting.\n(y)es/(n)o, default=(n)\n> ")
            if skip_tree.lower() == ("y" or "yes"):
                skip_tree = True
            else:
                skip_tree = False
            if index == 0:
                if not self.value_for_all("skip_tree", skip_tree, listlen):
                    self.config_dict[target].update({"skip_tree": skip_tree})
            else:
                self.config_dict[target].update({"skip_tree": skip_tree})
            print("skip_tree", skip_tree)

    def work_offline(self, target, index, listlen):
        if "offline" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            limit = input(
                "Work offline with local genome assemblies?"
                "\n(y)es/(n)o, default=(n)\n> ")
            if limit.lower() == ("y" or "yes"):
                offline = True
                skip_download = True
                assembly_list = ["offline"]

            else:
                offline = False
                skip_download = False
            if index == 0:
                if not self.value_for_all("offline", offline, listlen):
                    self.config_dict[target].update({"offline": offline})
                    if skip_download:
                        self.config_dict[target].update(
                            {"skip_download": skip_download})
                        self.config_dict[target].update(
                            {"assemblylevel": assembly_list})

                else:
                    if skip_download:
                        for target in self.target_list:
                            self.config_dict[target].update(
                                {"skip_download": skip_download})
                            self.config_dict[target].update(
                                {"assemblylevel": assembly_list})

            else:
                self.config_dict[target].update({"offline": offline})
                if skip_download:
                    self.config_dict[target].update(
                        {"skip_download": skip_download})
                    self.config_dict[target].update(
                        {"assemblylevel": assembly_list})

            print("offline", offline)

    def get_skip_download(self, target, index, listlen):
        if "skip_download" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            skip_download = input(
                "Skip the download of Genomes from NCBI?"
                "\n(y)es/(n)o, default=(n)\n> ")
            if skip_download.lower() == ("y" or "yes"):
                skip_download = True
                assembly_list = ["offline"]
            else:
                skip_download = False

            if index == 0:
                if not self.value_for_all(
                    "skip_download", skip_download, listlen
                ):
                    self.config_dict[target].update(
                        {"skip_download": skip_download})
                    if skip_download:
                        self.config_dict[target].update(
                            {"assemblylevel": assembly_list})
                else:
                    if skip_download:
                        for target in self.target_list:
                            self.config_dict[target].update(
                                {"assemblylevel": assembly_list})
            else:
                self.config_dict[target].update(
                    {"skip_download": skip_download})
                if skip_download:
                    self.config_dict[target].update(
                        {"assemblylevel": assembly_list})
            print("skip_download", skip_download)

    def assemblylevel(self, target, index, listlen):
        options = [
            'complete', 'chromosome', 'scaffold', 'contig', 'all', "offline"
        ]
        assembly_list = []
        valid = []
        if "assemblylevel" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            limit = input(
                "Limit downloads of Genomes to assembly status\n"
                "options: complete, chromosome, scaffold, contig, "
                "all (comma separated), default=all\n> ")
            if limit:
                items = list(limit.split(","))
                for i, item in enumerate(items):
                    x = item.strip()
                    assembly_list.insert(i, x)
                for status in assembly_list:
                    if status in options:
                        valid.append(status)
                    else:
                        print("\nNo valid assembly status:")
                        print(status)
                        return self.assemblylevel(target, index, listlen)
            else:
                valid = ['all']

            if index == 0:
                if not self.value_for_all("assemblylevel", valid, listlen):
                    self.config_dict[target].update({"assemblylevel": valid})
            else:
                self.config_dict[target].update({"assemblylevel": valid})
            print("assemblylevel:")
            print(valid)

    def get_customdb(self, target, index, listlen):
        if "customdb" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            customdb = input(
                "Do you want to use a custom database for blastn?\n"
                "Specifiy the absolute filepath of the custom database "
                "(e.g. '/home/blastdb/nontarget.fasta') or hit return"
                " to skip, default=None\n> ")
            if type(customdb) == str and len(customdb) > 0:
                if os.path.isfile(customdb + ".nsq"):
                    customdb = customdb
                else:
                    print("No BLAST DB was found " + customdb)
                    return self.get_customdb(target, index, listlen)
            else:
                customdb = None

            if index == 0:
                if not self.value_for_all("customdb", customdb, listlen):
                    self.config_dict[target].update(
                        {"customdb": customdb})
            else:
                self.config_dict[target].update({"customdb": customdb})
            print("customdb", customdb)

    def get_blastseqs(self, target, index, listlen):
        options = [100, 500, 1000, 2000, 5000]
        if "blastseqs" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            maxseqs = input(
                "Set the number of sequences per BLAST search. "
                "Decrease the number of sequences if BLAST slows down "
                "due to low memory."
                "\noptions = [100, 500, 1000, 2000, 5000]"
                ", default=1000\n> ")
            if maxseqs:
                if type(eval(maxseqs)) == int:
                    if int(maxseqs) in options:
                        valid = int(maxseqs)
                    else:
                        print("\nNo valid value for blastseqs:")
                        return self.get_blastseqs(target, index, listlen)
                else:
                    print("\nNo valid value for blastseqs:")
                    return self.get_blastseqs(target, index, listlen)
            else:
                valid = 1000

            if index == 0:
                if not self.value_for_all("blastseqs", valid, listlen):
                    self.config_dict[target].update({"blastseqs": valid})
            else:
                self.config_dict[target].update({"blastseqs": valid})

            print("blastseqs", valid)

    def get_qc_genes(self, target, index, listlen):
        options = ['rRNA', 'tuf', 'recA', 'dnaK', 'pheS']
        qc_list = []
        valid = []
        if "qc_gene" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            qc_genes = input(
                "Gene(s) (comma separated) for BLAST in the initial quality "
                "control step.\noptions: rRNA, tuf, recA, dnaK, pheS"
                ", default=rRNA\n> ")
            if qc_genes:
                items = list(qc_genes.split(","))
                for i, item in enumerate(items):
                    x = item.strip()
                    qc_list.insert(i, x)
                for gene in qc_list:
                    if gene in options:
                        valid.append(gene)
                    else:
                        print("\nNo valid gene:")
                        print(gene)
                        return self.get_qc_genes(target, index, listlen)
            else:
                valid = ["rRNA"]

            if index == 0:
                if not self.value_for_all("qc_gene", valid, listlen):
                    self.config_dict[target].update({"qc_gene": valid})
            else:
                self.config_dict[target].update({"qc_gene": valid})

            print("qc_gene:")
            print(valid)

    def get_exception(self, target, index, listlen):
        if "exception" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            exception = input(
                "Primer binding to this non-target species is tolerated.\n"
                "Provide a species name or hit return to skip:\n> ")
            if exception:
                exceptions = list(exception.split(","))
                if len(exceptions) > 1:
                    print("only one species allowed")
                    self.get_exception(target)
                else:
                    if " " in exception:
                        exception = "_".join(exception.split(" "))
            else:
                exception = None

            if index == 0:
                if not self.value_for_all("exception", exception, listlen):
                    self.config_dict[target].update({"exception": exception})
            else:
                self.config_dict[target].update({"exception": exception})
            print("exception", exception)

    def get_minsize(self, target, index, listlen):
        if "minsize" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            minsize = input("Minimal Amplicon size\ndefault=70\n> ")
            if minsize:
                if not type(eval(minsize)) == int:
                    print("No valid input")
                    return self.get_minsize(target, index, listlen)
            else:
                minsize = 70

            if index == 0:
                if not self.value_for_all("minsize", int(minsize), listlen):
                    self.config_dict[target].update({"minsize": int(minsize)})
            else:
                self.config_dict[target].update({"minsize": int(minsize)})
            print("minsize", minsize)

    def get_maxsize(self, target, index, listlen):
        if "maxsize" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            maxsize = input("Maximal amplicon size\ndefault=200\n> ")
            if maxsize:
                if not type(eval(maxsize)) == int:
                    print("No valid input")
                    return self.get_maxsize(target, index, listlen)
            else:
                maxsize = 200
            if index == 0:
                if not self.value_for_all("maxsize", int(maxsize), listlen):
                    self.config_dict[target].update({"maxsize": int(maxsize)})
            else:
                self.config_dict[target].update({"maxsize": int(maxsize)})

            print("maxsize", maxsize)

    def get_designprobe(self, target, index, listlen):
        if "probe" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            designprobe = input(
                "Do you want primer3 to design an internal probe?"
                "[Experimental!]"
                "\n(y)es/(n)o, default=(n)\n> ")
            if designprobe.lower() == ("y" or "yes"):
                designprobe = True
            else:
                designprobe = False
            if index == 0:
                if not self.value_for_all("probe", designprobe, listlen):
                    self.config_dict[target].update({"probe": designprobe})
            else:
                self.config_dict[target].update({"probe": designprobe})

            print("probe", designprobe)

    def get_mfold(self, target, index, listlen):
        if "mfold" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            mfold_th = input(
                "Delta G threshold for secondary structures in PCR products"
                " at 60 degree Celsius calculated by mfold\ndefault=-3.0)\n> ")
            if mfold_th:
                if not type(eval(mfold_th)) == float:
                    print("No valid input (float)")
                    return self.get_mfold(target, index, listlen)
            else:
                mfold_th = -3.0

            if index == 0:
                if not self.value_for_all("mfold", float(mfold_th), listlen):
                    self.config_dict[target].update({"mfold": float(mfold_th)})
            else:
                self.config_dict[target].update({"mfold": float(mfold_th)})
            print("mfold", mfold_th)

    def get_mpprimer(self, target, index, listlen):
        if "mpprimer" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            mpprimer_th = input(
                "Mpprimer threshold for delta G values calculated by "
                "MPprimer_dimer_check.pl\ndefault=-3.5\n> ")
            if mpprimer_th:
                if not type(eval(mpprimer_th)) == float:
                    print("No valid input (float)")
                    return self.get_mpprimer(target, index, listlen)
            else:
                mpprimer_th = -3.5
            if index == 0:
                if not self.value_for_all(
                    "mpprimer", float(mpprimer_th), listlen
                ):
                    self.config_dict[target].update(
                        {"mpprimer": float(mpprimer_th)})
            else:
                self.config_dict[target].update(
                    {"mpprimer": float(mpprimer_th)})
            print("mpprimer", mpprimer_th)

    def get_mfeprimer_threshold(self, target, index, listlen):
        if "mfethreshold" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            mfeprimer_th = input(
                "MFEprimer threshold for nontarget sequence PPC "
                "higher values mean more stringent selection.\ndefault=90\n> ")
            if mfeprimer_th:
                if not type(eval(mfeprimer_th)) == int:
                    print("No valid input (float)")
                    return self.mfeprimer_threshold(target, index, listlen)
            else:
                mfeprimer_th = 90
            if index == 0:
                if not self.value_for_all(
                    "mfethreshold", int(mfeprimer_th), listlen
                ):
                    self.config_dict[target].update(
                        {"mfethreshold": int(mfeprimer_th)})
            else:
                self.config_dict[target].update(
                    {"mfethreshold": int(mfeprimer_th)})
            print("mfeprimer threshold", mfeprimer_th)

    def no_singleton(self, target, index, listlen):
        if "singleton" not in self.config_dict[target].keys():
            self.config_dict[target].update({"singleton": False})

    def ignore_qc(self, target, index, listlen):
        if "ignore_qc" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            ignore_qc = input(
                "Do you want to include genomes that did not"
                " pass quality control?\ndefault=(n)\n> ")
            if ignore_qc.lower() == ("y" or "yes"):
                ignore_qc = True
            else:
                ignore_qc = False
            if index == 0:
                if not self.value_for_all("ignore_qc", ignore_qc, listlen):
                    self.config_dict[target].update({"ignore_qc": ignore_qc})
            else:
                self.config_dict[target].update({"ignore_qc": ignore_qc})

            print("ignore_qc", ignore_qc)

    def use_blastdbv5(self, target, index, listlen):
        if "blastdbv5" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            blastdbv5 = input(
                "Do you have the Version 5 of the BLAST DB? \ndefault=(n)\n> ")
            if blastdbv5.lower() == ("y" or "yes"):
                blastdbv5 = True
            else:
                blastdbv5 = False
            if index == 0:
                if not self.value_for_all("blastdbv5", blastdbv5, listlen):
                    self.config_dict[target].update({"blastdbv5": blastdbv5})
            else:
                self.config_dict[target].update({"blastdbv5": blastdbv5})

            print("blastdbv5", blastdbv5)

    def get_intermediate(self, target, index, listlen):
        if "intermediate" not in self.config_dict[target].keys():
            print("\n" + target + ":")
            intermediate = input(
                "Do you want to keep intermediate files?\ndefault=(n)\n> ")

            if intermediate.lower() == ("y" or "yes"):
                intermediate = True
            else:
                intermediate = False
            if index == 0:
                if not self.value_for_all(
                    "intermediate", intermediate, listlen
                ):
                    self.config_dict[target].update(
                        {"intermediate": intermediate})
            else:
                self.config_dict[target].update({"intermediate": intermediate})

            print("intermediate", intermediate)

    def get_nolist(self, target, index, listlen):
        if "nolist" not in self.config_dict[target].keys():
#            self.config_dict[target].update({"nolist": False})
            print("\n" + target + ":")
            nolist = input(
                "Do you want to perform specificity check without the "
                "(non-target) species list (for all sequences in the DB)?"
                "\nNot recommended for nt DB! May be used with a custom DB"
                "\ndefault=(n)\n> ")

            if nolist.lower() == ("y" or "yes"):
                nolist = True
            else:
                nolist = False
            if index == 0:
                if not self.value_for_all(
                    "nolist", nolist, listlen
                ):
                    self.config_dict[target].update(
                        {"nolist": nolist})
            else:
                self.config_dict[target].update({"nolist": nolist})

            print("nolist", nolist)

    def value_for_all(self, key, value, listlen):
        if listlen > 1:
            forall = input(
                "Use this value for all targets?\n(y)es/(n)o, default=(y)\n> ")
            if forall.lower() == "n":
                return False
            else:
                for target in self.target_list:
                    self.config_dict[target].update({key: value})
                return True
        else:
            return False

    def write_config_file(self, target):
        path = self.config_dict[target]["path"]
        target_dir = os.path.join(path, target)
        genomic_dir = os.path.join(target_dir, "genomic_fna")
        config_dir = os.path.join(target_dir, "config")
        config_file = os.path.join(config_dir, "config.json")
        G.create_directory(target_dir)
        G.create_directory(config_dir)
        G.create_directory(genomic_dir)
        self.config_paths.append(config_file)
        with open(config_file, "w") as f:
            f.write(json.dumps(self.config_dict[target]))

    def main(self):
        targets = self.targetinput()
        while (targets == "help" or targets == "" or targets is None):
            targets = self.helpmessage(targets)

        self.parse_targets(targets)
        listlen = len(self.target_list)
        # get input
        for i, target in enumerate(self.target_list):
            self.get_path(target, i, listlen)
            self.work_offline(target, i, listlen)
            self.get_customdb(target, i, listlen)
            self.get_blastseqs(target, i, listlen)
            self.use_blastdbv5(target, i, listlen)
            self.get_skip_download(target, i, listlen)
            self.assemblylevel(target, i, listlen)
            self.get_qc_genes(target, i, listlen)
            self.get_exception(target, i, listlen)
            self.get_skip_tree(target, i, listlen)
            self.get_minsize(target, i, listlen)
            self.get_maxsize(target, i, listlen)
            self.get_designprobe(target, i, listlen)
            self.get_mfold(target, i, listlen)
            self.get_mpprimer(target, i, listlen)
            self.get_mfeprimer_threshold(target, i, listlen)
            self.no_singleton(target, i, listlen)
            self.get_intermediate(target, i, listlen)
            self.ignore_qc(target, i, listlen)
            self.get_nolist(target, i, listlen)

        for target in self.target_list:
            self.write_config_file(target)

        return self.config_dict


class Output:
    def __init__(self):
        self.targets = []
        self.config_paths = []
        self.config_dict = {}
        self.default_dict = {
            "minsize": 70,
            "maxsize": 200,
            "mpprimer": -3.5,
            "exception": None,
            "path": os.getcwd(),
            "intermediate": False,
            "qc_gene": ["rRNA"],
            "mfold": -3.0,
            "skip_download": False,
            "assemblylevel": ["all"],
            "skip_tree": False,
            "nolist": False,
            "offline": False,
            "ignore_qc": False,
            "mfethreshold": 90,
            "customdb": None,
            "blastseqs": 1000,
            "probe": False,
            "blastdbv5": False}

    def get_path(self):
        path = input(
                "Please specify a path to use as the working directory "
                "or hit return to use the current working directory:\n")
        if path:
            pass
        else:
            path = os.getcwd()
        return path

    def search_configfiles(self, path):
        print("Search in " + path)
        for root, dirs, files in os.walk(path):
            for file_name in files:
                if "config.json" == file_name:
                    file_path = os.path.join(root, file_name)
                    print("found:", file_path)
                    G.logger("> Found config file: " + file_path)
                    target = "/".join(root.split("/")[-2:-1])
                    abbr = H.abbrev(target)
                    primer_csv = os.path.join(
                        path, "Summary", target, abbr + "_primer.csv")
                    if not os.path.isfile(primer_csv):
                        self.targets.append(target)
                        self.config_paths.append(file_path)
                    else:
                        G.logger("> Skip run found results for " + target)

        return self.targets, self.config_paths

    def read_config(self, target, config_path):
        self.config_dict.update({target: {}})
        with open(config_path, "r") as f:
            for line in f:
                self.config_dict[target] = json.loads(line)

        # new 20.11.2018
        # make old config files compatible with the updated pipeline
        for key in self.default_dict:
            try:
                self.config_dict[target][key]
            except KeyError:
                self.config_dict[target].update({key: self.default_dict[key]})
                info1 = (
                    "Warning! new option " + str(key)
                    + " was not found in config file")
                info2 = (
                    "Used default value: " + str(key) + " = "
                    + str(self.default_dict[key]))
                print("\n" + info1 + "\n")
                print("\n" + info2 + "\n")
                G.logger(info1)
                G.logger(info2)

    def run_output(self):
        targets = input(
            "Search for config files for (a)ll or (s)elect targets:\n")
        path = self.get_path()
        if targets.lower() == ("a" or "all"):
            self.search_configfiles(path)
        elif targets.lower() == ("s" or "select"):
            targets = Input().targetinput()
            while (targets == "help" or targets == "" or targets is None):
                targets = Input().helpmessage(targets)

            self.targets = Input().parse_targets(targets)
            for target in self.targets:
                configpath = os.path.join(
                        path, target, "config", "config.json")
                if os.path.isfile(configpath):
                    self.config_paths.append(configpath)
                else:
                    info = (
                        "No configuration files found for "
                        + target + " in specified path")
                    print(info)
                    G.logger(info)
                    return self.run_output()
        else:
            print("no valid input: choose (a)ll or (s)elect")
            sys.exit()
        for index, target in enumerate(self.targets):
            config_path = self.config_paths[index]
            self.read_config(target, config_path)

        return self.config_dict

    def run_gui_output(self, targets, path):
        if targets is None:
            self.search_configfiles(path)
        else:
            self.targets = Input().parse_targets(targets)
            for target in self.targets:
                configpath = os.path.join(
                        path, target, "config", "config.json")
                if os.path.isfile(configpath):
                    self.config_paths.append(configpath)
                else:
                    self.config_paths.append('None')
                    info = (
                        "No configuration files found for "
                        + target + " in specified path")
                    print(info)
                    G.logger(info)

        for index, target in enumerate(self.targets):
            config_path = self.config_paths[index]
            if not config_path == 'None':
                self.read_config(target, config_path)
        return self.config_dict


if __name__ == "__main__":
    Input().initiate()
