#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import time
import logging
import signal
import json
from basicfunctions import GeneralFunctions as G
from basicfunctions import HelperFunctions as H
from scripts.configuration import Config
from scripts.configuration import CLIconf
from scripts.configuration import PipelineStatsCollector
from scripts.datacollection import GenomeDownload
from scripts.datacollection import Annotation
from scripts.qualitycontrol import QualityControl
from scripts.coregenes import PangenomeAnalysis
from scripts.coregenes import CoreGenes
from scripts.coregenes import CoreGeneSequences
from scripts.primerdesign import PrimerDesign
from scripts.primerdesign import PrimerQualityControl
from scripts.summary import Summary


# paths
pipe_dir = os.path.dirname(os.path.abspath(__file__))
dict_path = os.path.join(pipe_dir, "dictionaries")
tmp_db_path = os.path.join(pipe_dir, 'tmp_config.json')
errors = []

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
        assemblylevel, skip_tree, nolist,
        offline, ignore_qc, mfethreshold, customdb,
        blastseqs, probe, virus, genbank,
        evalue, nuc_identity, runmode, strains, nontargetlist)

    return config


def run_pipeline_for_target(target, config):
    print("\nStart searching primer for " + target)
    G.logger("> Start searching primer for " + target)
    PipelineStatsCollector(config).write_stat(
        target + " pipeline statistics:")
    PipelineStatsCollector(config).write_stat(
        "Start: " + str(time.ctime()))
    newconfig = DataCollection(config).main()
    newconfig.save_config()
    if newconfig != 0:
        config = newconfig
    if config.virus is True:
        config.qc_gene = []

    QualityControl(config).main()
    PangenomeAnalysis(config).run_pangenome_analysis()

    if "strain" in config.runmode:
        import strainprimer
        msg = "Start searching for strain specific primers"
        print(msg)
        G.logger(msg)
        strainprimer.main(config)

    if "species" in config.runmode:
        G.comm_log("Start searching for species specific primers")
        CoreGenes(config).run_CoreGenes()
        CoreGeneSequences(config).run_coregeneanalysis()
        PrimerDesign(config).run_primerdesign()
        PrimerQualityControl(config).run_primer_qc()
        Summary(config).run_summary()


def get_configuration_from_args(target, args):
    if args.nolist:
        nontargetlist = []
    else:
        nontargetlist = H.create_non_target_list(target)

    config = CLIconf(
        args.minsize, args.maxsize, args.mpprimer, args.exception,
        target, args.path, args.intermediate,
        args.qc_gene, args.mfold, args.skip_download,
        args.assemblylevel,
        args.skip_tree, args.nolist, args.offline,
        args.ignore_qc, args.mfethreshold, args.customdb,
        args.blastseqs, args.probe, args.virus, args.genbank,
        args.evalue, args.nuc_identity, args.runmode, args.strains,
        nontargetlist)

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
            PipelineStatsCollector(config).write_stat(
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
