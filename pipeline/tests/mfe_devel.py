#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 12:35:52 2019

@author: ags-bioinfo
"""

import os
import csv
import sys
from Bio import SeqIO
import time
from datetime import timedelta
import shutil
import json

BASE_PATH = os.path.abspath(__file__).split("tests")[0]
new_path = os.path.join(BASE_PATH, "bin")
sys.path.append(new_path)
testfiles_dir = os.path.join(BASE_PATH, "tests", "testfiles")
from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G


confargs = {
    "ignore_qc": False, "mfethreshold": 80, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "home", "primerdesign", "test"),
    "probe": False, "exception": None, "minsize": 70, "skip_download": True,
    "customdb": None, "assemblylevel": ["all"], "qc_gene": ["rRNA"],
    "blastdbv5": True, "intermediate": True, "nontargetlist": []}

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def config():
    from speciesprimer import CLIconf
    args = AttrDict(confargs)
    nontargetlist = H.create_non_target_list(args.target)
    config = CLIconf(
            args.minsize, args.maxsize, args.mpprimer, args.exception,
            args.target, args.path, args.intermediate,
            args.qc_gene, args.mfold, args.skip_download,
            args.assemblylevel, nontargetlist,
            args.skip_tree, args.nolist, args.offline,
            args.ignore_qc, args.mfethreshold, args.customdb,
            args.blastseqs, args.probe, args.blastdbv5)

    config.save_config()

    return config

config = config()

def test_PrimerQualityControl_specificitycheck(config):
    from speciesprimer import PrimerQualityControl
    from speciesprimer import BlastPrep
    from speciesprimer import Blast

    tmpdir = os.path.join(BASE_PATH, "tests", "tmp")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    config.customdb = os.path.join(tmpdir, "primer_customdb.fas")
    config.blastdbv5 = False

#    def dbinputfiles():
#        filenames = [
#            "GCF_004088235v1_20191001.fna",
#            "GCF_002224565.1_ASM222456v1_genomic.fna"]
#        dbfile = os.path.join(testfiles_dir, "primer_customdb.fas")
#        with open(dbfile, "w") as f:
#            for filename in filenames:
#                filepath = os.path.join(testfiles_dir, filename)
#                records = SeqIO.parse(filepath, "fasta")
#                for record in records:
#                    if record.id == record.description:
#                        description = (
#                            record.id + " Lactobacillus curvatus strain SRCM103465")
#                        record.description = description
#                    SeqIO.write(record, f, "fasta")
#        return dbfile
#
#    def create_customblastdb(config, infile):
#        cmd = [
#            "makeblastdb", "-in", infile, "-parse_seqids", "-title",
#            "mockconservedDB", "-dbtype", "nucl", "-out", config.customdb]
#        G.run_subprocess(
#            cmd, printcmd=False, logcmd=False, log=False, printoption=False)
#
#    def test_collect_primer(config):
#        pqc = PrimerQualityControl(config, {})
#        exitstat = pqc.collect_primer() # returns 1 of len(primerlist) == 0, otherwise return 0
#        assert exitstat == 1
#        reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
#        with open(reffile) as f:
#            for line in f:
#                primer3dict = json.loads(line)
#        pqc = PrimerQualityControl(config, primer3dict)
#
#
#        exitstat = pqc.collect_primer()
#        item = ['comFA_5', 'Primer_pair_7', 11.374516]
#        pqc.get_blast_input(item) # test self.primerlist
#        assert pqc.primerlist[-2] == ['>Lb_curva_comFA_5_P7_F\n', 'ACAACGCTTATTATTATTTGTGCCA\n']
#        assert pqc.primerlist[-1] == ['>Lb_curva_comFA_5_P7_R\n', 'AAAGGCCGCTATCTTGTCTAAT\n']
#        del pqc.primerlist[-1]
#        del pqc.primerlist[-1]
#
#        return pqc
#
#    def primerBLAST(config):
#        dbfile = dbinputfiles()
#        create_customblastdb(config, dbfile)
#        if os.path.isfile(dbfile):
#            os.remove(dbfile)
#
#        G.create_directory(pqc.primerblast_dir)
#        prep = BlastPrep(
#            pqc.primerblast_dir, pqc.primerlist,
#            "primer", pqc.config.blastseqs)
#        use_cores, inputseqs = prep.run_blastprep()
#
##        Blast(pqc.config, pqc.primerblast_dir, "primer").run_blast(
##            "primer", use_cores)
#        G.create_directory(pqc.primer_qc_dir)
#        reffile = os.path.join(testfiles_dir, "primer_nontargethits.json")
#        tofile = os.path.join(pqc.primerblast_dir, "nontargethits.json")
#        shutil.copy(reffile, tofile)
#        reffile = os.path.join(testfiles_dir, "primerBLAST_DBIDS0.txt")
#        tofile = os.path.join(pqc.primer_qc_dir, "primerBLAST_DBIDS0.txt")
#        shutil.copy(reffile, tofile)
#        tofile = os.path.join(pqc.primer_qc_dir, "BLASTnontarget0.sequences")
#        reffile = os.path.join(testfiles_dir, "BLASTnontarget0.sequences")
#        shutil.copy(reffile, tofile)
#
#        pqc.call_blastparser.run_blastparser("primer")
#
#        if os.path.isdir(tmpdir):
#            shutil.rmtree(tmpdir)
#
#        return inputseqs
#
#    def test_get_primerinfo(inputseqs):
#        inputseqs.sort()
#        modes = ['mfeprimer', "mfold", "dimercheck"] #"results" after mfeprimer (PPC)
#        for mode in modes:
#            primerinfo = pqc.get_primerinfo([inputseqs[0]], mode)
#            if mode == modes[0]:
#                assert primerinfo == [[
#                    'Lb_curva_asnS_1_P0_F',
#                    'AAACCATGTCGATGAAGAAGTTAAA',
#                    'Lb_curva_asnS_1_P0_R',
#                    'TGCCATCACGTAATTGGAGGA',
#                    'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGG'
#                    + 'AAGATTTCATTCCTCCAATTACGTGATGGCACTGCCTTTTTCCAAGGGGTTGTCGTTAA'
#                    + 'AAG']]
#            elif mode == modes[1]:
#                assert primerinfo == [[
#                    'asnS_1', 'Primer_pair_0',
#                    'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGG'
#                    + 'AAGATTTCATTCCTCCAATTACGTGATGGCA', 'Lb_curva_asnS_1_P0']]
#            elif mode == modes[2]:
#                assert primerinfo == [[
#                    'Lb_curva_asnS_1_P0',
#                    'AAACCATGTCGATGAAGAAGTTAAA',
#                    'TGCCATCACGTAATTGGAGGA']]
#
#    def test_MFEprimer_DBs():
#        for files in os.listdir(pqc.primer_qc_dir):
#            if (
#                files.startswith("BLASTnontarget")
#                and files.endswith(".sequences")
#            ):
#                pqc.dbinputfiles.append(files)
#
#        pqc.make_nontargetDB(pqc.dbinputfiles[0])
#        qc_file = os.path.join(testfiles_dir, "Lb_curva_qc_sequences.csv")
#        if os.path.isdir(pqc.summ_dir):
#            shutil.rmtree(pqc.summ_dir)
#        G.create_directory(pqc.summ_dir)
#        shutil.copy(qc_file, pqc.summ_dir)
#        primerinfo = pqc.get_primerinfo(inputseqs, "mfeprimer")
#        pqc.prepare_MFEprimer_Dbs(primerinfo)
#        dbtype = [
#            "BLASTnontarget0.sequences",
#            "template.sequences",
#            "Lb_curva.genomic"]
#        dbfiles = [
#            ".fai",
#            ".json",
#            ".log",
#            ".primerqc",
#            ".primerqc.fai"]
#
#        for db in dbtype:
#            for end in dbfiles:
#                files =  db + end
#                print(files)
#                filepath = os.path.join(pqc.primer_qc_dir, files)
#                assert os.path.isfile(filepath) == True
#                assert os.stat(filepath).st_size > 0
#
#    def test_MFEprimer(inputseqs):
#        import itertools
#        primerinfo = pqc.get_primerinfo([inputseqs[0]], "mfeprimer")
#        os.chdir(pqc.primer_qc_dir)
#        # template
#        print("Test MFEprimer_template()")
#        ppcth = [100, 90, 80, 70]
#        outcome = [[None], [None], 'Lb_curva_asnS_1_P0_F', 'Lb_curva_asnS_1_P0_F']
#        results = []
#        for index, item in enumerate(ppcth):
#            pqc.mfethreshold = item
#            result = pqc.MFEprimer_template(primerinfo[0])
#            results.append(result)
#            if outcome[index] == [None]:
#                assert result[0] == outcome[index]
#            else:
#                assert result[0][0] == outcome[index]
#        check_nontarget = pqc.write_MFEprimer_results(results, "template")
#
#        # nontarget
#        print("Test MFEprimer_nontarget()")
#        ref = [[
#            'Lb_curva_asnS_1_P0_F',
#            'AAACCATGTCGATGAAGAAGTTAAA',
#            'Lb_curva_asnS_1_P0_R',
#            'TGCCATCACGTAATTGGAGGA',
#            'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGGAAGA'
#            + 'TTTCATTCCTCCAATTACGTGATGGCACTGCCTTTTTCCAAGGGGTTGTCGTTAAAAG'
#            , 17.700000000000003],
#            [
#            'AmpID\tFpID\tRpID\tHitID\tPPC\tSize\tAmpGC\tFpTm\tRpTm\tFpDg'
#             + '\tRpDg\tBindingStart\tBindingStop\tAmpSeq',
#             '1\tLb_curva_asnS_1_P0_F\tLb_curva_asnS_1_P0_F'
#             + '\tNZ_CP020459.1:1524567-1528567\t-16.4\t364\t37.09\t7.26'
#             + '\t15.22\t-4.51\t-6.20\t3075\t3409\t'
#             + 'aaaccatgtcgatgAAGAAGTTAAAaaaccagctgaattagacaaatacctaaagagtatt'
#             + 'gactaattcattacaaaaaaagatcccgtcaatgacgagatctttttttatgctcaattact'
#             + 'aaccgcgatggtcagcacccttcacttaaacgtgcttctcgaactgccttttccggttaacc'
#             + 'acaaaatgggaatcatggctaaacccgccataatccccattttcttgttaatcctcaaagcc'
#             + 'aacccgttcgagtcagcactttatttgttttctttcttttcgctttgctttaattcttcacc'
#             + 'caatttgatgatgtatttcttcaaatcatcTTTAACTTCTtcatcgacatggttt']]
#        results = []
#        nontarget_lists = []
#        for index, item in enumerate(ppcth):
#            pqc.mfethreshold = item
#            result = pqc.MFEprimer_nontarget(check_nontarget[1], pqc.dbinputfiles[0])
#            assert result == ref
#            results.append(result)
#        nontarget_lists = list(itertools.chain(nontarget_lists, results))
#        check_assembly = pqc.write_MFEprimer_results(nontarget_lists, "nontarget")
#
#        # assembly
#        print("Test MFEprimer_assembly()")
#        ref = [[
#            'Lb_curva_asnS_1_P0_F',
#            'AAACCATGTCGATGAAGAAGTTAAA',
#            'Lb_curva_asnS_1_P0_R',
#            'TGCCATCACGTAATTGGAGGA',
#            'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGGAAGA'
#            + 'TTTCATTCCTCCAATTACGTGATGGCACTGCCTTTTTCCAAGGGGTTGTCGTTAAAAG'
#            , 17.700000000000003]]
#        results = []
#        outcome = [[None], [None], [None], 'Lb_curva_asnS_1_P0_F']
#        for index, item in enumerate(ppcth):
#            pqc.mfethreshold = item
#            result = pqc.MFEprimer_assembly(check_assembly[0])
#            results.append(result)
#            if outcome[index] == [None]:
#                assert result[0] == outcome[index]
#            else:
#                assert result[0][0] == outcome[index]
#
#        check_final = pqc.write_MFEprimer_results(results, "assembly")
#        assert check_final == ref




    def test_fullMFEprimer_run():
        outcome = [15, 42, 54, 65]
        ppcth = [100, 90, 80, 70]
        primerinfos = pqc.get_primerinfo(inputseqs, "mfeprimer")
        for index, item in enumerate(ppcth):
            pqc.mfethreshold = item
            primer_name_list = pqc.MFEprimer_QC(primerinfos)
            assert len(primer_name_list) == outcome[index]

#    def create_primerlist(inputseqs):
#        primer_qc_list = pqc.get_primerinfo(inputseqs, "mfeprimer")
#        primerfile = os.path.join(pqc.primer_qc_dir, "primerfile.fa")
#        with open(primerfile, "w") as f:
#            for primerinfo in primer_qc_list:
#                [nameF, seqF, nameR, seqR, templ_seq] = primerinfo
#                f.write(
#                        ">" + nameF + "\n" + seqF + "\n>" + nameR + "\n" + seqR + "\n")

    ref_primer = [
        'Lb_curva_g1243_1_P0', 'Lb_curva_comFA_2_P0', 'Lb_curva_g1243_1_P2',
        'Lb_curva_gshAB_11_P0', 'Lb_curva_comFA_2_P1', 'Lb_curva_comFA_5_P0',
        'Lb_curva_g1243_1_P3', 'Lb_curva_g1243_2_P0', 'Lb_curva_g4295_1_P0',
        'Lb_curva_g1243_1_P4', 'Lb_curva_comFA_2_P2', 'Lb_curva_comFA_5_P1',
        'Lb_curva_comFA_4_P1', 'Lb_curva_g4430_1_P1', 'Lb_curva_g1243_1_P6',
        'Lb_curva_g4430_1_P0', 'Lb_curva_gshAB_3_P0', 'Lb_curva_gshAB_11_P1',
        'Lb_curva_comFA_5_P2', 'Lb_curva_comFA_4_P2', 'Lb_curva_gshAB_11_P2',
        'Lb_curva_g4295_1_P1', 'Lb_curva_g1243_1_P7', 'Lb_curva_g4430_1_P2',
        'Lb_curva_gshAB_2_P1', 'Lb_curva_comFA_4_P3', 'Lb_curva_gshAB_2_P2',
        'Lb_curva_gshAB_3_P1', 'Lb_curva_gshAB_1_P0', 'Lb_curva_g4295_1_P2',
        'Lb_curva_comFA_6_P0', 'Lb_curva_comFA_3_P1', 'Lb_curva_comFA_5_P3',
        'Lb_curva_g1243_1_P8', 'Lb_curva_comFA_4_P4', 'Lb_curva_comFA_6_P1',
        'Lb_curva_g1243_1_P9', 'Lb_curva_gshAB_2_P4', 'Lb_curva_g4295_1_P3',
        'Lb_curva_comFA_4_P5', 'Lb_curva_gshAB_2_P5', 'Lb_curva_gshAB_11_P3',
        'Lb_curva_g4295_1_P4', 'Lb_curva_g4430_1_P3', 'Lb_curva_gshAB_6_P0',
        'Lb_curva_comFA_5_P5', 'Lb_curva_comFA_5_P4', 'Lb_curva_comFA_5_P6',
        'Lb_curva_gshAB_2_P6', 'Lb_curva_g4295_1_P5', 'Lb_curva_g4295_1_P6',
        'Lb_curva_asnS_2_P0', 'Lb_curva_g4295_1_P7', 'Lb_curva_comFA_3_P2',
        'Lb_curva_gshAB_3_P2', 'Lb_curva_asnS_1_P0', 'Lb_curva_g4295_1_P8',
        'Lb_curva_gshAB_8_P0', 'Lb_curva_comFA_5_P7', 'Lb_curva_comFA_5_P8',
        'Lb_curva_g4295_1_P9']

    newmockprimer = {
        "mock_1": {
            "Primer_pairs": 2, "template_seq":
            "NNNGGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCCNNNCTTGTTTACGTGGCGGCGNNN",
            "Primer_pair_0": {
                "primer_P_penalty": 5.913076, "primer_L_TM": 58.15,
                "primer_R_sequence": "CGCCGCCACGTAAACAAG",
                "primer_L_penalty": 2.886877, "primer_R_penalty": 3.007342,
                "amplicon_seq": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCCNNNCTTGTTTACGTGGCGGCG",
                "primer_L_sequence": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCC",
                "product_TM": 80.1424, "primer_R_TM": 57.03, "product_size": 82},
            "Primer_pair_1": {
                "primer_P_penalty": 5.913076, "primer_L_TM": 58.15,
                "primer_R_sequence": "CGCCGCCACGTAAACAAG",
                "primer_L_penalty": 2.886877, "primer_R_penalty": 3.007342,
                "amplicon_seq": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCCNNNCTTGTTTACGTGGCGGCG",
                "primer_L_sequence": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCC",
                "product_TM": 80.1424, "primer_R_TM": 57.03, "product_size": 82}},
        "mock_9": {
            "Primer_pairs": 1, "template_seq": "NNNCGCAGCAACACACACACGNNNGGGGGGACACTCTTTCCCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCCNNN",
            "Primer_pair_3":{
                "primer_P_penalty": 2.427616, "primer_L_TM": 59.691,
                "primer_R_sequence": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCC",
                "primer_L_penalty": 0.345215, "primer_R_penalty": 2.064444,
                "amplicon_seq": "CGCAGCAACACACACACGNNNGGGGGGACACTCTTTCCCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCC",
                "primer_L_sequence": "CGCAGCAACACACACACG", "product_TM": 81.2095,
                "primer_R_TM": 60.029, "product_size": 88}},

        "dimer_1": {
                "Primer_pairs": 1, "template_seq": "NNNN",
                "Primer_pair_0":{
                "primer_P_penalty": 2.427616, "primer_L_TM": 59.691,
                "primer_R_sequence": "TTTTTTAAAAAA",
                "primer_L_penalty": 0.345215, "primer_R_penalty": 2.064444,
                "amplicon_seq": "NNNN",
                "primer_L_sequence": "TTTTTTAAAAAA", "product_TM": 81.2095,
                "primer_R_TM": 60.029, "product_size": 88}}}

    tmpdir = os.path.join(BASE_PATH, "tests", "tmp")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    config.customdb = os.path.join(tmpdir, "primer_customdb.fas")
    config.blastdbv5 = False

    def dbinputfiles():
        filenames = [
            "GCF_004088235v1_20191001.fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna"]
        dbfile = os.path.join(testfiles_dir, "primer_customdb.fas")
        with open(dbfile, "w") as f:
            for filename in filenames:
                filepath = os.path.join(testfiles_dir, filename)
                records = SeqIO.parse(filepath, "fasta")
                for record in records:
                    if record.id == record.description:
                        description = (
                            record.id + " Lactobacillus curvatus strain SRCM103465")
                        record.description = description
                    SeqIO.write(record, f, "fasta")
        return dbfile

    def create_customblastdb(config, infile):
        cmd = [
            "makeblastdb", "-in", infile, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", config.customdb]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, log=False, printoption=False)

    def get_primer3_dict():
        reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
        with open(reffile) as f:
            for line in f:
                primer3dict = json.loads(line)
        return primer3dict

    pqc = PrimerQualityControl(config, {})
    exitstat = pqc.collect_primer()
    assert exitstat == 1
    primer3dict = get_primer3_dict()
    primer3dict.update(newmockprimer)
    pqc = PrimerQualityControl(config, primer3dict)

    exitstat = pqc.collect_primer()
    assert exitstat == 0

    hairpins = pqc.hairpin_check()
    hairpins.sort()
    ref_hairpins = [
        'Lb_curva_mock_9_P3', 'Lb_curva_mock_1_P1', 'Lb_curva_mock_1_P0']
    ref_hairpins.sort()
    assert hairpins == ref_hairpins

    pqc.primerdimer_check(hairpins)
    primerlist, blastseqs = pqc.filter_primerdimer()
#    print(primer)
    primerlist.sort()
    ref_primer.sort()
    assert primerlist == ref_primer

    blastsum = os.path.join(pqc.primerblast_dir, "nontargethits.json")
    infile = dbinputfiles()
    create_customblastdb(config, infile)

    if not os.path.isfile(blastsum):
        G.create_directory(pqc.primerblast_dir)
        prep = BlastPrep(
            pqc.primerblast_dir, blastseqs,
            "primer", pqc.config.blastseqs)
        use_cores, inputseqs = prep.run_blastprep()
        bla = Blast(pqc.config, pqc.primerblast_dir, "primer")
        bla.run_blast("primer", use_cores)

    pqc.call_blastparser.run_blastparser("primer")

    primer_spec_list = pqc.get_primerinfo(primerlist, "mfeprimer")

    for files in os.listdir(pqc.primer_qc_dir):
        if (
            files.startswith("BLASTnontarget")
            and files.endswith(".sequences")
        ):
            pqc.dbinputfiles.append(files)

    pqc.prepare_MFEprimer_Dbs(primer_spec_list)
    os.chdir(pqc.primer_qc_dir)
#    for item in primer_spec_list:
#        pqc.MFEprimer_template(item)
#        break
#
    pqc.run_primer_qc()
#    resultlist = G.run_parallel(pqc.run_mfeprimerdimer, pqc.primerlist, hairpins)


def dev_new_QC(config):
    from speciesprimer import PrimerQualityControl
    newmockprimer = {
        "mock_1": {
            "Primer_pairs": 2, "template_seq":
            "NNNGGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCCNNNCTTGTTTACGTGGCGGCGNNN",
            "Primer_pair_0": {
                "primer_P_penalty": 5.913076, "primer_L_TM": 58.15,
                "primer_R_sequence": "CGCCGCCACGTAAACAAG",
                "primer_L_penalty": 2.886877, "primer_R_penalty": 3.007342,
                "template_seq": "NNNGGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCCNNNCTTGTTTACGTGGCGGCGNNN",
                "amplicon_seq": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCCNNNCTTGTTTACGTGGCGGCG",
                "primer_L_sequence": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCC",
                "product_TM": 80.1424, "primer_R_TM": 57.03, "product_size": 82},
            "Primer_pair_1": {
                "primer_P_penalty": 5.913076, "primer_L_TM": 58.15,
                "primer_R_sequence": "CGCCGCCACGTAAACAAG",
                "primer_L_penalty": 2.886877, "primer_R_penalty": 3.007342,
                "amplicon_seq": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCCNNNCTTGTTTACGTGGCGGCG",
                "primer_L_sequence": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCC",
                "product_TM": 80.1424, "primer_R_TM": 57.03, "product_size": 82}},
        "mock_9": {
            "Primer_pairs": 1, "template_seq": "NNNCGCAGCAACACACACACGNNNGGGGGGACACTCTTTCCCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCCNNN",
            "Primer_pair_3":{
                "primer_P_penalty": 2.427616, "primer_L_TM": 59.691,
                "primer_R_sequence": "GGACACTCTTTCCCTACACGACGCTCTTCCGATCTATGGGAAAGAGTGTCCCCCC",
                "primer_L_penalty": 0.345215, "primer_R_penalty": 2.064444,
                "amplicon_seq": "CGCAGCAACACACACACGNNNGGGGGGACACTCTTTCCCATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCC",
                "primer_L_sequence": "CGCAGCAACACACACACG", "product_TM": 81.2095,
                "primer_R_TM": 60.029, "product_size": 88}},

        "dimer_1": {
                "Primer_pairs": 1, "template_seq": "NNNN",
                "Primer_pair_0":{
                "primer_P_penalty": 2.427616, "primer_L_TM": 59.691,
                "primer_R_sequence": "TTTTTTAAAAAA",
                "primer_L_penalty": 0.345215, "primer_R_penalty": 2.064444,
                "amplicon_seq": "NNNN",
                "primer_L_sequence": "TTTTTTAAAAAA", "product_TM": 81.2095,
                "primer_R_TM": 60.029, "product_size": 88}}}

    tmpdir = os.path.join(BASE_PATH, "tests", "tmp")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    config.customdb = os.path.join(tmpdir, "primer_customdb.fas")
    config.blastdbv5 = False

    def dbinputfiles():
        filenames = [
            "GCF_004088235v1_20191001.fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna"]
        dbfile = os.path.join(testfiles_dir, "primer_customdb.fas")
        with open(dbfile, "w") as f:
            for filename in filenames:
                filepath = os.path.join(testfiles_dir, filename)
                records = SeqIO.parse(filepath, "fasta")
                for record in records:
                    if record.id == record.description:
                        description = (
                            record.id + " Lactobacillus curvatus strain SRCM103465")
                        record.description = description
                    SeqIO.write(record, f, "fasta")
        return dbfile

    def create_customblastdb(config, infile):
        cmd = [
            "makeblastdb", "-in", infile, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", config.customdb]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, log=False, printoption=False)

    def get_primer3_dict():
        reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
        with open(reffile) as f:
            for line in f:
                primer3dict = json.loads(line)
        return primer3dict

    def get_primer3_dict():
        reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
        with open(reffile) as f:
            for line in f:
                primer3dict = json.loads(line)
        return primer3dict

    primer3dict = get_primer3_dict()
    primer3dict.update(newmockprimer)

#    print(primer3dict)
    from speciesprimer import BlastPrep
    import multiprocessing
    pqc = PrimerQualityControl(config, primer3dict)
    dbfile = dbinputfiles()
    create_customblastdb(config, dbfile)

#    os.chdir(pqc.primer_qc_dir)
#    total_results = []
#    if pqc.collect_primer() == 0:
#        blastsum = os.path.join(pqc.primerblast_dir, "nontargethits.json")
#        hairpins = pqc.hairpin_check()
#        pqc.primerdimer_check(hairpins)
#
#    excluded = pqc.filter_primerdimer()
#    pqc.prepare_MFEprimer_Dbs(pqc.primerlist)
##    pqc.run_primer_qc()
#    nontarget_input = pqc.MFEprimer_QC_first(excluded)
#    nontarget_blastseqs = []
#    for primer in nontarget_input:
#        nontarget_blastseqs.append(primer[0:4])
#        ppc = primer[6]
#        pp_name = "_".join(primer[0].split("_")[0:-1])
#        target_id = "_".join(pp_name.split("_")[-3:-1])
#        primerpair = "Primer_pair_" + pp_name.split("_P")[-1]
#        pqc.primer3_dict[target_id][primerpair].update({"PPC": ppc})
#
#    blastsum = os.path.join(pqc.primerblast_dir, "nontargethits.json")
#
#    if not os.path.isfile(blastsum):
#        G.create_directory(pqc.primerblast_dir)
#        prep = BlastPrep(
#            pqc.primerblast_dir, nontarget_blastseqs,
#            "primer", pqc.config.blastseqs)
#        use_cores, inputseqs = prep.run_blastprep()
#        bla = Blast(pqc.config, pqc.primerblast_dir, "primer")
#        bla.run_blast("primer", use_cores)
#
#    pqc.call_blastparser.run_blastparser("primer")
##
##            primer_spec_list = self.get_primerinfo(primerlist, "mfeprimer")
##
#
#    for files in os.listdir(pqc.primer_qc_dir):
#        if (
#            files.startswith("BLASTnontarget")
#            and files.endswith(".sequences")
#        ):
#            pqc.dbinputfiles.append(files)
#
#    # parallelization try
#    pool = multiprocessing.Pool(processes=4)
#    results = [
#        pool.apply_async(
#            pqc.make_nontargetDB, args=(inputfiles,))
#        for inputfiles in pqc.dbinputfiles]
#    output = [p.get() for p in results]
#
#    for primerinfo in nontarget_input:
#        pqc.MFEprimer_nontarget(primerinfo)

    pqc.run_primer_qc()

#    hairpins = pqc.hairpin_check()
#
#    for i in range(0,2000):
#        pqc.primerlist.append(pqc.primerlist[0])

#    pqc.primerdimer_check(hairpins)
#    excluded = pqc.filter_primerdimer()
#    primerinfos = []
#    for primer in pqc.primerlist:
#        pp_name = "_".join(primer[0].split("_")[0:-1])
#        if pp_name not in excluded:
#            primerinfos.append(primer)
#    os.chdir(pqc.primer_qc_dir)
#
#    template_results = G.run_parallel(pqc.MFEprimer_template, primerinfos)
#    assembly_input = pqc.write_MFEprimer_results(
#            template_results, "template")
#    start = time.time()
#    assembly_results = G.run_parallel(pqc.MFEprimer_assembly, assembly_input)
#    nontarget_input = pqc.write_MFEprimer_results(assembly_results, "assembly")
#
#    print(len(nontarget_input))
#    duration = time.time() - start
#    print("parallel assembly: "
#                + str(timedelta(seconds=duration)).split(".")[0])


#dev_new_QC(config)

#test_PrimerQualityControl_specificitycheck(config)

def speed_up(config):
    from speciesprimer import BlastParser
    outputfiles = []
    bla = BlastParser(config, "primer")
    for files in os.listdir(bla.primer_qc_dir):

        if files.endswith(".sequences") or files.endswith(".txt"):
            outfile = os.path.join(bla.primer_qc_dir, files)
            outputfiles.append(outfile)

    for filename in outputfiles:
        if os.path.isfile(filename):
            print(filename)
            os.remove(filename)


    file_path = os.path.join(bla.primerblast_dir, "nontargethits.json")
    with open(file_path) as f:
        for line in f:
            nonred_dict = json.loads(line)
    print(len(nonred_dict))
#    additional_data = []
#    for key in nonred_dict.keys():
#        additional_data.append([key, nonred_dict[key]])
#    for item in additional_data:
#        for i in range(0, 50):
#            nonred_dict.update({item[0] + str(i): item[1]})
#
#    print(len(nonred_dict))
    print("start time")
    start = time.time()
    bla.highspeed_primerBLAST_DBIDS(nonred_dict)
    duration = time.time() - start
    print("new exract: "
                + str(timedelta(seconds=duration)).split(".")[0])

#    for filename in files:
#        if os.path.isfile(filename):
#            os.remove(filename)
#
#
#    start = time.time()
#    bla.create_primerBLAST_DBIDS(nonred_dict)
#    duration = time.time() - start
#    print("old exract: "
#                + str(timedelta(seconds=duration)).split(".")[0])


"""
Idee für redundanz zu reduizieren und weniger sequenzen zu extrahieren
neue dict accession == key
if key in newdict.keys update position
ergibt eine list von positionen diese können danach schön nach accession sortiert ausgewertet werden

"""
speed_up(config)

# 4x worker parallel assembly: 0:00:50 (normal)
# parallel exract: 0:00:09

def benchmark_primerblastparser(config):
# neue idee threshold für anzahl BLASTDB extraction für die extraktion der
# ganzen sequenz und anschliessende sequenz extraktion mittels indexierung

    # sort according accession
    # sub sort according position? possible?

    def write_DBIDS(self, prefix, part, suffix, info):
        filename = prefix + str(part) + suffix
        filepath = os.path.join(self.primer_qc_dir, filename)
        with open(filepath, "a+") as f:
            f.write(info)

    def create_primerBLAST_DBIDS(self, nonred_dict):
        print("\nGet sequence accessions of BLAST hits\n")
        G.logger("> Get sequence accessions of BLAST hits")
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
            written = []
            prefix = "primerBLAST_DBIDS"
            suffix = ".txt"

            for key in nonred_dict.keys():
                if not len(nonred_dict[key]) == 0:
                    for species in nonred_dict[key]:
                        db_id = nonred_dict[key][species]['main_id']
                        sbjct_start = (
                            nonred_dict[key][species]["subject_start"])
                        if not [db_id, sbjct_start] in written:
                            written.append([db_id, sbjct_start])
                            idcount = idcount + 1
                            part = idcount//self.maxgroupsize
                            info = str(db_id) + " " + str(sbjct_start) + "\n"
                            self.write_DBIDS(prefix, part, suffix, info)
                            filename = prefix + str(part) + suffix
                            if filename not in DBID_files:
                                DBID_files.append(filename)

        if idcount == 0:
            print("Error did not find any sequences for the non-target DB")
            G.logger("> Error did not find any sequences for non-target DB")
            return
        else:
            ntseq = str(idcount)
            info1 = ntseq + " sequences for non-target DB"
            info2 = "Start extraction of sequences from BLAST DB"
            print("\n" + info1)
            print(info2)
            G.logger("> " + info1)
            G.logger("> " + info2)
            pool = multiprocessing.Pool(processes=4)
            results = [
                pool.apply_async(
                    self.write_nontarget_sequences, args=(files,))
                for files in DBID_files]
            output = [p.get() for p in results]

    def run_blastparser(self, conserved_seq_dict):
        if self.mode == "primer":
            file_path = os.path.join(self.primerblast_dir, "nontargethits.json")
            G.logger("Run: run_blastparser(" + self.target + "), primer")
            align_dict = self.blast_parser(self.primerblast_dir)
            if not os.path.isfile(file_path):
                nonred_dict = self.remove_redundanthits(align_dict)
                self.write_nontargethits(
                        self.primerblast_dir, nonred_dict,
                        "json")
            else:
                nonred_dict = align_dict
            self.create_primerBLAST_DBIDS(nonred_dict)

            duration = time.time() - self.start
            G.logger(
                "> Primer blast parser time: "
                + str(timedelta(seconds=duration)).split(".")[0])

