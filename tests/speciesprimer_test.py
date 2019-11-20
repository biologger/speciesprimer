#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import pytest
import json
import time
import csv
from Bio import SeqIO
from Bio import Entrez
from basicfunctions import HelperFunctions as H
from basicfunctions import GeneralFunctions as G
from basicfunctions import ParallelFunctions as P


msg = (
    """
    Works only in the Docker container!!!
    - Start the container
        sudo docker start {Containername}
    - Start an interactive terminal in the container
        sudo docker exec -it {Containername} bash
    - Start the tests in the container terminal
        cd /
        pytest -vv --cov=pipeline /tests/
    """)


# /tests
BASE_PATH = os.path.dirname(os.path.abspath(__file__))
pipe_dir = os.path.join(BASE_PATH.split("tests")[0], "pipeline")
dict_path = os.path.join(pipe_dir, "dictionaries")
tmpdir = os.path.join("/", "primerdesign", "tmp")
testfiles_dir = os.path.join(BASE_PATH, "testfiles")
ref_data = os.path.join(BASE_PATH, "testfiles", "ref")
confargs = {
    "ignore_qc": False, "mfethreshold": 90, "maxsize": 200,
    "target": "Lactobacillus_curvatus", "nolist": False, "skip_tree": False,
    "blastseqs": 1000, "mfold": -3.0, "mpprimer": -3.5,
    "offline": False,
    "path": os.path.join("/", "primerdesign", "test"),
    "probe": False, "exception": [], "minsize": 70, "skip_download": True,
    "customdb": None, "assemblylevel": ["all"], "qc_gene": ["rRNA"],
    "blastdbv5": False, "intermediate": True,
    "nontargetlist": ["Lactobacillus sakei"]}


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


@pytest.fixture
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


def compare_ref_files(results_dir, ref_dir):

    resfileslist = []
    reffileslist = []
    resrecords = []
    refrecords = []

    def compare_files(resfiles, reffiles):
        resrecords = []
        refrecords = []
        records = SeqIO.parse(resfiles, "fasta")
        for record in records:
            resrecords.append([str(record.id), str(record.seq)])
        records = SeqIO.parse(reffiles, "fasta")
        for record in records:
            refrecords.append([str(record.id), str(record.seq)])

        resrecords.sort()
        refrecords.sort()

        for index, item in enumerate(resrecords):
            assert item[0] == refrecords[index][0]
            assert item[1] == refrecords[index][1]

    if os.path.isdir(results_dir):
        for files in os.listdir(results_dir):
            if not files.endswith(".txt"):
                resfiles = os.path.join(results_dir, files)
                resfileslist.append(resfiles)
        resfileslist.sort()

        for files in resfileslist:
            resrecords.append(files)

        for files in os.listdir(ref_dir):
            if not files.endswith(".txt"):
                reffiles = os.path.join(ref_dir, files)
                reffileslist.append(reffiles)
        reffileslist.sort()

        for files in reffileslist:
            refrecords.append(files)

        for index, item in enumerate(resrecords):
            compare_files(item, refrecords[index])

    else:
        resfiles = os.path.join(results_dir)
        reffiles = os.path.join(ref_dir)
        compare_files(resfiles, reffiles)


def test_commandline():
    os.chdir("/")
    import speciesprimer
    parser = speciesprimer.commandline()
    args = parser.parse_args([])
    assert args.assemblylevel == ['all']
    assert args.blastdbv5 is False
    assert args.blastseqs == 1000
    assert args.customdb is None
    assert args.email is None
    assert args.exception == []
    assert args.ignore_qc is False
    assert args.intermediate is False
    assert args.maxsize == 200
    assert args.mfethreshold == 90
    assert args.mfold == -3.0
    assert args.minsize == 70
    assert args.mpprimer == -3.5
    assert args.nolist is False
    assert args.offline is False
    assert args.path == '/'
    assert args.probe is False
    assert args.qc_gene == ['rRNA']
    assert args.skip_download is False
    assert args.skip_tree is False
    assert args.target is None


def test_CLIconf(config):
    assert config.minsize == confargs['minsize']
    assert config.maxsize == confargs['maxsize']
    assert config.ignore_qc == confargs['ignore_qc']
    assert config.mfethreshold == confargs['mfethreshold']
    assert config.target == confargs['target']
    assert config.nolist == confargs['nolist']
    assert config.skip_tree == confargs['skip_tree']
    assert config.blastseqs == confargs['blastseqs']
    assert config.mfold == confargs['mfold']
    assert config.mpprimer == confargs['mpprimer']
    assert config.offline == confargs['offline']
    assert config.probe == confargs['probe']
    assert config.exception == confargs['exception']
    assert config.customdb == confargs['customdb']
    assert config.skip_download == confargs['skip_download']
    assert config.assemblylevel == confargs['assemblylevel']
    assert config.qc_gene == confargs['qc_gene']
    assert config.blastdbv5 == confargs['blastdbv5']


def test_auto_run_config():

    def prepare_tmp_db():
        t = os.path.join(BASE_PATH, "testfiles", "tmp_config.json")
        # Docker only
        tmp_path = os.path.join(pipe_dir, "tmp_config.json")
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        shutil.copy(t, tmp_path)

    def change_tmp_db():
        tmp_path = os.path.join(pipe_dir, "tmp_config.json")
        with open(tmp_path) as f:
            for line in f:
                tmp_dict = json.loads(line)
        tmp_dict["new_run"].update(
                {'modus': "continue", "targets": "Lactobacillus_curvatus"})
        with open(tmp_path, "w") as f:
            f.write(json.dumps(tmp_dict))

    def run_autorun():
        from speciesprimer import auto_run
        from speciesprimer import CLIconf
        targets, conf_from_file, use_configfile = auto_run()
        nontargetlist = []
        assert targets == ["Lactobacillus_curvatus"]
        assert use_configfile is True
        for target in targets:
            target = target.capitalize()
            if use_configfile:
                (
                    minsize, maxsize, mpprimer, exception, target, path,
                    intermediate, qc_gene, mfold, skip_download,
                    assemblylevel, skip_tree, nolist,
                    offline, ignore_qc, mfethreshold, customdb,
                    blastseqs, probe, blastdbv5
                ) = conf_from_file.get_config(target)

        assert minsize == confargs['minsize']
        assert maxsize == confargs['maxsize']
        assert ignore_qc == confargs['ignore_qc']
        assert mfethreshold == confargs['mfethreshold']
        assert target == confargs['target']
        assert nolist == confargs['nolist']
        assert skip_tree == confargs['skip_tree']
        assert blastseqs == confargs['blastseqs']
        assert mfold == confargs['mfold']
        assert mpprimer == confargs['mpprimer']
        assert offline == confargs['offline']
        assert probe == confargs['probe']
        assert exception == confargs['exception']
        assert customdb == confargs['customdb']
        assert skip_download == confargs['skip_download']
        assert assemblylevel == confargs['assemblylevel']
        assert qc_gene == confargs['qc_gene']
        assert blastdbv5 == confargs['blastdbv5']

        tmpconfig = CLIconf(
            minsize, maxsize, mpprimer, exception, target, path,
            intermediate, qc_gene, mfold, skip_download,
            assemblylevel, nontargetlist, skip_tree,
            nolist, offline, ignore_qc, mfethreshold, customdb,
            blastseqs, probe, blastdbv5)

        tmpconfig.save_config()

    prepare_tmp_db()
    run_autorun()
    change_tmp_db()
    run_autorun()


def test_DataCollection(config, monkeypatch):
    dtd_path = os.path.join(
            "/", "root", ".config", "biopython", "Bio", "Entrez", "DTDs")
    dtd1 = os.path.join(dtd_path, "esummary_assembly.dtd")
    dtd2 = os.path.join(dtd_path, "efetch.dtd")
    mockdtd1 = os.path.join(
                    testfiles_dir, "entrezmocks", "esummary_assembly.dtd")
    mockdtd2 = os.path.join(
                    testfiles_dir, "entrezmocks", "efetch.dtd")
    mocks = [mockdtd1, mockdtd2]
    if not os.path.isdir(dtd_path):
        G.create_directory(dtd_path)
    for i, filepath in enumerate([dtd1, dtd2]):
        if not os.path.isfile(filepath):
            shutil.copy(mocks[i], filepath)


    from speciesprimer import DataCollection
    shutil.rmtree(os.path.join(config.path, config.target))
    downdir = os.path.join(
            tmpdir, "download", "GCF_902362325.1_MGYG-HGUT-00020")
    G.create_directory(downdir)
    genomefile = os.path.join(
        testfiles_dir, "GCF_902362325.1_MGYG-HGUT-00020_genomic.fna.gz")
    downfile = os.path.join(
            downdir, "GCF_902362325.1_MGYG-HGUT-00020_genomic.fna.gz")
    shutil.copy(genomefile, downfile)
    assert os.path.isfile(downfile)

    DC = DataCollection(config)

    def test_get_email_from_config(config):
        email = H.get_email_for_Entrez()
        assert email == "biologger@protonmail.com"

    def test_get_taxid(target, monkeypatch):
        def mock_taxid(db, term):
            mockfile = os.path.join(
                    testfiles_dir, "entrezmocks", "esearchmock01.xml")
            f = open(mockfile)
            return f

        def mock_syn(db, id):
            mockfile = os.path.join(
                    testfiles_dir, "entrezmocks", "efetchmock01.xml")
            f = open(mockfile)
            return f

        monkeypatch.setattr(Entrez, "esearch", mock_taxid)
        monkeypatch.setattr(Entrez, "efetch", mock_syn)
        syn, taxid = DC.get_taxid(target)
        assert taxid == str(28038)
        assert syn == [
            'Bacterium curvatum', 'Lactobacillus curvatus subsp. curvatus',
            'Lactobacillus sp. N55', 'Lactobacillus sp. N61']

        # get_taxid fails
        def mock_taxid(db, term):
            mockfile = os.path.join(
                    testfiles_dir, "entrezmocks", "esearchmock02.xml")
            f = open(mockfile)
            return f

        monkeypatch.setattr(Entrez, "esearch", mock_taxid)
        syn, taxid = DC.get_taxid(target)
        assert taxid is None

    def test_syn_exceptions(config):
        # standard case
        confdict = config.__dict__
        config.exception = ["Lactobacillus curvatus"]
        DC = DataCollection(config)
        G.create_directory(DC.config_dir)
        conffile = os.path.join(DC.config_dir, "config.json")
        with open(conffile, "w") as f:
            f.write(json.dumps(confdict))
        syn = ['Bacterium_curvatum']
        exceptions = DC.add_synonym_exceptions(syn)
        assert exceptions == ['Lactobacillus curvatus', 'Bacterium_curvatum']
        config.exception = ['Lactobacillus curvatus', 'Bacterium_curvatum']
        DC = DataCollection(config)
        syn = ['Lactobacillus curvatus', 'Bacterium_curvatum']
        exceptions = DC.add_synonym_exceptions(syn)
        assert exceptions == ['Lactobacillus curvatus', 'Bacterium_curvatum']

    def test_create_GI_list():
        filepath = os.path.join(DC.config_dir, "no_blast.gi")
        indictpath = os.path.join(dict_path, "no_blast.gi")
        with open(indictpath, "w") as f:
            f.write("1231231231")
        DC.create_GI_list()
        gi_list = []
        with open(filepath) as f:
            for line in f:
                gi_list.append(line.strip())

        assert gi_list == ["1231231231"]

        defaultpath = os.path.join(dict_path, "default", "no_blast.gi")
        if os.path.isfile(filepath):
            shutil.copy(defaultpath, indictpath)

        if os.path.isfile(filepath):
            os.remove(filepath)

    def mocked_urlopen_fail(monkeypatch):
        import urllib.error
        import urllib.request

        def mocked(url, data):
            raise urllib.error.HTTPError(url, 400, "Bad request", None, None)

        monkeypatch.setattr(urllib.request, "urlopen", mocked)

        with pytest.raises(urllib.error.HTTPError):
            DC.ncbi_download()

    def test_ncbi_download(taxid, monkeypatch):
        def mock_getsummary(db, term, retmax):
            mockfile = os.path.join(
                    testfiles_dir, "entrezmocks", "getsummarymock.xml")
            f = open(mockfile)
            return f

        def mock_getlinks(db, id, rettype, retmode):
            mockfile = os.path.join(
                    testfiles_dir, "entrezmocks", "getlinksmock.xml")
            f = open(mockfile)
            return f

        monkeypatch.setattr(Entrez, "esearch", mock_getsummary)
        monkeypatch.setattr(Entrez, "efetch", mock_getlinks)

        DC.get_ncbi_links(taxid, 1)
        DC.ncbi_download()
        filepath = os.path.join(
            DC.target_dir, "genomic_fna",
            "GCF_902362325.1_MGYG-HGUT-00020_genomic.fna")
        # will not download the file a second time
        # because it is already extracted
        time.sleep(2)
        DC.ncbi_download()
        # prep annotated files
        os.remove(filepath)
        dirs = [DC.gff_dir, DC.ffn_dir, DC.fna_dir]
        for direct in dirs:
            G.create_directory(direct)
        mockname = "GCF_902362325v1_2019118"
        endings = ["fna", "gff", "ffn"]
        for end in endings:
            filepath = os.path.join(
                    DC.target_dir, end + "_files", mockname + "." + end)
            with open(filepath, "w") as f:
                f.write("Mock " + mockname + "." + end)
        # No download since annotated files are present
        DC.ncbi_download()
        annotation_dirs, annotated = DC.run_prokka()
        assert annotation_dirs == []
        assert annotated == []
        assert os.path.isfile(filepath + ".gz") is False
        for direct in dirs:
            shutil.rmtree(direct)
        for direct in dirs:
            G.create_directory(direct)
        # test excluded files
        excluded_dir = os.path.join(
                DC.config.path, "excludedassemblies", DC.config.target)
        excludedfile = os.path.join(excluded_dir, "excluded_list.txt")
        G.create_directory(excluded_dir)
        with open(excludedfile, "w") as f:
            f.write("GCF_902362325v1\n")
        annotation_dirs, annotated = DC.run_prokka()
        assert annotation_dirs == []
        # will not download the file because it is in excluded list
        time.sleep(2)
        DC.ncbi_download()
        assert os.path.isfile(filepath) is False
        shutil.rmtree(excluded_dir)
        mocked_urlopen_fail(monkeypatch)
        fna = os.path.join(config.path, config.target, "genomic_fna")
        shutil.rmtree(fna)
        G.create_directory(fna)
        os.chdir(config.path)

    def prepare_prokka(config):
        targetdir = os.path.join(config.path, config.target)
        fileformat = ["fna", "gff", "ffn"]
        for fo in fileformat:
            for files in os.listdir(testfiles_dir):
                if files.endswith("." + fo):
                    fromfile = os.path.join(testfiles_dir, files)
                    tofile = os.path.join(targetdir, fo + "_files", files)
                    shutil.copy(fromfile, tofile)
                    if fo == "fna":
                        tofile = os.path.join(targetdir, "genomic_fna", files)
                        shutil.copy(fromfile, tofile)

    def test_run_prokka():
        annotation_dirs, annotated = DC.run_prokka()
        assert annotated == ["GCF_004088235v1"]
        DC.copy_genome_files()

    def remove_prokka_testfiles():
        fileformat = ["fna", "gff", "ffn"]
        targetdir = os.path.join(config.path, config.target)
        for form in fileformat:
            dirpath = os.path.join(targetdir, form + "_files")
            if os.path.isdir(dirpath):
                shutil.rmtree(dirpath)
        dirpath = os.path.join(targetdir, "fna_genomic")
        if os.path.isdir(dirpath):
            shutil.rmtree(dirpath)
            G.create_directory(dirpath)

    test_get_email_from_config(config)
    DC.prepare_dirs()
    test_get_taxid(config.target, monkeypatch)
    test_ncbi_download("28038", monkeypatch)
    test_syn_exceptions(config)
    test_create_GI_list()
    G.create_directory(DC.gff_dir)
    G.create_directory(DC.ffn_dir)
    G.create_directory(DC.fna_dir)
    prepare_prokka(config)
    test_run_prokka()
    remove_prokka_testfiles()


def test_QualityControl(config):
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    targetdir = os.path.join(config.path, config.target)
    config.blastdbv5 = False
    qc_gene = config.qc_gene[0]
    G.create_directory(tmpdir)
    config.customdb = os.path.join(tmpdir, "customdb.fas")

    def create_customblastdb(infile=None):
        if infile is None:
            infile = os.path.join(testfiles_dir, "customdb.fas")
        cmd = [
            "makeblastdb", "-in", infile, "-parse_seqids", "-title",
            "mock16SDB", "-dbtype", "nucl", "-out", config.customdb]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    from speciesprimer import QualityControl
    QC = QualityControl(config)

    def make_duplicate(filename):
        lines = []
        with open(filename, "r") as f:
            for line in f:
                if "v1" in line:
                    line = "v2".join(line.split("v1"))
                    lines.append(line)
                else:
                    lines.append(line)

        with open(filename, "w") as f:
            for line in lines:
                f.write(line)

    def prepare_QC_testfiles(config):
        QC = QualityControl(config)
        G.create_directory(QC.gff_dir)
        G.create_directory(QC.ffn_dir)
        G.create_directory(QC.fna_dir)
        excluded_dir = os.path.join(
                QC.config.path, "excludedassemblies", QC.config.target)
        G.create_directory(excluded_dir)
        # GCF_004088235v1_20191001
        testcases = [
            "004088235v1_20191001", "maxcontigs_date",
            "noseq_date", "duplicate"]
        fileformat = ["fna", "gff", "ffn"]
        for testcase in testcases:
            for fo in fileformat:
                for files in os.listdir(testfiles_dir):
                    if files.endswith(testcases[0] + "."+fo):
                        fromfile = os.path.join(testfiles_dir, files)
                        if testcase == "duplicate":
                            files = "v2".join(files.split("v1"))
                            tofile = os.path.join(
                                targetdir, fo + "_files", files)
                            shutil.copy(fromfile, tofile)
                            make_duplicate(tofile)
                        else:
                            tofile = os.path.join(
                                targetdir, fo + "_files",
                                "GCF_" + testcase + "." + fo)
                            shutil.copy(fromfile, tofile)

        def make_maxcontigs():
            name = "GCF_maxcontigs_date.fna"
            filepath = os.path.join(targetdir, "fna_files", name)
            nonsense = ">nonsense_contig\nATTAG\n"
            with open(filepath, "w") as f:
                for i in range(0, 500):
                    f.write(nonsense)

        def remove_qc_seq():
            keeplines = []
            name = "GCF_noseq_date.ffn"
            filepath = os.path.join(targetdir, "ffn_files", name)
            with open(filepath, "r") as f:
                for line in f:
                    keeplines.append(line)

            for gene in QC.searchdict.values():
                for index, line in enumerate(keeplines):
                    if gene in line:
                        del keeplines[index]

            with open(filepath, "w") as f:
                for line in keeplines:
                    f.write(line)

        make_maxcontigs()
        remove_qc_seq()

    prepare_QC_testfiles(config)
    create_customblastdb()

    def prepare_GI_list():
        filepath = os.path.join(QC.config_dir, "no_blast.gi")
        with open(filepath, "w") as f:
            f.write("1231231231")

    def remove_GI_list():
        filepath = os.path.join(QC.config_dir, "no_blast.gi")
        if os.path.isfile(filepath):
            os.remove(filepath)

    # test without GI list QC should fail
    def test_get_excluded_gis():
        remove_GI_list()
        excluded_gis = QC.get_excluded_gis()
        assert excluded_gis == []
        prepare_GI_list()
        excluded_gis = QC.get_excluded_gis()
        assert excluded_gis == ["1231231231"]
        remove_GI_list()
        excluded_gis = QC.get_excluded_gis()
        assert excluded_gis == []

    def test_search_qc_gene():
        gff = []
        for files in os.listdir(QC.gff_dir):
            if files not in gff:
                gff.append(files)
        for file_name in gff:
            QC.search_qc_gene(file_name, qc_gene)
        assert len(QC.qc_gene_search) == 12
        return gff

    def test_count_contigs(gff_list, contiglimit):
        gff_list = QC.count_contigs(gff_list, contiglimit)
        gff_list.sort()
        assert gff_list == [
                'GCF_004088235v1_20191001.gff',
                'GCF_004088235v2_20191001.gff',
                'GCF_noseq_date.gff']

        return gff_list

    def test_identify_duplicates(gff_list):
        gff_list = QC.identify_duplicates(gff_list)
        gff_list.sort()
        assert gff_list == [
            'GCF_004088235v2_20191001.gff',
            'GCF_noseq_date.gff']
        return gff_list

    def test_check_no_sequence(qc_gene, gff):
        ffn = QC.check_no_sequence(qc_gene, gff)
        assert ffn == ['GCF_004088235v2_20191001.ffn']
        qc_gene = "tuf"
        ffn = QC.check_no_sequence(qc_gene, gff)
        assert ffn == ['GCF_004088235v2_20191001.ffn']

        for files in os.listdir(QC.ffn_dir):
            if files in ffn:
                if files not in QC.ffn_list:
                    QC.ffn_list.append(files)
        assert len(QC.ffn_list) == 1

    def test_choose_sequence(qc_gene):
        qc_dir = os.path.join(QC.target_dir, qc_gene + "_QC")
        G.create_directory(qc_dir)
        qc_seqs = QC.choose_sequence(qc_gene)
        assert len(qc_seqs) == 1
        assert len(qc_seqs[0]) == 2
        # these values changed due to removal of fasta style ">" and new line
        assert len(qc_seqs[0][0]) == 21
        assert len(qc_seqs[0][1]) == 1568
        return qc_seqs

    def qc_blast(qc_gene):
        from speciesprimer import Blast
        from speciesprimer import BlastPrep
        qc_dir = os.path.join(QC.target_dir, qc_gene + "_QC")
        use_cores, inputseqs = BlastPrep(
                        qc_dir, qc_seqs, qc_gene,
                        QC.config.blastseqs).run_blastprep()
        Blast(
            QC.config, qc_dir, "quality_control"
            ).run_blast(qc_gene, use_cores)

    def test_qc_blast_parser(gi_list=False):
        from speciesprimer import errors
        passed = QC.qc_blast_parser(qc_gene)
        if gi_list is True:
            error_msg = "Error: Less than two genomes survived QC"
            assert passed == [[
                'GCF_004088235v2_00210', '343201711',
                'NR_042437', 'Lactobacillus curvatus',
                'Lactobacillus curvatus', 'passed QC']]
            exstat = QC.check_passed_list(passed, qc_gene)
            assert exstat == 1
            assert errors[-1] == [QC.config.target, error_msg]


        else:
            error_msg = "Error: No genomes survived QC"
            exstat = QC.check_passed_list(passed, qc_gene)
            assert passed == []
            assert exstat == 1
            assert errors[-1] == [QC.config.target, error_msg]

        qc_dir = os.path.join(QC.target_dir, qc_gene + "_QC")
        corrfile = os.path.join(testfiles_dir, "rRNA_0_results_err.xml")
        errfile = os.path.join(qc_dir, "rRNA_0_results.xml")
        if os.path.isfile(errfile):
            os.remove(errfile)
        shutil.copy(corrfile, errfile)
        with pytest.raises(Exception):
            passed = QC.qc_blast_parser(qc_gene)
        assert os.path.isfile(errfile) is False

    def qc_blast_fail(qc_gene):
        qc_dir = os.path.join(QC.target_dir, qc_gene + "_QC")
        if os.path.isdir(qc_dir):
            shutil.rmtree(qc_dir)
        G.create_directory(qc_dir)
        QC.config.customdb = "/tmp/not_a_db.fas"
        from speciesprimer import Blast
        from speciesprimer import BlastPrep
        qc_dir = os.path.join(QC.target_dir, qc_gene + "_QC")
        use_cores, inputseqs = BlastPrep(
                        qc_dir, qc_seqs, qc_gene,
                        QC.config.blastseqs).run_blastprep()
        Blast(
            QC.config, qc_dir, "quality_control"
            ).run_blast(qc_gene, use_cores)

        QC.config.customdb = os.path.join(tmpdir, "customdb.fas")

    def test_qc_blast_parser_fail():
        with pytest.raises(Exception):
            QC.qc_blast_parser(qc_gene)
        filename = (
            '/primerdesign/test/Lactobacillus_curvatus/rRNA_QC/'
            'rRNA_0_results.xml')
        error_msg = (
                "A problem with the BLAST results file " + filename +
                " was detected. Please check"
                " if the file was removed and start the run again")
        from speciesprimer import errors
        assert errors[-1] == ['Lactobacillus_curvatus', error_msg]
        qc_dir = os.path.join(QC.target_dir, qc_gene + "_QC")
        if os.path.isdir(qc_dir):
            shutil.rmtree(qc_dir)
        G.create_directory(qc_dir)

    def test_DBError():
        from speciesprimer import Blast
        from speciesprimer import BlastPrep
        infile = os.path.join(testfiles_dir, "customdb_err.fas")
        create_customblastdb(infile)
        qc_dir = os.path.join(QC.target_dir, qc_gene + "_QC")
        use_cores, inputseqs = BlastPrep(
                        qc_dir, qc_seqs, qc_gene,
                        QC.config.blastseqs).run_blastprep()
        Blast(
            QC.config, qc_dir, "quality_control"
            ).run_blast(qc_gene, use_cores)

        with pytest.raises(Exception):
            QC.qc_blast_parser(qc_gene)

        from speciesprimer import errors
        errstart = errors[-1][1]
        assert errstart[0:29] == "Error: No definition line in "

    def test_remove_qc_failures():
        gen_dir = os.path.join(
                QC.config.path, "excludedassemblies", "Lactobacillus_curvatus",
                "genomic_fna")
        G.create_directory(gen_dir)
        infile = os.path.join(
            testfiles_dir, "GCF_002224565.1_ASM222456v1_genomic.fna")
        to_dir = os.path.join(
            QC.ex_dir, "genomic_fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna")
        shutil.copy(infile, to_dir)
        delete = QC.remove_qc_failures(qc_gene)
        fna_files = os.path.join(QC.ex_dir, "fna_files")
        genomic_files = os.path.join(QC.ex_dir, "genomic_fna")
        sumfile = os.path.join(QC.ex_dir, "excluded_list.txt")
        fna = []
        genomic = []
        sumf = []
        for files in os.listdir(fna_files):
            fna.append(files)
        for files in os.listdir(genomic_files):
            genomic.append(files)
        with open(sumfile) as f:
            for line in f:
                sumf.append(line.strip())
        delete.sort()
        assert delete == ['GCF_004088235v1', 'GCF_maxcontigs', 'GCF_noseq']
        fna.sort()
        assert fna == [
            'GCF_004088235v1_20191001.fna',
            'GCF_maxcontigs_date.fna',
            'GCF_noseq_date.fna']
        assert genomic == ["GCF_002224565.1_ASM222456v1_genomic.fna"]
        assert sumf == ['GCF_noseq', 'GCF_maxcontigs', 'GCF_004088235v1']

    test_get_excluded_gis()
    gff_list = test_search_qc_gene()
    gff_list = test_count_contigs(gff_list, QC.contiglimit)
    gff = test_identify_duplicates(gff_list)
    ffn_list = test_check_no_sequence(qc_gene, gff)
    qc_seqs = test_choose_sequence(qc_gene)
    qc_blast(qc_gene)
    test_qc_blast_parser()
    if os.path.isdir(QC.ex_dir):
        shutil.rmtree(QC.ex_dir)

    qc_blast_fail(qc_gene)
    test_qc_blast_parser_fail()

    # Remove Fake species by GI
    QC = QualityControl(config)
    prepare_QC_testfiles(config)
    prepare_GI_list()
    gff_list = test_search_qc_gene()
    gff_list = test_count_contigs(gff_list, QC.contiglimit)
    gff = test_identify_duplicates(gff_list)
    ffn_list = test_check_no_sequence(qc_gene, gff)
    qc_seqs = test_choose_sequence(qc_gene)
    qc_blast(qc_gene)
    test_qc_blast_parser(gi_list=True)
    if os.path.isdir(QC.ex_dir):
        shutil.rmtree(QC.ex_dir)
    prepare_QC_testfiles(config)
    test_remove_qc_failures()

    test_DBError()

    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    if os.path.isdir(QC.ex_dir):
        shutil.rmtree(QC.ex_dir)


def test_skip_pangenome_analysis(config):
    from speciesprimer import PangenomeAnalysis
    PA = PangenomeAnalysis(config)
    G.create_directory(PA.pangenome_dir)
    fromfile = os.path.join(testfiles_dir, "gene_presence_absence.csv")
    tofile = os.path.join(PA.pangenome_dir, "gene_presence_absence.csv")
    if os.path.isfile(tofile):
        os.remove(tofile)
    shutil.copy(fromfile, tofile)
    exitstat = PA.run_pangenome_analysis()
    assert exitstat == 2


def test_CoreGenes(config):
    from speciesprimer import CoreGenes
    CG = CoreGenes(config)

    def prepare_tests():
        if os.path.isdir(CG.ffn_dir):
            shutil.rmtree(CG.ffn_dir)
        new_ffn_dir = os.path.join(testfiles_dir, "ffn_files")
        shutil.copytree(new_ffn_dir, CG.ffn_dir)

    def test_get_singlecopy_genes(config):
        coregenesummary = CG.get_singlecopy_genes
        assert coregenesummary == [8, 14, 2, 8]

    def test_coregene_extract(config):
        G.create_directory(CG.results_dir)
        fasta_dir = os.path.join(CG.results_dir, "fasta")
        G.create_directory(fasta_dir)
        CG.run_CoreGenes()
        ref_dir = os.path.join(ref_data, "fasta")
        fasta_dir = os.path.join(CG.results_dir, "fasta")
        compare_ref_files(fasta_dir, ref_dir)

    prepare_tests()
    test_coregene_extract(config)

def test_CoreGeneSequences(config):
    from speciesprimer import CoreGeneSequences
    from speciesprimer import BlastParser
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    config.customdb = os.path.join(tmpdir, "customdb.fas")
    config.blastdbv5 = False
    CGS = CoreGeneSequences(config)

    def test_seq_alignments():
        # skip this test because of variation in alignments
        # CGS.seq_alignments()
        # results_dir = CGS.alignments_dir
        ref_dir = os.path.join(ref_data, "alignments")
        # compare_ref_files(results_dir, ref_dir)
        shutil.copytree(ref_dir, CGS.alignments_dir)

    def test_seq_consenus():
        CGS.seq_consensus()
        ref_dir = os.path.join(ref_data, "consensus")
        results_dir = CGS.consensus_dir
        compare_ref_files(results_dir, ref_dir)

    def test_conserved_seqs():
        cons_summary = os.path.join(
            CGS.consensus_dir, "consensus_summary.txt")
        with open(cons_summary, "w") as f:
            f.write("")
        conserved_sequences = CGS.conserved_seqs()
        ref_file = os.path.join(ref_data, "conserved_seqs.fas")
        assert conserved_sequences == 1
        os.remove(cons_summary)
        CGS.conserved_seqs()
        result_file = os.path.join(CGS.blast_dir, "Lb_curva_conserved")
        compare_ref_files(result_file, ref_file)

    def create_mockdict():
        ref_file = os.path.join(ref_data, "conserved_seqs.fas")
        ref_dataset = []
        with open(ref_file, "r") as f:
            datapoint = []
            for index, line in enumerate(f):
                if index % 2 == 0:
                    datapoint.append(line.split(">")[1].strip())
                else:
                    datapoint.append(line.strip())
                    ref_dataset.append(datapoint)
                    datapoint = []
        tmp_dict = {}
        for item in ref_dataset:
            tmp_dict.update({item[0]: item[1]})

        return tmp_dict

    def test_noconservedseqs():
        sakei = []
        with open(os.path.join(CGS.blast_dir, "nontargethits.json")) as f:
            for line in f:
                non_dict = json.loads(line)
        for key in non_dict.keys():
            try:
                newval = non_dict[key]["Lactobacillus sakei"]["main_id"]
                if newval not in sakei:
                    sakei.append(newval)
            except KeyError:
                pass
        for key in non_dict.keys():
            for i, item in enumerate(sakei):
                try:
                    if not (
                        item == non_dict[key]["Lactobacillus sakei"]["main_id"]
                    ):
                        non_dict[key].update(
                            {
                                "Lactobacillus sakei_" + str(i + 2):
                                    {"main_id": item, "gi_ids": []}})
                except KeyError:
                    non_dict[key].update(
                        {
                            "Lactobacillus sakei_" + str(i + 2):
                            {"main_id": item, "gi_ids": []}})
        with open(os.path.join(CGS.blast_dir, "nontargethits.json"), "w") as f:
            f.write(json.dumps(non_dict))

    def test_run_coregeneanalysis(config):
        G.create_directory(tmpdir)

        def create_customblastdb(config):
            infile = os.path.join(testfiles_dir, "conserved_customdb.fas")
            cmd = [
                "makeblastdb", "-in", infile, "-parse_seqids", "-title",
                "mockconservedDB", "-dbtype", "nucl", "-out", config.customdb]
            G.run_subprocess(
                cmd, printcmd=False, logcmd=False, printoption=True)

        create_customblastdb(config)
        conserved_seq_dict = CGS.run_coregeneanalysis()
        shutil.rmtree(tmpdir)
        return conserved_seq_dict

    test_seq_alignments()
    test_seq_consenus()
    test_conserved_seqs()
    tmp_dict = create_mockdict()
    conserved_seq_dict = test_run_coregeneanalysis(config)
    assert conserved_seq_dict == tmp_dict

    conserved = BlastParser(
            config).run_blastparser(conserved_seq_dict)
    assert conserved == 0
    test_noconservedseqs()
    conserved = BlastParser(
            config).run_blastparser(conserved_seq_dict)
    assert conserved == 1
    nt_file = os.path.join(CGS.blast_dir, "nontargethits.json")
    if os.path.isfile(nt_file):
        os.remove(nt_file)
    conserved_seq_dict = test_run_coregeneanalysis(config)
    conserved = BlastParser(
            config).run_blastparser(conserved_seq_dict)

def test_PrimerDesign(config):
    reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
    from speciesprimer import PrimerDesign
    pd = PrimerDesign(config)
    G.create_directory(pd.primer_dir)
    p3_output = os.path.join(pd.primer_dir, "primer3_output")
    with pytest.raises(Exception):
        settings_file = os.path.join(BASE_PATH, "p3parameters")
        errorsettings = os.path.join(BASE_PATH, "p3parameters_none")
        try:
            os.rename(settings_file, errorsettings)
            pd.run_primerdesign()
        finally:
            os.rename(errorsettings, settings_file)

    if os.path.isfile(p3_output):
        os.remove(p3_output)

    pd.run_primer3()

    assert os.path.isfile(p3_output) is True
    pd.run_primerdesign()

    with open(reffile) as f:
        for line in f:
            refdict = json.loads(line)

    assert refdict == pd.p3dict

    # test primer3 error
    p3_error = os.path.join(testfiles_dir, "primer3_output_err")
    pd.p3dict = {}
    pd.parse_Primer3_output(p3_error)
    errorreport = os.path.join(pd.primer_dir, "primer3_errors.csv")
    assert os.path.isfile(errorreport) is True
    with open(errorreport) as f:
        reader = csv.reader(f)
        for row in reader:
            assert row == [
                'SEQUENCE_ID=yfnB_2',
                'SEQUENCE_TEMPLATE='
                'GCCAANACGCAATATCGGCGGTTACAAGATTCAGGATTAATCACATATT',
                'PRIMER_PRODUCT_SIZE_RANGE=70-200',
                'PRIMER_ERROR=SEQUENCE_INCLUDED_REGION length'
                ' < min PRIMER_PRODUCT_SIZE_RANGE']


def test_BlastDBError(config):
    from speciesprimer import CoreGeneSequences
    from speciesprimer import BlastParser
    from speciesprimer import errors
    from basicfunctions import BlastDBError
    dbpath = os.path.join(tmpdir, "customdb.fas")
    blapa = BlastParser(config)
    nontargetfile = os.path.join(blapa.blast_dir, "nontargethits.json")
    blastresult = os.path.join(blapa.blast_dir, "conserved_0_results.xml")
    primer3in = os.path.join(blapa.results_dir, "primer3_input")
    consfile = os.path.join(testfiles_dir, "consensus_summary.txt")
    consens = os.path.join(
            blapa.results_dir, "consensus", "consensus_summary.txt")
    if os.path.isfile(consens):
        os.remove(consens)
    shutil.copy(consfile, consens)

    config.customdb = dbpath
    config.blastdbv5 = False

    def dbinputfiles(dbpath, nodesc=False):
        filenames = [
            "GCF_004088235v1_20191001.fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna"]
        with open(dbpath, "w") as f:
            for filename in filenames:
                filepath = os.path.join(testfiles_dir, filename)
                records = SeqIO.parse(filepath, "fasta")
                for record in records:
                    if nodesc:
                        record.description = ""
                    else:
                        if record.id == record.description:
                            description = (
                                record.id +
                                " Lactobacillus curvatus strain SRCM103465")
                            record.description = description
                    SeqIO.write(record, f, "fasta")
            mockseqs = os.path.join(testfiles_dir, "mocktemplate.seqs")
            records = list(SeqIO.parse(mockseqs, "fasta"))
            for record in records:
                if nodesc:
                    record.description = ""
                SeqIO.write(record, f, "fasta")

        return dbpath

    def create_customblastdb(dbpath, noseq=False):
        if noseq:
            cmd = [
                "makeblastdb", "-in", dbpath, "-title",
                "mockconservedDB", "-dbtype", "nucl", "-out", dbpath]
        else:
            cmd = [
                "makeblastdb", "-in", dbpath, "-parse_seqids", "-title",
                "mockconservedDB", "-dbtype", "nucl", "-out", dbpath]

        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    def remove_blastresults():
        if os.path.isfile(nontargetfile):
            os.remove(nontargetfile)
        if os.path.isfile(blastresult):
            os.remove(blastresult)
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)

    print("\n>>> Start normal customDB:\n")
    remove_blastresults()
    G.create_directory(os.path.dirname(dbpath))
    dbinputfiles(dbpath)
    create_customblastdb(dbpath)
    conserved_seq_dict = CoreGeneSequences(config).run_coregeneanalysis()
    blapa.run_blastparser(conserved_seq_dict)
    assert os.path.isfile(primer3in) is True

    print("\n>>> Start custom DB no parse_seqids option:\n")
    remove_blastresults()
    G.create_directory(os.path.dirname(dbpath))
    dbinputfiles(dbpath)
    create_customblastdb(dbpath, noseq=True)
    conserved_seq_dict = CoreGeneSequences(config).run_coregeneanalysis()
    with pytest.raises(BlastDBError):
        blapa.run_blastparser(conserved_seq_dict)

    error_msg = (
        "Problem with custom DB, Please use the '-parse_seqids'"
        " option for the makeblastdb command")
    assert errors[-1] == [config.target, error_msg]

    print("\n>>> Start missing description in BLASTDB:\n")
    remove_blastresults()
    G.create_directory(os.path.dirname(dbpath))
    dbinputfiles(dbpath, nodesc=True)
    create_customblastdb(dbpath)
    conserved_seq_dict = CoreGeneSequences(config).run_coregeneanalysis()
    with pytest.raises(BlastDBError):
        blapa.run_blastparser(conserved_seq_dict)
    errstart = errors[-1][1]
    assert errstart[0:29] == "Error: No definition line in "

    remove_blastresults()


def test_blastprep(config):
    from speciesprimer import BlastPrep
    testconditions = [[100, 50], [500, 10], [1000, 5], [2000, 3], [5000, 1]]
    listlengths = [50*[100], 10*[500], 5*[1000], [1667, 1667, 1666], [5000]]
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    G.create_directory(tmpdir)
    input_list = []
    for i in range(0, 5000):
        i = i + 1
        elem = ["seq_" + str(i), i * "A"]
        input_list.append(elem)
    name = "test"
    directory = tmpdir
    for index, [maxpart, dicts] in enumerate(testconditions):
        maxpart, dicts = [maxpart, dicts]
        bp = BlastPrep(directory, input_list, name, maxpartsize=maxpart)
        bp.create_listdict()
        assert len(bp.list_dict) == dicts
        bp.get_equalgroups()
        print("inputlist maxpartsize dicts")
        print(len(input_list), maxpart, dicts,)
        for i, value in enumerate(bp.list_dict.values()):
            assert len(value) == listlengths[index][i]

    bp = BlastPrep(directory, input_list, name, maxpartsize=1000)
    bp.create_listdict()
    bp.get_equalgroups()
    inlist = bp.write_blastinput()
    reflist = []
    rangelist = [5001, 5000, 4999, 4998, 4997]
    for i, r in enumerate(rangelist):
        for x in reversed(range(5-i, r, 5)):
            if x not in reflist:
                reflist.append("seq_" + str(x))
    assert reflist == inlist


def test_BLASTsettings(config):
    import multiprocessing
    from speciesprimer import Blast
    bl = Blast(config, tmpdir, "test")
    blastfiles = bl.search_blastfiles(tmpdir)
    assert len(blastfiles) == 5
    blastfile = "test.part-0"
    modes = ["quality_control", "conserved", "primer"]

    db_settings = [
        [False, None],
        [True, None],
        [False, os.path.join(tmpdir, "customdb.fas")],
        [True, os.path.join(tmpdir, "customdb_v5.fas")]]
    db_outcome = [
            "nt", "nt_v5",
            os.path.join(tmpdir, "customdb.fas"),
            os.path.join(tmpdir, "customdb_v5.fas")]
    for i, settings in enumerate(db_settings):
        config.blastdbv5 = settings[0]
        config.customdb = settings[1]
        for mode in modes:
            bl = Blast(config, tmpdir, mode)
            cores = cores = multiprocessing.cpu_count()
            cmd = bl.get_blast_cmd(
                    blastfile, blastfile + "_results.xml", cores)
            if mode == "quality_control":
                assert cmd[2] == "megablast"
            elif mode == "conserved":
                assert cmd[2] == "dc-megablast"
            elif mode == "primer":
                assert cmd[2] == "blastn-short"
            else:
                # there should not be any other modes
                assert cmd[2] == "error"
            assert cmd[-1] == db_outcome[i]
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)


def test_blastparser(config):
    corrfile = os.path.join(testfiles_dir, "conserved_0_results_err.xml")
    from speciesprimer import BlastParser
    from speciesprimer import CoreGeneSequences
    blapa = BlastParser(config)
    nontarg = os.path.join(blapa.blast_dir, "nontargethits.json")
    errfile = os.path.join(blapa.blast_dir, "conserved_0_results.xml")
    if os.path.isfile(nontarg):
        os.remove(nontarg)
    if os.path.isfile(errfile):
        os.remove(errfile)
    shutil.copy(corrfile, errfile)
    conserved_seq_dict = CoreGeneSequences(config).run_coregeneanalysis()
    with pytest.raises(Exception):
        blapa.run_blastparser(conserved_seq_dict)
    assert os.path.isfile(errfile) is False


def prepare_test(config):
    testdir = os.path.join("/", "primerdesign", "test")
    target_dir = os.path.join(testdir, "Lactobacillus_curvatus")
    pangenome_dir = os.path.join(target_dir, "Pangenome")
    results_dir = os.path.join(pangenome_dir, "results")
    primer_dir = os.path.join(results_dir, "primer")
    primerblast_dir = os.path.join(primer_dir, "primerblast")
    primer_qc_dir = os.path.join(primer_dir, "primer_QC")
    dbpath = os.path.join(tmpdir, "customdb.fas")

    def dbinputfiles():
        filenames = [
            "GCF_004088235v1_20191001.fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna"]
        with open(dbpath, "w") as f:
            for filename in filenames:
                filepath = os.path.join(testfiles_dir, filename)
                records = SeqIO.parse(filepath, "fasta")
                for record in records:
                    if record.id == record.description:
                        description = (
                            record.id
                            + " Lactobacillus curvatus strain SRCM103465")
                        record.description = description
                    SeqIO.write(record, f, "fasta")
        return dbpath

    def create_customblastdb():
        cmd = [
            "makeblastdb", "-in", dbpath, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", dbpath]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    G.create_directory(os.path.dirname(dbpath))
    dbinputfiles()
    create_customblastdb()

    from speciesprimer import PrimerQualityControl
    from speciesprimer import Blast
    from speciesprimer import BlastPrep
    G.create_directory(primer_dir)

    config.customdb = dbpath
    config.blastdbv5 = False
    reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
    with open(reffile) as f:
        for line in f:
            primer3dict = json.loads(line)
    pqc = PrimerQualityControl(config, primer3dict)
    pqc.collect_primer()
    G.create_directory(pqc.primerblast_dir)
    file_path = os.path.join(primerblast_dir, "nontargethits.json")
    if os.path.isfile(file_path):
        os.remove(file_path)
    prep = BlastPrep(
        pqc.primerblast_dir, pqc.primerlist,
        "primer", pqc.config.blastseqs)
    use_cores, inputseqs = prep.run_blastprep()

    Blast(pqc.config, pqc.primerblast_dir, "primer").run_blast(
        "primer", use_cores)

    return primer_dir, primerblast_dir, primer_qc_dir


def test_primerblastparser(config):
    from speciesprimer import BlastParser
    primer_dir, primerblast_dir, primer_qc_dir = prepare_test(config)
    if os.path.isdir(primer_qc_dir):
        shutil.rmtree(primer_qc_dir)
    file_path = os.path.join(primerblast_dir, "nontargethits.json")
    if os.path.isfile(file_path):
        os.remove(file_path)
    dbpath = os.path.join(tmpdir, "customdb.fas")
    config.customdb = dbpath
    config.blastdbv5 = False
    blapa = BlastParser(config, "primer")
    align_dict = blapa.blast_parser(blapa.primerblast_dir)

    if not os.path.isfile(file_path):
        nonred_dict = blapa.remove_redundanthits(align_dict)
        blapa.write_nontargethits(
                            blapa.primerblast_dir, nonred_dict,
                            "json")
    else:
        nonred_dict = align_dict

    G.create_directory(blapa.primer_qc_dir)
    nonreddata = blapa.sort_nontarget_sequences(nonred_dict)
    assert len(nonreddata) == 97
    blapa.write_nontarget_sequences(nonreddata)
    nontarget_seqs = os.path.join(primer_qc_dir, "BLASTnontarget0.sequences")
    DBIDS = os.path.join(primer_qc_dir, "primerBLAST_DBIDS.csv")
    assert os.path.isfile(nontarget_seqs) is True
    assert os.path.isfile(DBIDS) is True
    count = 0
    with open(nontarget_seqs) as f:
        for line in f:
            if ">" in line:
                count += 1
    assert count == 97
    # skip
    exitstatus = blapa.get_primerBLAST_DBIDS(nonred_dict)
    assert exitstatus == 0
    assert os.path.isfile(nontarget_seqs) is True
    assert os.path.isfile(DBIDS) is True


def test_primerblastparser_exceptions(config):
    from speciesprimer import BlastParser
    primer_dir, primerblast_dir, primer_qc_dir = prepare_test(config)
    if os.path.isdir(primer_qc_dir):
        shutil.rmtree(primer_qc_dir)
    file_path = os.path.join(primerblast_dir, "nontargethits.json")
    if os.path.isfile(file_path):
        os.remove(file_path)
    dbpath = os.path.join(tmpdir, "customdb.fas")
    config.customdb = dbpath
    config.blastdbv5 = False
    config.exception = [
        'Lactobacillus_sakei', 'Lactobacillus_sakei']
    blapa = BlastParser(config, "primer")
    align_dict = blapa.blast_parser(blapa.primerblast_dir)

    if not os.path.isfile(file_path):
        nonred_dict = blapa.remove_redundanthits(align_dict)
        blapa.write_nontargethits(
                            blapa.primerblast_dir, nonred_dict,
                            "json")
    else:
        nonred_dict = align_dict

    G.create_directory(blapa.primer_qc_dir)
    exitstatus = blapa.get_primerBLAST_DBIDS(nonred_dict)
    assert exitstatus == 1
    nonreddata = blapa.sort_nontarget_sequences(nonred_dict)
    assert len(nonreddata) == 0
    nontarget_seqs = os.path.join(primer_qc_dir, "BLASTnontarget0.sequences")
    DBIDS = os.path.join(primer_qc_dir, "primerBLAST_DBIDS.csv")
    assert os.path.isfile(nontarget_seqs) is False
    assert os.path.isfile(DBIDS) is False

    if os.path.isdir(os.path.dirname(dbpath)):
        shutil.rmtree(os.path.dirname(dbpath))
    shutil.rmtree(primer_dir)
#    write_primer3_input(self, selected_seqs, conserved_seq_dict)
#    get_alignmentdata(self, alignment)


def test_PrimerQualityControl_specificitycheck(config):
    from speciesprimer import PrimerQualityControl
    from speciesprimer import BlastPrep
    from speciesprimer import Blast

    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    G.create_directory(tmpdir)
    config.customdb = os.path.join(tmpdir, "primer_customdb.fas")
    config.blastdbv5 = False

    def dbinputfiles():
        filenames = [
            "GCF_004088235v1_20191001.fna",
            "GCF_002224565.1_ASM222456v1_genomic.fna"]
        dbfile = os.path.join(tmpdir, "primer_customdb.fas")
        with open(dbfile, "w") as f:
            for filename in filenames:
                filepath = os.path.join(testfiles_dir, filename)
                records = SeqIO.parse(filepath, "fasta")
                for record in records:
                    if record.id == record.description:
                        description = (
                            record.id
                            + " Lactobacillus curvatus strain SRCM103465")
                        record.description = description
                    SeqIO.write(record, f, "fasta")
        return dbfile

    def create_customblastdb(config, infile):
        cmd = [
            "makeblastdb", "-in", infile, "-parse_seqids", "-title",
            "mockconservedDB", "-dbtype", "nucl", "-out", config.customdb]
        G.run_subprocess(
            cmd, printcmd=False, logcmd=False, printoption=False)

    def test_collect_primer(config):
        pqc = PrimerQualityControl(config, {})
        # returns 1 if len(primerlist) == 0
        exitstat = pqc.collect_primer()
        assert exitstat == 1
        reffile = os.path.join(testfiles_dir, "ref_primer3_summary.json")
        with open(reffile) as f:
            for line in f:
                primer3dict = json.loads(line)
        config.exception = [
            'Bacterium curvatum', 'Lactobacillus curvatus subsp. curvatus',
            'Lactobacillus sp. N55', 'Lactobacillus sp. N61']
        pqc = PrimerQualityControl(config, primer3dict)

        exitstat = pqc.collect_primer()
        item = ['comFA_5', 'Primer_pair_7', 11.374516]
        # test self.primerlist
        pqc.get_blast_input(item)
        assert pqc.primerlist[-2] == [
                'Lb_curva_comFA_5_P7_F', 'ACAACGCTTATTATTATTTGTGCCA']
        assert pqc.primerlist[-1] == [
                'Lb_curva_comFA_5_P7_R', 'AAAGGCCGCTATCTTGTCTAAT']
        del pqc.primerlist[-1]
        del pqc.primerlist[-1]

        return pqc

    def primerBLAST(config):
        dbfile = dbinputfiles()
        create_customblastdb(config, dbfile)
        if os.path.isfile(dbfile):
            os.remove(dbfile)

        G.create_directory(pqc.primerblast_dir)
        prep = BlastPrep(
            pqc.primerblast_dir, pqc.primerlist,
            "primer", pqc.config.blastseqs)
        use_cores, inputseqs = prep.run_blastprep()
#        Blast(pqc.config, pqc.primerblast_dir, "primer").run_blast(
#            "primer", use_cores)
        G.create_directory(pqc.primer_qc_dir)
        reffile = os.path.join(testfiles_dir, "primer_nontargethits.json")
        tofile = os.path.join(pqc.primerblast_dir, "nontargethits.json")
        shutil.copy(reffile, tofile)
        pqc.call_blastparser.run_blastparser("primer")

        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)

        return inputseqs

    def test_get_primerinfo(inputseqs):
        ref1 = [[
            'Lb_curva_asnS_1_P0_F',
            'AAACCATGTCGATGAAGAAGTTAAA',
            'Lb_curva_asnS_1_P0_R',
            'TGCCATCACGTAATTGGAGGA',
            'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGG'
            + 'AAGATTTCATTCCTCCAATTACGTGATGGCACTGCCTTTTTCCAAGGGGTTGTCGTTAA'
            + 'AAG']]
        ref2 = [[
            'asnS_1', 'Primer_pair_0',
            'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGG'
            + 'AAGATTTCATTCCTCCAATTACGTGATGGCA', 'Lb_curva_asnS_1_P0']]
        inputseqs.sort()
        modes = ['mfeprimer', "mfold", "dimercheck"]
        for mode in modes:
            primerinfo = pqc.get_primerinfo([inputseqs[0]], mode)
            if mode == modes[0]:
                assert primerinfo == ref1
            elif mode == modes[1]:
                assert primerinfo == ref2
            elif mode == modes[2]:
                assert primerinfo == [[
                    'Lb_curva_asnS_1_P0',
                    'AAACCATGTCGATGAAGAAGTTAAA',
                    'TGCCATCACGTAATTGGAGGA']]

    def test_MFEprimer_DBs():
        for files in os.listdir(pqc.primer_qc_dir):
            if (
                files.startswith("BLASTnontarget")
                and files.endswith(".sequences")
            ):
                pqc.dbinputfiles.append(os.path.join(pqc.primer_qc_dir, files))
        qc_file = os.path.join(testfiles_dir, "Lb_curva_qc_sequences.csv")
        if os.path.isdir(pqc.summ_dir):
            shutil.rmtree(pqc.summ_dir)
        G.create_directory(pqc.summ_dir)
        shutil.copy(qc_file, pqc.summ_dir)
        primerinfo = pqc.get_primerinfo(inputseqs, "mfeprimer")
        pqc.prepare_MFEprimer_Dbs(primerinfo)
        dbtype = [
            "BLASTnontarget0.sequences",
            "template.sequences",
            "Lb_curva.genomic"]
        dbfiles = [
            ".2bit",
            ".sqlite3.db",
            ".uni"]
        for db in dbtype:
            for end in dbfiles:
                files = db + end
                print(files)
                filepath = os.path.join(pqc.primer_qc_dir, files)
                assert os.path.isfile(filepath) is True
                assert os.stat(filepath).st_size > 0

    def test_MFEprimer(inputseqs):
        import itertools
        primerinfo = pqc.get_primerinfo([inputseqs[0]], "mfeprimer")
        os.chdir(pqc.primer_qc_dir)
        # template
        print("Test MFEprimer_template()")
        ppcth = [100, 90, 80, 70]
        outcome = [
            [None], [None], 'Lb_curva_asnS_1_P0_F', 'Lb_curva_asnS_1_P0_F']
        results = []
        for index, item in enumerate(ppcth):
            pqc.mfethreshold = item
            result = P.MFEprimer_template(
                    primerinfo[0], [pqc.primer_qc_dir, pqc.mfethreshold])
            results.append(result)
            if outcome[index] == [None]:
                assert result[0] == outcome[index]
            else:
                assert result[0][0] == outcome[index]
        check_nontarget = pqc.write_MFEprimer_results(results, "template")

        # nontarget
        print("Test MFEprimer_nontarget()")
        ref = [[
            'Lb_curva_asnS_1_P0_F',
            'AAACCATGTCGATGAAGAAGTTAAA',
            'Lb_curva_asnS_1_P0_R',
            'TGCCATCACGTAATTGGAGGA',
            'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGGAAGA'
            + 'TTTCATTCCTCCAATTACGTGATGGCACTGCCTTTTTCCAAGGGGTTGTCGTTAAAAG',
            17.700000000000003],
            [
                'AmpID\tFpID\tRpID\tHitID\tPPC\tSize\tAmpGC\tFpTm\tRpTm\tFpDg'
                '\tRpDg\tBindingStart\tBindingStop\tAmpSeq',
                '1\tLb_curva_asnS_1_P0_F\tLb_curva_asnS_1_P0_F'
                '\tNZ_CP020459.1_1524567_1528567\t-16.4\t364\t37.09\t7.26'
                '\t15.22\t-4.51\t-6.20\t3075\t3409\t'
                'aaaccatgtcgatgAAGAAGTTAAAaaaccagctgaattagacaaatacctaaagagtat'
                'tgactaattcattacaaaaaaagatcccgtcaatgacgagatctttttttatgctcaatt'
                'actaaccgcgatggtcagcacccttcacttaaacgtgcttctcgaactgccttttccggt'
                'taaccacaaaatgggaatcatggctaaacccgccataatccccattttcttgttaatcct'
                'caaagccaacccgttcgagtcagcactttatttgttttctttcttttcgctttgctttaa'
                'ttcttcacccaatttgatgatgtatttcttcaaatcatcTTTAACTTCTtcatcgacatg'
                'gttt']]
        results = []
        nontarget_lists = []
        for index, item in enumerate(ppcth):
            pqc.mfethreshold = item
            dbfile = pqc.dbinputfiles[0]
            result = P.MFEprimer_nontarget(
                    check_nontarget[1], [dbfile, pqc.primer_qc_dir])
            assert result == ref
            results.append(result)
        nontarget_lists = list(itertools.chain(nontarget_lists, results))
        check_assembly = pqc.write_MFEprimer_results(
                                                nontarget_lists, "nontarget")

        # assembly
        print("Test MFEprimer_assembly()")
        ref = [[
            'Lb_curva_asnS_1_P0_F',
            'AAACCATGTCGATGAAGAAGTTAAA',
            'Lb_curva_asnS_1_P0_R',
            'TGCCATCACGTAATTGGAGGA',
            'AAACCATGTCGATGAAGAAGTTAAAATTGGCGTTTGGTTAACCGACAAACGNTCAAGTGGGAAGA'
            'TTTCATTCCTCCAATTACGTGATGGCACTGCCTTTTTCCAAGGGGTTGTCGTTAAAAG',
            17.700000000000003]]
        results = []
        outcome = [[None], [None], [None], 'Lb_curva_asnS_1_P0_F']
        dbfile = dbfile = H.abbrev(pqc.target) + ".genomic"
        for index, item in enumerate(ppcth):
            result = P.MFEprimer_assembly(
                    check_assembly[0], [pqc.primer_qc_dir, dbfile, item])

            results.append(result)
            if outcome[index] == [None]:
                assert result[0] == outcome[index]
            else:
                assert result[0][0] == outcome[index]

        check_final = pqc.write_MFEprimer_results(results, "assembly")
        assert check_final == ref

    def test_fullMFEprimer_run():
        outcome = [15, 42, 54, 65]
        ppcth = [100, 90, 80, 70]
        primerinfos = pqc.get_primerinfo(inputseqs, "mfeprimer")
        for index, item in enumerate(ppcth):
            pqc.mfethreshold = item
            primer_name_list = pqc.MFEprimer_QC(primerinfos)
            assert len(primer_name_list) == outcome[index]

    pqc = test_collect_primer(config)
    inputseqs = primerBLAST(config)
    test_get_primerinfo(inputseqs)
    test_MFEprimer_DBs()
    test_MFEprimer(inputseqs)
    test_fullMFEprimer_run()

    tmpdict = os.path.join(pqc.primer_qc_dir, "tmp_p3dict.json")
    with open(tmpdict, "w") as f:
        f.write(json.dumps(pqc.primer3_dict))


def test_PrimerQualityControl_structures(config):
    from speciesprimer import PrimerQualityControl
    pqc = PrimerQualityControl(config, {})
    tmpdict = os.path.join(pqc.primer_qc_dir, "tmp_p3dict.json")
    with open(tmpdict) as f:
        for line in f:
            primer3dict = json.loads(line)
    pqc = PrimerQualityControl(config, primer3dict)

    if os.path.isfile(tmpdict):
        os.remove(tmpdict)

    def test_mfold():
        spec_output = [
            'Lb_curva_g1243_1_P2', 'Lb_curva_g4430_1_P0',
            'Lb_curva_gshAB_2_P0', 'Lb_curva_comFA_2_P0']
        ref_output = [
             'Lb_curva_gshAB_2_P0', 'Lb_curva_g4430_1_P0',
             'Lb_curva_g1243_1_P2']
        ref_output.sort()
        primerinfos = pqc.get_primerinfo(spec_output, "mfold")
        pqc.mfold_analysis(primerinfos)
        pqc.config.mfold = -1.0
        selected_primer, excluded_primer = pqc.mfold_parser()
        assert len(selected_primer) == 3
        assert len(excluded_primer) == 1

        pqc.config.mfold = -12.0
        selected_primer, excluded_primer = pqc.mfold_parser()
        assert len(selected_primer) == 4
        assert len(excluded_primer) == 0
        # default
        pqc.config.mfold = -3.0
        selected_primer, excluded_primer = pqc.mfold_parser()
        selected_primer.sort()
        assert selected_primer == ref_output
        assert excluded_primer == ['Lb_curva_comFA_2_P0']

    def test_dimercheck():
        selected_primer = [
            'Lb_curva_g1243_2_P0', 'Lb_curva_comFA_5_P1',
            'Lb_curva_g4295_1_P0', 'Lb_curva_g4295_1_P4',
            'Lb_curva_gshAB_2_P0', 'Lb_curva_gshAB_2_P4',
            'Lb_curva_comFA_4_P2', 'Lb_curva_comFA_4_P5',
            'Lb_curva_g4430_1_P0', 'Lb_curva_g1243_1_P0',
            'Lb_curva_g1243_1_P2', 'Lb_curva_comFA_6_P1',
            'Lb_curva_comFA_2_P1', 'Lb_curva_asnS_2_P0']
        excluded_primer = ['Lb_curva_comFA_2_P0']
        ref_choice = [
            'Lb_curva_g1243_2_P0', 'Lb_curva_comFA_5_P1',
            'Lb_curva_g4295_1_P0', 'Lb_curva_g4295_1_P4',
            'Lb_curva_gshAB_2_P4', 'Lb_curva_comFA_4_P2',
            'Lb_curva_comFA_4_P5', 'Lb_curva_g1243_1_P2',
            'Lb_curva_comFA_6_P1', 'Lb_curva_comFA_2_P1']
        ref_choice.sort()
        dimercheck = pqc.dimercheck_primer(selected_primer, excluded_primer)
        choice = pqc.check_primerdimer(dimercheck)
        choice.sort()
        assert choice == ref_choice
        return choice

    test_mfold()
    choice = test_dimercheck()
    total_results = pqc.write_results(choice)
    total_results.sort()
    assert total_results[0] == [
        'Lb_curva_comFA_2_P1', 100.0, 1.64, 'comFA_2', 'TACCAAGCAACAACGCCATG',
        59.4, 0.63, 'ACACACACGCTGCCCATTAG', 60.95, 0.99,
        'None', 'None', 'None', 106, 82.53,
        'TACCAAGCAACAACGCCATGTCTTAACTGCAGTGACCGGTGCTGGTAAAACTGAGATGTTATTTCAAGG'
        + 'CATTGCGACGGCTTTNGCTAATGGGCAGCGTGTGTGT',
        'GCGGTTAGTGGAAGGGGATCAGTTGNTTTACGTGCCTGAATGCAATCTGTTTGAGTCGATNGTAGAACC'
        + 'GCTAACTTGGTCGGGGACACTAACCCCTTTTCAAGCACAAGCCGCAAAGCAGATTGTTGCGGTCATT'
        + 'ACCAAGCAACAACGCCATGTCTTAACTGCAGTGACCGGTGCTGGTAAAACTGAGATGTTATTTCAAG'
        + 'GCATTGCGACGGCTTTNGCTAATGGGCAGCGTGTGTGTGTTGCTGCGCCGCGGGTGGCGGTTTGTTT'
        + 'AGAACTCTATCCGCGCTTGCAAGCAGCGTTTGCTAACACACCAAT']


def test_summary(config):
    total_results = [
        [
            'Lb_curva_comFA_2_P1', 100.0, 1.64, 'comFA_2',
            'TACCAAGCAACAACGCCATG', 59.4, 0.63,
            'ACACACACGCTGCCCATTAG', 60.95, 0.99,
            'None', 'None', 'None', 106, 82.53,
            'TACCAAGCAAC...',
            'TACCAAGCAACAACGCCATGTC...']]

    from speciesprimer import Summary
    su = Summary(config, total_results)
    su.run_summary("normal")
    files = []
    for filename in os.listdir(su.summ_dir):
        files.append(filename)

    ref_files = [
        "Lb_curva_primer.csv", "Lb_curva_qc_sequences_details.csv",
        "mostcommonhits.csv", "Lb_curva_qc_sequences.csv"]
    files.sort()
    ref_files.sort()
    assert files == ref_files

    shutil.rmtree(su.summ_dir)

    su.run_summary("last")
    today = time.strftime("%Y_%m_%d", time.localtime())
    configfile = "config_" + today + ".json"
    stats = "Lb_curva_pipeline_stats_" + today + ".txt"

    ref_files = [
        "Lb_curva_primer.csv", "Lb_curva_qc_sequences_details.csv",
        "mostcommonhits.csv", "Lb_curva_qc_sequences.csv", configfile,
        stats]

    files = []
    for filename in os.listdir(su.summ_dir):
        files.append(filename)

    files.sort()
    ref_files.sort()
    assert files == ref_files


def test_end(config):
    def remove_test_files(config):
        test = config.path
        shutil.rmtree(test)
        tmp_path = os.path.join("/", "pipeline", "tmp_config.json")
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
        os.chdir(BASE_PATH)
        assert os.path.isdir(test) is False

    remove_test_files(config)


if __name__ == "__main__":
    print(msg)
