
# SpeciesPrimer

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://img.shields.io/badge/License-GPLv3-blue.svg)
[![Build Status](https://travis-ci.com/biologger/speciesprimer.svg?branch=master)](https://travis-ci.com/biologger/speciesprimer)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/biologger/speciesprimer)](https://img.shields.io/docker/cloud/build/biologger/speciesprimer)
[![codecov](https://codecov.io/gh/biologger/speciesprimer/branch/master/graph/badge.svg)](https://codecov.io/gh/biologger/speciesprimer)
[![CodeFactor](https://www.codefactor.io/repository/github/biologger/speciesprimer/badge)](https://www.codefactor.io/repository/github/biologger/speciesprimer)


New in SpeciesPrimer v2.1
* Configfile option for pipeline setup (v2.1.1)
* Custom Blast DB support
* Email option for command line
* Increased speed
* Species synonyms are added to exceptions
* Bugfixes and KeyboardInterrupt rollback
* Simpler directory structure


## Contents
* [Hardware recommendations](https://github.com/biologger/speciesprimer/tree/master#hardware-recommendations)
* [quick start (Ubuntu 16.04)](https://github.com/biologger/speciesprimer/tree/master#quick-start-ubuntu-1604)
* [Introduction](https://github.com/biologger/speciesprimer/tree/master#introduction)
	* [Pipeline workflow and tools](https://github.com/biologger/speciesprimer/tree/master#pipeline-workflow-and-tools)
	* [Command line options](https://github.com/biologger/speciesprimer/tree/master#run-settings)

## Docs
* [Pipeline and Docker tutorial](https://github.com/biologger/speciesprimer/tree/master/docs/tutorial.md)
* [Advanced command line usage](https://github.com/biologger/speciesprimer/blob/master/docs/cmdlineonly.md)
* [Pipeline setup](https://github.com/biologger/speciesprimer/tree/master/docs/pipelinesetup.md)
* [Primerdesign](https://github.com/biologger/speciesprimer/tree/master/docs/primerdesign.md)
* [Troubleshooting](https://github.com/biologger/speciesprimer/tree/master/docs/troubleshooting.md)
* [Custom BLAST DB tutorial](https://github.com/biologger/speciesprimer/tree/master/docs/customdbtutorial.md)
* [More troubleshooting (Docker)](https://github.com/biologger/speciesprimer/tree/master/docs/dockertroubleshooting.md)
* [Docker and proxy settings](https://github.com/biologger/speciesprimer/tree/master/docs/dockerproxy.md)

## Minimum system requirements

* Quad core processor
* 16 GB RAM
* SSD / fast hard disk (recommended)
* 60 GB free space for nt database
* 4.5 GB for the docker image
* 5 - 20 GB for each analysis

## quick start (Ubuntu 16.04)

* [Download](https://www.docker.com/get-docker) and [install](https://docs.docker.com/install/) docker

		$ sudo docker pull biologger/speciesprimer
		$ mkdir $HOME/primerdesign
		$ mkdir $HOME/blastdb
		$ sudo docker run \
		-v $HOME/blastdb:/blastdb \
		-v $HOME/primerdesign:/primerdesign \
		-p 5000:5000 -p 9001:9001 \
		--name speciesprimer biologger/speciesprimer


* Open the address [http://localhost:5000](http://localhost:5000) or [http://127.0.0.1:5000](http://127.0.0.1:5000) in your favorite webbrowser
* Enter your E-mail address (required for the biopython NCBI Entrez module)
* Download the nt_v5 BLAST DB (>60 GB) or the ref_prok_rep_genomes DB (~6.5 GB). [BLAST DB](http://localhost:5000/blastdb)
* Customize the species list and other parameters if required. [SpeciesPrimer settings](http://localhost:5000/pipelineconfig)
* Navigate to Primer design and start primer design for new targets. [Primer design](http://localhost:5000/primerdesign)
* If you want to use the ref_prok_rep_genomes DB provide the path (/blastdb/ref_prok_rep_genomes) in the customdb settings field
* The results can be found in the Summary directory e.g. /primerdesign/Summary (container) or $HOME/primerdesign/Summary (host)

#### Use the pipeline with the command line

* After the docker run command open a new terminal

		# open an interactive terminal in the docker container
		$ sudo docker exec -it speciesprimer bash

* Download the nt BLAST DB (>60 GB):

		$ getblastdb.py -dbpath /blastdb --delete

* or download the ref_prok_rep_genomes DB (~6.5 GB):

		$ getblastdb.py -db ref_prok_rep_genomes -dbpath /blastdb --delete

* or alternatively

		$ cd /blastdb

		$ update_blastdb.pl --passive --decompress --blastdb_version 5 nt_v5
		# or
		$ update_blastdb.pl --passive --decompress ref_prok_rep_genomes

		$ cd /primerdesign

* Customize the species list and other parameters if required (see docs/pipelinesetup.md for more info):

		$ nano /pipeline/dictionaries/species_list.txt
		$ nano /pipeline/dictionaries/p3parameters
		$ nano /pipeline/dictionaries/no_blast.gi

* Start primer design

		$ speciesprimer.py

* Starting the script will start an assistant for the configuration of a new run

For more information and advanced settings see [Advanced command line usage](https://github.com/biologger/speciesprimer/blob/master/docs/cmdlineonly.md)

#### If you want to use the ref_prok_rep_genomes DB use the customdb option with the path

		/blastdb/ref_prok_rep_genomes

# Introduction
The SpeciesPrimer pipeline is intended to help researchers finding specific primer pairs for the detection and quantification of bacterial species in complex ecosystems. The pipeline uses genome assemblies of the target species to identify core genes (genes which are present in all assemblies) and checks the specificity for the target species using BLAST. Primer design is performed by primer3, followed by a stringent primer quality control. To make the evaluation of primer specificity faster and simpler, not all sequences of all bacterial species in the BLAST database are considered, the user has to provide a list of organisms which are expected to be present in the investigated ecosystem and should not be detected by the primer pair. The output of the pipeline is a comma separated file with possible primer pairs for the target species, which can be further tested and evaluated by the user.


## Pipeline workflow and tools

|Pipeline workflow|Tools|Reference|
|--|--|--|
|<tr><td colspan=3>__Input genome assemblies__</td></tr>|
|-  download|NCBI Entrez (Biopython)|[Cock et al. 2009](https://doi.org/10.1093/bioinformatics/btp163); [Sayers 2009](https://www.ncbi.nlm.nih.gov/books/NBK25499/)|
|- annotation|Prokka|[Seemann 2014](https://doi.org/10.1093/bioinformatics/btu153)|
|- quality control|BLAST+|[Altschul et al. 1990](https://doi.org/10.1016/s0022-2836%2805%2980360-2)|
|<tr><td colspan=3>__Core gene sequences__</td></tr>|
|- identification|Roary|[Page et al. 2015](https://doi.org/10.1093/bioinformatics/btv421)|
|- phylogeny|FastTree 2|[Price et al. 2010](https://doi.org/10.1371/journal.pone.0009490)|
|- selection of conserved sequences|Prank, consambig (EMBOSS),GNU parallel|[Löytynoja 2014](https://doi.org/10.1007/978-1-62703-646-7_10); [Rice et al. 2000](https://doi.org/10.1016/S0168-9525%2800%2902024-2); [Tange 2011](https://www.usenix.org/publications/login/february-2011-volume-36-number-1/gnu-parallel-command-line-power-tool)|
|- evaluation of specificity|BLAST+|[Altschul et al. 1990](https://doi.org/10.1016/s0022-2836%2805%2980360-2)|
||||
|<tr><td colspan=3>__Primer__</td></tr>|
|- design|Primer3|[Untergasser et al. 2012](https://doi.org/10.1093/nar/gks596)|
|- quality control|BLAST+, Mfold, MFEPrimer 2.0, MPprimer|[Altschul et al. 1990](https://doi.org/10.1016/s0022-2836%2805%2980360-2); [Zuker et al. 1999](https://doi.org/10.1007/978-94-011-4485-8_2); [Qu et al. 2012](https://doi.org/10.1093/nar/gks552); [Shen et al. 2010](https://doi.org/10.1186/1471-2105-11-143)|

The DBGenerator.py script from [Microbial Genomics Lab at CBIB](https://github.com/microgenomics/tutorials) and SQlite3 was used in an earlier version to create an SQL database from the Roary output.

Python modules and software used for the GUI:

[flask](https://github.com/pallets/flask)

[flask-wtf](https://github.com/lepture/flask-wtf)

[gunicorn](https://gunicorn.org/)

[MyDaemon](https://github.com/emrahcom/MyDaemon)


## Run settings
|Section|Command line option [Input]|Description|Default|
|--|--|--|--|
|General| target [str]|Name of the target species|None (required)|
|	|exception [str]|Name of a non-target bacterial species for which primer binding is tolerated|None|
|	|path [str]|Absolute path of the working directory|Current working directory|
|	|offline|Work offline with local genome assemblies|False|
|	|skip\_download|Skips download of genome assemblies from NCBI RefSeq FTP server|False|
|	|assemblylevel [all, complete, chromosome, scaffold, contig]| Only genome assemblies with the selected assembly status will be downloaded from the NCBI RefSeq FTP server|['all']|
|	|customdb [str]|Use the NCBI ref_prok_rep_genomes database or any other BLAST DB|None|
|	|blastseqs [100, 500, 1000, 2000, 5000]|Set the number of sequences per BLAST search. Decreasing the number of sequences requires less memory|1000|
|	|blastdbv5 |Uses the nt_v5 database or a v5 custom DB and limits all BLAST searches to taxid:2 (bacteria). May increase speed.|False|
|	|email [str]|Provide your email in the command line to access NCBI. No input required during the run.|None|
|	|intermediate|Select this option to keep intermediate files.|False|
|	|nolist|Do not use the (non-target) species list, only sequences without Blast hits are selected for primer design. May be used with a custom Blast DB|False|
|	|configfile [str]|Path to configuration file (json) to use custom species_list.txt, p3parameters, genus_abbrev.csv and no_blast.gi files|None|
|Quality control|qc\_gene  [rRNA, recA, dnaK, pheS, tuf]|Selection of housekeeping genes for BLAST search to determine the species of input genome assemblies|['rRNA']
|	 |ignore\_qc|Keep genome assemblies, which fail to meet the criteria of the quality control step|False|
|Pan-genome analysis|skip_tree|Skips core gene alignment (Roary) and core gene phylogeny (FastTree)|False|
|Primer design|minsize [int] | Minimal accepted amplicon size of PCR primer pairs|70|
|	|maxsize [int] | Maximal accepted amplicon size of PCR primer pairs|200|
|Primer quality control|mfold [float] | Set the deltaG threshold (max. deltaG) for the secondary structures at 60 °C in the PCR product, calculated by Mfold|-3.5|
|	|mpprimer [float] |Set the deltaG threshold (max. deltaG)  for the primer-primer 3’-end binding, calculated by MPprimer|-3.0|
|	|mfethreshold [int] | Threshold for MFEprimer primer pair coverage (PPC) score. Higher values: select for better coverage for target and lower coverage for for non-target sequences  (recommended range 80 - 100).|90|
