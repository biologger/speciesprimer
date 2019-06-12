
# SpeciesPrimer

## Contents
* [Hardware recommendations](https://github.com/biologger/speciesprimer/blob/master/README.md#hardware-recommendations)
* [quick start (Ubuntu 16.04)](https://github.com/biologger/speciesprimer#quick-start-ubuntu-1604)
* [Introduction](https://github.com/biologger/speciesprimer#introduction)
	* [Pipeline workflow and tools](https://github.com/biologger/speciesprimer#pipeline-workflow-and-tools)
	* [Command line options](https://github.com/biologger/speciesprimer#command-line-options)
* [Tutorial (Ubuntu 16.04)](https://github.com/biologger/speciesprimer#tutorial-ubuntu-1604)
	* [Docker Setup](https://github.com/biologger/speciesprimer#docker-setup)
	* [BLAST DB setup](https://github.com/biologger/speciesprimer#blast-db-setup)
	* [Pipeline setup](https://github.com/biologger/speciesprimer#pipeline-setup)
	* [Primerdesign](https://github.com/biologger/speciesprimer#primerdesign)
	* [Troubleshooting](https://github.com/biologger/speciesprimer#troubleshooting)
	* [More troubleshooting (Docker)](https://github.com/biologger/speciesprimer#More-troubleshooting-docker)

## Hardware recommendations

* Quad core processor
* 16 GB RAM
* SSD / fast hard disk (recommended)
* 60 GB free space for nt database
* 4.5 GB for the docker image
* 5 - 20 GB for each analysis

## quick start (Ubuntu 16.04)


		$ sudo docker pull biologger/speciesprimer
		$ mkdir $HOME/primerdesign
		$ mkdir $HOME/blastdb
		$ sudo docker run \
		-v $HOME/blastdb:/home/blastdb \
		-v $HOME/primerdesign:/home/primerdesign \
		--name speciesprimer -it biologger/speciesprimer
		$ cd $HOME/blastdb
		$ update_blastdb.pl --passive --decompress nt

			
*  create the species_list.txt file: a list of bacterial species abundant in the ecosystem you want to investigate and copy the list into the container or use the one provided (dairy bacteria).

		$ cp path_to_host_dir/species_list.txt \
		/home/primerpipeline/dictionaries/species_list.txt

*  Start the pipeline

		$ speciesprimer.py

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
|- selection of conserved sequences|SQlite3, Prank, consambig (EMBOSS),GNU parallel|[Löytynoja 2014](https://doi.org/10.1007/978-1-62703-646-7_10); [Rice et al. 2000](https://doi.org/10.1016/S0168-9525%2800%2902024-2); [Tange 2011](https://www.usenix.org/publications/login/february-2011-volume-36-number-1/gnu-parallel-command-line-power-tool)|
|- evaluation of specificity|BLAST+|[Altschul et al. 1990](https://doi.org/10.1016/s0022-2836%2805%2980360-2)|
||||
|<tr><td colspan=3>__Primer__</td></tr>|
|- design|Primer3|[Untergasser et al. 2012](https://doi.org/10.1093/nar/gks596)|
|- quality control|BLAST+, Mfold, MFEPrimer 2.0, MPprimer|[Altschul et al. 1990](https://doi.org/10.1016/s0022-2836%2805%2980360-2); [Zuker et al. 1999](https://doi.org/10.1007/978-94-011-4485-8_2); [Qu et al. 2012](https://doi.org/10.1093/nar/gks552); [Shen et al. 2010](https://doi.org/10.1186/1471-2105-11-143)|


## Command line options
|Section|Command line option [Input]|Description|Default|
|--|--|--|--|
|General| target [str]|Name of the target species|None (required)|
|	|exception [str]|Name of a non-target bacterial species for which primer binding is tolerated|None|
|	|path [str]|Absolute path of the working directory|Current working directory|
|	|offline|Work offline with local genome assemblies|False|
|	|skip\_download | Skips download of genome assemblies from NCBI RefSeq FTP server|False|
|	|assemblylevel [all, complete, chromosome, scaffold, contig]| Only genome assemblies with the selected assembly status will be downloaded from the NCBI RefSeq FTP server| ['all']|
|	|remote|Use the BLAST+ remote flag for BLAST searches|False|
|	|blastseqs [100, 500, 1000, 2000, 5000]|Set the number of sequences per BLAST search. Decreasing the number of sequences requires less memory|1000|
|Quality control|qc\_gene  [rRNA, recA, dnaK, pheS, tuf]|Selection of housekeeping genes for BLAST search to determine the species of input genome assemblies|['rRNA']
|	 |ignore\_qc|Keep genome assemblies, which fail to meet the criteria of the quality control step|False|
|Pan-genome analysis|skip_tree|Skips core gene alignment (Roary) and core gene phylogeny (FastTree)|False|
|Primer design|minsize [int] | Minimal accepted amplicon size of PCR primer pairs|70|
|	|maxsize [int] | Maximal accepted amplicon size of PCR primer pairs|200|
|Primer quality control|mfold [float] | Set the deltaG threshold (max. deltaG) for the secondary structures at 60 °C in the PCR product, calculated by Mfold|-3.5|
|	|mpprimer [float] |Set the deltaG threshold (max. deltaG)  for the primer-primer 3’-end binding, calculated by MPprimer|-3.0|
|	|mfethreshold [int] | Threshold for MFEprimer primer pair coverage (PPC) score. Higher values: selct for better coverage for target and lower coverage for for non-target sequences  (recommended range 80 - 100).|90|

# Tutorial (Ubuntu 16.04):

## Docker setup

### Download and install docker
Download from <https://www.docker.com/get-docker> and
see the docs for installation instructions <https://docs.docker.com/>

### Download images from docker hub
1. Open a terminal:

		
	* __HOST:__
  
			$ sudo docker pull biologger/speciesprimer

2. Now you have the image, you can display the image with
	* __HOST:__
  
			$ sudo docker images

3. If there is more than one image from the repository __biologger/speciesprimer__, you can remove the image with the <none\> Tag
 	* __HOST:__
 
			$ sudo docker rmi {image_id}
 		
### Choose directories

1. Decide which directories (on the host) should be used by the container

* If the pre-formatted nucleotide (nt) database from NCBI is already downloaded and unpacked on your computer, just add the path to the directory in the docker run command (-v path\_to\_host\_blastdb_dir:/home/blastdb) 

* Create a directory for primerdesign and one for the BLAST database

* __Example:__

	* Create two new directories in the home directory
 
 	* __HOST:__
	
 			# one for the primer design files
			$ mkdir /home/biologger/primerdesign
  			# one for the nucleotide blast database
 			$ mkdir /home/biologger/blastdb

### Run a container instance
Create the container instance using the host directories as volumes for the docker container. In the container these directories are then located in /home/blastdb and /home/primerdesign. The name of the container can be changed (--name), the -it command immediately starts an interactive terminal in the container.

1. __HOST:__

		$ sudo docker run \
		-v path_to_host_blastdb_dir:/home/blastdb \
		-v path_to_host_primerdesign_dir:/home/primerdesign \
		--name speciesprimer_pipeline -it biologger/speciesprimer

2. After a while you should now see something like this in the terminal
__root@{containerID}:/home/primerdesign#__

__Example:__

* __HOST:__

		$ sudo docker run \
		-v /home/biologger/blastdb:/home/blastdb \
		-v /home/biologger/primerdesign:/home/primerdesign \
		--name speciesprimer_pipeline -it biologger/speciesprimer

### Test the container
If you see __root@{containerID}:/home/primerdesign#__ in the terminal, you have now access to the terminal of the container.

Test if you have mounted the volumes correctly 

* __CONTAINER:__
  
		$ echo test > test.txt

	* Check if you find test.txt
  
			$ ls -l
		
	
  	 
* __HOST:__

	* Check if you find test.txt on the host
  
			$ ls -l /home/{linux_username}/primerdesign

If you want to delete this test.txt file there are two options

1. Do it in the container 
	* __CONTAINER:__
 
			$ rm test.txt

2. Do it on the host
 	* Change the owner of the files in the primerdesign directory on the host (recursively).

 	* __HOST:__
 
			$ sudo chown -R {linux_username} {path_to_primerdesign_dir}
 
	 * Now you can move and delete the files and directories.

	__Example__:

	* __HOST:__
 
			$ sudo chown -R biologger /home/biologger/primerdesign

### Stop the container

* Exit the container by typing 

	* __CONTAINER:__
  
			$ exit
  
* The container is still running and you could return to a container terminal with:
  
  	* __HOST:__

			$ sudo docker exec -it {ContainerID/ContainerName} bash

* Stop the container

	* __HOST:__
 
			$ sudo docker stop {ContainerID/ContainerName}

### Start the container

1. Start the container

	* __HOST:__
 
			$ sudo docker start {ContainerID/ContainerName}
		
2. Start a terminal in the container

	* __HOST:__
 
			$ sudo docker exec -it {ContainerID/ContainerName} bash

-----

## BLAST DB setup

### Download nt BLASTDB (~60 GB) (recommended)
1. In the container terminal change to the blastdb directory in the container 
 	* __CONTAINER:__
 
			$ cd /home/blastdb

2. Start the perl script from the BLAST+ program 
 	* __CONTAINER:__
 
			$ update_blastdb --passive --decompress nt
 
			# or

			$  /usr/bin/update_blastdb --passive --decompress nt

	* Grab a coffee... This may take a couple of hours.
	
3. Change back to the primerdesign directory after the script has finished
	
  	* __CONTAINER:__
  
			$ cd /home/primerdesign

### Download nt BLASTDB (~60 GB) (alternative)
Python script for the download of the database (Useful for Networks with a Proxyserver).

* Go to the [NCBI Server](https://ftp.ncbi.nlm.nih.gov/blast/db/) and check in how many parts the nt database is splitted. This is the {number\_of\_parts} you have to give in the command line in order to download all parts of the database
* __CONTAINER:__
 
		$ get_blastdb.py -dbpath /home/blastdb -delete -parts {number_of_parts}

* The download of the tar archives and the md5 files starts and the archives are unpacked and deleted if the -delete option is used. 

* Grab a coffee... This may take a couple of hours.

__Example__:

* Last part of the nt database:
  * nt.57.tar.gz	174 MB	Apr 29 15:32

* __CONTAINER:__

		$ get_blastdb.py -dbpath /home/blastdb -delete -parts 57

__Troubleshooting__

If a problem occurs, the connection is lost and restarting the script shows an error message like this:

		FileNotFoundError: [Errno 2] No such file or directory: '/home/blastdb/nt.00.tar.gz.md5'

* Change to the blastdb directory and delete the last unpacked tar archive and start the script again.

__Example__:

* __CONTAINER:__

		$ cd /home/blastdb
		$ rm nt.00.tar.gz
		$ cd /home/primerdesign

-----

## Pipeline setup

### Create a species_list.txt file
The species list is used to exclude unspecific primer pairs it should contain the typical bacterial species you would expect in the ecosystem you want to investigate with your PCR assay. 

 * Add and remove single species names

  * use nano text editor in the terminal
  * __CONTAINER:__
  
		$ nano /home/primerpipeline/dictionaries/species_list.txt

 * Create a species list with an text editor on the host computer and copy it to the container.
 __If you create a new container instance from the image this step has to be repeated!!!__

1. Change to the primerdesign directory on the host 

2. Create the file with a text editor

3. Copy the file to the container

__Example:__

* __HOST:__

		$ cd /home/biologger/primerdesign
		$ gedit species_list.txt
		
  * File example: (one species name per line)
   
	  Microbacterium oxydans 
	  Micrococcus luteus 
	  Moraxella osloensis 
	  Morganella morganii 
	  Morganella psychrotolerans 
	  Mycobacterium vaccae 
	  Ornithinibacillus bavariensis 
	  Paenibacillus lactis 
	  Paludibacter propionicigenes 
	  Pantoea agglomerans 
	  Pediococcus acidilactici 
	  ...

* __CONTAINER:__

		$ cp /home/primerdesign/species_list.txt \
		/home/primerpipeline/dictionaries/species_list.txt

### Optional parameters:

### Abbreviations of genus names
If you wish to use abbreviations of genus names for the pipeline output files,
add the genus name and the desired abbreviation to the "genus\_abbrev.csv" file
in the "/home/primerpipeline/dictionaries" directory.

1. Use the nano text editor in the container 
	* __CONTAINER:__

			$ nano /home/primerpipeline/dictionaries/genus_abbrev.csv

2. or copy the file on your host, 
	* __CONTAINER:__
  
			$ cp /home/primerpipeline/dictionaries/genus_abbrev.csv \
			/home/primerdesign/genus_abbrev.csv
 
 	* Change the file on the host 

 	* Copy to container
 	__The file will only be in this container instance, if another container is created from the image this step has to be repeated!!!__.

  	* __CONTAINER:__
  	
  			$ cp /home/primerpipeline/dictionaries/genus_abbrev.csv \
  			/home/primerdesign/genus_abbrev.csv

### Primer3 settings
Parameters in this file can be changed with the same methods as for the genus_abbrev.csv file.

* Primer3 settings are changed in the "/home/primerpipeline/p3parameters" file
* see primer3 documentation for an explanation of the different parameters
<http://primer3.sourceforge.net/primer3_manual.htm>
* minimal and maximal amplicon size can be changed in the command line using 
--minsize [int] --maxsize [int]  commands



## Primerdesign

### First run
A good starting point is to use speciesprimer.py
choose a target (e.g. __Lactobacillus curvatus__) and use the default values.  
Quality control gene(s) and assembly level have to be specified 
(for assembly level e.g. "__all__" and "__tuf__" for quality control genes).  
30 - 200 genome assemblies are good as a starting point.
Fewer or more genome assemblies have to be of good quality to give good results.

* Just start the pipeline script


__Example:__

* __CONTAINER:__

		$ speciesprimer.py
		$ n # new / (hit return)

Answer the questions, primer design starts immediately after the last question is answered.

__If the pipeline was stopped you can continue primer design like this__

* __CONTAINER:__

		$ speciesprimer.py
		$ s # start / (hit return)

The log file of the pipeline run is saved in the primerdesign directory.
After the run has finished you will find a Summary directory in the primerdesign 
directory with one directory for each species. The selected configuration for 
the pipeline run is saved in a json file in the config directory and is also 
copied to the summary directory after the run has finished.

* __Example:__
__/home/biologger/primerdesign/Summary/Lactobacillus_curvatus__
In there you will find:
	* a list with designed primers 	
		* __Lb_curva_primer.csv__
	
	* statistics and a summary of the initial quality control of the genome assemblies 
		
		* __Lb_curva_pipeline_stats_2018_04_29.txt__
		* __Lb_curva_qc_sequences.csv__
		
	* the phylogeny files (if skip_tree = False)
		* __core_gene_alignment.aln__
		* __Lb_curva_tree.newick__
		
	* the configuration file
	
		* __config_2018_04_29.json__
		
The newick tree can be opened for example with seaview (http://doua.prabi.fr/software/seaview).

The configuration selected for the pipeline run is also saved in the  
__/home/biologger/primerdesign/Lactobacillus_curvatus/config/config.json__ file.

If no primer pairs were found or problems occur during the run, check your input genomes (their quality is crucial) and the troubleshooting section below.

### Directory tree
![directory_tree](directory_tree.png  "Directory tree example")

In the container the directory "primerdesign" is the starting directory, on the host it is the directory you specified in the initial docker run command. In there are directories for each target species, a Summary directory and an excludedassemblies directory. In the summary directory the results of the pipeline runs are saved for each target species separately. In the excludedassemblies directory the prokka annotation files of assemblies which did not pass the quality control are saved for each species separately. You can delete the directories with the annotation files if you want to save hard disk space. The information is still saved in the text file in the target species directory.

### Using user provided genomes
* Names of input fasta files should not contain dots and not more than one underscore.
1. Start the batchassist.py script and step through the configuration. 
2. In the primerdesign directory on the host you now see directories with the species names containing a config and a genomic_fna directory.
3. Copy your genome fasta files into the genomic_fna directory

Using your own prokka annotated files is not recommended. Depending on the version of prokka some gene names could have changed, this can cause problems during the quality control and the pan-genome analysis. The name convention has to be exactly the same for all input genomes, errors there can lead to a variety of strange errors during the pipeline run.

### batchassist.py
The batchassist.py script creates configuration files for future primer search. 
It offers individual or global selection of the command line options for multiple target species.  

To start the pipeline with this configuration files, use speciesprimer.py
and choose "start" (s) and then (a)ll or (s)elect, to start primer search for all species or to select certain target species.

__Example:__ search primer pairs for all species (without a primer.csv file) with a configuration file 

* __CONTAINER:__

		$ speciesprimer.py
		$ s # (hit return)
		$ a # (hit return)

## Troubleshooting

#### Delete files created by the docker container from the host.

* __HOST:__

		$ sudo chown -R {username} {path_to_primerdesign_dir}

__Example:__

* __HOST:__

		$ sudo chown -R biologger /home/biologger/primerdesign

Now you are the owner and can move and delete the files and directories.

### Troubleshooting pipeline runs

#### "fatal error while working on {genus_species} check logfile "
The log file can be found in the primerdesign directory. It is named speciesprimer\_yyyy\_mm\_dd.log.
* __Example:__

		speciesprimer_2018_06_08.log

##### Naming of user provided genome assemblies
* Do not use dots in the name (use filename.fna instead of file.name.fna) 

##### Errors during download of genome assemblies

		<urlopen error ftp error: TimeotError(110, 'Connection timed out')> 

Sometimes it can happen (especially on Virtual Machines) that the internet connection is interrupted for some seconds. Just restart the pipeline and it should work again. The configuration selected for the pipeline run has been saved in the  __{genus_species}/config/config.json__ file. The pipeline will check the files already downloaded and will not try to download them again. 

__Example:__ Restart the pipeline with the same configuration

		$ speciesprimer.py
		# select start
		$ s (Enter/Return)
		# select specific
		$ s (Enter/Return)
		# select the working directory (primerdesign = current working directory)
		$ (Enter/Return)
		# Type the species name for which you would like to repeat the run
		$ Lactobacillus curvatus
		
##### Error during quality control of genome assemblies

		Error: No genomes survived QC
		
		# or
		
		Error: Less than two genomes survived QC

You can check the results of the quality control in the __{genus_species}\_qc_sequences.csv__ file in the __"primerdesign/Summary"__ directory	

__Example:__

* __primerdesign/Summary/Lactobacillus\_curvatus/Lb\_curva\_qc\_sequences.csv__

* There you will see why the quality control for the genomes failed.

|qc gene		|BLAST result                   |Explanation							|
|-----------------------|:------------------------------|:--------------------------------------------------------------|
|passed QC   		|Lactobacillus curvatus   	|correct species matched best in BLAST search			|
|failed QC      	|Lactobacillus casei		|wrong species matched best in BLAST search			|
|Duplicate     		|				|NCBI genomes only, a newer assembly was selected		|
|Max contigs		|				|Too many small contigs >500					|
|QC gene missing 	|				|The sequence of the QC gene was not identified in the assembly	|


It could be that the nt database contains a sequence with a wrong species assigned to it.

1. select another QC gene and try again

2. Troubleshoot:
	* check the “{genus_species}\_qc\_sequences\_details.csv” file
	* check the Hit DB_id (NCBI accession of the sequence of the BLAST hit)
	* Go to NCBI and check e.g. by BLAST if the sequence is assigned to the correct taxonomic species
	* Sequence can be “banned” from future quality control BLAST searches 

	__Example:__ Bann a sequence from future QC BLAST searches

	* Identify the Hit GI (Gene identifier) in the “{genus_species}\_qc\_sequences\_details.csv” file

	* __CONTAINER:__

			$ nano /home/primerpipeline/NO_Blast/NO_BLAST.gi

	* Copy the Hit GI in the NO_BLAST.gi file and save it

3. Restart the run
	* Go to the excludedassemblies directory and copy the removed assemblies back to the target species directory
	
	* Delete "excluded_list.txt"
	
	* __CONTAINER:__
			
			$ rm /home/primerdesign/excludedassemblies/Lactobacillus_curvatus/excluded_list.txt

	* Restart the pipeline run
	
			$ speciesprimer.py

## Troubleshooting Docker
#### Conflict with the container name 
	"docker: Error response from daemon: Conflict. The container name "/speciesprimer_pipeline" is already in use by container 	"e9d0de003ce8eff06b34f8f46e4934797052e16dcdbd7e60214d05ea3828a70", You have to remove (or rename) that container to be able to reuse that name"
1. Display your containers
	* __HOST:__
 
			$ sudo docker ps -a

2. Stop the container with the container ID or the container name
	* 
 			
			$ sudo docker stop {ContainerID/containername}

3. Delete the container
	* 
 			
			$ sudo docker rm {ContainerID/containername}
 		
 

Now you can try again to create a container with the sudo docker run  command

* __HOST:__

		$ sudo docker run \
		-v /home/{linux_username}/blastdb:/home/blastdb \
		-v /home/{linux_username}/primerdesign:/home/primerdesign \
		--name speciesprimer_pipeline -it biologger/speciesprimer

__Example:__

* __HOST:__

		$ sudo docker run \
		-v /home/biologger/blastdb:/home/blastdb \
		-v /home/biologger/primerdesign:/home/primerdesign \
		--name speciesprimer_pipeline -it biologger/speciesprimer__

### Configuration of Docker with a DNS server and/or proxy server

#### DNS server (Docker daemon)
1. Open /etc/docker/daemon.json using a text editor
	* __Example:__
		
		$ sudo gedit /etc/docker/daemon.json
		
2. Copy the text below with the brackets and your DNS server adress into the deamon.json file

		{
		  "dns": ["DNSadress1", "DNSadress2"]
		}
		
	* __Example:__ (Google Public DNS IP addresses)

			{
			    "dns": ["8.8.8.8, "8.8.4.4"]
			}

3. Restart docker
$ sudo service docker restart

#### Proxy server (Docker daemon)

1. Create the docker.service.d directory if it does not exist
	* __Example:__ 

			$ sudo mkdir /etc/systemd/system/docker.service.d

2. Open / create  /etc/systemd/system/docker.service.d/http-proxy.conf to add the http proxy configuration
	* __Example:__ 
	
			$ sudo gedit /etc/systemd/system/docker.service.d/http-proxy.conf

3. Copy the text below with your {proxy_address:port} for the http proxy in the http-proxy.conf file

		[Service] 
		Environment="HTTP_PROXY=http://{proxyadress:port}/" 
	
	* __Example:__

			[Service] 
			Environment="HTTP_PROXY=http://myproxy.proxy.com:8080/" 

4. Open / create  /etc/systemd/system/docker.service.d/https-proxy.conf to add the https proxy configuration
	* __Example:__

			$ sudo gedit /etc/systemd/docker.service.d/https-proxy.conf

5. Copy the text below with your {proxy_address:port} for the https proxy in the https-proxy.conf file

		[Service] 
		Environment="HTTPS_PROXY=http://{proxyadress:port}/" 

	* __Example:__

			[Service] 
			Environment="HTTPS_PROXY=http://myproxy.proxy.com:8080/" 
		
6. Reload the docker daemon

		$ sudo systemctl daemon-reload

7. Check if the environment variables are set

		$ sudo systemctl show --property Environment docker

8. Restart docker

		$ sudo systemctl restart docker

#### DNS server (Docker container)
1. Open /etc/default/docker
	* __Example:__

			$ sudo gedit /etc/default/docker

2. Remove the hash (#) from the line  with DOCKER_OPTS= and add the  DNS server adress

		DOCKER_OPTS="--dns DNSadress1  --dns DNSadress2"
		
	* __Example:__ (Google Public DNS IP addresses)

			DOCKER_OPTS="--dns 8.8.8.8  --dns 8.8.4.4"

#### Proxy server (Docker container)
1. Create a .docker directory in your $HOME directory

		mkdir .docker

2. Create a config.json file

		$ sudo gedit ~/.docker/config.json

3. Add your proxy configuration as in the example below
	* __Example:__

			{
				"proxies": {
					"default": {	
						 "httpProxy": "http://myproxy.proxy.com:8080/",
						 "httpsProxy": "http://myproxy.proxy.com:8080/",
						 "ftpProxy": "http://myproxy.proxy.com:8080/"
					}
				 }
			}


