
# SpeciesPrimer

## Contents
* [Hardware recommendations](https://github.com/biologger/speciesprimer/tree/GUI#hardware-recommendations)
* [quick start (Ubuntu 16.04)](https://github.com/biologger/speciesprimer/tree/GUI#quick-start-ubuntu-1604)
* [Introduction](https://github.com/biologger/speciesprimer/tree/GUI#introduction)
	* [Pipeline workflow and tools](https://github.com/biologger/speciesprimer/tree/GUI#pipeline-workflow-and-tools)
	* [Command line options](https://github.com/biologger/speciesprimer/tree/GUI#run-settings)
* [Tutorial Docker Setup](https://github.com/biologger/speciesprimer/tree/GUI#Tutorial-ubuntu-1604)

## Docs
* [Pipeline setup](https://github.com/biologger/speciesprimer/tree/GUI/docs/pipelinesetup.md)
* [Primerdesign](https://github.com/biologger/speciesprimer/tree/GUI/docs/primerdesign.md)
* [Troubleshooting](https://github.com/biologger/speciesprimer/tree/GUI/docs/troubleshooting.md)
* [More troubleshooting (Docker)](https://github.com/biologger/speciesprimer/tree/GUI/docs/dockertroubleshooting.md)
* [Docker and proxy settings](https://github.com/biologger/speciesprimer/tree/GUI/docs/dockerproxy.md)

## Hardware recommendations

* Quad core processor
* 16 GB RAM
* SSD / fast hard disk (recommended)
* 60 GB free space for nt database
* 4.5 GB for the docker image
* 5 - 20 GB for each analysis

## quick start (Ubuntu 16.04)

* [Download](https://www.docker.com/get-docker) and install docker
		
		$ sudo docker pull biologger/speciesprimergui
		$ mkdir $HOME/primerdesign
		$ mkdir $HOME/blastdb
		$ sudo docker run \
		-v $HOME/blastdb:/home/blastdb \
		-v $HOME/primerdesign:/home/primerdesign \
		-p 5000:5000 -p 9001:9001 \
		--name speciesprimergui biologger/speciesprimergui


* Open the address [http://localhost:5000] or [http://127.0.0.1:5000] in your favorite webbrowser
* Enter your E-mail address (required for the biopython NCBI Entrez module)
* Navigate to SpeciesPrimer settings [http://localhost:5000/pipelineconfig]
* Download the nt BLAST DB
* Customize the species list and other parameters if required
* Navigate to Primer design [http://localhost:5000/primerdesign] and start primer design for new targets


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


## Run settings
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
|	|blastdbv5 | Uses the nt_v5 database and limits all BLAST searches to taxid:2 (bacteria). Increases speed.|False|
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
  
			$ sudo docker pull biologger/speciesprimergui

2. Now you have the image, you can display the image with
	* __HOST:__
  
			$ sudo docker images

3. If there is more than one image from the repository __biologger/speciesprimergui__, you can remove the image with the <none\> Tag
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
Create the container instance using the host directories as volumes for the docker container. In the container these directories are then located in /home/blastdb and /home/primerdesign. The name of the container can be changed (--name).
The -p option defines the ports which are open for the container so you can access the container app http://127.0.0.1:{hostport1/2} / http://localhost:{hostport1/2}}.
On the left side the host port is given and on the right side the container port. The container port is fixed and cannot be changed, if the host port is already used another port can be selected.
The link on the page where you can control the runs is however fixed to port 9001, but you can open the log file stream by opening http://localhost:{hostport2} in your browser.

#### docker run

1. __HOST:__

		$ sudo docker run \
		-v path_to_host_blastdb_dir:/home/blastdb \
		-v path_to_host_primerdesign_dir:/home/primerdesign \
		-p {hostport1}:5000 -p {hostport2}:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimergui

__Example:__

* __HOST:__

		$ sudo docker run \
		-v /home/biologger/blastdb:/home/blastdb \
		-v /home/biologger/primerdesign:/home/primerdesign \
		-p 5000:5000 -p 9001:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimergui

In the terminal you see that the server in the container was started.
Afterwards you can open the address [http://localhost:5000] or what port you have choosen for {hostport1} in your webbrowser.

#### docker stop		
You can shutdown the container by opening a terminal and the command 

        $ sudo docker stop {containername/id}

__Example:__

* __HOST:__

        $ sudo docker stop speciesprimer_pipeline

#### docker start        
The next time you do not have to repeat the docker run command (this would create a new container, without your modified settings)
Instead you simply start the container with the command 

        $ sudo docker start {containername/id}

__Example:__

* __HOST:__

        $ sudo docker start speciesprimer_pipeline
        
Afterwards you can open the address [http://localhost:5000] or what port you have choosen for {hostport1} in your webbrowser.

#### docker attach
If you want to see the status of the webserver in the container in your terminal (like after the docker run command)

        $ sudo docker attach {containername/id}

__Example:__

* __HOST:__

        $ sudo docker attach speciesprimer_pipeline
        
#### docker exec
If you w to access the container with the terminal you can use (the -it option is for the interactive terminal)

        $ sudo docker exec -it {containername/id} bash
        
__Example:__

* __HOST:__

        $ sudo docker exec -it speciesprimer_pipeline bash

#### Leave the container terminal
You can leave the docker container by typing exit        

__Example:__

* __CONTAINER:__

        $ exit
        

### Test the container

If not already started 

$ sudo docker start {containername/id} 

$ sudo docker exec -it {containername/id} bash

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
 		
 

Now you can try again to create a container with the sudo docker run command

* __HOST:__

		$ sudo docker run \
		-v /home/{linux_username}/blastdb:/home/blastdb \
		-v /home/{linux_username}/primerdesign:/home/primerdesign \
		-p {hostport1}:5000 -p {hostport2}:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimer

__Example:__

* __HOST:__

		$ sudo docker run \
		-v /home/biologger/blastdb:/home/blastdb \
		-v /home/biologger/primerdesign:/home/primerdesign \
		-p 5000:5000 -p 9001:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimer__
