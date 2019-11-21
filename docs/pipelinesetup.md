## Pipeline setup

### Download the pre-formatted nucleotide collection (nt) BLAST database from NCBI

If you do not yet have a local copy of the nucleotide collection (nt) BLAST database:

In your internetbrowser [http://localhost:5000] Navigate to the BLAST DB link.

To download the database with the update_blastdb perl script just click on the Start download button below. On the next page Click Start Blast DB download. This will download the pre-formatted nt BLAST database from NCBI to your blastdb directory (specified with the docker run command). You can follow the progress by clicking stream log in new tab.

#### Troubleshooting:
The script changes to the "/blastdb" directory in the container and runs the command below

	$ update_blastdb.pl --passive --decompress --blastdb_version 5 nt_v5
	
To get possible hints or error messages you can run this command from a docker terminal.

__Example__

* __HOST__
	
		$ sudo docker exec -it {containername} bash

* __Container__
	
		$ cd /blastdb
		$ update_blastdb.pl --passive --decompress --blastdb_version 5 nt_v5

As an alternative the python script __getblastdb.py__ can be used.

The database has currently a size of about 65 GB, the archives are also quite large, therefore select the delete option if harddisk space is limited on your computer. If selected the script will automatically delete the archive files after extraction.

In your internetbrowser [http://localhost:5000] Navigate to the BLAST DB link.

Decide if you want to delete the archives after extraction and click on the second Start download button, 
the download starts after you clicked on "Start BLAST DB download" button on the next page. You can follow the progress by clicking stream log in new tab.

The click on the button is equivalent to the following command in the container.

__Example:__

* __Container__

		$ getblastdb.py -dbpath /blastdb --delete	

--------------------------------------------------
### Create a species_list.txt file
The species list is used to evaluate the specificity of the target sequences and to exclude unspecific primer pairs. 
It should contain the typical bacterial species you expect to be in the ecosystem you want to investigate with your PCR assay. (Including your target species of course).

In the browser [http://localhost:5000] Navigate to Pipeline settings.

Create a species_list.txt file or download the default species_list.txt file and adapt it according your needs.

Then upload it with the corresponding form.

If something was going wrong you can always reset the list to the default values with the Reset to default button.
Your adapted file in the container will be overwritten.

__If you create a new container with the docker run command the default settings are applied.
Keep a copy of your adapted files on your host computer in case you want to use it in a different container later on.__

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

------------------------

### Optional:

### Update the taxid list for the BLAST search with the nt BLAST database Version 5
The taxid list contains all NCBI taxonomic ids for bacteria (NCBI:txid2). From time to time NCBI may add new taxids if new species are added.

As an alternative if the __Entrez EDirect Utility__ and __the BLAST+ command line application__ are installed on the host computer, the __get_species_taxids.sh__ script from the BLAST+ program can be used to create this list.

__Example__

		$ get_species_taxids.sh -t 2 > /home/pipeline/dictionaries/2.txids

------------------------------------------------------
### Abbreviations of genus names
If you wish to use abbreviations of genus names for the pipeline output files,
add the genus name and the desired abbreviation to the "__genus_abbrev.csv__" file.

You can download, change and then upload it on the Pipeline settings page.

__If you create a new container with the docker run command the default settings are applied.
Keep a copy of your adapted files on your host computer in case you want to use it in a different container later on.__

-----------------------------------
### Primer3 settings
Parameters in this file can be changed with the same method as for the genus_abbrev.csv file. 
All settings except the amplicon minimal and maximal size are saved in the container and affect the primer design for all targets.

* Primer3 settings are changed in the "p3parameters" file
* see primer3 documentation for an explanation of the different parameters
<http://primer3.sourceforge.net/primer3_manual.htm>
* minimal and maximal amplicon size can be changed in the command line using the
--minsize [int] --maxsize [int] commands or on the settings page.

You can download, change and then upload it on the Pipeline settings page.

__If you create a new container with the docker run command the default settings are applied.
Keep a copy of your adapted files on your host computer in case you want to use it in a different container later on.__

-----------------------------------------------------------
### Exclude geneidentifier from BLAST

If certain gene identifier for sequences in the BLAST nt database should be ignored, because they have a wrong taxonomic classification for example. The gene identifier (GI) can be added to this list, one GI per line. For more information see the troubleshooting help.

You can download, change and then upload it on the Pipeline settings page.

__If you create a new container with the docker run command the default settings are applied.
Keep a copy of your adapted files on your host computer in case you want to use it in a different container later on.__
