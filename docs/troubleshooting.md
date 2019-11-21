## Troubleshooting

#### nt BLAST DB  / nt BLAST DB version 5

Because a new version of BLAST databases was released recently [(BLAST DB V5)](https://ftp.ncbi.nlm.nih.gov/blast/db/v5/blastdbv5.pdf), there are two versions of the nt database available. Therefore, if you have the BLAST DB nt version 5 then you should activate this option in the pipeline settings.

If the blastdbv5 option is selected, the pipeline uses the __nt_v5.nal__ (alias) file to identify the different parts of the database. Otherwise, it looks for the __nt.nal__ file. For the blastn commands this changes the database command as follows.

		$blastn -db nt_v5 -query ... 	# new database
		
		# or
		
		$blastn -db nt -query ... 		# old database

However, if you do want to have the "old" behaviour you can copy the __nt_v5.nal__ file in the blastdb directory and rename it to __nt.nal__. Then the BLAST searches will use the BLAST nt v5 without using the taxonomy awareness features.

-------------------------------------------------------------------------------------------------------------
#### GUI trouble

* If pages are not loaded and an Error message appears.

1. Check if the container is running

	__HOST:__

		$ sudo docker ps

Something like this should show up
CONTAINER ID	IMAGE			COMMAND		CREATED		STATUS		PORTS						NAMES
ae681870458a	biologger/speciesprimer   "/boot.sh"    2 weeks ago     Up 5 hours      0.0.0.0:5000->5000/tcp, 0.0.0.0:9001->9001/tcp  specprime2.1

2. Try to navigate to the login page and enter (again) your e-mail address.

3. If this does not help...There could be a problem with the code. Please write to biologger@protonmail.com or on GitHub https://github.com/biologger/speciesprimer

If you want to help...

* Open a new terminal and type:

	__HOST:__

		$ sudo docker attach {containerid/name}

	__example__:

		$ sudo docker attach specprime2.1

* Repeat the action you did before the error page was shown and provide the error message in the attached terminal.

	__example__:

		[2019-11-21 14:55:53,018] ERROR in app: Exception on /change_settings [POST]
		Traceback (most recent call last):
		  File "/usr/local/lib/python3.5/dist-packages/flask/app.py", line 2446, in wsgi_app
		    response = self.full_dispatch_request()
		  File "/usr/local/lib/python3.5/dist-packages/flask/app.py", line 1951, in full_dispatch_request
		    rv = self.handle_user_exception(e)
		  File "/usr/local/lib/python3.5/dist-packages/flask/app.py", line 1820, in handle_user_exception
		    reraise(exc_type, exc_value, tb)
		  File "/usr/local/lib/python3.5/dist-packages/flask/_compat.py", line 39, in reraise
		    raise value
		  File "/usr/local/lib/python3.5/dist-packages/flask/app.py", line 1949, in full_dispatch_request
		    rv = self.dispatch_request()
		  File "/usr/local/lib/python3.5/dist-packages/flask/app.py", line 1935, in dispatch_request
		    return self.view_functions[rule.endpoint](**req.view_args)
		  File "/pipeline/gui/app/routes.py", line 168, in change_settings
		    path = old_settings['path']
		UnboundLocalError: local variable 'old_settings' referenced before assignment


-------------------------------------------------------------------------------------------------------------
#### Delete files created by the docker container from the host.

* __HOST:__

		$ sudo chown -R {username} {path_to_primerdesign_dir}

__Example:__

* __HOST:__

		$ sudo chown -R biologger /home/biologger/primerdesign

Now you are the owner and can move and delete the files and directories.

-------------------------------------------------------------------------------------------------------------
#### Setting the timezone for the container

* __Host:__

		$ sudo docker exec -it {container name} bash

* __Container:__

		$ dpkg-reconfigure tzdata

### Troubleshooting pipeline runs
------------------------------------------
#### Primer quality control
If the primer quality control gives many index errors the creation of the databases did probablly not work. Then the files ending with sqlite3.db are all of the size 2.0 kb, which means empty. The MFEprimer IndexDb.sh script uses the python2 psutil module. Try to install the psutil module, delete the primer_QC directory and try again. 

---------------------------------------------------------------------------------------------------
#### "fatal error while working on {genus_species} check logfile "

This error is shown if the speciesprimer.py script was started via the terminal.

With the graphic user interface the logfile can be streamed during pipeline runs http://localhost:9001/ including errors and exceptions.

The log file can be found in the primerdesign directory. It is named speciesprimer\_yyyy\_mm\_dd.log.

* __Example:__

		speciesprimer_2018_06_08.log

------------------------------------------------------------------------------------------------
#### Naming of user provided genome assemblies
* Do not use dots in the name (use filename.fna instead of file.name.fna) 

----------------------------------------------------------------------------------
#### Errors during download of genome assemblies

		<urlopen error ftp error: TimeotError(110, 'Connection timed out')> 

Sometimes it can happen (especially on Virtual Machines) that the internet connection is interrupted for some seconds. Just restart the pipeline and it should work again. The configuration selected for the pipeline run has been saved in the  __{genus_species}/config/config.json__ file. The pipeline will check the files already downloaded and will not try to download them again. 

---------------------------------------------------------------------------
#### Restart / continue stopped pipeline runs

* With the graphic user interface navigate to the [primerdesign page](http://localhost:5000/primerdesign) and select the "Search" button. 

* Type the path to the primerdesign directory

* Type the species name

* Hit the "Search configuration files" button

* Hit the "Start Primerdesign" button to start

The pipeline searches in the config directory for config.json files and continues the run 

 __Example:__

 
* Path to search for configuration files: __/primerdesign__
* Search only configuration files for the selected species: __Lactobacillus curvatus__
* --> SpeciesPrimer will search for config files in: __/primerdesign/Lactobacillus_curvatus/config__


__Example:__ 

Restart the pipeline with the same configuration (commandline)

		$ speciesprimer.py
		# select start
		$ s (Enter)
		# select specific
		$ s (Enter)
		# select the working directory (primerdesign = current working directory)
		$ (Enter)
		# Type the species name for which you would like to repeat the run
		$ Lactobacillus curvatus
		
#### Error during quality control of genome assemblies

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
	* Go to NCBI and check e.g. by BLAST if the sequence is assigned to the correct species
    

3. Restart the run
	* Go to the excludedassemblies directory and copy the removed assemblies you want to include back to the target species directory	
	* Delete "excluded_list.txt" or remove the entries of the assemblies you want to keep.
	* Add the GI of sequences in the nt database with wrong taxonomic classification to the NO_BLAST.gi file (see below)
	* Alternatively, change the ignore_qc option to true (change settings (GUI) or in the config.json file)
	* Restart the pipeline run

	* __Example__

		* __CONTAINER:__

				# remove the entries of the assemblies you want to keep.
				$ nano /primerdesign/excludedassemblies/Lactobacillus_curvatus/excluded_list.txt
            	
		* __Host__:            

				$ gedit $HOME/primerdesign/excludedassemblies/Lactobacillus_curvatus/excluded_list.txt

-----------------------------------------------------------------
#### Only a few core genes are reported

In the summary directory and in the primerdesign/{species} directory, there is  a file with statistics of the pipeline run ("pipeline_stats_yyyy_mm_dd.txt"). If there are only few core genes reported you may have some problems with the input genome assemblies. Roary which is used to create the pan-genome works only with quite similar genomes (species level). Check the quality of your input genomes and that there are no contaminations and genomes from differnt species.

* Check the quality of the input genome assemblies
* Check the phylogeny tree made from the core_gene_alignment file using SeaView or another graphical tree viewer. If there are genome assemblies which build distinct clusters maybe try to repeat the analysis without them.

-----------------------------------------------------------------
#### Thousands of primer are reported and primer QC takes forever 
Probably you had only a few input genomes (2 are the minimum), however to identify relevant core genes it is better to use at least 10 input genomes. If the genome assemblies are very similar more conserved sequences are identified. 

* Try to increase the number of input genomes
* Select more stringent primer3 settings
* Check if you have a proper selection of non-target species in the species list (some species closely related to the target and many others). (If not, to many sequences are identified as species specific)

--------------------------------------------------------------------------------	
#### Exclude gene identifier (GI) from BLAST results

__Example__

We assume that we are searching for primer for Lactobacillus curvatus and we receive only few species specific conserved sequences. If the input assemblies are of good quality and many core genes are identified, we can continue to troubleshoot by looking at the __mostcommonhits.csv__ file in the "Summary" directory.
If we see there a GI and a species name different from Lactobacillus curvatus with hits for (almost) every core gene (query), we can try to check if the taxonomy of the genome sequence is correct and if not (if it is actually a Lactobacillus curvatus genome) we can exclude the corresponding GI by adding it to the list.

During parsing of BLAST results or the BLAST summary file the GIs of the nontarget hits are searched for the excluded GIs. If the sequence with the excluded GI is the only nontarget hit for the query the corresponing core gene is considered specific.

* navigate to SpeciesPrimer settings [http://localhost:5000/pipelineconfig](http://localhost:5000/pipelineconfig)
* Download the current NO_BLAST.gi file
* Copy the Hit GI in the NO_BLAST.gi file and upload it
	
#### Taxonomic vagueness
If many sequences in the nt database are historically or because of some other reasons wrongly assigned you can add the species name as an exception. The pipeline handles this species like the actual target species. 

