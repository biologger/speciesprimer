# Advanced command line usage

There are some cases a user cannot or does not want to run the GUI. 
The command line version of SpeciesPrimer is always available by accessing 
the container via a interactive terminal. 
However, it is not necessary to run a local server at the same time and to have an extra 
terminal with the server status after the docker run command.
Therefore there is a way to run a container without the local server and the GUI.

* This docker run command will override the start of the GUI server (boot.sh)

		$ sudo docker run -v $HOME/blastdb:/blastdb \
		-v $HOME/primerdesign:/primerdesign \
		--name cmdline -it biologger/speciesprimer /bin/bash

* This docker run command will search for primers (using the default species list and primer settings) and removes the container afterwards (--rm).

		$ sudo docker run --rm -v $HOME/blastdb:/blastdb -v $HOME/primerdesign:/primerdesign \
		biologger/speciesprimer speciesprimer.py --target Lactobacillus_curvatus \
		--email your@email.com --assemblylevel complete

* Including the docker -it (interactive terminal) command will give you the opportunity to see the progress and potential error messages

		$ sudo docker run --rm -it -v $HOME/blastdb:/blastdb -v $HOME/primerdesign:/primerdesign \
		biologger/speciesprimer speciesprimer.py --target Lactobacillus_curvatus \
		--email your@email.com --assemblylevel complete

* Instead of using the --rm option (removing the container after the run) you can also just exit and stop a container and reuse it later. In this case giving the container a name (e.g. --name cmdline) is advantageous.

		$ exit

		# or

		Ctrl + D

	__HOST:__
 
		$ sudo docker stop cmdline

* And start it again later for a new run

	__HOST:__

		$ sudo docker start cmdline

		$ sudo docker exec -it cmdline /bin/bash

		# or

		$ sudo docker exec -it cmdline speciesprimer.py --target Lactobacillus_curvatus \
		--email your@email.com --assemblylevel complete


* Advanced settings

There is the option to provide a configfile in json format speciesprimer.py --configfile {path_to_file} e.g. to provide a custom species list without first accessing the container.

The filepath given by the user has to be absolute and located in one of the two volumes attached to the Docker container. 
As for example, if the /primerdesign directory (-v {HOSTPATH}:/primerdesign) was attached to $HOME/primerdesign (-v $HOME/primerdesign:/primerdesign) and the configuration file is in $HOME/primerdesign/customsettings/ the path for the config file (myconfigfile.json) is /primerdesign/customsettings/myconfigfile.json.

		$ sudo docker run --rm -it -v $HOME/blastdb:/blastdb -v $HOME/primerdesign:/primerdesign \
		biologger/speciesprimer speciesprimer.py --target Lactobacillus_curvatus \
		--email your@email.com --assemblylevel complete --configfile /primerdesign/customsettings/myconfigfile.json

The configuration file allows to provide a custom species list, Primer3 settings, genus abbreviations file and a list of excluded GIs. See [Pipline Setup](https://github.com/biologger/speciesprimer/blob/master/docs/pipelinesetup.md) for more information. Further, there is the option to provide a custom certificate for working with a TLS termination proxy.

| Setting | Key | supported file extension | default filename |
| Species list | species\_list | .txt | species_list.txt |
| Primer3 settings | p3settings |  | p3parameters |
| Genus abbreviations | genus\_abbrev | .csv | genus_abbrev.csv |
| Excluded GIs | excludedgis | .gi | no_blast.gi |
| TLS proxy certificate | certificate | .crt |  |

The simplest way to provide the settings is to save the custom files in a directory attached to the container (e.g. $HOME/primerdesign/customsettings)
examples for the file structure and format can be found in the [dictionaries directory](https://github.com/biologger/speciesprimer/tree/master/pipeline/dictionaries/default).

Examples for configuration files in json format (one line only and .json file extension)

* Configuration file for changing the species list {key: pathtofile}

		{"species_list": "/primerdesign/customsettings/myspecieslist.txt"}

* Configuration file for changing two or more files {key: pathtofile, key: pathtofile, ...}

		{"species_list": "/primerdesign/customsettings/myspecieslist.txt", "genus_abbrev": "/primerdesign/customsettings/mygenusabbrev.csv"}

* Configuration file to provide a TLS termination proxy certificate

		{"certificate": "/primerdesign/customsettings/my_TLS_Proxy_Cert.crt"}

* Alternatively provide the settings as list (works with species list and Excluded GIs)

		{"species_list": ["Acidipropionibacterium thoenii", "Streptococcus thermophilus", "Lactobacillus curvatus"], "excludedgis": ['1231231','1231232','1231233']}

* or as list of lists for Genus abbreviations

		{"genus_abbrev": [["Lactococcus", "Lc"], ["Propionibacterium", "Pb"], ["Enterococcus", "Ec"]]}

To reset to default settings just create a new container (docker run ...) or get the files from GitHub [dictionaries directory](https://github.com/biologger/speciesprimer/tree/master/pipeline/dictionaries/default) and include them in the configuration file.


