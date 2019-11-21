# Tutorial for command line only

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
