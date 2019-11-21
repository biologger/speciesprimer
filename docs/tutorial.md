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

* If the pre-formatted nucleotide (nt) database from NCBI is already downloaded and unpacked on your computer, just add the path to the directory in the docker run command (-v path\_to\_host\_blastdb_dir:/blastdb) 

* Create a directory for primerdesign and one for the BLAST database

* __Example:__

	* Create two new directories in the home directory
 
 	* __HOST:__
	
 			# one for the primer design files
			$ mkdir $HOME/primerdesign
  			# one for the nucleotide blast database
 			$ mkdir $HOME/blastdb

### Run a container instance
Create the container instance using the host directories as volumes for the docker container. In the container these directories are then located in /blastdb and /primerdesign. The name of the container can be changed (--name).
The -p option defines the ports which are open for the container so you can access the container app http://127.0.0.1:{hostport1/2} / http://localhost:{hostport1/2}}.
On the left side the host port is given and on the right side the container port. The container port is fixed and cannot be changed, if the host port is already used another port can be selected.
The link on the page where you can control the runs is however fixed to port 9001, but you can open the log file stream by opening http://localhost:{hostport2} in your browser.

#### docker run

1. __HOST:__

		$ sudo docker run \
		-v path_to_host_blastdb_dir:/blastdb \
		-v path_to_host_primerdesign_dir:/primerdesign \
		-p {hostport1}:5000 -p {hostport2}:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimer

__Example:__

* __HOST:__

		$ sudo docker run \
		-v /home/biologger/blastdb:/blastdb \
		-v /home/biologger/primerdesign:/primerdesign \
		-p 5000:5000 -p 9001:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimer

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

If you see __root@{containerID}:/primerdesign#__ in the terminal, you have now access to the terminal of the container.

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
		-v /home/{linux_username}/blastdb:/blastdb \
		-v /home/{linux_username}/primerdesign:/primerdesign \
		-p {hostport1}:5000 -p {hostport2}:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimer

__Example:__

* __HOST:__

		$ sudo docker run \
		-v /home/biologger/blastdb:/blastdb \
		-v /home/biologger/primerdesign:/primerdesign \
		-p 5000:5000 -p 9001:9001 \
		--name speciesprimer_pipeline -it biologger/speciesprimer__
