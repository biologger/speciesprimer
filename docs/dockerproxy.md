### Configuration of Docker with DNS and proxy server

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

2. Open or create the __http-proxy.conf__ file to add the http proxy configuration
	* __Example:__ 
	
			$ sudo gedit /etc/systemd/system/docker.service.d/http-proxy.conf

3. Copy the text below with your {proxy_address:port} for the http proxy in the __http-proxy.conf__ file

		[Service] 
		Environment="HTTP_PROXY=http://{proxyadress:port}/" 
	
	* __Example:__

			[Service] 
			Environment="HTTP_PROXY=http://myproxy.proxy.com:8080/" 

4. Open or create the __https-proxy.conf__ file to add the https proxy configuration
	* __Example:__

			$ sudo gedit /etc/systemd/system/docker.service.d/http-proxy.conf/https-proxy.conf

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

#### TLS termination (Docker container)
1. Copy the certificate file (cer / crt) to a directory which is mounted in the container 
2. Copy the certificate file to /usr/local/share/ca-certificates/
3. Run update-ca-certificates
4. Export requests variable (add it to .bashrc to make it permanent)

    *__Example:__
    
        * __HOST:__
        
            $ mkdir $HOME/primerdesign/certificates
            $ cp $HOME/Proxy_CA.cer $HOME/primerdesign/certificates
         
         *__CONTAINER:__
         
            $ cp /home/primerdesign/certificates/Proxy_CA.cer /usr/local/share/ca-certificates/Proxy_CA.crt
            $ update-ca-certificates 
            $ export REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
            $ echo export REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt >> /root/.bashrc
