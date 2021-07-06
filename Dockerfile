
FROM biologger/speciesprimerdeps:v2.2

LABEL maintainer="biologger@protonmail.com"

# Copy the directory contents into the docker directory
COPY . /
COPY boot.sh /

# Set env variables and change mod
ENV FLASK_APP="/pipeline/gui/speciesprimergui.py"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV SSL_CERT_DIR="/etc/ssl/certs"
ENV PATH="/pipeline/:${PATH}"
ENV PATH="/pipeline/gui/daemon/:${PATH}"
ENV PATH="/pipeline/ext-scripts/:${PATH}"
RUN chmod +x /pipeline/*.py \
&& chmod +x /pipeline/gui/daemon/*.py \
&& chmod +x /boot.sh

# directories for tests
RUN chown -R primer /tests \
&& chown -R primer /pipeline \
&& chown -R root /pipeline/dictionaries/default \
&& chown -R primer /programs \
&& chown -R primer /usr/local/share/ca-certificates \
&& chown -R primer /opt/conda/lib/python3.7/site-packages/Bio/Entrez/DTDs

CMD ["/boot.sh"]
WORKDIR /primerdesign
USER primer
