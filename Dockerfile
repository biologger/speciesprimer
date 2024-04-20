FROM biologger/speciesprimerdeps:V2.1

LABEL maintainer="biologger@protonmail.com"

# Copy the directory contents into the docker directory
COPY . /
COPY boot.sh /

# include ncbi datasets and dataformats for genome download to keep speciesprimer alive for now
RUN curl -o /programs/datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' && \
    curl -o /programs/dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat' && \
    chmod +x /programs/datasets /programs/dataformat

RUN pip3 install --no-cache-dir jsonlines pandas

# Set env variables and change mod
ENV FLASK_APP="/pipeline/gui/speciesprimergui.py"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH="/pipeline/:${PATH}"
ENV PATH="/pipeline/gui/daemon/:${PATH}"
ENV PATH="/pipeline/ext-scripts/:${PATH}"
RUN chmod +x /pipeline/*.py
RUN chmod +x /pipeline/gui/daemon/*.py
RUN chmod +x /boot.sh

CMD ["/boot.sh"]
WORKDIR /primerdesign
