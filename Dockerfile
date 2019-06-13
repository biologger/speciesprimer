FROM biologger/speciesprimerdeps:V1

MAINTAINER biologger

# Copy the directory contents into the docker directory
COPY pipeline /home/pipeline
ENV PATH="/home/pipeline/bin/:${PATH}"
ENV PATH="/home/pipeline/ext-scripts/:${PATH}"
RUN chmod +x /home/pipeline/bin/*.py
RUN chmod +x /home/pipeline/ext-scripts/*.py

# Workdir
WORKDIR /home/primerdesign
