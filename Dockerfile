FROM biologger/speciesprimerdeps

MAINTAINER biologger

# Copy the directory contents into the docker directory
COPY pipeline /home/pipeline
COPY boot.sh /
# Set env variables and change mod
ENV FLASK_APP="/home/pipeline/bin/speciesprimergui.py"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH="/home/pipeline/bin/:${PATH}"
ENV PATH="/home/pipeline/ext-scripts/:${PATH}"
RUN chmod +x /home/pipeline/bin/*.py
RUN chmod +x /boot.sh
RUN chmod +x /home/pipeline/ext-scripts/*.py

CMD ["/boot.sh"]
WORKDIR /home/primerdesign
