FROM biologger/speciesprimerdeps:V3

MAINTAINER biologger

# Copy the directory contents into the docker directory
COPY . /home/
COPY boot.sh /
# Set env variables and change mod
ENV FLASK_APP="/home/pipeline/bin/speciesprimergui.py"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH="/home/pipeline/bin/:${PATH}"
RUN chmod +x /home/pipeline/bin/*.py
RUN chmod +x /boot.sh

CMD ["/boot.sh"]
WORKDIR /home/primerdesign
