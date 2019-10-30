FROM biologger/speciesprimerdeps:simple

MAINTAINER biologger

# Copy the directory contents into the docker directory
COPY . /
COPY boot.sh /
# Set env variables and change mod
ENV FLASK_APP="/pipeline/app/speciesprimergui.py"
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH="/pipeline/:${PATH}"
ENV PATH="/pipeline/daemon/:${PATH}"
ENV PATH="/pipeline/ext-scripts/:${PATH}"
RUN chmod +x /pipeline/*.py
RUN chmod +x /pipeline/daemon/*.py
RUN chmod +x /boot.sh
RUN chmod +x /pipeline/ext-scripts/*.py

CMD ["/boot.sh"]
WORKDIR /home/primerdesign
