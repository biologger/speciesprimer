FROM specdeps:test

LABEL maintainer="biologger@protonmail.com"

# Copy the directory contents into the docker directory
COPY . /

# Set env variables and change mod
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH="/pipeline/:${PATH}"
ENV PATH="/pipeline/gui/daemon/:${PATH}"
ENV PATH="/pipeline/ext-scripts/:${PATH}"
RUN chmod +x /pipeline/*.py

SHELL [ "/bin/bash", "--login", "-c" ]
RUN conda activate base

WORKDIR /pipeline
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
