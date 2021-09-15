# ================================== BUILDER ===================================
FROM continuumio/miniconda3 AS builder

ENV http_proxy ${HTTP_PROXY}
ENV https_proxy ${HTTPS_PROXY}
ENV no_proxy localhost
ENV GIT_SSL_NO_VERIFY: "True"

MAINTAINER Michael Hoffmann <michael.hoffmann@ifam.fraunhofer.de>

RUN apt-get update -y
RUN apt-get install -y less
RUN apt-get install -y nano
RUN apt-get install -y build-essential

WORKDIR /app

COPY environment.yml /app/

RUN git clone -b unger_add_correlation --recurse-submodules https://github.com/BAMresearch/ModelCalibration.git

RUN conda update --update-all
RUN conda env update --name root --file environment.yml

WORKDIR /app/ModelCalibration

RUN git pull --recurse-submodules
RUN pip install -e BayesianInference

# ================================= PRODUCTION =================================
FROM builder as production

WORKDIR /app

RUN useradd -m sid
RUN chown -R sid:sid /app
USER sid
ENV PATH="/app:/home/sid/.local/bin:${PATH}"
RUN echo "alias ll='ls -alF'" >> ~/.bashrc

CMD [ "/bin/bash" ]
