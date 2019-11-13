FROM continuumio/miniconda3

LABEL maintainer "Rick Tankard <rick.tankard@murdoch.edu.au>" \
      version "0.1" \
      description "circexplorer Conda environment"

ADD circexplorer_env.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml \
    && conda clean --all \
    && echo "source activate $(grep '^name: ' /tmp/environment.yml | cut -d' ' -f2)" >> ~/.bashrc

# Pull the environment name out of the environment.yml
ENV PATH /opt/conda/envs/$(grep '^name: ' /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH

