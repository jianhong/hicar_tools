FROM nfcore/base:dev
LABEL authors="Jianhong Ou, Yu Xiang, Yarui Diao" \
      description="Docker image containing all software requirements for the nf-core/hicar pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a
RUN apt-get update --fix-missing && \
  apt-get install --yes gcc && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*
RUN /opt/conda/envs/nf-core-hicar-1.0.1/bin/pip install cooler

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-hicar-1.0.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-hicar-1.0.1 > nf-core-hicar-1.0.1.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
