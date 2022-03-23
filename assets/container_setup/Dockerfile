FROM zilinli/staarpipeline:0.9.6
RUN echo -ne '\n' | /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" && eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
ENV PATH=/home/linuxbrew/.linuxbrew/bin:$PATH

RUN brew install xsv
RUN eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
RUN which xsv
RUN which Rscript
RUN Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'
# RUN brew install r
# RUN brew install wget

# RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && bash ~/miniconda.sh -b -p $HOME/miniconda && export PATH="/home/linuxbrew/miniconda/bin:$PATH" && which conda
# ENV PATH=/home/linuxbrew/miniconda/bin:$PATH
# COPY environment.yml /tmp/environment.yml
# RUN brew install allureofthestars
# RUN conda config --set always_yes yes --set changeps1 no && \
#     conda update conda --yes && \
#     conda env update -n root --file /tmp/environment.yml 

# # ENV export TAR="/bin/tar"
# # RUN tar --version
# RUN export TAR="/bin/tar" && Rscript -e 'devtools::install_github("xihaoli/STAAR"); devtools::install_github("xihaoli/STAARpipeline",ref="main");'