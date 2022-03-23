FROM zilinli/staarpipeline:0.9.6
RUN echo -ne '\n' | /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" && eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
ENV PATH=/home/linuxbrew/.linuxbrew/bin:$PATH

RUN brew install xsv
RUN eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
RUN which xsv
RUN which Rscript
RUN Rscript -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'
