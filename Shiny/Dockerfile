FROM rocker/shiny:3.6.1

# # system libraries of general use
 RUN apt-get update && apt-get install -y \
     libssl-dev \
     libmagick++-dev

RUN R -e "install.packages(c('DT', 'htmlwidgets', 'ggplot2', 'gridExtra', 'reticulate', 'networkDynamic', 'RcolorBrewer', 'ndtv', 'tsna', 'shinydashboard', 'shinyBS', 'shinyjs', 'zip'), repos='http://cran.rstudio.com/')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE))  install.packages('BiocManager')"
RUN R -e "BiocManager::install('topGO')"

# Download and install miniconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Export path for conda
ENV PATH="/opt/conda/bin:${PATH}"

# Create Inter-ViSTA env
RUN conda create -n intervista biopython requests pandas

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY app.R /srv/shiny-server/app.R
COPY AppFiles /srv/shiny-server/AppFiles

ENV PATH_WITH_SPACE "Sample Datasets"
COPY ${PATH_WITH_SPACE} /srv/shiny-server/${PATH_WITH_SPACE}

COPY server /srv/shiny-server/server
COPY ui /srv/shiny-server/ui
COPY www /srv/shiny-server/www

EXPOSE 3838

COPY shiny-server.sh /usr/bin/shiny-server.sh

# Start the shiny server
CMD ["/usr/bin/shiny-server.sh"]
