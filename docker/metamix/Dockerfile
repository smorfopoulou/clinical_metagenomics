FROM ubuntu
ENV DEBIAN_FRONTEND=noninteractive
RUN apt -y  update && apt -y upgrade && apt -y install r-base openmpi-bin libopenmpi-dev r-cran-rmpi
COPY metaMix_0.4.tar.gz .
RUN R -e "install.packages('data.table', repos='https://cloud.r-project.org/')"\
    && R -e "install.packages('gtools', repos='https://cloud.r-project.org/')"\
    && R -e "install.packages('ggplot2', repos='https://cloud.r-project.org/')"\
    && R -e "install.packages('metaMix_0.4.tar.gz', type='source')"
