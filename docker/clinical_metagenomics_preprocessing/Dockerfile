From centos:7
RUN yum -y upgrade
RUN yum -y install wget


############################## installing blast ##########################
RUN yum -y install https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-1.x86_64.rpm # install the current latest version

############################## installing diamond ##########################
# version speicfied in email was v0.9.30.131
RUN pushd /tmp && wget https://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-linux64.tar.gz && tar -xf diamond-linux64.tar.gz && rm diamond-linux64.tar.gz && mv diamond /usr/local/bin/ && popd

############################## installing trimg-galore ##########################
CMD [ "/bin/bash" ]
ENV PATH /opt/conda/bin:$PATH
RUN pushd /tmp && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && popd # installs conda in /root/miniconda3
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc &&  echo "conda activate base" >> ~/.bashrc
RUN /opt/conda/bin/conda config --add channels conda-forge # needed to build python-isal
RUN /opt/conda/bin/conda install -y -c bioconda trim-galore
RUN yum -y install pigz # pigz is needed to run trim-galore in parallel

############################## installing bowtie2 #########################
RUN /opt/conda/bin/conda install -y -c bioconda bowtie2


############################## installing samptools ##########################
RUN yum -y install epel-release && yum -y upgrade && yum -y install -y samtools # install samtools also implicetly install perl. samtools version is 1.9 but script uses 0.1.19 so might be incompatable

############################## installing seqtk ##########################
# seqtk Version specified in email 1.0-r45
RUN yum -y install make gcc zlib-devel # ToDo move this to build container
RUN pushd /tmp && wget https://github.com/lh3/seqtk/archive/refs/tags/v1.3.tar.gz && tar -xf v1.3.tar.gz && cd seqtk-1.3 && make && make install && cd .. && rm -rf v1.3.tar.gz seqtk-1.3 && popd # install latest version (version not specified in script)

############################## installing bedtools ##########################
RUN yum -y install BEDTools # needs epel installed too. installed version is 2.26, script specifies 2.25 thus might be incompatable

############################## installing prinseq ##########################
# potentually should run upgrade to c++ version
RUN pushd /tmp && wget -O prinseq-lite-0.20.4.tar.gz https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download && tar -xf prinseq-lite-0.20.4.tar.gz && mv prinseq-lite-0.20.4 /opt/ && chmod +x /opt/prinseq-lite-0.20.4/prinseq-lite.pl && ln -s /opt/prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq-lite.pl && popd # I took the last version ever published which is 1 minor verson greater that specified in the script


############################## dev tools ##########################
RUN yum -y install less man which htop # for development convienence, remove later
RUN yum -y install nodejs
