# Source Image
 FROM ubuntu:latest

  # Set noninterative mode
 ENV DEBIAN_FRONTEND noninteractive

  # apt update and install global requirements
 RUN apt-get clean all && \
     apt-get update && \
     apt-get upgrade -y && \
     apt-get install -y  \
         autoconf \
         build-essential \
         cmake \
         git \
         libbz2-dev \
         libcurl4-openssl-dev \
         libssl-dev \
         zlib1g-dev \
         liblzma-dev \
         tabix

  # apt clean and remove cached source lists
 RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

 RUN git clone https://github.com/samtools/htslib.git
 RUN cd htslib && \
     autoheader && \
     autoconf -Wno-syntax&& \
     ./configure --prefix=/usr/local/ && \
     make && \
     make install

 RUN git clone https://github.com/samtools/samtools.git
 RUN cd samtools && \
     autoheader && \
     autoconf -Wno-syntax&& \
     ./configure --prefix=/usr/local/ --without-curses&& \
     make && \
     make install
  # Install FASTQuick

 RUN git clone https://github.com/Griffan/FASTQuick.git
 RUN cd FASTQuick && \
     mkdir build && \
     cd build && \
     cmake .. && \
     make && \
     make test
 RUN cp /FASTQuick/bin/FASTQuick /usr/local/bin

  # Define default command
 CMD ["cd //FASTQuick/example/ && FASTQuick index --siteVCF hapmap.test.vcf.gz --dbsnpVCF dbsnp.test.vcf.gz --ref test.fa --out_prefix test_ref && FASTQuick align --fq_list test.fq.list --index_prefix test_ref --out_prefix test_out"]
