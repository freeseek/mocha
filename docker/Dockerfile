# use the ubuntu base image
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
    	    	 wget \
		 autoconf \
		 zlib1g-dev \
		 bwa \
		 gzip \
		 unzip \
		 samtools \
		 bcftools \
		 r-cran-optparse \
		 r-cran-ggplot2 \
		 r-cran-data.table \
		 r-cran-gridextra \
		 git \ 
		 gcc \
		 g++ \
		 make && \
    wget https://personal.broadinstitute.org/giulio/gtc2vcf/gtc2vcf_1.0-1.deb && \
    dpkg -i gtc2vcf_1.0-1.deb && \
    rm gtc2vcf_1.0-1.deb && \
    wget https://personal.broadinstitute.org/giulio/mocha/mochatools_1.0-1.deb && \
    dpkg -i mochatools_1.0-1.deb && \
    rm mochatools_1.0-1.deb && \
    wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz && \
    mkdir /opt/iaap-cli && \
    tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz -C /opt/iaap-cli iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1 && \
    ln -s /opt/iaap-cli/iaap-cli/iaap-cli /usr/local/bin/ && \
    rm iaap-cli-linux-x64-1.1.0.tar.gz && \
    wget https://downloads.thermofisher.com/APT/APT2.11.0/apt_2.11.0_linux_64_bit_x86_binaries.zip && \
    unzip -ojd /usr/local/bin apt_2.11.0_linux_64_bit_x86_binaries.zip apt_2.11.0_linux_64_bitx86_binaries/bin/apt-probeset-genotype && \
    chmod a+x /usr/local/bin/apt-probeset-genotype && \
    rm apt_2.11.0_linux_64_bit_x86_binaries.zip && \
    wget -O /usr/local/bin/eagle https://data.broadinstitute.org/alkesgroup/Eagle/downloads/dev/eagle_v2.4.1 && \
    chmod a+x /usr/local/bin/eagle && \
    cd /tmp/ && \
    git clone git://github.com/statgen/Minimac3.git && \
    sed -i 's/USER_WARNINGS ?= -Werror/USER_WARNINGS ?= -Wno-format-truncation/' Minimac3/Library/libStatGenForMinimac3/general/Makefile && \
    sed -i 's/bool legacy_count = 0/int legacy_count = 0/' Minimac3/Library/libStatGenForMinimac3/general/Parameters.cpp && \
    sed -i 's/"\([0-9][0-9]*\)"/"\1","chr\1"/g;s/,"X"/,"X","chrX"/;s/finChromosome=="X"/(finChromosome=="X" || finChromosome=="chrX")/;s/finChromosome!="X"/(finChromosome!="X" \&\& finChromosome!="chrX")/' Minimac3/src/HaplotypeSet.cpp && \
    sed -i 's/rHap.finChromosome!="X"/rHap.finChromosome!="X" \&\& rHap.finChromosome!="chrX"/' Minimac3/src/Imputation.cpp && \
    cd Minimac3 && make && cd /tmp && \
    cp /tmp/Minimac3/bin/Minimac3 /usr/local/bin/ && \
    cp /tmp/Minimac3/bin/Minimac3-omp /usr/local/bin/ && \
    apt-get -qqy purge autoconf automake gcc g++ make && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* \
	   /tmp/Minimac3 