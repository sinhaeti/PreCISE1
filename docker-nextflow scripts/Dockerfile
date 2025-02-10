FROM ubuntu:23.04 as base_image
#===============================#
# Docker Image Configuration	#
#===============================#
LABEL org.opencontainers.image.source='https://github.com/eipm/PreCise' \
	vendor="Englander Institute for Precision Medicine" \
	description="Data processing image for the PreCISE and PanHeme assays. It includes FASTQ trimming, alignment, and variant calling." \
	maintainer="zli4004@med.cornell.edu" \
	base_image="ubuntu" \
	base_image_version="23.04"

ENV APP_NAME="PreCISE" \
	TZ='US/Eastern' \
	PROGRAMS="opt"

RUN apt-get update -y && apt-get upgrade -y && \
	apt-get install -y wget \
		build-essential \
		zlib1g-dev \
		libvcflib-dev \
		libvcflib-tools \
  		vcftools \
		bcftools \
  		tabix \
		bwa \
  		samtools \
                git \
		libbz2-dev \
		liblzma-dev \
		autoconf \
		default-jre \
	    && rm -rf /var/lib/apt/lists/*

#===========================#
# Install FASTP #
#===========================#
LABEL	tool=fastp version=0.20.0
ENV 	FASTP_VERSION 0.20.0
RUN     cd /tmp \
&&      wget https://github.com/OpenGene/fastp/archive/v${FASTP_VERSION}.tar.gz\
&&      tar xf v${FASTP_VERSION}.tar.gz \
&&      cd fastp-${FASTP_VERSION} \
&&      make \
&&      make install \
&&      rm -rf /tmp/*

#===========================#
# Install htslib #
#===========================#
ENV HTSLIB_VERSION 1.9 
RUN  cd /tmp \
&& wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
&&  tar -vxjf htslib-${HTSLIB_VERSION}.tar.bz2 \
&&  cd htslib-${HTSLIB_VERSION} \
&&  ./configure \
&&  make \
&&  make install

#===========================#
# Install mosdepth #
#===========================#
LABEL	tool=mosdepth version=0.2.9
ENV 	MOSDEPTH_VERSION 0.2.9
RUN     cd /tmp \
&&      wget https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth -P /usr/local/bin/ \
&& 		chmod 755 /usr/local/bin/mosdepth


#===========================#
# Install VarDictJava #
#===========================#
LABEL	tool=vardict version=1.8.3
ENV 	VARDICT_VERSION 1.8.3
ADD     https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.8.3/VarDict-${VARDICT_VERSION}.tar /VarDict/
RUN     cd /VarDict \
&&	   	tar xf VarDict-${VARDICT_VERSION}.tar \
&&      rm -rf /VarDict/VarDict-${VARDICT_VERSION}.tar

#===========================#
# Install VEP #
#===========================#
## I'd suggest using the docker image already available at https://hub.docker.com/r/ensemblorg/ensembl-vep

# RUN		apt-get update -y && apt-get install -y libdbi-perl libdbd-mysql-perl
# LABEL	tool=VEP
# ENV 	VEP_VERSION 109.3
# ADD		https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/${VEP_VERSION}.tar.gz /tmp/
# RUN     cd /tmp \
# &&      tar xzf ${VEP_VERSION}.tar.gz \
# &&      cd ensembl-vep-release-${VEP_VERSION} \
# &&      perl INSTALL.pl --NO_HTSLIB \
# &&	    ls -la 
# # &&	  	mv VEP /usr/local/bin/ \
# # &&      rm -rf /tmp/*

#===========================#
# Install gencore #
#===========================#
LABEL	tool=gencore
## The wget command downloads an executable. The chmod simply gives execute permissions to everyone. You don't need to compile or install. However, you'd need to make it availalbe to the second stage with COPY --from=base_image
RUN     cd /tmp \
&&      wget http://opengene.org/gencore/gencore -P /usr/local/bin/ \
&&      chmod 755 /usr/local/bin/gencore

#===========================#
# Install vawk #
#===========================#
LABEL   tool=vawk version=0.0.1
ENV     VAWK_VERSION 0.0.1
RUN     cd /tmp \
&&      wget https://github.com/cc2qe/vawk/archive/refs/tags/0.0.1.tar.gz \
&&      tar xf 0.0.1.tar.gz \
&&      cd vawk-0.0.1 \
&&      cp vawk /usr/local/bin/ \
&&      chmod +x /usr/local/bin/vawk \
&&      rm -rf /tmp/*

#===========================#
# Second stage of build     #
#===========================#
FROM ubuntu:23.04
ENV TZ='US/Eastern'
ENV 	VARDICT_VERSION 1.8.3
RUN     apt-get update -y && apt-get upgrade -y && apt-get install -y default-jre && rm -rf /var/lib/apt/lists/*
COPY --from=base_image /usr/local/bin/ /usr/local/bin/
COPY --from=base_image /usr/bin/ /usr/bin/
COPY --from=base_image /usr/lib/x86_64-linux-gnu/ /usr/lib/x86_64-linux-gnu/
COPY --from=base_image /VarDict/VarDict-${VARDICT_VERSION}/ /VarDict-${VARDICT_VERSION}/
ENV PATH $PATH:/VarDict-${VARDICT_VERSION}/bin/
WORKDIR	/data
