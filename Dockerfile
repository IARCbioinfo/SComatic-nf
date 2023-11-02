################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="SComatic-nf"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for SComatic-nf"
LABEL about.home="http://github.com/IARCbioinfo/SComatic-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/SComatic-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/SComatic-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **nalcala** <**alcalan@iarc.who.int**>

################## INSTALLATION ######################
COPY environment.yml /
COPY requirements.txt /
COPY r_requirements_install.R /
RUN conda env update -n root -f /environment.yml && conda clean -a /
RUN pip install -r requirements.txt /
RUN Rscript r_requirements_install.R
