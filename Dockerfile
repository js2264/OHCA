ARG BIOC_VERSION
ARG VERSION

# Pull the base image from bioconductor/bioconductor_docker
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}

LABEL org.opencontainers.image.source=https://github.com/js2264/OHCA
LABEL org.opencontainers.image.documentation=https://js2264.github.io/OHCA
LABEL org.opencontainers.image.authors="OHCA authors"
LABEL org.opencontainers.image.description="Orchestrating Hi-C analysis with Bioconductor"
LABEL org.opencontainers.image.licenses=MIT
LABEL org.opencontainers.image.version ${VERSION}

# Install quarto
# Install OHCA package suite
# Install OHCA package itself
# Print info 
# Cleanup files
WORKDIR /opt/OHCA
COPY . /opt/OHCA
RUN make 
