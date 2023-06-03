# Pull the base image from bioconductor/bioconductor_docker
FROM bioconductor/bioconductor_docker:devel

ARG VERSION

LABEL org.opencontainers.image.source=https://github.com/js2264/OHCA
LABEL org.opencontainers.image.documentation=https://js2264.github.io/OHCA
LABEL org.opencontainers.image.authors="OHCA authors"
LABEL org.opencontainers.image.description="Orchestrating Hi-C analysis with Bioconductor"
LABEL org.opencontainers.image.licenses=MIT
LABEL org.opencontainers.image.version ${VERSION}

# Copy everything (minus anything specified in .dockerignore) into the image
COPY . /opt/OHCA
WORKDIR /opt/OHCA

RUN make quarto
RUN make info
RUN make setup
RUN make info
RUN make install
RUN make info
RUN make render
