# Pull a base image from js2264/hicexperiment (should already contain most needed deps)
FROM js2264/hicexperiment

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

RUN make setup
RUN make install
RUN make info
RUN make render
