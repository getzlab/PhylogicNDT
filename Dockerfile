FROM --platform=linux/amd64 python:2.7-stretch

# normal debian mirrors no longer host packages
RUN sed -i 's|deb.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i 's|security.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i '/stretch-updates/d' /etc/apt/sources.list && \
    echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99no-check-valid

RUN apt-get update && apt-get install -y \
    build-essential \
    python-dev \
    r-base \
    r-base-dev \
    git \
    graphviz \
    libgraphviz-dev \
    pkg-config \
    python-tk && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade \
    "pip<21" \
    "setuptools<45" \
    "wheel<0.35"

RUN pip install \
    "numpy<1.17" \
    "scipy<1.3" \
    "matplotlib<3" \
    "pandas<0.25"

COPY req /tmp/req
RUN pip install -r /tmp/req
RUN pip install git+https://github.com/rmcgibbo/logsumexp.git#egg=sselogsumexp
RUN mkdir /phylogicndt/
COPY PhylogicSim /phylogicndt/PhylogicSim
COPY GrowthKinetics /phylogicndt/GrowthKinetics
COPY BuildTree /phylogicndt/BuildTree
COPY Cluster /phylogicndt/Cluster
COPY SinglePatientTiming /phylogicndt/SinglePatientTiming
COPY LeagueModel /phylogicndt/LeagueModel
COPY data /phylogicndt/data
COPY ExampleData /phylogicndt/ExampleData
COPY ExampleRuns /phylogicndt/ExampleRuns
COPY output /phylogicndt/output
COPY utils /phylogicndt/utils
COPY PhylogicNDT.py /phylogicndt/PhylogicNDT.py
COPY LICENSE /phylogicndt/LICENSE
COPY req /phylogicndt/req
COPY README.md /phylogicndt/README.md