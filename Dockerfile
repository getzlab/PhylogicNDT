FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.5

# adapted from Phylogic's Dockerfile
RUN apt-get update && apt-get install -y r-base r-base-dev graphviz python-tk libgraphviz-dev python2-dev libxml2-dev libxslt-dev
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py && python2 get-pip.py && rm get-pip.py
RUN python2 -m pip install setuptools wheel
RUN python2 -m pip install intervaltree==2.1.0 scikit-learn==0.18.1 networkx==1.11 seaborn numpy scipy matplotlib pandas
RUN pip install -e git+https://github.com/rmcgibbo/logsumexp.git#egg=sselogsumexp

COPY BuildTree /build/PhylogicNDT/
COPY Cluster /build/PhylogicNDT/
COPY data /build/PhylogicNDT/
COPY ExampleData /build/PhylogicNDT/
COPY ExampleRuns /build/PhylogicNDT/
COPY GrowthKinetics /build/PhylogicNDT/
COPY LeagueModel /build/PhylogicNDT/
COPY output /build/PhylogicNDT/
COPY PhylogicSim /build/PhylogicNDT/
COPY SinglePatientTiming /build/PhylogicNDT/
COPY utils /build/PhylogicNDT/
COPY LICENSE /build/PhylogicNDT/
COPY PhylogicNDT.py /build/PhylogicNDT/

WORKDIR /build/PhylogicNDT/

ENV PATH=$PATH:/build/PhylogicNDT/