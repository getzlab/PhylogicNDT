FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.5

# adapted from Phylogic's Dockerfile
RUN apt-get update && apt-get install -y r-base r-base-dev graphviz python-tk libgraphviz-dev python2-dev libxml2-dev libxslt-dev python-is-python2
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py && python2 get-pip.py && rm get-pip.py
RUN python2 -m pip install setuptools wheel
RUN python2 -m pip install intervaltree==2.1.0 scikit-learn==0.18.1 networkx==1.11 seaborn numpy scipy matplotlib pandas
RUN pip install -e git+https://github.com/rmcgibbo/logsumexp.git#egg=sselogsumexp

WORKDIR /build/PhylogicNDT/

COPY BuildTree BuildTree/
COPY Cluster Cluster/
COPY data data/
COPY ExampleData ExampleData/
COPY ExampleRuns ExampleRuns/
COPY GrowthKinetics GrowthKinetics/
COPY LeagueModel LeagueModel/
COPY output output/
COPY PhylogicSim PhylogicSim/
COPY SinglePatientTiming SinglePatientTiming/
COPY utils utils/
COPY LICENSE LICENSE/
COPY PhylogicNDT.py .

ENV PATH=$PATH:/build/PhylogicNDT/