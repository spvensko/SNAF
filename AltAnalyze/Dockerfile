# syntax=docker/dockerfile:1

# AltAnalyze container with both human (Hs) and mouse (Mm) database support.
#
# Build:
#   docker build -t spvensko/altanalyze:0.7.0.1-fix-nowd-mm .
#
# Push:
#   docker push spvensko/altanalyze:0.7.0.1-fix-nowd-mm
#
# The Mm (mouse) EnsMart31 annotations are generated at build time from the
# GENCODE vM31 GTF by build_mm39_reference_db.py, then processed into
# AltAnalyze's RNASeq database format by build_mm_db.py. This avoids
# committing ~200MB of annotation files to git.

FROM python:2

WORKDIR /usr/src/app

# Debian Buster is archived; point apt to the archive mirror
RUN sed -i 's|http://deb.debian.org/debian|http://archive.debian.org/debian|g' /etc/apt/sources.list \
    && sed -i 's|http://security.debian.org/debian-security|http://archive.debian.org/debian-security|g' /etc/apt/sources.list \
    && sed -i '/buster-updates/d' /etc/apt/sources.list

RUN apt-get clean \
    && apt-get update \
    && apt-get install -y parallel libxml2-dev zlib1g-dev build-essential pkg-config python3 \
    && rm -rf /var/lib/apt/lists/*

# Build igraph C core from source (libigraph0-dev not available in Buster archive)
RUN curl -L https://github.com/igraph/igraph/releases/download/0.7.1/igraph-0.7.1.tar.gz -o /tmp/igraph.tar.gz \
    && cd /tmp && tar xzf igraph.tar.gz \
    && cd igraph-0.7.1 \
    && ./configure \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && rm -rf /tmp/igraph*

# Install AltAnalyze and its dependencies
RUN cd /usr/src/app \
    && curl -o requirements.txt https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/requirements_slim_test.txt \
    && pip install python-igraph==0.7.1.post6 \
    && pip install --no-cache-dir -r requirements.txt \
    && git clone https://github.com/nsalomonis/altanalyze.git \
    && cd altanalyze \
    && python AltAnalyze.py --species Hs --update Official --version EnsMart91 --additional all \
    && cd .. \
    && curl -o Hs.bed https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/Hs.bed \
    && mv ./Hs.bed ./altanalyze/AltDatabase/EnsMart91/ensembl/Hs

# Generate Mm (mouse) EnsMart31 annotations from GENCODE vM31 GTF at build time.
# This downloads the GTF (~50MB), parses it, and writes the annotation files
# that build_mm_db.py expects in AltDatabase/EnsMart31/ensembl/Mm/.
# --skip-fasta avoids downloading the ~800MB genome FASTA (not needed by AltAnalyze).
COPY build_mm39_reference_db.py /usr/src/app/build_mm39_reference_db.py
RUN python3 /usr/src/app/build_mm39_reference_db.py --output-dir /usr/src/app/altanalyze --skip-fasta \
    && rm -f /usr/src/app/altanalyze/gencode.vM31.primary_assembly.annotation.gtf.gz \
    && rm -f /usr/src/app/altanalyze/Mm_Alt31_db

# Build the processed Mm RNASeq database from the generated annotations.
# This creates AltDatabase/EnsMart31/Mm/RNASeq/ with:
#   - Mm_Ensembl_exons.txt (reformatted from Mm_Ensembl_exon.txt)
#   - Mm_Ensembl_junctions.txt (reformatted from Mm_Ensembl_junction.txt)
#   - chr/Mm_Ensembl_exons.bed (BED file for BAM-to-BED conversion)
#   - probeset-domain-annotations-exoncomp.txt (stub, passes checkForLocalArraySupport)
#   - platform.txt (identifies as RNASeq platform)
#   - Mm_junction_comps.txt and Mm_junction_comps_updated.txt (stubs)
# Without this step, AltAnalyze tries to download Mm_RNASeq.zip at runtime
# and fails because Singularity containers have no internet access.
COPY build_mm_db.py /usr/src/app/build_mm_db.py
RUN python /usr/src/app/build_mm_db.py /usr/src/app/altanalyze/AltDatabase

# Symlink GO-Elite data from EnsMart91 to EnsMart31.
# The goelite directory contains species-agnostic GO ontology annotations that
# AltAnalyze loads during expression building (even with --runGOElite no).
RUN cd /usr/src/app/altanalyze/AltDatabase/EnsMart31 \
    && ln -s ../EnsMart91/goelite goelite

# Patch AltAnalyze's version list to recognize EnsMart31
# This adds: "Mus musculus\tEnsMart31\tRNASeq"
RUN cd /usr/src/app/altanalyze \
    && printf 'Mus musculus\tEnsMart31\tRNASeq\n' >> Config/array_versions.txt

# Set Config/version.txt to EnsMart31 so AltAnalyze resolves Mm database
# paths correctly. The --version flag at runtime overwrites this, but having
# a valid version here ensures correct path resolution during startup checks.
RUN cd /usr/src/app/altanalyze \
    && printf 'EnsMart31\t01/01/2025\n' > Config/version.txt

# Copy local scripts (includes species-parameterized AltAnalyze.sh and prune.py)
COPY AltAnalyze.sh /usr/src/app/AltAnalyze.sh
COPY prune.py /usr/src/app/prune.py

RUN chmod 777 /usr/src/app/prune.py \
    && chmod 777 /usr/src/app/AltAnalyze.sh \
    && chmod -R 777 /usr/src/app/altanalyze

ENTRYPOINT ["/usr/src/app/AltAnalyze.sh"]