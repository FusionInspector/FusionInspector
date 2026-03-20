FROM trinityctat/fusioninspector:2.11.1

MAINTAINER bhaas@broadinstitute.org

# Install pysam for BAM processing
RUN pip install --break-system-packages pysam

# Copy the CUDLL BAM to FASTQ conversion script
COPY cudll_bams_to_fastq.py /usr/local/bin/
RUN chmod +x /usr/local/bin/cudll_bams_to_fastq.py

# Set working directory
WORKDIR /data
