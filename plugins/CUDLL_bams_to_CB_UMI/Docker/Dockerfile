FROM trinityctat/fusioninspector:latest

MAINTAINER bhaas@broadinstitute.org

# Install dependencies for BAM processing and plotting
RUN pip install --break-system-packages pysam matplotlib

# Copy the CUDLL BAM to CB/UMI conversion script
COPY cudll_bams_to_cb_umi.py /usr/local/bin/
RUN chmod +x /usr/local/bin/cudll_bams_to_cb_umi.py

# Set working directory
WORKDIR /data
