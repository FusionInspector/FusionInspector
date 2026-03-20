# FusionInspector Docker Image

## Current Version
- **Image**: `trinityctat/fusioninspector:2.11.1`
- **Base**: Ubuntu 24.04 LTS
- **Build Date**: March 2026

## Key Components

### System Environment
- **Ubuntu**: 24.04 LTS
- **Python**: 3.12 (with PEP 668 externally-managed-environment)
- **R**: 4.4.2
- **gcc**: 13.x

### R Packages
- BiocManager
- argparse
- **tidyverse** (full installation with system dependencies)
- cowplot
- ranger

### Bioinformatics Tools
- **STAR**: 2.7.11b
- **Trinity**: 2.15.2 (commit 4be803497fd22ce8461a9637eff46bc3b75a594a)
- **Samtools**: 1.20
- **Bowtie2**: 2.5.4
- **Salmon**: 1.10.0
- **Picard**: 3.2.0
- **minimap2**: 2.26 (with k8 and paftools.js)
- **Jellyfish**: 2.2.7

### Perl Modules
- PerlIO::gzip
- Set::IntervalTree
- DB_File
- URI::Escape
- Carp::Assert
- JSON::XS

### Python Packages
- requests
- igv-reports 1.8.0
- numpy

## Building the Docker Image

```bash
cd Docker
./build_docker.sh
```

This builds both versioned and latest tags:
- `trinityctat/fusioninspector:2.11.1`
- `trinityctat/fusioninspector:latest`

## Pushing to Docker Hub

```bash
cd Docker
./push_docker.sh
```

## Important Build Notes

### Trinity Compilation
- Requires **autoconf 2.69** specifically (not 2.71+) for bamsifter htslib compatibility
- Must use **git clone --recursive** (not tarball) for proper submodule initialization
- Trinity kept in `/usr/local/src/trinityrnaseq` with `TRINITY_HOME` and PATH settings

### FusionInspector Installation
- FusionInspector kept in `/usr/local/src/FusionInspector`
- **Do NOT copy PerlLib or PyLib** to avoid conflicts with Trinity's versions
- PATH includes both `/usr/local/src/FusionInspector/util` and `/usr/local/src/FusionInspector`

### PATH Ordering
Critical for using patched versions of utilities:
```
/usr/local/src/FusionInspector/util:/usr/local/src/FusionInspector:/usr/local/src/trinityrnaseq:...:/usr/local/bin
```

This ensures:
- FusionInspector's **patched paftools.js** (handles fusion transcript CDS coordinates) takes precedence over minimap2's version at `/usr/local/bin/paftools.js`
- Programs find their own PerlLib/PyLib relative to their location

### Python 3.12 / PEP 668
Ubuntu 24.04 uses Python 3.12 with externally-managed-environment protection. All pip installs require `--break-system-packages` flag (safe in Docker containers).

### Tidyverse Dependencies
Full tidyverse requires these system packages:
- Arrow C++ libraries (cmake, libutf8proc-dev, liblz4-dev, libzstd-dev, libsnappy-dev, libbrotli-dev)
- Font/graphics libraries (libfontconfig1-dev, libwebp-dev, libharfbuzz-dev, libfribidi-dev)
- Standard build tools (libgit2-dev, libxml2-dev, libssl-dev, libcurl4-openssl-dev)

## Troubleshooting

### Trinity bamsifter fails to compile
**Issue**: `HTSCODECS_VERSION_TEXT undeclared` or similar errors
**Solution**: Ensure autoconf 2.69 is installed before Trinity build

### paftools.js error on fusion transcripts
**Issue**: `Error: inconsistent thick start or end for transcript`
**Solution**: Verify `/usr/local/src/FusionInspector/util` is in PATH before `/usr/local/bin`

### SAM_entry.pm method not found
**Issue**: `Can't locate object method "get_read_span"`
**Solution**: Don't copy PerlLib between Trinity and FusionInspector; programs must use their own versions

### Salmon not found in PATH
**Issue**: `which salmon` fails
**Solution**: Verify symlink at `/usr/local/bin/salmon` points to `/usr/local/src/salmon-latest_linux_x86_64/bin/salmon`

## Version History

See `../Changelog.txt` for detailed version history and feature changes.
