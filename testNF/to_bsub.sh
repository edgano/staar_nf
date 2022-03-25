#!/usr/bin/env bash

# temporary cache dir used by singularity when pulling images from dockerhub
# this is different than the nextflow singularity cache dir when it caches images (defined in NF conf)
export SINGULARITY_CACHEDIR="$PWD/singularity_cache"
mkdir -p $SINGULARITY_CACHEDIR

# tmp dir used by singularity when pulling images from dockerhub
export TMPDIR="$PWD/tmpdir"
mkdir -p $TMPDIR

nextflow run main.nf \
     --singularity_use_pre_cached_images \
	 -profile lsf \
	 -resume

