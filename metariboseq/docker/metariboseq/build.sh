#!/bin/bash

docker build --platform=linux/amd64 \
    -t ghcr.io/bhattlab/bhattlab-metariboseq .
