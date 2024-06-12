#!/bin/bash

if [ ! -d taxonomy ]; then
  mkdir -p taxonomy
  wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
  tar -xvf taxdump.tar.gz -C taxonomy  && rm -f taxdump.tar.gz
fi

