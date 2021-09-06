#!/bin/bash

set -e
set -u
set -o pipefail

find -L ../raw_data/ -name '*R1*' | perl -npe 's/\.\.\/raw_data\///' | sort | perl -ne 'BEGIN{print "sample\tcontrol\tfq1\tfq2\n"}; chomp; /([^_\/]+)_[^\/\s]+$/; print "$1\tNA\t"; print "$_\t"; s/_R1_/_R2_/; print "$_\n"' > samples_template.tsv
