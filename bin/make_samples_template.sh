#!/bin/bash

set -e
set -u
set -o pipefail

ls -1 ../raw_data/*R1* | perl -npe 's/\.\.\/raw_data\///' | perl -ne 'BEGIN{print "sample\tcontrol\tfq1\tfq2\n"}; chomp; /^([^_]+)/; print "$1\tNA\t"; print "$_\t"; s/_R1_/_R2_/; print "$_\n"' > samples_template.tsv
