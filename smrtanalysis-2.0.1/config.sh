#!/bin/bash

. /net/eichler/vol20/projects/pacbio/opt/smrtanalysis-2.0.1/etc/setup.sh

. /etc/profile.d/modules.sh
if test ! -z $MODULESHOME; then
   module load modules modules-init modules-gs/prod modules-eichler/prod
fi

module load cross_match/latest
module load bedtools/latest
export MAKEDIR=/net/eichler/vol4/home/jlhudd/projects/pacbio/smrtanalysis-2.1.0
export VECTOR=/net/eichler/vol4/home/jlhudd/projects/pacbio/vector/vector.fasta
