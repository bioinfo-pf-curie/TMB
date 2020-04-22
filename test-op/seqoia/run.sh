#!/bin/bash

time python /data/users/nservant/GitLab/tmb/bin/pyTMB.py \
-i /data/tmp/tgutman/SeqOIA/TMB/MAP586_cancer_04022020_wes_only.vcf \
--bed /data/tmp/tgutman/SeqOIA/TMB/Merge_CORE-SPIKE_GRCh38.merged.bed \
--dbConfig /data/users/nservant/GitLab/tmb/config/snpeff.yml \
--varConfig /data/users/nservant/GitLab/tmb/config/mutect2.yml \
--filterNonCoding \
--filterSplice \
--filterSyn \
--filterPolym --polymDb 1k,gnomad \
--filterLowQual \
--export

