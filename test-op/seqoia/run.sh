#!/bin/bash

python ./bin/pyTMB.py -i /data/tmp/tgutman/SeqOIA/TMB/MAP586_cancer_04022020_wes_only.vcf \
--filterNonCoding \
--filterSplice \
--filterSyn \
--filterPolym --polymDb 1k,gnomad\
--filterLowQual \
--effGenomeSize 50000000 \
--caller mutect --annot snpeff \
--dbConfig /data/users/nservant/GitLab/tmb/config/databases.yml --varConfig /data/users/nservant/GitLab/tmb/config/calling.yml --debug --verbose
