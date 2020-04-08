# TMB1 : Basic SNVs only 
# TMB2 : SNVs+Indels PASS - splicing - Saved = Basic
# TMB3 : SNVs PASS - splicing - Saved + 10% allelic ratio
# TMB4 : SNVs+Indels PASS - splicing - Saved + 10% allelic ratio

## Filters: nonCoding, splicing, syn, polym, recurrent, low_depth
configs="--caller varscan --annot annovar"
filters="--filterNonCoding --filterSplice --filterSyn --minMAF 0.001 --filterPolym --filterLowQual --effGenomeSize 1590000 --debug"


VCF=/data/tmp/egirard/dragon_tests/D320/ANALYSIS/D320R01/CALLING/D320R01.hg19_multianno.vcf

OFILE=$(basename $VCF | sed -e 's/.hg19_multianno.vcf//')

# TMB1: 
cmd="time python bin/pyTMB.py -i ${VCF} --filterIndels --minVAF 5 ${filters} ${configs}"
echo $cmd > ${OFILE}_TMB1.txt
eval $cmd >> ${OFILE}_TMB1.txt

# TMB2:
cmd="time python bin/pyTMB.py -i ${VCF} --minVAF 5 ${filters} ${configs}"
#echo $cmd > ${OFILE}_TMB2.txt
#eval $cmd >> ${OFILE}_TMB2.txt

# TMB3
cmd="time python bin/pyTMB.py -i ${VCF} --filterIndels --minVAF 10 ${filters} ${configs}"
#echo $cmd > ${OFILE}_TMB3.txt
#eval $cmd >> ${OFILE}_TMB3.txt

# TMB4
cmd="time python bin/pyTMB.py -i ${VCF} --minVAF 10 ${filters} ${configs}"
#echo $cmd > ${OFILE}_TMB4.txt
#eval $cmd >> ${OFILE}_TMB4.txt
