# TMB1 : Basic SNVs splicing
# TMB2 : SNVs+Indels PASS - splicing - Saved = Basic
# TMB3 : SNVs PASS - splicing - Saved + 10% allelic ratio
# TMB4 : SNVs+Indels PASS - splicing - Saved + 10% allelic ratio

## Filters: nonCoding, splicing, syn, polym, recurrent, low_depth
configs="--dbConfig /data/users/nservant/GitLab/tmb/config/annovar.yml --varConfig /data/users/nservant/GitLab/tmb/config/varscan2.yml --export --debug"
filters="--minDepth 100 --filterNonCoding --filterSplice --filterSyn --filterPolym --filterLowQual --filterRecurrence --maf 0.001 --polymDb 1k,gnomad,esp,exac --effGenomeSize 1590000"

RUN=D326
IDIR=/data/kdi_prod/dataset_all/2011015/backup/
echo -e "Barcode\tNb_Mut1\tTMB1\tNb_Mut2\tTMB2\tNb_Mut3\tTMB3\tNb_Mut4\tTMB4" > TMB_${RUN}_results.tsv
for i in $(find ${IDIR} -maxdepth 1 -type d -name "${RUN}*")
do
    SAMPLE=$(basename $i)
    VCF=${IDIR}/${SAMPLE}/CALLING/${SAMPLE}.hg19_multianno.vcf
    REC=${IDIR}/RESULTS/VARIANTS/${RUN}_table_report_tagged_tmb_final.tsv

    echo ${SAMPLE}

    awk -v sample=${SAMPLE} '$1==sample{print}' ${REC} > ${SAMPLE}_table_report_rec.tsv
    cmd="time python addRec.py -i ${VCF} -r ${SAMPLE}_table_report_rec.tsv -o ${SAMPLE}_rec.vcf"
    echo $cmd
    eval $cmd

    IVCF="${SAMPLE}_rec.vcf"
    OFILE=$(basename $VCF | sed -e 's/.hg19_multianno.vcf//')

    # TMB1: 
    cmd="time python ../../bin/pyTMB.py -i ${IVCF} --sample ${SAMPLE} --filterIndels --vaf 5 ${filters} ${configs}"
    echo $cmd
    eval $cmd >> ${OFILE}_TMB1.txt
    nbvar1=$(grep "Variants after filters=" ${OFILE}_TMB1.txt | awk -F"=" '{print $2}')
    tmb1=$(grep "TMB=" ${OFILE}_TMB1.txt | awk -F"=" '{print $2}')
    
    # TMB2:
    cmd="time python ../../bin/pyTMB.py -i ${IVCF} --sample ${SAMPLE} --vaf 5 ${filters} ${configs}"
    echo $cmd
    eval $cmd >> ${OFILE}_TMB2.txt
    nbvar2=$(grep "Variants after filters=" ${OFILE}_TMB2.txt | awk -F"=" '{print $2}')
    tmb2=$(grep "TMB=" ${OFILE}_TMB2.txt | awk -F"=" '{print $2}')
    
    # TMB3
    cmd="time python ../../bin/pyTMB.py -i ${IVCF} --sample ${SAMPLE} --filterIndels --vaf 10 ${filters} ${configs}"
    echo $cmd
    eval $cmd >> ${OFILE}_TMB3.txt
    nbvar3=$(grep "Variants after filters=" ${OFILE}_TMB3.txt | awk -F"=" '{print $2}')
    tmb3=$(grep "TMB=" ${OFILE}_TMB3.txt | awk -F"=" '{print $2}')
    
    # TMB4
    cmd="time python ../../bin/pyTMB.py -i ${IVCF} --sample ${SAMPLE} --vaf 10 ${filters} ${configs}"
    echo $cmd
    eval $cmd >> ${OFILE}_TMB4.txt
    nbvar4=$(grep "Variants after filters=" ${OFILE}_TMB4.txt | awk -F"=" '{print $2}')
    tmb4=$(grep "TMB=" ${OFILE}_TMB4.txt | awk -F"=" '{print $2}')
    
    cmd="rm ${SAMPLE}_table_report_rec.tsv ${SAMPLE}_rec.vcf *.txt"
    eval $cmd
    
    echo -e "${SAMPLE}\t${nbvar1}\t${tmb1}\t${nbvar2}\t${tmb2}\t${nbvar3}\t${tmb3}\t${nbvar3}\t${tmb3}" >> TMB_${RUN}_results.tsv
done

sort -k1,1 TMB_${RUN}_results.tsv > tmp
mv tmp TMB_${RUN}_results.tsv
