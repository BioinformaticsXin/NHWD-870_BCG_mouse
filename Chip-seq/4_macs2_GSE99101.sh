##GSE99101
##call peak
wkd=/HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101 #设置工作目录
mkdir $wkd/macs2
cd $wkd/align

cat $wkd/macs2/macs2 |while read id
do
    arr=(${id})
    treatment=${arr[0]}
    control=${arr[2]}
    outname=${arr[1]}
    # echo ${arr[1]}
    macs2 callpeak -c ${control}.sort.bam -t ${treatment}.sort.bam -q 0.05 -f BAM -g mm --keep-dup auto --bdg --outdir $wkd/macs2 -n ${outname} 2> $wkd/macs2/${outname}.log &
done

# cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101
# nohup sh 4_macs2_GSE99101.sh > macs2_out &


cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101/align
macs2 callpeak -c SRR5578785.sort.bam -t SRR5578767.sort.bam -q 0.05 -f BAM -g mm --keep-dup auto --bdg --outdir $wkd/macs2 -n Brd4_ChIPSeq_D7_Brd4_KO 2> $wkd/macs2/Brd4_ChIPSeq_D7_Brd4_KO.log &
