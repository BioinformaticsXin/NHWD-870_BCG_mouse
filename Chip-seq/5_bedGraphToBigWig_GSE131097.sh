##GSE131097
cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE131097
cd macs2
ls *_treat_pileup.bdg >config
cat config |while read str
do
    arr=(${str//_treat_pileup./ })  
    id=${arr[0]}
    echo $id
    bedGraphToBigWig ${id}_treat_pileup.bdg /hengya/apps/references/Homo_sapiens.GRCh38.chrom.sizes ${id}_treat_pileup.bw
done


# cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE131097
# nohup sh 5_bedGraphToBigWig_GSE131097.sh > bedGraphToBigWig_out &
