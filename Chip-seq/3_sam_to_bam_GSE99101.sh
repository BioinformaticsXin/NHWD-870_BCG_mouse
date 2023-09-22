##GSE99101
wkd=/HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101 #设置工作目录
cd $wkd/align
ls *.sam|cut -d"." -f 1 |while read id ;do
    samtools view -@ 10 -S $id.sam -1b -o $id.bam
    samtools sort -@ 10 -l 5 -o $id.sort.bam $id.bam
done

ls *.sort.bam |xargs -i samtools index {}

# cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101
# nohup sh 3_sam_to_bam_GSE99101.sh > sam_to_bam &

cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101/align
samtools sort -@ 10 -l 5 -o SRR5578767.sort.bam SRR5578767.bam
samtools index SRR5578767.sort.bam

