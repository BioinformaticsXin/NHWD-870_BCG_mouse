##GSE131097
wkd=/HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE131097 #设置工作目录
cd $wkd/align
ls *.sam|cut -d"." -f 1 |while read id ;do
    samtools view -@ 10 -S $id.sam -1b -o $id.bam
    samtools sort -@ 10 -l 5 -o $id.sort.bam $id.bam
done

ls *.sort.bam |xargs -i samtools index {}

# cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE131097
# nohup sh 3_sam_to_bam_GSE131097.sh > sam_to_bam &


