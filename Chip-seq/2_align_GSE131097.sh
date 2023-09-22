##GSE131097
cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE131097
mkdir "align"
cd cleandata
ls *gz|cut -d"_" -f 1 |sort -u |while read id;do
ls -lh ${id}_trimmed.fq.gz
bowtie2 -q -p 10 -x /hengya/apps/references/bowtie2-build/Homo_sapiens/Homo_sapiens -U ${id}_trimmed.fq.gz -S /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE131097/align/${id}.sam
done

# cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE131097
# nohup sh 3_align_GSE131097.sh > align_out &
