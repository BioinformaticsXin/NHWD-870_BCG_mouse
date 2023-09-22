##GSE99101
cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101
mkdir "align"
cd cleandata
ls *gz|cut -d"_" -f 1 |sort -u |while read id;do
ls -lh ${id}_trimmed.fq.gz
bowtie2 -q -p 10 -x /hengya/apps/references/bowtie2-build/Mus_musculus/Mus_musculus -U ${id}_trimmed.fq.gz -S /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101/align/${id}.sam
done

# cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101
# nohup sh 3_align_GSE99101.sh > align_out &
