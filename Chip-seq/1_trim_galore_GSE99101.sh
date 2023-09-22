##GSE99101
cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101
mkdir "cleandata"
cd cleandata
ls /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101/*/*.fastq.gz >config
##开始质控
wkd=/HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101/cleandata #设置工作目录
cd $wkd
cat config |while read id
do
    arr=(${id})
    fq1=${arr[0]}
    # fq2=${arr[1]} 
    echo " trim_galore cut adapters started at $(date)"
    trim_galore -q 25 --phred33 --stringency 3 --gzip -o $wkd  $fq1
    echo "trim_galore cut adapters finished at $(date)"
done

# cd /HSCR/hengya_work/lixin/Project/870_Mouse/Data/Chip-seq/GSE99101
# nohup sh 2_trim_galore_GSE99101.sh > trim_galore_out &
