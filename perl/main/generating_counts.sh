# Generating counts from the BAM files:
#base_filesystem=/Volumes/tankard
base_filesystem=/home/users/allstaff/tankard

mkdir -p output

for bamfile in $base_filesystem/SNPchipdata/MPS_samples/MCRI/AGIP/kinghorn_2017_01_17/repexp_kinghorn_2017_01_17/devel/repexp_2017_01_17_pipeline/bam_recal/*.bam
do
    echo Counting file $bamfile 1>&2
    echo -en "$bamfile\\t"; samtools view -c -F 0x404 -q 1 "$bamfile"
done > output/read_counts_rediscovery_totals_01.txt

