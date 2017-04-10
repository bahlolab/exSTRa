# STR expansion on WGS score round 2
version="01"
#base_filesystem=/Volumes/tankard
base_filesystem=/home/users/allstaff/tankard
bams="$base_filesystem/SNPchipdata/MPS_samples/MCRI/AGIP/kinghorn_2017_01_17/repexp_kinghorn_2017_01_17/devel/repexp_2017_01_17_pipeline/bam_recal/*.bam"
repeat_database=$base_filesystem/projects/research/STRs/disorders/repeat_disorders.xlsx
mkdir -p output

perl extract_real_data_read_score.pl $bahlolab_db/hg19/standard_gatk/hg19.fa "$repeat_database" $bams > output/repeat_rediscovery_wgs_2017_01_17_score_${version}.txt 2> output/repeat_rediscovery_wgs_2017_01_17_score_${version}.txt.stderr.txt

