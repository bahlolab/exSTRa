# Extract that data from those BAMs
version=01
#base_filesystem=/Volumes/tankard
base_filesystem=/home/users/allstaff/tankard
bams="$base_filesystem/SNPchipdata/MPS_samples/MCRI/AGIP/kinghorn_2017_01_17/repexp_kinghorn_2017_01_17/devel/repexp_2017_01_17_pipeline/bam_recal/*.bam"
repeat_database=$base_filesystem/projects/research/STRs/disorders/repeat_disorders.xlsx
mkdir -p output_by_alignment

perl extract_real_data_info_with_module.pl --by_alignment $bahlolab_db/hg19/standard_gatk/hg19.fa "$repeat_database" $bams \
   > output_by_alignment/repeat_rediscovery_wgs_2017_01_17_byalignment_${version}.txt \
   2> output_by_alignment/repeat_rediscovery_wgs_2017_01_17_byalignment_${version}.txt.stderr.txt

