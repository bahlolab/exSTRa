# Running exSTRa

bam_glob="path/to/bams/*.bam"
output=output/exSTRa_scores.txt

repeat_database=path/to/repeat_disorders.xlsx
reference=path/to/hg19.fa

mkdir -p $(dirname $output)

perl exSTRa_score.pl 
    "$reference" \
    "$repeat_database" 
    $bams \
    > "$output 

