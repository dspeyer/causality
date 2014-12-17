#!/bin/sh

ROOT=`echo $FN | sed 's/.sra//'`

./bin/fastq-dump  --fasta 0  -O ${ROOT}.fasta --split-spot ${ROOT}.sra 2>${ROOT}.errs1

./bin/blastn -db 16SMicrobial -query ${ROOT}.fasta/${ROOT}.fasta -outfmt "6 qseqid stitle evalue" -max_target_seqs 1 > ${ROOT}.blasout 2>${ROOT}.errs2

cat ${ROOT}.blastout  | awk '{print $2 " " $3}' | tr -d '][' | sort | uniq -c > ${ROOT}.results

