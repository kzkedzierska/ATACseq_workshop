#!/bin/bash

mkdir -p ./peaks
mkdir -p ./motifs
while read line; do
	accession_number=$(echo $line | cut -f1 -d' ');
	shiftedBamFile=$(echo $line | cut -f6 -d' ');
	macs2 callpeak \
		--verbose 3 \
		--treatment ${shiftedBamFile} \
		-g mm \
		-B \
		-q 0.05 \
		--extsize 200 \
		--nomodel \
		--shift -100 \
		--nolambda \
		--keep-dup all \
		-f BAM \
		--outdir ./peaks/${accesssion_number} \
		--call-summits
	cat ./peaks/${accession_number}/NA_peaks.narrowPeak | cut -f1-3 > ./peaks/${accession_number}.bed
	bedtools intersect -abam ${shiftedBamFile} -b ./peaks/${accession_number}.bed -bed > ./peaks/${accession_number}_frip.bed
	#findMotifsGenome.pl ./peaks/${accesion_number}.bed mm10 ./motifs/${accession_number} -size 200 -p 16 -S 15 -len 8
done < sampleSheet_up.csv
