# ATACseq workshop

Katarzyna Kedzierska
September 16, 2017
Jachranka, Poland

## Table of Content
	* [ATACseq workshop](#atacseq-workshop)
      * [Table of Content](#table-of-content)
      * [Plan](#plan)
      * [Data](#data)
      * [Before we start](#before-we-start)
      * [Quality metrics](#quality-metrics)
      * [Reads shifting](#reads-shifting)
      * [Peak calling](#peak-calling)
      * [Consensus peakset](#consensus-peakset)
         * [Enrichr](#enrichr)
      * [Differential analysis](#differential-analysis)
      * [Footprinting](#footprinting)

## Plan

1) Quality metrics, fragment size distribution and reads shifting with R package ATACseqQC.
2) Calling peaks with MACS2 
3) Creating consensus peakset with bedtools 
4) Motif search with HOMER 
5) Footprinting with R package ATACseqQC

## Data 

We will be working on data from the Enocde project. It'll be embryonic liver in 0 and 12.5 days post fertilization. This can gives insights into chromatin dynamics in mouse development.

> Systematic mapping of chromatin state landscapes during mouse development
David Gorkin, Iros Barozzi, Yanxiao Zhang, Ah Young Lee, Bin Lee, Yuan Zhao, Andre Wildberg, Bo Ding, Bo Zhang, Mengchi Wang, J. Seth Strattan, Jean M. Davidson, Yunjiang Qiu, Veena Afzal, Jennifer A. Akiyama, Ingrid Plajzer-Frick, Catherine S. Pickle, Momoe Kato, Tyler H. Garvin, Quan T. Pham, Anne N. Harrington, Brandon J. Mannion, Elizabeth A. Lee, Yoko Fukuda-Yuzawa, Yupeng He, Sebastian Preissl, Sora Chee, Brian A. Williams, Diane Trout, Henry Amrhein, Hongbo Yang, J. Michael Cherry, Yin Shen, Joseph R. Ecker, Wei Wang, Diane E. Dickel, Axel Visel, Len A. Pennacchio, Bing Ren
bioRxiv 166652; doi: https://doi.org/10.1101/166652

I downloaded the small piece data from Encode servers:

accession number | days | tissue | replicate
---------------- | ----------- | ------ | ---------
ENCFF109LQF | 12.5 | liver | 2
ENCFF146ZCO | 12.5 | liver | 1
ENCFF848NLJ | 0 | liver | 2
ENCFF929LOH | 0 | liver | 1

The reads were processed and aligned according to he following pipeline (paired end).

[Kundaje lab ATAC-seq pipeline specifications - updated July 11th 2017](./ATACSeqPipeline.pdf)

I prepared the data for this workshop by executing following script. 

```bash
while read line; do 
	access_number=$(echo $line | cut -f1 -d' '); 
	echo ${access_number}; 
	bash workshopPreparation.sh -a ${access_number} -o mouse --directory ./mouse --threads 16 -c -p 0.02;
done < sampleSheet.csv 2>&1 | tee workshopPreparation.log
```

## Before we start

We need to do two things - one, copy bam files and two, check if every package in R is already installed. 

```
rsync -avz --exclude="*.git/" USERNAME@192.168.1.111:/ngschool/2017/ATACseq_workshop ~/ngschool/ATACseq_workshop
```
Open R and run the following code. 

```R
###  checks for the installed packages 
packages.needed <- c("ATACseqQC", "Diffbind", "MotifDb", "BSgenome.Mmusculus.UCSC.mm10")
which.ones <- packaged.needed %in% installed.packages()
if (sum(which.ones) > 0) {
	source("https://bioconductor.org/biocLite.R")
}
for (f in pakgaes.needed[]) {
	biocLite(f)
}

### check again
packages.needed <- c("ATACseqQC", "Diffbind", "MotifDb", "BSgenome.Mmusculus.UCSC.mm10)
which.ones <- packaged.needed %in% installed.packages()
if (sum(which.ones) > 0) {
	print("Houston we have a problem!")
	print(paste0(paste(packages.needed[which.ones], sep=', '), " were not installed."))
}

```

## Quality metrics

```R
require(ATACseqQC)
sampleSheet <- read.csv("sampleSheet.tsv", d="\t")
bamfiles <- sampleSheet$V6
bamfiles.labels <- gsub(".bam", "", basename(bamfiles))
fragSize <- fragSizeDist(bamfiles, bamfiles.labels)

```

Fragment size distribution plot.



## Reads shifting

```R
require("BSgenome.Mmusculus.UCSC.mm10")
require(ATACseqQC)

## bamfile tags
#tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
## files will be output into outPath

outPath <- "splited"
dir.create(outPath)

## shift the bam file by the 5'ends
require("BSgenome.Mmusculus.UCSC.mm10")

gal <- readBamFile(bamfiles[1], tag=tags, asMates=TRUE)
gal1 <- shiftGAlignmentsList(gal)

shiftedBamfile <- file.path(outPath, paste0(bamfiles.labels, "_shifted.bam"))
export(gal1, shiftedBamfile)

```

## Peak calling

```bash
mkdir ./peaks
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
	bedtools intersect -abam ./shifted/${accesssion_number} -b ./peaks/${accesssion_number}/NA_peaks.narrowPeak -bed > ./peaks/${accesssion_number}_frip.bed
	cat ./peaks/${accesssion_number}/NA_peaks.narrowPeak | cut -f1-3 > ${accesion_number}.bed
	findMotifsGenome.pl $f mm10 ./motifs/${accession_number} -size 200 -p 2 -S 15 -len 8
done < sampleSheet.csv
```
FRiP

bedtools intersect -abam ./shifted/${accesssion_number} -b ./peaks/${accesssion_number}/NA_peaks.narrowPeak -bed > ./peaks/${accesssion_number}_frip.bed




```R
require(ChIPseeker)
require(rtracklayer)

# To import narrowPeak files
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peakFiles <- list.files(pattern = "./peaks/*/NA_peaks.narrowPeak")
peakGRs <- lapply(peakFiles, import, format = "BED", extraCols = extraCols_narrowPeak)

saveRDS <- (peaksGRs, "peaksGRs.rds")
peakGRs <- readRDS("peaksGRs.rds")

peakAnnoList <- lapply(peaksGRs, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)

```

## Consensus peakset

```bash
wc -l ./peaks/*/NA_peaks.narrowPeak
bedtools intersect -wa -wb -q 0.5 r -a ./peaks/ENCFF848NLJ/NA_peaks.narrowPeak -b ./peaks/ENCFF929LOH/NA_peaks.narrowPeak | sort -k1.1 -k2,2n | bedtools merge > embryoLiver_12.5.bed
bedtools intersect -wa -wb -q 0.5 r -a ./peaks/ENCFF848NLJ/NA_peaks.narrowPeak -b ./peaks/ENCFF929LOH/NA_peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge> embryoLiver_0.bed
wc -l ./peaks/*.bed
```
### Enrichr

Use the bed files to run Enrichr on 2k closest genes.

[embryoLiver_12.5days]()
[embryoLiver_0days]()

## Differential analysis

```R

require(Diffbind)
exp <- dba(sampleSheet = "./sampleSheet.csv")
plot(exp)
dba.ap(exp, mode=DBA_OLAP_RATE)
exp <- dba.count(exp)
exp <- dba.contrast(exp, minMembers=2)
exp <- dba.analyze(exp)
plot(exp)
dba.plotPCA(exp, DBA_FACTOR, label=DBA_TISSUE)
saveRDS(object = exp, file = "exp.rds")
dba.plotMA(exp)
dba.plotMA(exp, contrast=2)
dba.plotHeatmap(exp, contrast=2)
dba.plotHeatmap(exp, contrast=2,correlations=FALSE)

```

## Footprinting

```R
library(MotifDb)
require(ATACseqQC)
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
print(CTCF[[1]], digits=2)

factorFootprints(shiftedBamfile, pfm=CTCF[[1]], 
                 genome=genome,
                 min.score="95%", seqlev=seqlev,
                 upstream=100, downstream=100)

```