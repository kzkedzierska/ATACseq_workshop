#!bin/bash 
### Author: Katarzyna Kedzierska
### Jachranka, September 2017
### for ATAC-seq  workshop purposes

error()
{
    echo -e "$@" 1>&2
    usage_and_exit 1
}

usage() 
{     
    echo -e "Usage: 
    bash filter.sh -a ACCESSION_NUMBER -o ORGANISM [-b BASENAME] [-f FLAGS_OUT] [-i FLAGS_IN] [-q QUALITY] [-ckr]

Filters the alignment based on quality and samtools flags, removes duplicates and blacklisted regions if specified. On default produces only the last file, if keep option specified produces all the intermediate files.

1) Sorts the alignment and removes duplicates: _rmdup.bam
2) Removes blacklisted regions: _no.blacklist.bam
3) Filter the regions based on quality and flags: _filtered.bam
4) Index the filtered file and filter only the canonical chromosomes: _cleaned.bam
5) Subsample with a given percentage of a file: _part.bam 

Options required
	-a / --accession ACCESSION_NUMBER 	an accession number or sample name
    -o / --organism ORGANISM			organism, either human or mouse
Optional
	-b / --basename BASENAME	if provided, output files will be named with BASENAME, default: ACCESSION_NUMBER
	-e / --directory INPUT_DIR	directory in which the script looks for bams and were it saves outpu [./]
	-f / --flags FLAGS_OUT 		flags used to filtered out alignments (-F samtools view option) [3844]
	-i / --flagsin FLAGS_IN		flags used to keep alignments (-f samtools view option), default [2]
	-q / --quality QUALITY		quality used to filtered out the reads [30]
	-k / --keep					keep the intermediate files [keeps only the _cleaned and _part files]
	-d / --duplicates			removes duplicates , default [skips]
	-s / --skip-sorting			only valid if removing duplicates, skips sorting, if not specified alignment will be sorted based on coordinates
	-c / --chromosome			should be used if alignment has chr1, chr2 etc. names for chromosomes, otherwise 1,2 etc. will be used.	
	-p / --part	PART			if specified, the bam will be subsampled for a fraction (0 < PART < 1) or a specific chromosome (1 <= PART <= max for an organism).
	-t / --threads THREADS 		number of threads used for sorting and saving bams"	
}  

usage_and_exit()
{
    usage
    exit $1
}

# Sets defaults
ACCESSION_NUMBER=
QUALITY=30
ORGANISM=
KEEP=
CHR=
FLAGS_OUT=3844
FLAGS_IN=2
THREADS=1
INPUT_DIR="./"

# Option
while test $# -gt 0
do 
    case $1 in
    --accession | -a )
        shift  
        ACCESSION_NUMBER="$1"
        ;;
    --duplicates | -d )
		DUPLICATES=true
		;;
	--directory | -e )
		shift
		INPUT_DIR=$1
		;;
	--skip-sorting | -s )
		SKIP=true
		;;
    --base | -b )  
        shift
        BASENAME="$1"
        ;;
	--chromosome | -c )
		CHR=true
		;;
	--part | -p )
		shift
		PART=$1
		;;
	--flags | -f )
		shift
		FLAGS_OUT=$1
		;;
	--flagsin | -i )
		shift
		FLAGS_IN=$1
		;;
	--keep | -k )
		KEEP=true
		;;
	--threads | -t )
		shift
		THREADS=$1
		;;
	--quality | -q )
		shift
		QUALITY=$1
		;;
    --organism | -o )
        shift
        ORGANISM="$1"
        ##works for human and mouse only
        ;;
    --help | -h | '-?' | '--?' ) 
         usage_and_exit 0
         ;;
    -*) 
        error "Unrecognized option: $1"
        ;;
    *)
        break
        ;;
    esac
    shift
done

#check for the necessary arguments
if [ -z "${ACCESSION_NUMBER}" -o -z "${ORGANISM}" ]; then
    error "Accession number or organism not specified."
fi 

### Assigns values based (chromosomes and blacklisted regions) on organism.
if [ "X${ORGANISM}" == "Xhuman" ]; then
	BLACKLISTED_REGION="/home/kzk5f/scratch/ref/hg19.blacklisted.ENCFF001TDO.bed"
	if [ "X${CHR}" == "Xtrue" ]; then
		CHROMOSOMER=$(echo $(for i in {1..22}; do echo "chr"$i; done | xargs) "chrX" "chrY");
	else 
		CHROMOSOMES=$(echo $(seq 22) X Y)
	fi
elif [ "X${ORGANISM}" == "Xmouse" ]; then
	BLACKLISTED_REGION="/home/kzk5f/scratch/ref/mm10-blacklist.bed"
	if [ "X${CHR}" == "Xtrue" ]; then
		CHROMOSOMER=$(echo $(for i in {1..19}; do echo "chr"$i; done | xargs) "chrX" "chrY");
	else 
		CHROMOSOMES=$(echo $(seq 19) X Y);
	fi
else 
	error "Only human and mouse supported!"
fi

### Removes duplicates
TMP=${INPUT_DIR}/${ACCESSION_NUMBER}.bam
if [ "X${DUPLICATES}" == "Xtrue" ]; then
	if [ "X${SKIP}" == "Xtrue" ]; then
		samtools rmdup ${TMP} ${INPUT_DIR}/${ACCESSION_NUMBER}_rmdup.bam;
	fi

	echo "samtools sort -n -@ ${THREADS} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_sortedName.bam ${TMP}";
	samtools sort -n -@ ${THREADS} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_sortedName.bam ${TMP}
	TMP=${INPUT_DIR}/${ACCESSION_NUMBER}_sortedName.bam;
	echo $TMP;
	echo "samtools rmdup ${TMP} ${INPUT_DIR}/${ACCESSION_NUMBER}_rmdup.bam"
	samtools rmdup ${TMP} ${INPUT_DIR}/${ACCESSION_NUMBER}_rmdup.bam
	TMP=${INPUT_DIR}/${ACCESSION_NUMBER}_rmdup.bam;
	echo $TMP
	if [ "X${KEEP}" != "Xtrue" -a -f ${INPUT_DIR}/${ACCESSION_NUMBER}_sortedName.bam ]; then
		rm -f ${INPUT_DIR}/${ACCESSION_NUMBER}_sortedName.bam;
	fi
fi

### Removes blacklisted region
echo "bedtools intersect -v -abam ${TMP} -b ${BLACKLISTED_REGION} > ${INPUT_DIR}/${ACCESSION_NUMBER}_no.blacklist.bam"
bedtools intersect -v -abam ${TMP} -b ${BLACKLISTED_REGION} > ${INPUT_DIR}/${ACCESSION_NUMBER}_no.blacklist.bam;
TMP=${INPUT_DIR}/${ACCESSION_NUMBER}_no.blacklist.bam
if [ "X${KEEP}" != "Xtrue" -a -f ${INPUT_DIR}/${ACCESSION_NUMBER}_rmdup.bam ]; then
	rm -f ${INPUT_DIR}/${ACCESSION_NUMBER}_rmdup.bam;
fi

### Filters based on flags and quality
echo "samtools view -h -b -@ ${THREADS} -f ${FLAGS_IN} -F ${FLAGS_OUT} -q ${QUALITY} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_filtered.bam ${TMP}";
samtools view -h -b -@ ${THREADS} -f ${FLAGS_IN} -F ${FLAGS_OUT} -q ${QUALITY} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_filtered.bam ${TMP};
if [ "X${KEEP}" != "Xtrue"  -a -f ${INPUT_DIR}/${ACCESSION_NUMBER}_no.blacklist.bam ]; then
	rm -f ${INPUT_DIR}/${ACCESSION_NUMBER}_no.blacklist.bam;
fi
TMP=${INPUT_DIR}/${ACCESSION_NUMBER}_filtered.bam
### Outputs only valid chromosomes
echo "samtools sort ${TMP} -@ ${THREADS} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_sorted.bam";
samtools sort ${TMP} -@ ${THREADS} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_sorted.bam
TMP=${INPUT_DIR}/${ACCESSION_NUMBER}_sorted.bam
echo "samtools index ${TMP}";
samtools index ${TMP}
echo "samtools view -h -b -@ ${THREADS} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_cleaned.bam ${TMP} ${CHROMOSOMES};"
samtools view -h -b -@ ${THREADS} -o ${INPUT_DIR}/${ACCESSION_NUMBER}_cleaned.bam ${TMP} ${CHROMOSOMES};
if [ "X${KEEP}" != "Xtrue" -a -f ${INPUT_DIR}/${ACCESSION_NUMBER}_filtered.bam ]; then
	rm -f ${INPUT_DIR}/${ACCESSION_NUMBER}_filtered.bam;
fi
TMP=${INPUT_DIR}/${ACCESSION_NUMBER}_cleaned.bam

if [ "X${PART}" != "X" ]; then
	if [ ${PART} -ge 1 ]; then ### /TODO/ it would be good to check for maximum for a given organism
		echo "samtools view -b -h -@ ${THREADS} -o {ACCESSION_NUMBER}_part.bam ${TMP} chr${PART}"
		samtools view -b -h -@ ${THREADS} -o {ACCESSION_NUMBER}_part.bam ${TMP} "chr${PART}"
	else
		### /TODO/ should check if float! maybe use awk or bc
		echo "samtools view -b -h -@ ${THREADS} -s ${PART} -o {ACCESSION_NUMBER}_part.bam ${TMP}"
		samtools view -b -h -@ ${THREADS} -s ${PART} -o {ACCESSION_NUMBER}_part.bam ${TMP}
	fi
fi