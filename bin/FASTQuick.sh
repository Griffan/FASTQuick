#!/bin/bash
#
# FASTQuick: a ultra-fast full-alignment-free quality control toolkit
#
# Example ./FASTQuick.sh  -t 4 -b wgEncodeDacMapabilityConsensusExcludable.bed -r ../hg19.fa -w out -o out/FASTQuick.full.chr12.1527326.DEL1024.vcf -a out/FASTQuick.full.chr12.1527326.DEL1024.assembly.bam -j ../target/FASTQuick-2.6.2-FASTQuick-jar-with-dependencies.jar --jvmheap 8g chr12.1527326.DEL1024.bam
set -o errexit -o pipefail -o noclobber -o nounset
last_command=""
current_command=""
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
trap 'echo "\"${last_command}\" command completed with exit code $?."' EXIT
! getopt --test > /dev/null
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo '`getopt --test` failed in this environment.'
    exit 1
fi
USAGE_MESSAGE="
Usage: FASTQuick.sh --candidateVCF <1000g.phase3.site.vcf> --reference <reference.fa> --output <output.prefix> --dbSNP <dbSNP.vcf.gz> [--workingdir <directory>] [--callableRegion <callableRegion.bed>] [--targetRegion <targetRegion.bed>] [--steps All|AllButIndex|Index|Align|ContaminationAndAncestry] [--fastqList <one_pair_of_fq_or_single_fq_per_line>]

	-l/--candidateVCF:  VCF format candidate variant list to choose from.
	-r/--reference: reference genome fasta file to use.
	-t/--targetRegion: target region to focus on.
	-o/--output: output VCF.
	-d/--dbSNP: location of the dbSNP vcf file.
	-w/--workingdir: directory to place FASTQuick intermediate and temporary files. .FASTQuick.working subdirectories will be created. Defaults to the current directory.
	-c/--callableRegion: BED file containing regions to consider
	-s/--steps: processing steps to run. Defaults to AllButIndex steps. Multiple steps are specified using comma separators
	-f/--fastqList: tab-separated list of fastq files, one pair of fq files or single fq files per line,
	"


OPTIONS=l:r:o:d:t:w:c:s:f:
LONGOPTS=candidateVCF:,reference:,output:,dbSNP:,targetRegion:,workingdir:,callableRegion:,steps:,fastqList:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
	echo "$USAGE_MESSAGE" 1>&2
    exit 2
fi
eval set -- "$PARSED"
workingdir="."
candidateVCF=""
reference=""
output_prefix=""
dbSNP=""
threads=$(nproc)
callableRegion=""
targetRegion=""
fastqList=""
metricsrecords=10000000
steps="AllButIndex"

while true; do
  case "$1" in
    -l|--candidateVCF)
            candidateVCF="$2"
            shift 2
            ;;
    -r|--reference)
            reference="$2"
            shift 2
            ;;
    -s|--steps)
            steps="${2,,}"
            shift 2
            ;;
    -w|--workingdir)
            workingdir="$2"
            shift 2
            ;;
    -o|--output)
            output_prefix="$2"
            shift 2
            ;;
    -d|--dbSNP)
            dbSNP="$2"
            shift 2
            ;;
    -c|--callableRegion)
            callableRegion="$2"
            shift 2
            ;;
    -f|--fastqList)
            fastqList="$2"
            shift 2
            ;;
    -t|--targetRegion)
            targetRegion="$2"
            shift 2
            ;;
       --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

do_index=false
do_align=false
do_cont_anc=false
echo "--steps:$steps"
if [[ "$steps" == *"allbutindex"* ]] ; then
	do_align=true
	do_cont_anc=true
elif [[ "$steps" == *"all"* ]] ; then
  echo "Running All steps at once is not suggested because Index step only requires to be run once!"
	do_index=true
	do_align=true
	do_cont_anc=true
elif [[ "$steps" == *"index"* ]] ; then
	do_index=true
elif [[ "$steps" == *"align"* ]] ; then
	do_align=true
elif [[ "$steps" == "contaminationandancestry" ]] ; then
	do_cont_anc=true
fi


##### --workingdir
echo "Using working directory \"$workingdir\"" 1>&2
if [[ "$workingdir" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Working directory must be specified. Specify using the --workingdir command line argument" 1>&2
	exit 3
fi
if [[ "$(tr -d ' 	\n' <<< "$workingdir")" != "$workingdir" ]] ; then
		echo "workingdir cannot contain whitespace" 1>&2
		exit 16
	fi
if [[ ! -d $workingdir ]] ; then
	if ! mkdir -p $workingdir ; then
		echo Unable to create working directory $workingdir 1>&2
		exit 2
	fi
fi
workingdir=$(dirname $workingdir/placeholder)
##### --reference
echo "Using reference genome \"$reference\"" 1>&2
if [[ "$reference" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Reference genome must be specified. Specify using the --reference command line argument" 1>&2
	exit 6
fi
if [ ! -f $reference ] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing reference genome $reference. Specify reference location using the --reference command line argument" 1>&2
	exit 6
fi


##### --dbSNP
if [[ $do_index == "true" ]] ; then
	if [[ "$dbSNP" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Specify dbSNP location using the --dbSNP command line argument." 1>&2
		exit 4
	fi
fi

##### --output
	if [[ "$output_prefix" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Output VCF not specified. Use --output to specify output prefix." 1>&2
		exit 9
	fi
	mkdir -p $(dirname $output_prefix) || echo "Unable to create directory $(dirname $output_prefix) for output files." 1>&2
	echo "Using output prefix $output_prefix" 1>&2

##### --threads
if [[ "$threads" -lt 1 ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument" 1>&2
	exit 10
fi
if [[ "$threads" -gt 4 ]] ; then
	echo "WARNING: FASTQuick scales sub-linearly at high thread count. Up to 4 threads is the recommended level of parallelism." 1>&2
fi
echo "Using $threads worker threads." 1>&2
if [[ "$callableRegion" == "" ]] ; then
	callableRegion_arg=""
	echo "Using no callableRegion bed. The 20141020.strict_mask.whole_genome.bed callableRegion is recommended for hg19." 1>&2
elif [[ ! -f $callableRegion ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing callableRegion file $callableRegion" 1>&2
	exit 11
fi

if [[ ! -f $fastqList ]] ; then
	echo "Input file $f does not exist"  1>&2
	exit 12
fi


# Validate tools exist on path
for tool in tabix sort bcftools; do
	if ! which $tool >/dev/null; then echo "Error: unable to find $tool on \$PATH" 1>&2 ; exit 2; fi
	echo "Found $(which $tool)" 1>&2
done


timestamp=$(date +%Y%m%d_%H%M%S)
mkdir -p $workingdir
logfile=$workingdir/FASTQuick.full.$timestamp.$HOSTNAME.$$.log
timinglogfile=$workingdir/FASTQuick.timing.$timestamp.$HOSTNAME.$$.log

ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles
unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device
echo "Max file handles: $(ulimit -n)" 1>&2

echo "$(date)	Running FASTQuick. The full log is in $logfile"

#TODO:start from here tomorrow
FASTQuick_BIN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
FASTQuick_PROGRAM="$FASTQuick_BIN_DIR/FASTQuick"

if [[ $do_index == true ]] ; then
	  if [[ "$targetRegion" == "" ]] ; then
		echo "$(date)	Start indexing on whole genome..."| tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			$FASTQuick_PROGRAM index \
			--siteVCF $candidateVCF \
			--dbsnpVCF $dbSNP \
			--ref $reference  \
			--mask $callableRegion \
			--out_prefix $output_prefix \
		; } 1>&2 2>> $logfile
		else
		echo "$(date)	Start indexing on target region..."| tee -a $timinglogfile
		{ /usr/bin/time -a -o $timinglogfile \
			$FASTQuick_PROGRAM index \
			--siteVCF $candidateVCF \
			--dbsnpVCF $dbSNP \
			--ref $reference  \
			--mask $callableRegion \
			--out_prefix $output_prefix \
			--regionList $targetRegion \
		; } 1>&2 2>> $logfile
		fi
		echo "$(date)	Finished selecting markers..."| tee -a $timinglogfile
    if [[ -f "$output_prefix.FASTQuick.fa" ]] ; then
      echo "$(date)	Start preparing eigen space in 1000g phase3 genotype matrix..."| tee -a $timinglogfile
      #TODO:consider hg19 and hg38
      for chr in {1..22};do
        bcftools view -v snps -O z -R  ${output_prefix}.FASTQuick.fa.bed http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > ${output_prefix}.FASTQuick.fa.bed.chr${chr}.vcf.gz &
        pids[${chr}]=$!
      done

      for pid in ${pids[*]}; do
        if wait $pid; then
           echo "Process $pid finished"
        else
           echo "Process $pid failed"
        fi
      done

      echo "ls ${output_prefix}.FASTQuick.fa.bed.chr*.vcf.gz |bcftools concat -O z -f - -o ${output_prefix}.FASTQuick.fa.bed.phase3.vcf.gz"|bash

      if [[ -f "$output_prefix.FASTQuick.fa.bed.phase3.vcf.gz" ]] ; then
        echo "$(date) Extract subset of genotype matrix finished"| tee -a $timinglogfile
        echo "
        { /usr/bin/time -a -o $timinglogfile \
          $FASTQuick_PROGRAM pop+con \
          --RefVCF ${output_prefix}.FASTQuick.fa.bed.phase3.vcf.gz \
          --Reference $reference \
        ; } 1>&2 2>> $logfile"
      else
        echo "$(date) Extract subset of genotype matrix failed"| tee -a $timinglogfile
      fi
    else
      echo "$(date)	Marker selction failed..."| tee -a $timinglogfile
      exit 13
    fi
else
	echo "$(date)	Skipping indexing."| tee -a $timinglogfile
fi

if [[ $do_align == true ]] ; then
	echo "$(date)	Start analyzing	$fastqList" | tee -a $timinglogfile
	if [[ -f $fastqList ]] ; then
		{ /usr/bin/time -a -o $timinglogfile \
			$FASTQuick_PROGRAM align \
			--index_prefix $output_prefix \
			--fq_list $fastqList \
			--out_prefix ${output_prefix}_post_alignment \
		; } 1>&2 2>> $logfile
	else
		echo "$(date)	Abort analyzing as $fastqList does not exist." | tee -a $timinglogfile
	fi
	echo "$(date)	Complete analyzing	$fastqList" | tee -a $timinglogfile
else
	echo "$(date)	Skipping analyzing	$fastqList" | tee -a $timinglogfile
fi
if [[ $do_cont_anc == true ]] ; then
	echo "$(date)	Start estimating contamination and genetic ancestry..." | tee -a $timinglogfile
	if [[ -f "${output_prefix}_post_alignment.UD" ]] ; then
		{ /usr/bin/time -a -o $timinglogfile \
			$FASTQuick_PROGRAM pop+con \
			--PileupFile ${output_prefix}_post_alignment.Pileup \
			--Reference $reference \
			--SVDPrefix ${output_prefix}.FASTQuick.fa.bed.phase3.vcf.gz \
			--Output ${output_prefix}_post_alignment \
		; } 1>&2 2>> $logfile
	else
		echo "$(date)	Failed to load eigen space files:${output_prefix}_post_alignment.UD"
	fi
	echo "$(date)	Complete estimating contamination and genetic ancestry" | tee -a $timinglogfile
fi
echo "$(date)	Run complete with $(grep WARNING $logfile | wc -l) warnings and $(grep ERROR $logfile | wc -l) errors."