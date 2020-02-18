#
# FASTQuick: a ultra-fast full-alignment-free quality control toolkit
#
set -o errexit -o pipefail -o noclobber -o nounset
last_command=""
current_command=""
#trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
#trap 'echo "\"${last_command}\" command completed with exit code $?."' EXIT
! getopt --test > /dev/null
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo '`getopt --test` failed in this environment.'
    echo 'If you are under Mac OS, please consider to follow these steps:'
    brew -v > /dev/null
    if [ "$?" != "0" ];then
        echo 'Please install brew for Mac OS: ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"'
    fi
    echo 'Please install gnu-getopt for Mac OS: brew install gnu-getopt'
    echo 'Please add `export PATH=$(find /usr/local/Cellar/gnu-getopt/*/bin/ -type d):$PATH` at the begining of FASTQuick.sh'
    exit 1
fi
USAGE_MESSAGE="
Usage: FASTQuick.sh [--steps All|AllButIndex|Index|Align|Contamination|Ancestry|Visualize] --candidateVCF <1000g.phase3.site.vcf> --reference <reference.fa> --output <output.prefix> --index <index.prefix> --dbSNP <dbSNP.vcf.gz> --fastqList <one_pair_of_fq_or_single_fq_per_line> [--workingdir <directory>] [--callableRegion <callableRegion.bed>] [--targetRegion <targetRegion.bed>]

	-l/--candidateVCF:  VCF format candidate variant list to choose from.
	-r/--reference: reference genome fasta file to use.
	-t/--targetRegion: target region to focus on.
	-o/--output: prefix of output files.
	-i/--index: prefix of index files.
	-d/--dbSNP: location of the dbSNP vcf file.
	-w/--workingdir: directory to place FASTQuick intermediate and temporary files(.FASTQuick.working subdirectories will be created).
	-c/--callableRegion: BED file to specify callable region. For example: 20141020.strict_mask.whole_genome.bed
	-s/--steps: processing steps to run. Defaults to AllButIndex steps. Multiple steps are specified using comma separators.
	-f/--fastqList: tab-separated list of fastq files, one pair of fq files or single fq files per line.
	"


OPTIONS=l:r:o:i:d:t:w:c:s:f:
LONGOPTS=candidateVCF:,reference:,output:,index:,dbSNP:,targetRegion:,workingdir:,callableRegion:,steps:,fastqList:
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
outputPrefix="default.output"
indexPrefix="default.index"
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
            steps="$(echo $2 | tr '[:upper:]' '[:lower:]')"
            shift 2
            ;;
    -w|--workingdir)
            workingdir="$2"
            shift 2
            ;;
    -o|--output)
            outputPrefix="$2"
            shift 2
            ;;
    -i|--index)
            indexPrefix="$2"
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
do_viz=false
#echo "--steps:$steps"
if [[ "$steps" == *"allbutindex"* ]] ; then
	do_align=true
	do_cont_anc=true
elif [[ "$steps" == *"all"* ]] ; then
  echo "Running All steps at once is not suggested because indices can be reused!"
	do_index=true
	do_align=true
	do_cont_anc=true
	do_viz=true
elif [[ "$steps" == *"index"* ]] ; then
	do_index=true
elif [[ "$steps" == *"align"* ]] ; then
	do_align=true
	do_viz=true
elif [[ "$steps" == *"contamination"* ]] ; then
	do_cont_anc=true
	do_viz=true
elif [[ "$steps" == *"ancestry"* ]] ; then
	do_cont_anc=true
	do_viz=true
elif [[ "$steps" == *"visualize"* ]] ; then
	do_viz=true
fi


##### --workingdir
echo "Using working directory \"$workingdir\"" 1>&2
if [[ "$workingdir" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Working directory must be specified. Specify using the --workingdir command line argument" 1>&2
	exit 4
fi
if [[ "$(tr -d ' 	\n' <<< "$workingdir")" != "$workingdir" ]] ; then
		echo "workingdir cannot contain whitespace" 1>&2
		exit 4
	fi
if [[ ! -d $workingdir ]] ; then
	if ! mkdir -p $workingdir ; then
		echo Unable to create working directory $workingdir 1>&2
		exit 4
	fi
fi
workingdir=$(dirname $workingdir/placeholder)
##### --reference
echo "Using reference genome \"$reference\"" 1>&2
if [[ "$reference" == "" ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Reference genome must be specified. Specify using the --reference command line argument" 1>&2
	exit 5
fi
if [ ! -f $reference ] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing reference genome $reference. Specify reference location using the --reference command line argument" 1>&2
	exit 5
fi


##### --dbSNP
if [[ $do_index == "true" ]] ; then
	if [[ "$dbSNP" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Specify dbSNP location using the --dbSNP command line argument." 1>&2
		exit 6
	fi
fi

##### --output
	if [[ "$outputPrefix" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "WARNING: Output prefix not specified. Use --output to specify output prefix." 1>&2
		exit 9
	fi
	#mkdir -p $(dirname $outputPrefix) || echo "Unable to create directory $(dirname $outputPrefix) for output files." 1>&2
	if [[ "${outputPrefix:0:1}" != "/" && "${outputPrefix:0:2}" != "~[/a-z]" ]] ; then
    echo "Absolute path is required for --outputPrefix." 1>&2
    exit 9
  fi
	echo "Using output prefix $outputPrefix" 1>&2

##### --index
	if [[ "$indexPrefix" == "" ]] ; then
		echo "$USAGE_MESSAGE"  1>&2
		echo "Index prefix not specified. Use --index to specify index prefix." 1>&2
		exit 10
	fi
	#mkdir -p $(dirname $outputPrefix) || echo "Unable to create directory $(dirname $outputPrefix) for output files." 1>&2
	if [[ "${indexPrefix:0:1}" != "/" && "${indexPrefix:0:2}" != "~[/a-z]" ]] ; then
    echo "Absolute path is required for --indexPrefix." 1>&2
    exit 10
  fi
	echo "Using index prefix $indexPrefix" 1>&2

##### --threads
if [[ "$threads" -lt 1 ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Illegal thread count: $threads. Specify an integer thread count using the --threads command line argument" 1>&2
	exit 11
fi
if [[ "$threads" -gt 4 ]] ; then
	echo "WARNING: FASTQuick scales sub-linearly at high thread count. Up to 4 threads is the recommended level of parallelism." 1>&2
fi
echo "Using $threads worker threads." 1>&2
if [[ "$callableRegion" == "" ]] ; then
	callableRegion_arg=""
	echo "Skip input of callableRegion. The callableRegion(ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed) is recommended for hg19." 1>&2
elif [[ ! -f $callableRegion ]] ; then
	echo "$USAGE_MESSAGE"  1>&2
	echo "Missing callableRegion file $callableRegion" 1>&2
	exit 12
fi

if [[ ! -f $fastqList ]] ; then
	echo "Input file $f does not exist"  1>&2
	exit 13
fi


# Validate tools exist on path
for tool in sort bcftools pandoc; do
	if ! which $tool >/dev/null; then
	  echo "Error: unable to find $tool on \$PATH, please install before continue" 1>&2 ; exit 2; fi
	echo "Found $(which $tool)" 1>&2
done


timestamp=$(date +%Y%m%d_%H%M%S)
mkdir -p $workingdir
logfile=$workingdir/FASTQuick.full.$timestamp.$HOSTNAME.$$.log
timinglogfile=$workingdir/FASTQuick.timing.$timestamp.$HOSTNAME.$$.log

#ulimit -n $(ulimit -Hn) # Reduce likelihood of running out of open file handles
#unset DISPLAY # Prevents errors attempting to connecting to an X server when starting the R plotting device
echo "Max file handles: $(ulimit -n)" 1>&2

echo "$(date)	Running FASTQuick. The full log is in $logfile"

FASTQuick_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
FASTQuick_BIN_DIR="${FASTQuick_SRC_DIR}/bin/"
#FASTQuick_RESOURCE_DIR="${FASTQuick_SRC_DIR}/resource/"

FASTQuick_PROGRAM="${FASTQuick_SRC_DIR}/bin/FASTQuick"
if [[ ! -f "${FASTQuick_PROGRAM}" ]] ; then
  echo "${FASTQuick_PROGRAM} not found, please make sure execute FASTQuick.sh with full installation path..."
  exit 14
fi

###analysis start
if [[ $do_index == true ]] ; then
	  if [[ "$targetRegion" == "" ]] ; then
		echo "$(date)	Start indexing on whole genome..."| tee -a $timinglogfile
		{ /usr/bin/time \
			$FASTQuick_PROGRAM index \
			--siteVCF $candidateVCF \
			--dbsnpVCF $dbSNP \
			--ref $reference  \
			--mask $callableRegion \
			--out_prefix $indexPrefix \
		1>&2 2>> $logfile; } 1>&2 2>>$timinglogfile
		else
		echo "$(date)	Start indexing on target region..."| tee -a $timinglogfile
		{ /usr/bin/time \
			$FASTQuick_PROGRAM index \
			--siteVCF $candidateVCF \
			--dbsnpVCF $dbSNP \
			--ref $reference  \
			--mask $callableRegion \
			--out_prefix $indexPrefix \
			--regionList $targetRegion \
		1>&2 2>> $logfile; } 1>&2 2>>$timinglogfile
		fi
		echo "$(date)	Finished selecting markers..."| tee -a $timinglogfile
    if [[ -f "$indexPrefix.FASTQuick.fa" ]] ; then
      echo "$(date)	Start preparing eigen space in 1000g phase3 genotype matrix..."| tee -a $timinglogfile
      #TODO:consider hg19 and hg38
      genoPrefix=""
      genoTestFile="ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      ebiPrefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
      ucscPrefix="http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/"
      if curl --head --fail --silent "$ebiPrefix$genoTestFile" >/dev/null; then
          genoPrefix=$ebiPrefix
      elif curl --head --fail --silent "$ucscPrefix$genoTestFile" >/dev/null; then
          genoPrefix=$ucscPrefix
      fi
      echo "Extracting genotype matrices from ${genoPrefix}..."
      for chr in {1..22};do
        bcftools view -v snps -O z -R  ${indexPrefix}.FASTQuick.fa.SelectedSite.vcf "${genoPrefix}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" > ${indexPrefix}.FASTQuick.fa.bed.chr${chr}.vcf.gz &
        pids[${chr}]=$!
        sleep 1s
      done

      for pid in ${pids[*]}; do
        if wait $pid; then
           echo "Process $pid finished"
        else
           echo "Process $pid failed"
        fi
      done

      echo "ls ${indexPrefix}.FASTQuick.fa.bed.chr*.vcf.gz |bcftools concat -O z -f - -o ${indexPrefix}.FASTQuick.fa.bed.phase3.vcf.gz"|bash
      rm ${indexPrefix}.FASTQuick.fa.bed.chr*.vcf.gz

      if [[ -f "$indexPrefix.FASTQuick.fa.bed.phase3.vcf.gz" ]] ; then
        echo "$(date)	Extract subset of genotype matrix finished"| tee -a $timinglogfile
        { /usr/bin/time \
          $FASTQuick_PROGRAM pop+con \
          --RefVCF ${indexPrefix}.FASTQuick.fa.bed.phase3.vcf.gz \
          --Reference $reference \
        1>&2 2>> $logfile; } 1>&2 2>>$timinglogfile
        echo "$(date)	Build SVD files finished" | tee -a $timinglogfile
      else
        echo "$(date)	Extract subset of genotype matrix failed"| tee -a $timinglogfile
      fi
    else
      echo "$(date)	Marker selction failed..."| tee -a $timinglogfile
      exit 15
    fi
else
	echo "$(date)	Skipping indexing."| tee -a $timinglogfile
fi

if [[ $do_align == true ]] ; then
  ###check overwrite
  if [[ -f "${outputPrefix}.Summary" ]] ; then
    echo "${outputPrefix}.Summary exists, please manually remove existing files before restart..."
    exit 16
  fi
	echo "$(date)	Start analyzing	$fastqList" | tee -a $timinglogfile
	if [[ -f $fastqList ]] ; then
	  echo "$(date)	Align fastq file in: $fastqList" | tee -a $timinglogfile
		{ /usr/bin/time \
			$FASTQuick_PROGRAM align \
			--index_prefix $indexPrefix \
			--fq_list $fastqList \
			--out_prefix ${outputPrefix} \
			--q 15 \
		1>&2 2>> $logfile; } 1>&2 2>>$timinglogfile
		if [[ -f "${outputPrefix}.bam" ]]; then
		  echo "$(date)	Start sorting bam file..." | tee -a $timinglogfile
		  samtools sort -O BAM ${outputPrefix}.bam > ${outputPrefix}.sorted.bam
		  samtools index ${outputPrefix}.sorted.bam
		  rm ${outputPrefix}.bam
		else
		  echo "$(date)	Analyzing $fastqList failed!" | tee -a $timinglogfile
		  exit 16
		fi
	else
		echo "$(date)	Abort analyzing as $fastqList does not exist." | tee -a $timinglogfile
	fi
	echo "$(date)	Complete analyzing	$fastqList" | tee -a $timinglogfile
else
	echo "$(date)	Skipping analyzing	$fastqList" | tee -a $timinglogfile
fi
if [[ $do_cont_anc == true ]] ; then
	echo "$(date)	Start estimating contamination and genetic ancestry..." | tee -a $timinglogfile
	if [[ -f "${indexPrefix}.FASTQuick.fa.bed.phase3.vcf.gz.UD" ]] ; then
		{ /usr/bin/time \
			$FASTQuick_PROGRAM pop+con \
			--PileupFile ${outputPrefix}.Pileup \
			--Reference $reference \
			--SVDPrefix ${indexPrefix}.FASTQuick.fa.bed.phase3.vcf.gz \
			--Output ${outputPrefix} \
		1>&2 2>> $logfile; } 1>&2 2>>$timinglogfile
	else
		echo "$(date)	Failed to load eigen space files:${indexPrefix}.FASTQuick.fa.bed.phase3.vcf.gz.UD" | tee -a $timinglogfile
	fi
	echo "$(date)	Complete estimating contamination and genetic ancestry" | tee -a $timinglogfile
else
	echo "$(date)	Skipping estimating contamination and genetic ancestry..."| tee -a $timinglogfile
fi

if [[ $do_viz == true ]] ; then
echo "$(date)	Visualize QC statistics..." | tee -a $timinglogfile
{ /usr/bin/time \
	Rscript ${FASTQuick_BIN_DIR}/RPlotScript.R \
	${outputPrefix} \
	${indexPrefix}.FASTQuick.fa.bed.phase3.vcf.gz \
	${FASTQuick_SRC_DIR} \
1>&2 2>> $logfile; } 1>&2 2>>$timinglogfile
  outputDir=$(dirname ${outputPrefix})
{ /usr/bin/time \
	Rscript -e "rmarkdown::render('${FASTQuick_BIN_DIR}/FinalReport.rmd', params=list(input = '${outputPrefix}', SVDPrefix = '${indexPrefix}.FASTQuick.fa.bed.phase3.vcf.gz', FASTQuickInstallDir = '${FASTQuick_SRC_DIR}'), output_dir = '$outputDir')"
	1>&2 2>> $logfile;} 1>&2 2>>$timinglogfile

fi
echo "$(date)	Run complete with $(grep WARNING $logfile | wc -l) warnings and $(grep ERROR $logfile | wc -l) errors."