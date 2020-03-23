[![Build Status](https://travis-ci.org/Griffan/FASTQuick.png?branch=master)](https://travis-ci.org/Griffan/FASTQuick)
[![GitHub Downloads](https://img.shields.io/github/downloads/Griffan/FASTQuick/total.svg?style=flat)](https://github.com/Griffan/FASTQuick/releases)

### OVERVIEW
   FASTQuick is an **ultra-fast** QC tool for NGS sequencing fastq files. It generates a comprehensive list of QC statistics, including **ancestry** and **contamination estimation**, at ~50x faster turnaround time than alignment-based QC tools.
   
### CONTENTS

- [GETTING STARTED](#getting-started)
- [INSTALL](#install)
- [WIKI PAGE](#wiki-page)
- [FAQ](#faq)
- [AUTHOR](#author)
- [COPYRIGHT](#copyright)

### GETTING STARTED

Follow the procedures below to quickly get started using FASTQuick.

#### Clone and Install FASTQuick

First, to start using FASTQuick, clone and install the repository.

```
git clone https://github.com/Griffan/FASTQuick.git
cd FASTQuick
mkdir build
cd build
cmake ..
make   
make test
```

Please refer to [INSTALL](#install) for more comprehensive guide on how to download and install FASTQuick.

#### Perform a Test Run

To perform a test run to make sure that FASTQuick runs as expected with a very small-sized example (assuming that you are still inside `build` directory, run

```
cd ../example/
bash example.sh
```

for more example scripts to test whether the software tool works as expected or not, see [EXAMPLES](#examples).

**Note** that you need `samtools` installed in your system and included in the `${PATH}` directory to run the test successfully

#### Download Resource Files

To run FASTQuick on real human sequence data, you need to download resource files using the following commands. (Before downloading, you may want to change your current directory.)

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/dbsnp132_20101103.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz
```

Please refer to [RESOURCE FILES](#resource-files) for more details. Note that you do not need to run FASTQuick on GRCh38 reference because the alignment of FASTQuick will be used only for internal purpose and will provide very similar results regardless of which reference will be used.


#### Run FASTQuick for Your Own FASTQ Files

For simplicity, we prepared an all-in-one script to process the whole FASTQuick pipeline or choose any start point of the pipeline (All | AllButIndex | Index | Align | Contamination | Ancestry | Visualize) in one command line.

```
${FASTQUICK_HOME}/bin/FASTQuick.sh 
--steps All \
--reference /path/to/hs37d5.fa \
--dbSNP /path/to/dbsnp132_20101103.vcf.gz \
--callableRegion /path/to/20141020.strict_mask.whole_genome.bed \
--output <output.prefix> \
--fastq1 <input.R1.fastq.gz> \
--fastq2 <input.R2.fastq.gz> \
--candidateVCF /path/to/hapmap_3.3.b37.sites.vcf.gz \
[--SVDPrefix ${FASTQUICK_HOME}/resource/1000g.phase3.10k.b37.vcf.gz] \
[--targetRegion <targetRegion.bed>] 
```

Please replace `/path/to/` the directory that contains the downloaded reference files (or use `.` if everything happened in the same directory). You will need to specify the input and output file names denoted as `<...>`.

**Note** that you only need to build indices once, hence `--steps AllButIndex` should be the preferred option once indices are ready.

#### Resource Files

These resource files can be shared and reused by different samples:

**reference genome**(**--reference**) [hs37d5.fa](http://tinyurl.com/jvflzg3)

**dbSNP VCF**(**--dbSNP**) [dbsnp132_20101103.vcf.gz](http://tinyurl.com/sl2kgof)

**1000 strict masked region**(**--callableRegion**) [20141020.strict_mask.whole_genome.bed](http://tinyurl.com/sjhb5nn)

**candidate variant list**(**--candidateVCF**) [hapmap_3.3.b37.vcf.gz](https://tinyurl.com/u69z6ts) (It's recommended to "shuffle" candidate variant list with [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html) before using it.)

For a quick start, it's recommended that you can feed **--candidateVCF** with a set of predefined markers, and **--SVDPrefix** with ready SVD files(to avoid SVD calculation on the fly). For example: 
```
--candidateVCF ${FASTQUICK_HOME}/resource/1000g.phase3.10k.b37.vcf.gz
--SVDPrefix ${FASTQUICK_HOME}/resource/1000g.phase3.10k.b37.vcf.gz

```

**Optionally**, you can enable **_target region_** mode by specifying **--targetRegion** with a bed format file.

**Note** that if other build version of reference genomes are needed, all these resource files are required to be the same build version with the reference genome.

#### Input Files

**--fastq_1** and **--fastq_2** expect pair-end fastq files(omit --fastq_2 for single-end dataset).

You can download fastq files of HG00553 from 1000 genome to reproduce the low-coverage [FinalReport.html](https://www.dropbox.com/s/7fbtpq82zduk4la/FinalReport.html?dl=1) in our example:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR013/ERR013170/ERR013170_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR013/ERR013170/ERR013170_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR015/ERR015764/ERR015764_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR015/ERR015764/ERR015764_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR018/ERR018525/ERR018525_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR018/ERR018525/ERR018525_2.fastq.gz
```
You can also use **--fastqList** to provide fastq files in format described in [FAQ](#faq)

#### Output Files

Once the process finished, you'll find summary statistics in various files starting with the same prefix(provided by **--output**). 

You also will find a similar [FinalReport.html](https://www.dropbox.com/s/7fbtpq82zduk4la/FinalReport.html?dl=1) in your output directory(base directory of prefix provided by **--output**). 

### INSTALL

To install FASTQuick, run the following series of commands.

```
git clone https://github.com/Griffan/FASTQuick.git   
mkdir build
cd build
cmake ..
make   
make test
```

Installation is complete if all tests finish successfully.

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

```

For libhts:
  - cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/  -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.a ..

For bzip2:
  - cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
```

**Note** that if you use docker to deploy, the minimal memory requirement is 4GB.

### EXAMPLES

You can find example scripts for each single step in example directory as template for customized usage.

For example:
 * the script **example.sh** is the simplified template for one-stop analysis.(our bin/FASTQuick.sh is more comprehensive and hence recommended)
 * the script **example.index.sh** is the template for selection new marker set and indexing reference data structures.
 * the script **example.align.sh** is the template for primary analysis.
 * the script **example.pop+con.sh** is the template to estimate contamination level and genetic ancestry of the intended sample.
 * the script **example.predefine.marker.index.sh** is the template to use pre-defined marker set to build indices.

### WIKI PAGE

We encourage users to refer to FASTQuick wiki page for more detailed description. [https://github.com/Griffan/FASTQuick/wiki]
   
### FAQ

1. **--fastqList** expects tab-delimited format as follows:

```
read.group.A.read_1.fq.gz   read.group.A.read_2.fq.gz
read.group.A.single.end.fq.gz
read.group.B.read_1.fq.gz   read.group.B.read_2.fq.gz
read.group.C.read_1.fq.gz   read.group.C.raed_2.fq.gz
read.group.C.single.end.fq.gz
```
2. A full list of required libraries an packages that are required to run the pipeline:
```
binary libraries:
zlib
libbzip2
libcurl
libssl

R libraries:
ggplot2
scales
knitr
rmarkdown
```

### AUTHOR
Fan Zhang (email:fanzhang@umich.edu)

### COPYRIGHT
   The full FASTQuick package is distributed under MIT License.


