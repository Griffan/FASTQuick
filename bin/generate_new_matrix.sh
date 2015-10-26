#$1 is your own *.SelectedSite.vcf.gz 
#$2 is resource directory
#$3 is the src/bin directory which contains auxilary tools
#$4 is exact path of hapmap_3.3.b37.vcf.gz, e.g. GATK-resources/hapmap_3.3.b37.vcf.gz
#$5 is exact path of 1000g/release/20110521/ e.g. /net/1000g/1000g/release/20110521/

PRG=$0
usage()
{
    echo "usage:$PRG [your own *.SelectedSite.vcf.gz] [src/resource directory] [ src/bin directory] [path of hapmap_3.3.b37.vcf.gz] [path of 1000g/release/20110521/]...." >&2
    exit 1
}

[ "$#" -ne 5 ] && usage

# parse commandline
#while [ $# -gt 0 ]
#do
#      arg="$1"
#      case "$arg" in
#     -d|--debug) echo "debug" ;;
#     -h|--help) usage ;;
     #--) shift; break;;  # no more options
#     -*) usage ;; 
#      *) break;; # not option, its some argument
#    esac
#    shift
#done



##prepare hapmap_3.3.b37.dat file and bed file
echo "selected vcf file: $1
resource directory: $2
auxilary tools directory $3"
gzip -dc $1 |grep -P -v "^#|^Y|^X"|awk '{print $1"\t"$2-1"\t"$2"\t"$4"\t"$5}' >$2/choose.bed
$3/vcfast convert --bedf $2/choose.bed --vcf $4 --out $2/hapmap_3.3.b37.dat --minMAF 0.01
tail -n +2 $2/hapmap_3.3.b37.dat |cut -f1|cut -f1 -d "_"|awk -F ":" '{print $1"\t"$2-1"\t"$2}' >$2/choose.bed.post.bed
tail -n +2 $2/hapmap_3.3.b37.dat |cut -f1|cut -d "_" -f2-3|cut -f1 -d "_"|awk -F "/" '{print $1"\t"$2}' >$2/choose.bed.allele
paste $2/choose.bed.post.bed $2/choose.bed.allele >$2/choose.bed.post.bed.allele.bed
rm $2/choose.bed.post.bed $2/choose.bed.allele
echo "Finish updating choose.bed.post.bed.allele.bed\n"
##prepare 1kg.phase1.selected.GLs.dat file
echo "perl $3/a03-extract-GLs-from-phase1.pl $2 $5 "|sh
echo "Finish updating 1kg.phase1.selected.GLs.dat file\n"
##generate new UD V mu matrix and 1kg PCs
Rscript $3/run.SVD.r $2/hapmap_3.3.b37.dat $2
echo "Finish updating SVD matrices\n"
echo "Update all the files successfully! Now you can proceed using FastPopCon pop function to predict population identification now."
