cat $1|tail -n +2|cut -f1|cut -f1 -d "_"|awk -F ":" '{print $1"\t"$2-1"\t"$2}' >choose.bed.post.bed
cat $1|tail -n +2 |cut -f1|cut -d "_" -f2-3|cut -f1 -d "_"|awk -F "/" '{print $1"\t"$2}' >choose.bed.allele
paste choose.bed.post.bed choose.bed.allele >choose.bed.post.bed.allele.bed
