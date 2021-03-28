wc -l */*txt | grep -v -e scripts -e total | awk '{split($2,a,"/"); split(a[2],b,"QC_"); split(b[2],c,"."); print a[1],c[1],c[3],$1}' > accessory/chromfilter.stats
wc -l */*dupvar | awk '{split($2,a,"/"); split(a[2],b,"QC_"); split(b[1],c,"."); print a[1],c[1],c[2],$1}' | grep -v total >> accessory/chromfilter.stats
grep variants */*raw*log | grep QC | awk '{split($1,b,"/"); n=split(b[2],c,"_"); print b[1],c[n]}' | tr ":" " " | sed 's/raw/chr/' | sed 's/.log/ raw/' >> accessory/chromfilter.stats
wc -l */*-QC-Ref_chr*.bim | awk '{split($2,a,"/"); n=split(a[2],b,"_"); printf("%s\t%s\tpostQC\t%s\n",a[1],b[n],$1)}' | sed 's/.bim//' | grep -v total >> accessory/chromfilter.stats
grep people *QC_chr*log | grep -v males | awk '{n=split($1,a,"_"); print a[n]}' | sed 's/.log:/\tMAF\t/' >> accessory/merge.stats
wc -l *mergeSNPs.txt | grep -v total | awk '{n=split($2,b,"_"); printf("%s\tovlp\t%s\n",b[n-1],$1)}' >> accessory/merge.stats
