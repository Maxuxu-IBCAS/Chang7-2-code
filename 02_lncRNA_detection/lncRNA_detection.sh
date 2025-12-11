#lncRNA detection
python CPC2.py -i all.collapsed.rep.availabel.all.candidate -o cpc2_predicted_184091
awk '$3<100' cpc2_predicted_184091.non_coding >cpc2_predicted_184091.non_coding.100aa

#genomic feature of lncRNA
FEELnc_classifier.pl -i cpc2_predicted_184091.non_coding.100aa.gtf -a ZM11D.genome.gff --maxwindow=100000 -l lncRNA_classes.C72_T2T.log > lncRNA_classes.C72_T2T.txt

