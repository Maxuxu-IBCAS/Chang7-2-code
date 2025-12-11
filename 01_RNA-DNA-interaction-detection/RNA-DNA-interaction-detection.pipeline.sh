/mnt2/maxuxu/software/FLASH-1.2.11-Linux-x86_64/flash --threads 20 -M 145 -O --output-prefix ALL-PET-maize-30min.flash --output-directory ./ ../ALL-PET-maize-30min_1.fq.gz ../ALL-PET-maize-30min_2.fq.gz

cutadapt -b file:DNA_linker -n 14 --no-indels -o ALL-PET-maize-30min.combine.noLinker --info-file ALL-PET-maize-30min.combine.Linker_info --discard -O 18 ALL-PET-maize-30min.flash.extendedFrags.fastq >ALL-PET-maize-30min.combine.stat
cutadapt -b file:DNA_linker -n 8 --no-indels -o ALL-PET-maize-30min.notCombined.1.noLinker --info-file ALL-PET-maize-30min.notCombined.1.Linker_info --discard -O 18 ALL-PET-maize-30min.flash.notCombined_1.fastq >ALL-PET-maize-30min.flash.notCombined_1.stat
cutadapt -b file:DNA_linker -n 8 --no-indels -o ALL-PET-maize-30min.notCombined.2.noLinker --info-file ALL-PET-maize-30min.notCombined.2.Linker_info --discard -O 18 ALL-PET-maize-30min.flash.notCombined_2.fastq >ALL-PET-maize-30min.flash.notCombined_2.stat

perl ../script/ts_cutadapt2oneline.pl ALL-PET-maize-30min.combine.Linker_info ALL-PET-maize-30min.combine.cut.info >ALL-PET-maize-30min.combine.cut.info.log
perl ../script/ts_cutadapt2oneline.pl ALL-PET-maize-30min.notCombined.1.Linker_info ALL-PET-maize-30min.notCombined.1.cut.info > ALL-PET-maize-30min.notCombined.1.cut.info.log
perl ../script/ts_cutadapt2oneline.pl ALL-PET-maize-30min.notCombined.2.Linker_info ALL-PET-maize-30min.notCombined.2.cut.info > ALL-PET-maize-30min.notCombined.2.cut.info.log

perl ../script/split.combine.pl ALL-PET-maize-30min.combine.cut.info ALL-PET-maize-30min.combine 1 >ALL-PET-maize-30min.combine.fq.stat 2>ALL-PET-maize-30min.combine.fq.stat.log
perl ../script/split.only.pl ALL-PET-maize-30min.notCombined.1.cut.info ALL-PET-maize-30min.notCombined.2.cut.info ALL-PET-maize-30min.notCombined.only 1> ALL-PET-maize-30min.notCombined.only.fq.stat 2>ALL-PET-maize-30min.notCombined.only.fq.stat.log

cat /mnt/maxuxu/data/minor/36_Chang72_RNA_interaction/01_RNA-DNA/ALL-PET-maize-30min.combine.combine.RNA.fq /mnt/maxuxu/data/minor/36_Chang72_RNA_interaction/01_RNA-DNA/ALL-PET-maize-30min.notCombined.only.RNA.fq ALL-PET-maize-30min.combine.combine.RNA.fq ALL-PET-maize-30min.notCombined.only.RNA.fq >ALL-PET-maize-30min.RNA.fastq
cat /mnt/maxuxu/data/minor/36_Chang72_RNA_interaction/01_RNA-DNA/ALL-PET-maize-30min.combine.combine.DNA.fq /mnt/maxuxu/data/minor/36_Chang72_RNA_interaction/01_RNA-DNA/ALL-PET-maize-30min.notCombined.only.DNA.fq ALL-PET-maize-30min.combine.combine.DNA.fq ALL-PET-maize-30min.notCombined.only.DNA.fq >ALL-PET-maize-30min.DNA.fastq

#map
perl ../script/split.pl ALL-PET-maize-30min.DNA.fastq ALL-PET-maize-30min.DNA 70
bwa aln -t 10 -f ALL-PET-maize-30min.DNA.DNA.short.sai /mnt/maxuxu/data/genome/maize_chang7-2_T2T/ZM1-1D.fasta ALL-PET-maize-30min.DNA.DNA.short.fastq
bwa samse -f ALL-PET-maize-30min.DNA.DNA.short.sam /mnt/maxuxu/data/genome/maize_chang7-2_T2T/ZM1-1D.fasta ALL-PET-maize-30min.DNA.DNA.short.sai ALL-PET-maize-30min.DNA.DNA.short.fastq
bwa mem -t 10 /mnt/maxuxu/data/genome/maize_chang7-2_T2T/ZM1-1D.fasta ALL-PET-maize-30min.DNA.DNA.long.fastq >ALL-PET-maize-30min.DNA.DNA.long.sam

hisat2 -p 10 -x /mnt/maxuxu/data/genome/maize_chang7-2_T2T/ZM1-1D.fasta.hisat2 --rna-strandness F -U ALL-PET-maize-30min.RNA.fastq -S ALL-PET-maize-30min.RNA.sam && python ../script/flt_bam.py ALL-PET-maize-30min.RNA.sam ALL-PET-maize-30min.RNA && samtools sort -@ 10 -o ALL-PET-maize-30min.RNA.hisat.bam ALL-PET-maize-30min.RNA.uniqmap.bam
samtools view -Sb -F 256 ALL-PET-maize-30min.RNA.multimap.bam | bamToFastq -i - -fq ALL-PET-maize-30min.RNA.tmp.multi.fq
bamToFastq -i ALL-PET-maize-30min.RNA.unmap.bam -fq ALL-PET-maize-30min.RNA.tmp.unmap.fq
cat ALL-PET-maize-30min.RNA.tmp.multi.fq  ALL-PET-maize-30min.RNA.tmp.unmap.fq >ALL-PET-maize-30min.RNA.remap.fq
bowtie2 -p 10 --local --very-sensitive-local -x /mnt/maxuxu/data/genome/maize_chang7-2_T2T/ZM1-1D.fasta.bt2 -U ALL-PET-maize-30min.RNA.remap.fq -S ALL-PET-maize-30min.RNA.remap.sam && samtools view -Sb -q 2 ALL-PET-maize-30min.RNA.remap.sam | samtools sort -@ 10 -o ALL-PET-maize-30min.RNA.remap.bam -
samtools merge -f ALL-PET-maize-30min.RNA.bam ALL-PET-maize-30min.RNA.hisat.bam ALL-PET-maize-30min.RNA.remap.bam
samtools index ALL-PET-maize-30min.RNA.bam


java -cp ../script/LGL.jar LGL.util.UniqueSam ALL-PET-maize-30min.DNA.DNA.long.sam ALL-PET-maize-30min.DNA.DNA.long.rmdup.sam
python ../script/combine_DNA_file.py ALL-PET-maize-30min.DNA.DNA.short.sam ALL-PET-maize-30min.DNA.DNA.long.rmdup.sam ALL-PET-maize-30min.DNA.sam
samtools view -Sb ALL-PET-maize-30min.DNA.sam | samtools sort -@ 10 -o ALL-PET-maize-30min.DNA.bam -
samtools index ALL-PET-maize-30min.DNA.bam

htseq-count -o ALL-PET-maize-30min.RNA.yes.gene.sam -f bam -r name -s no -t gene -i gene_id -m intersection-nonempty ALL-PET-maize-30min.RNA.bam /mnt/maxuxu/data/genome/maize_chang7-2_T2T/Chang72.5Kb.window.gff >ALL-PET-maize-30min.RNA.yes.gene.count

samtools view -H ALL-PET-maize-30min.RNA.bam | cat - ALL-PET-maize-30min.RNA.yes.gene.sam >ALL-PET-maize-30min.RNA.yes.gene.cg.sam

htseq-count -o ALL-PET-maize-30min.RNA.yes.exon.sam -f sam -r name -s no -t exon -i gene_id -m intersection-nonempty ALL-PET-maize-30min.RNA.yes.gene.cg.sam /mnt/maxuxu/data/genome/maize_chang7-2_T2T/Chang72.5Kb.window.gff >ALL-PET-maize-30min.RNA.yes.exon.count

samtools view -H ALL-PET-maize-30min.RNA.bam | cat - ALL-PET-maize-30min.RNA.yes.exon.sam | perl -lane 'if(/XF:Z/){$_=~s/XF:Z/GE:Z/;print $_}else{print $_}' > ALL-PET-maize-30min.RNA.yes.sam 

python ../script/sam2bed2.py ALL-PET-maize-30min.RNA.yes.sam 1> ALL-PET-maize-30min.RNA.yes.bed 2>ALL-PET-maize-30min.RNA.yes.bed.log

bamToBed -cigar -i ALL-PET-maize-30min.DNA.bam | awk -v OFS="\t" '{if($6=="+"){print $1,$2,$2+500,$4,$5,$6,$7}else{start=$3-500;if(start<0){start=0;}print $1,start,$3,$4,$5,$6,$7}}' | intersectBed -a - -b /mnt/maxuxu/data/genome/maize_chang7-2_T2T/Chang72.5Kb.window.anchor -wao | python ../script/uniq_DNA.py >ALL-PET-maize-30min.DNA.bed

python ../script/combineDNAandRNA_withanchor.py ALL-PET-maize-30min.DNA.bed ALL-PET-maize-30min.RNA.yes.bed ALL-PET-maize-30min.DNARNA.bedpe

awk -v qu=20 '$8>=qu' ALL-PET-maize-30min.DNARNA.bedpe | sort --parallel=10 -k14,14 -k13,13 -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n -k8,8n -k9,9 -k10,10 -k11,11 -k12,12 >ALL-PET-maize-30min.DNARNA.sorted.bedpe

python ../script/uniqread.py ALL-PET-maize-30min.DNARNA.sorted.bedpe >ALL-PET-maize-30min.DNARNA.uniq.bedpe

python ../script/PetClusterWithGivenAnchors.py /mnt/maxuxu/data/genome/maize_chang7-2_T2T/Chang72.5Kb.window.anchor /mnt/maxuxu/data/genome/maize_chang7-2_T2T/Chang72.5Kb.window.2.bed ALL-PET-maize-30min.DNARNA.uniq.bedpe | sort --parallel=20 -k1,1 -k4,4 -k2,2n -k5,5n -k3,3n -k6,6n > ALL-PET-maize-30min.DNARNA.givenanchor.cluster

awk '$9>1' ALL-PET-maize-30min.DNARNA.givenanchor.cluster >ALL-PET-maize-30min.DNARNA.givenanchor.gt1.cluster

tmp=`wc -l ALL-PET-maize-30min.DNARNA.uniq.bedpe | awk '{print $1}'`

Rscript ../script/hypergeometric5.r $tmp ALL-PET-maize-30min.DNARNA.givenanchor.gt1.cluster ALL-PET-maize-30min.DNARNA.givenanchor.gt1.cluster.withpvalue

awk '$13+0.0<=0.05' ALL-PET-maize-30min.DNARNA.givenanchor.gt1.cluster.withpvalue >ALL-PET-maize-30min.DNARNA.givenanchor.FDRfiltered.txt
awk '$9>2' ALL-PET-maize-30min.DNARNA.givenanchor.FDRfiltered.txt >ALL-PET-maize-30min.DNARNA.givenanchor.FDRfiltered.gt2.txt

