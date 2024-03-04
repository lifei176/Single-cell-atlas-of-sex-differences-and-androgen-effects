#Code script for cellranger pipeline which is executed on Linux System. 
#Pipeline can be built accodrding to the tutorials provided by 10xgenomics at https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials.
#Example is taken on the basis of adrenal sample 1 under FD condition  (FD_1_adrenal) which can be accessed and downloaded from GSA under CRA006610. 
#/../bin/cellranger: cellranger installation location.
#/../mm10: reference file location.
#/../raw/adrenal: raw fastq files location.
---------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------
################# 1. Download mm10 genome file and genome annotation file
#################
#################
#################
#################
#################
wget ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gunzip Mus_musculus.GRCm38.93.gtf.gz

################# 2. Build mm10 reference for cellranger downstream analysis
#################
#################
#################
#################
#################
/../bin/cellranger mkgtf Mus_musculus.GRCm38.93.gtf Mus_musculus.GRCm38.93.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

/../bin/cellranger mkref --genome=mm10 \
                 --fasta=Mus_musculus.GRCm38.dna.primary_assembly.fa \
                 --genes=Mus_musculus.GRCm38.93.filtered.gtf \
                 --ref-version=3.0.0
                 
################# 3. Aligns sequencing reads in FASTQ files to generate expression matrix files for further analysis.
#################
#################
#################
#################
#################
/../bin/cellranger count --id=FD_1_adrenal \
                   --transcriptome=/../mm10 \
                   --fastqs=/../raw/adrenal \
                   --sample=FD-1-adrenal \
                   --expect-cells=10000 \
                   --localcores=8 \
                   --localmem=64

