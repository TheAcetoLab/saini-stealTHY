#!/bin/bash

# Setup variables and directories
# #+name: setup_env

FASTQLIST=$PWD/../../rawdata/fastq.tsv
SAMPLE_ANNOT=$PWD/../../rawdata/samples.csv
DATASET=$PWD/../../rawdata/dataset.csv
PIPELINE=rnaseq_star_featurecounts_cas9

# Setup modules
# #+name: setup_modules

#FASTQSCREEN_ML=fastq-screen/0.11.2
BOWTIE2_ML=bowtie2/2.4.4
FASTQC_ML=fastqc/0.11.9
TRIMGALORE_ML=trimgalore/0.6.6
KRAKEN_ML=kraken2/2.1.2
STAR_VS=2.7.9a
STAR_ML=star/$STAR_VS
SAMTOOLS_ML=samtools/1.12
SUBREAD_ML=subread/2.0.3
R_ML=r/4.1.3
PYTHON_ML=python/3.10.4
MULTIQC_ML=multiqc/1.9
RSeQC=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Apps/RSeQC/4.0.0/scripts
FASTQSCREEN=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Apps/fastq-screen/0.15.2/fastq_screen

# Pipeline configuration
# #+name: setup_modules

FASTQSCREEN_SUBSET=100000
FASTQSCREEN_CONF=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Pipelines/configuration/fastq_screen/variousSpecies_rRNA_silva138.1_DNA_mm_hs.conf
ORGANISM=Mus_musculus
GENOMEBUILD=GRCm39_CAS9
READ_LENGTH=100
STAR_INDEX=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Genomes/$ORGANISM/$GENOMEBUILD/Sequence/Indexes/STAR/$STAR_VS
GENCODE_VS=release_M29
GTF=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Genomes/$ORGANISM/$GENOMEBUILD/Annotation/gencode/$GENCODE_VS/genes/genes.gtf
GTF_BED=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Genomes/$ORGANISM/$GENOMEBUILD/Annotation/gencode/$GENCODE_VS/genes/genes.bed
MULTIQC_CONF=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Pipelines/configuration/multiqc/multiqc_config.yaml
KRAKEN2_DB_PATH=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Databases/kraken2/minikraken2_v2_8GB_201904_UPDATE
PIPELINE_UTILS=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Pipelines/code/utils

# Prepare folder structure
# #+name: setup_dir_strc

mkdir -p fastq_screen
mkdir -p kraken2
mkdir -p fastqc_raw
mkdir -p trim_galore
mkdir -p fastqc_trimmed
mkdir -p star
mkdir -p rseqc
mkdir -p featureCounts
mkdir -p featureCounts_unique
mkdir -p rseqc
mkdir -p multiqc/fastqc_raw
mkdir -p multiqc/fastq_processed
mkdir -p multiqc/featureCounts_unique
mkdir -p summarizedExperiment
mkdir -p cas9

# Run main pipeline
# #+name: run_pipeline

unset multiqc_dep
unset multiqc_raw_dep
unset multiqc_featureCounts_unique_dep

for FQ in $(<$FASTQLIST); do
  SM=$(basename $FQ | sed -r 's/^[0-9]+.[A-Z]-(.*)_R..fastq.gz/\1/' | sed -r 's/^(.*)_S[0-9]+_R._[0-9]+.fastq.gz/\1/')
  BASEFQ=$(basename $FQ | sed 's/.fastq.gz//')

  if ! grep -q $SM $SAMPLE_ANNOT; then
    printf "Sample $SM not foundin sample annotation file ($SAMPLE_ANNOT)"
  else
    printf "\n\n####### Submit jobs for $SM\n"

    # FastQ Screen
    CWDIR=fastq_screen
    pushd $CWDIR
    NREADS=100000

    NCPU=8
    let MEM_MB=1024*5

    echo "
    module purge
    module load gcc/8.2.0 $BOWTIE2_ML
    $FASTQSCREEN --threads $NCPU --aligner bowtie2 --subset $FASTQSCREEN_SUBSET -conf $FASTQSCREEN_CONF --threads $NCPU --bowtie2 '-k 10 --very-sensitive --trim5 4 --trim3 4' $FQ
  " > ${SM}.sh
    JOBID_fastq_screen=$(bsub -W 24:00 -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    multiqc_dep+=" && done($JOBID_fastq_screen)"
    popd

    # Kraken2
    CWDIR=kraken2
    pushd $CWDIR

    NCPU=8
    let MEM_MB=1024*2

    echo "
    module load gcc/8.2.0 $KRAKEN_ML
    kraken2 --db $KRAKEN2_DB_PATH \
      --gzip-compressed \
      --classified-out ${SM}.classified-out.txt \
      --unclassified-out ${SM}.unclassified-out.txt \
      --unclassified-out ${SM}.unclassified-out.txt \
      --output ${SM}.kraken2 \
      --report ${SM}.k2report \
      --report-minimizer-data \
      --minimum-hit-groups 3 \
      --threads $NCPU \
      $FQ
    " > ${SM}.sh
    JOBID_kraken2=$(bsub -W 24:00 -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    multiqc_dep+=" && done($JOBID_kraken2)"
    popd

    # FastQC rawdata
    CWDIR=fastqc_raw
    pushd $CWDIR

    NCPU=2
    let MEM_MB=1024*2

    echo "
    module purge
    module load gcc/8.2.0 $FASTQC_ML
    fastqc -q -o ./  -t $NCPU $FQ
    " > ${SM}.sh
    JOBID_fastqc_raw=$(bsub -W 24:00 -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    multiqc_raw_dep+=" && done($JOBID_fastqc_raw)"
    popd


    # Trim_galore
    CWDIR=trim_galore
    QCDIR=$PWD/fastqc_trimmed
    pushd $CWDIR

    NCPU=1
    let MEM_MB=1024*5

    echo "
    module purge
    module load gcc/8.2.0 $TRIMGALORE_ML
    trim_galore -q 20 --fastqc --length 20 --fastqc_args \"--outdir $QCDIR\" $FQ
    mv ${BASEFQ}_trimmed.fq.gz ${BASEFQ}.fastq.gz
    " > ${SM}.sh
    JOBID_trim_galore=$(bsub -W 24:00 -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    popd

    # STAR-2pass alignment
    CWDIR=star
    FQTRIMMED=$PWD/trim_galore/${BASEFQ}.fastq.gz
    let OVERHANG=${READ_LENGTH}-1

    pushd $CWDIR

    NCPU=8
    let MEM_MB=(1024*6)

    echo "
    module purge
    module load gcc/8.2.0 $STAR_ML $SAMTOOLS_ML
    STAR \
      --genomeDir $STAR_INDEX \
      --sjdbGTFfile $GTF \
      --runThreadN $NCPU \
      --readFilesIn $FQTRIMMED \
      --readFilesCommand gunzip -c \
      --sjdbOverhang $OVERHANG \
      --outFilterType BySJout \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMultimapNmax 50 \
      --outFilterMatchNmin 30 \
      --outFilterMismatchNmax 10 \
      --outFilterMismatchNoverLmax 0.05 \
      --alignIntronMax 100000 \
      --alignMatesGapMax 100000 \
      --outMultimapperOrder Random \
      --outSAMmapqUnique 60 \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 30000000000 \
      --twopassMode Basic \
      --outTmpDir ${SM}.2pass \
      --outSAMattrRGline ID:$SM SM:$SM PL:ILLUMINA \
      --outFileNamePrefix ${SM}_
    mv ${SM}_Aligned.sortedByCoord.out.bam ${SM}.bam
    samtools index ${SM}.bam
    " > ${SM}.sh
    JOBID_star=$(bsub -W 24:00 -w "$JOBID_trim_galore" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    popd

    # CAS9 extraction
    CWDIR=cas9
    BAM=$PWD/star/${SM}.bam

    pushd $CWDIR

    NCPU=1
    let MEM_MB=1024*10

    echo "
    module purge
    module load gcc/8.2.0 $SAMTOOLS_ML
    samtools view -b $BAM "CAS9" > ${SM}.bam
    samtools index ${SM}.bam
    " > ${SM}.sh
    JOBID_cas9=$(bsub -W 01:00 -w "$JOBID_cas9" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    popd

    # RSeQC
    CWDIR=rseqc
    BAM=$PWD/star/${SM}.bam

    pushd $CWDIR

    NCPU=1
    let MEM_MB=1024*10

    echo "
    module purge
    module load gcc/8.2.0 $PYTHON_ML $SAMTOOLS_ML $R_ML
    infer_experiment.py -r $GTF_BED -i $BAM > ${SM}.strandeness
    junction_saturation.py -i $BAM -o $SM -r $GTF_BED
    bam_stat.py -i $BAM -q 20 > ${SM}.bam_stat.out
    read_distribution.py -i $BAM -r $GTF_BED > ${SM}.read_distribution.out
    read_duplication.py -i $BAM -q 20 -o $SM
    read_GC.py -i $BAM -o $SM
    #geneBody_coverage.py -i $BAM -r $GTF_BED -o $SM
    " > ${SM}.sh
    JOBID_rseqc=$(bsub -W 24:00 -w "$JOBID_star" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    multiqc_dep+=" && done($JOBID_rseqc)"
    popd


    # featurecounts
    CWDIR=featureCounts
    BAM=$PWD/star/${SM}.bam
    INF_EXP=$PWD/rseqc/${SM}.strandeness

    pushd $CWDIR
    mkdir -p gene exon transcript

    NCPU=1
    let MEM_MB=1024*5

    echo "
    module purge
    module load gcc/8.2.0 $SUBREAD_ML
    source $PIPELINE_UTILS/*
    STD=\$(strandedness $INF_EXP)
    featureCounts -t exon -g gene_id -a $GTF -F GTF -M -O --minOverlap 10 -Q 10 -s \$STD -o gene/${SM}.counts.txt $BAM
    featureCounts -t exon -f -J -a $GTF -F GTF -M -O --minOverlap 10 -Q 10 -s \$STD -o exon/${SM}.counts.txt $BAM
    featureCounts -t exon -g transcript_id -J -a $GTF -F GTF -M -O --minOverlap 10 -Q 10 -s \$STD -o transcript/${SM}.counts.txt $BAM
    " > ${SM}.sh
    JOBID_featureCounts=$(bsub -W 24:00 -w "$JOBID_rseqc" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    multiqc_dep+=" && done($JOBID_featureCounts)"
    popd

    # featurecounts removing multi-mapping reads
    CWDIR=featureCounts_unique
    BAM=$PWD/star/${SM}.bam
    INF_EXP=$PWD/rseqc/${SM}.strandeness

    pushd $CWDIR
    mkdir -p gene exon transcript

    NCPU=1
    let MEM_MB=1024*5

    echo "
    module purge
    module load gcc/8.2.0 $SUBREAD_ML
    source $PIPELINE_UTILS/*
    STD=\$(strandedness $INF_EXP)
    featureCounts -t exon -g gene_id -a $GTF -F GTF -O --minOverlap 10 -Q 10 -s \$STD -o gene/${SM}.counts.txt $BAM
    featureCounts -t exon -f -J -a $GTF -F GTF --minOverlap 10 -Q 10 -s \$STD -o exon/${SM}.counts.txt $BAM
    featureCounts -t exon -g transcript_id -J -a $GTF -F GTF -O --minOverlap 10 -Q 10 -s \$STD -o transcript/${SM}.counts.txt $BAM
    " > ${SM}.sh
    JOBID_featureCounts_unique=$(bsub -W 24:00 -w "$JOBID_rseqc" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J ${SM}-${CWDIR} -o ${SM}.out -e ${SM}.err < ${SM}.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
    multiqc_featureCounts_unique_dep+=" && done($JOBID_featureCounts_unique)"
    popd

  fi

done
multiqc_dep=${multiqc_dep:4} # remove leading " && "
multiqc_raw_dep=${multiqc_raw_dep:4} # remove leading " && "
multiqc_featureCounts_unique_dep=${multiqc_featureCounts_unique_dep:4} # remove leading " && "

# Run MultiQC
# #+name: run_multiqc

NCPU=8
MEM_MB=8192
let MEM_MB=($NCPU*$MEM_MB)-2

# MultiQC raw fastqc
CWDIR=multiqc/fastqc_raw
PDIR=$PWD
pushd $CWDIR
NCPU=1
let MEM_MB=1024*5
echo "
module purge
module load gcc/6.3.0 $MULTIQC_ML
multiqc \
  -c ${MULTIQC_CONF} \
  -fo ./ \
  $PDIR/fastqc_raw
" > run.sh
JOBID_multiqc_fastqc_raw=$(bsub -W 05:00 -w "$multiqc_raw_dep" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J mq_multiqc_fastqc_raw -o run.out -e run.err < run.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
popd

# MultiQC processed
CWDIR=multiqc/fastq_processed
PDIR=$PWD
pushd $CWDIR
NCPU=1
let MEM_MB=1024*5
echo "
module purge
module load gcc/8.2.0 python/3.10.4
multiqc \
  -c ${MULTIQC_CONF} \
  -fo ./ \
  $PDIR/fastq_screen \
  $PDIR/kraken2 \
  $PDIR/trim_galore \
  $PDIR/fastqc_trimmed \
  $PDIR/star \
  $PDIR/rseqc \
  $PDIR/featureCounts/gene
" > run.sh
JOBID_multiqc_fastq_processed=$(bsub -W 05:00 -w "$multiqc_dep" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J mq_fastq_processed -o run.out -e run.err < run.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
popd

# MultiQC featureCounts unique
CWDIR=multiqc/featureCounts_unique
PDIR=$PWD
pushd $CWDIR
NCPU=1
let MEM_MB=1024*5
echo "
module purge
module load gcc/8.2.0 python/3.10.4

multiqc \
  -c ${MULTIQC_CONF} \
  -fo ./ \
  $PDIR/featureCounts_unique/gene
" > run.sh
JOBID_multiqc_featureCounts_unique=$(bsub -W 05:00 -w "$multiqc_featureCounts_unique_dep" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J mq_featureCounts_unique -o run.out -e run.err < run.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
popd

# Generate SCE objects
# #+name: run_generate_se

CWDIR=summarizedExperiment
PDIR=$PWD
GEN_SE=$PDIR/generate_se.R
pushd $CWDIR
NCPU=8
let MEM_MB=1024*2

echo "
module purge
module load gcc/8.2.0 $R_ML
Rscript $GEN_SE $PDIR/featureCounts/gene gene $SAMPLE_ANNOT $GTF $PDIR/multiqc/fastq_processed $NCPU ./se_gene.rds
" > generate_se_gene.sh
JOBID_generate_se_gene=$(bsub -W 05:00 -w "$JOBID_multiqc_fastq_processed" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J generate_se_gene -o generate_se_gene.out -e generate_se_gene.err < generate_se_gene.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

echo "
module purge
module load gcc/8.2.0 $R_ML
Rscript $GEN_SE $PDIR/featureCounts_unique/gene gene $SAMPLE_ANNOT $GTF $PDIR/multiqc/featureCounts_unique $NCPU ./se_gene_fc_unique.rds
" > generate_se_gene_unique.sh
JOBID_generate_se_gene_unique=$(bsub -W 05:00 -w "$JOBID_multiqc_featureCounts_unique" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J generate_se_gene_unique -o generate_se_gene_unique.out -e generate_se_gene_unique.err < generate_se_gene_unique.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

echo "
module purge
module load gcc/8.2.0 $R_ML
Rscript $GEN_SE $PDIR/featureCounts/transcript transcript $SAMPLE_ANNOT $GTF $PDIR/multiqc/fastq_processed $NCPU ./se_tx.rds
" > generate_se_tx.sh
JOBID_generate_se_tx=$(bsub -W 05:00 -w "$JOBID_multiqc_fastq_processed" -n $NCPU -R  "rusage[mem=$MEM_MB]" -J generate_se_tx -o generate_se_tx.out -e generate_se_tx.err < generate_se_tx.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')

popd

# Transfer to NAS
# Run only after succesful execution of the pipeline
# #+name: clean_data_and_transfer

PROJECT_ID=$(sed 1d $SAMPLE_ANNOT | head -1 | cut -d "," -f 1)
ORDER_ID=$(sed 1d $SAMPLE_ANNOT | head -1 | cut -d "," -f 2)
PO_ID=p${PROJECT_ID}_o${ORDER_ID}
PDIR=/nfs/nas22.ethz.ch/fs2202/biol_imhs_aceto/Bioinformatics/Analysis/${PO_ID}/pipelines/$PIPELINE
mkdir -p $PDIR
echo "
for i in run.sh generate_se.R multiqc summarizedExperiment featureCounts featureCounts_unique fastq_screen kraken2 cas9; do
  rsync -tr --exclude 'kraken2/*classified*txt' --exclude 'kraken2/.kraken2' \$i $PDIR/.
done
" > rsync.sh
#JOBID_rsync=$(bsub -W 05:00 -w "$rsync_dep" -J rsync -o rsync.out -e rsync.err < rsync.sh | awk '/is submitted/{print substr($2, 2, length($2)-2);}')
