PID=$(basename $PWD)

PIPELINE=rrnaseq_star_featurecounts_pe_trim3prime
mkdir -p pipelines/$PIPELINE
rsync -ztatrvhPe ssh fcastro@euler.ethz.ch:~/bioinformatics/Analysis/${PID}/rawdata .
rsync -ztatrvhPe ssh fcastro@euler.ethz.ch:~/bioinformatics/Analysis/${PID}/pipelines/$PIPELINE/run.sh pipelines/$PIPELINE/.
