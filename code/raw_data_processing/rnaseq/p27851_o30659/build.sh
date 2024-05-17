PID=$(basename $PWD)

PIPELINE=rnaseq_star_featurecounts_cas9
mkdir -p pipelines/$PIPELINE
rsync -ztatrvhPe ssh fcastro@euler.ethz.ch:~/bioinformatics/Analysis/${PID}/rawdata .
rsync -ztatrvhPe ssh fcastro@euler.ethz.ch:~/bioinformatics/Analysis/${PID}/pipelines/$PIPELINE/run.sh pipelines/$PIPELINE/.
    