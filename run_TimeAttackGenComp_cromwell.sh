#!/bin/bash
#SLURM --mem=4G
#SLURM --time=20:00:00

java -Xmx4G \
    -Dconfig.file=./application.slurm.singularity.conf \
    -DLOG_MODE=standard \
    -jar cromwell.jar \
    run \
    TimeAttackGenComp.wdl \
    --inputs TimeAttackGenComp_inputs.json
