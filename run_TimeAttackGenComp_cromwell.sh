#!/bin/bash
#PBS -V
#PBS -l mem=4G,walltime=20:00:00

module load java/1.8.0_101

java -Xmx4G \
    -Dconfig.file=./application.conf.cic \
    -DLOG_MODE=standard \
    -jar cromwell.jar \
    run \
    TimeAttackGenComp.wdl \
    --inputs TimeAttackGenComp_inputs.json
