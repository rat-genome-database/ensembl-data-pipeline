#!/usr/bin/env bash

APPDIR=/home/rgddata/pipelines/ensembl-data-pipeline

# by default all species and assemblies configured in AppConfigure.xml are processed;
# optional pass-through args: -species <species> and/or -mapKey <k1[,k2,...]> restrict the run,
# f.e.:  run_ensembl.sh -species rat -mapKey 301,302,303
$APPDIR/run.sh -useGff3Loader "$@"
