#!/usr/bin/env bash

APPDIR=/home/rgddata/pipelines/ensembl-data-pipeline

# load PolyPhen-2 predictions for human variants (GRCh37 and GRCh38) from Ensembl variation database dumps;
# the model (humdiv/humvar) and the source files are configured in AppConfigure.xml, bean 'polyphenLoader'
$APPDIR/run.sh -polyphen
