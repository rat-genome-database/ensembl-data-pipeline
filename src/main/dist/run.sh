#!/usr/bin/env bash
. /etc/profile
APPNAME="ensembl-data-pipeline"
APPDIR=/home/rgddata/pipelines/$APPNAME
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

EMAILLIST=mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
  EMAILLIST=mtutaj@mcw.edu,jrsmith@mcw.edu,jdepons@mcw.edu
fi

cd $APPDIR
java -jar -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/$APPNAME.jar "$@" 2>&1 > run.log

mailx -s "[$SERVER] EnsemblDataPipeline OK!" $EMAILLIST < $APPDIR/logs/summary.log
mailx -s "[$SERVER] EnsemblDataPipeline conflicts" $EMAILLIST < $APPDIR/logs/conflicts.log


