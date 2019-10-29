#!/usr/bin/env bash
. /etc/profile

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
APPNAME=EnsemblDataPipeline
APPDIR=/home/rgddata/pipelines/$APPNAME
EMAILLIST=selanki@mcw.edu,mtutaj@mcw.edu

if [ "$SERVER" == "REED" ]; then
  EMAILLIST=selanki@mcw.edu
fi

cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/$APPNAME.jar "$@" 2>&1

mailx -s "[$SERVER] Ensembl Data Pipeline" $EMAILLIST < $APPDIR/logs/status.log
