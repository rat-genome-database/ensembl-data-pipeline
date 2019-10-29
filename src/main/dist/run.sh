#!/usr/bin/env bash
. /etc/profile

SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
APPNAME=EnsemblDataPipeline
APPDIR=/home/rgddata/pipelines/$APPNAME
EMAILLIST=selanki@mcw.edu

if [ "$SERVER" == "REED" ]; then
  EMAILLIST=selanki@mcw.edu
fi

cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/$APPNAME.jar "$@" 2 > status.log

mailx -s "[$SERVER] EnsemblDataPipeline OK!" $EMAILLIST < status.log
