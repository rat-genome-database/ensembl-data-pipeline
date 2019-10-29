#!/usr/bin/env bash
. /etc/profile

APPNAME=EnsemblDataPipeline
APPDIR=/home/rgddata/pipelines/$APPNAME
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`
EMAIL_LIST=selanki@mcw.edu

if [ "$SERVER" == "REED" ]; then
  EMAIL_LIST=selanki@mcw.edu
fi

cd $APPDIR

java -jar -Dspring.config=$APPDIR/../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/$APPNAME.jar "$@" 2>&1

mailx -s "[$SERVER] EnsemblDataPipeline OK!" $EMAIL_LIST < $APPDIR/logs/status.log
