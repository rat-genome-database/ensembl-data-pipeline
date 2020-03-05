#!/usr/bin/env bash
. /etc/profile
APPNAME=EnsemblDataPipeline
APPDIR=/home/rgddata/pipelines/$APPNAME
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

EMAILLIST=hsnalabolu@mcw.edu,mtutaj@mcw.edu
if [ "$SERVER" == "REED" ]; then
  EMAILLIST=hsnalabolu@mcw.edu,mtutaj@mcw.edu,jrsmith@mcw.edu,jdepons@mcw.edu
fi

cd $APPDIR
java -jar -Dspring.config=$APPDIR/../properties/default_db.xml \
    -Dlog4j.configuration=file://$APPDIR/properties/log4j.properties \
    -jar lib/$APPNAME.jar "$@" run.sh 2>&1

mailx -s "[$SERVER] EnsemblDataPipeline OK!" $EMAILLIST < $APPDIR/logs/status.log

