APPDIR=/home/rgddata/pipelines/EnsemblDataPipeline
# set variable HOST to uppercase hostname (f.e. KIRWAN, REED)
HOST=`hostname -s | tr '[a-z]' '[A-Z]'`

EMAIL_LIST=sellanki@mcw.edu
if [ "$HOST" == "REED" ]; then
  EMAIL_LIST="sellanki@mcw.edu"
fi

$APPDIR/run.sh -species Rat
$APPDIR/run.sh -species Mouse
$APPDIR/run.sh -species Human

mailx -s "[$HOST] Ensembl pipeline" $EMAIL_LIST < $APPDIR/logs/status.log
