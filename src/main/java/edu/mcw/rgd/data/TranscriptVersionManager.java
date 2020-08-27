package edu.mcw.rgd.data;

import edu.mcw.rgd.process.CounterPool;
import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Map;

/**
 * manages transcript versions in STABLE_TRANSCRIPTS table
 */
public class TranscriptVersionManager {

    // singleton
    private static TranscriptVersionManager _instance = new TranscriptVersionManager();

    public static TranscriptVersionManager getInstance() {
        return _instance;
    }

    private TranscriptVersionManager() {
    }

    // hashmap of tr accession to Info
    private Map<String, Info> map = new HashMap<String, Info>();

    Logger log = Logger.getLogger("transcriptVersions");

    public void init() {
        map.clear();
    }

    public void addVersion(String acc, String version) {

        Info info = map.get(acc);
        if( info==null ) {
            info = new Info();
            info.version = version;
            map.put(acc, info);
        } else {
            if( info.version==null ) {
                info.version = version;
            }
            else if( !info.version.equals(version) ) {
                throw new RuntimeException("transcript acc ver mismatch: "+info.version+" vs "+version);
            }
        }
    }

    public void addRgdId(String acc, int rgdId) {

        Info info = map.get(acc);
        if( info==null ) {
            info = new Info();
            info.rgdId = rgdId;
            map.put(acc, info);
        } else {
            if( info.rgdId==0 ) {
                info.rgdId = rgdId;
            }
            else if( info.rgdId != rgdId ) {
                log.warn("WARNING: transcript rgd id mismatch: "+acc+" for "+info.rgdId+" vs "+rgdId);
                info.conflict = true;
            }
        }
    }

    public void qcAndLoad(CounterPool counters) throws Exception {

        EnsemblDAO dao = new EnsemblDAO();

        for( Map.Entry<String, Info> entry: map.entrySet() ) {
            String acc = entry.getKey();
            Info info = entry.getValue();
            counters.increment("TRANSCRIPT_VERSIONS__PROCESSED");

            if( info.conflict ) {
                counters.increment("TRANSCRIPT_VERSIONS_WITH_CONFLICT");
                continue;
            }

            String accVer = info.version;

            String accVerInDb = dao.getTranscriptVersionInfo(acc);
            if( accVerInDb==null ) {
                log.debug("INSERTED "+acc+" RGD:"+info.rgdId+" "+accVer);
                dao.insertTranscriptVersionInfo(acc, accVer, info.rgdId);
                counters.increment("TRANSCRIPT_VERSIONS_INSERTED");
            } else {
                if( accVerInDb.equals(accVer) ) {
                    counters.increment("TRANSCRIPT_VERSIONS_UP_TO_DATE");
                } else {
                    if( accVer.length() > accVerInDb.length() || accVerInDb.compareTo(accVer)<0 ) {
                        log.debug("UPDATED "+acc+" RGD:"+info.rgdId+" "+accVerInDb+" ==> "+accVer);
                        dao.updateTranscriptVersionInfo(acc, accVer);
                        counters.increment("TRANSCRIPT_VERSIONS_MODIFIED");
                    }
                }
            }
        }
    }

    class Info {
        public String version;
        public int rgdId;
        public boolean conflict = false;

        public String toString() {
            return "VER="+version+", RGD:"+rgdId;
        }
    }
}
