package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.MapData;
import org.apache.logging.log4j.Logger;

import java.util.*;

public class AssemblyPositions {

    private Map<String, MapData> positionsInRgd;
    private Map<String, MapData> positionsIncoming;

    public int loadPositionsInRgd(int mapKey, EnsemblDAO dao) throws Exception {

        positionsInRgd = new HashMap<>();
        positionsIncoming = new HashMap<>();

        List<MapData> mds = dao.getGenePositions(mapKey);
        for( MapData md: mds ) {
            String key = computeKey(md);
            positionsInRgd.put(key, md);
        }
        return positionsInRgd.size();
    }

    /**
     *
     * @param md
     * @return true if this position was already among incoming positions (uncommon)
     */
    public boolean addIncomingPosition(MapData md) {
        String key = computeKey(md);
        return positionsIncoming.put(key, md) != null;
    }

    public void qcAndLoad(Logger log, EnsemblDAO dao) throws Exception {

        log.info("gene positions in RGD: "+positionsInRgd.size());
        log.info("gene positions incoming: "+positionsIncoming.size());

        // determine positions for insert
        Set<String> keysForInsert = new HashSet<>(positionsIncoming.keySet());
        keysForInsert.removeAll(positionsInRgd.keySet());

        for( String keyForInsert: keysForInsert ) {
            MapData posForInsert = positionsIncoming.get(keyForInsert);
            dao.insertMapData(posForInsert, "GENE ");
        }
        log.info("gene positions inserted: "+keysForInsert.size());


        // determine positions for delete
        //
        // be conservative: limit positionsInRgd to RGD IDS that are present
        //  (that means if for any reason a gene that exists in RGD is not on the processed list of incoming rgd ids
        //   then we do not delete any positions for these genes)
        Set<Integer> processedRgdIds = new HashSet<>();
        for( Map.Entry<String,MapData> entry: positionsIncoming.entrySet() ) {
            processedRgdIds.add(entry.getValue().getRgdId());
        }
        log.info("genes processed: "+processedRgdIds.size());

        Set<String> keysForDelete = new HashSet<>(positionsInRgd.keySet());
        keysForDelete.removeAll(positionsIncoming.keySet());
        log.info("gene positions for delete (prospective): "+keysForDelete.size());

        List<MapData> posListForDelete = new ArrayList<>();
        for( String keyForDelete: keysForDelete ) {
            MapData posForDelete = positionsInRgd.get(keyForDelete);
            if( processedRgdIds.contains(posForDelete.getRgdId()) ) {
                posListForDelete.add(posForDelete);
            }
        }
        dao.deleteMapData(posListForDelete);
        log.info("gene positions actually deleted (for processed genes): "+posListForDelete.size());
    }

    String computeKey(MapData md) {
        return md.getRgdId()+"|"+md.getStartPos()+"|"+md.getStopPos()+"|"+md.getStrand()+"|"+md.getChromosome()+"|"+md.getMapKey()+"|"+md.getSrcPipeline();
    }
}
