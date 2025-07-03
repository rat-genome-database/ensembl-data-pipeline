package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.AbstractDAO;
import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.spring.MapDataQuery;
import edu.mcw.rgd.dao.spring.StringListQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.List;


/**
 * Created by sellanki on 8/6/2019.
 */
public class EnsemblDAO extends AbstractDAO {

    Logger logInsertedXdbs = LogManager.getLogger("inserted_xdbs");
    Logger logInsertedGenePos = LogManager.getLogger("inserted_gene_pos");
    Logger logDeletedGenePos = LogManager.getLogger("deleted_gene_pos");
    Logger logExonsInserted = LogManager.getLogger("exons_inserted");
    Logger logExonsRemoved = LogManager.getLogger("exons_removed");
    Logger logUtrsRemoved = LogManager.getLogger("utrs_removed");
    Logger logUpdatedPos = LogManager.getLogger("pos_updated");
    Logger logGenesInserted = LogManager.getLogger("genes_inserted");

    AliasDAO aliasDAO = new AliasDAO();
    GeneDAO geneDAO = new GeneDAO();
    TranscriptDAO transcriptDAO = new TranscriptDAO();
    RGDManagementDAO managementDAO = new RGDManagementDAO();
    NomenclatureDAO nomenclatureDAO = new NomenclatureDAO();
    MapDAO mapDAO = new MapDAO();
    XdbIdDAO xdbDAO = new XdbIdDAO();

    RgdId getRgdId(int rgdId) throws Exception {
        return managementDAO.getRgdId2(rgdId);
    }

    RgdId createRgdId(int objectKey, int speciesTypeKey) throws Exception {
        return managementDAO.createRgdId(objectKey, "ACTIVE", "created by Ensembl pipeline", speciesTypeKey);
    }

    public boolean existsGeneType(String geneType) throws Exception {
        return geneDAO.existsGeneType(geneType);
    }

    public void createGeneType(String geneType) throws Exception {
        geneDAO.createGeneType(geneType, geneType, geneType);
    }

    public void insertGene(Gene gene) throws Exception {
        logGenesInserted.debug(gene.dump("|"));
        geneDAO.insertGene(gene);
    }

    public void updateGene(Gene gene) throws Exception {
        geneDAO.updateGene(gene);
    }

    /// TODO: gene module -- processing bottleneck
    public Gene getGene(int rgdId) throws Exception {
        return geneDAO.getGene(rgdId);
    }

    public void insertTranscript(Transcript transcript,int speciesTypeKey) throws Exception {
        transcriptDAO.createTranscript(transcript,speciesTypeKey);
    }
    public void updateTranscript(Transcript transcript) throws Exception {
        transcriptDAO.updateTranscript(transcript);
    }
    public Transcript getTranscript(int rgdId) throws Exception {
        return transcriptDAO.getTranscript(rgdId);
    }
    List<Transcript> getTranscriptsForGene(int geneRgdId) throws Exception {
        return transcriptDAO.getTranscriptsForGene(geneRgdId);
    }
    List<Transcript> getTranscriptsByAccId(String accId) throws Exception {
        return transcriptDAO.getTranscriptsByAccId(accId);
    }

    public void insertTranscriptFeature(TranscriptFeature tf, int speciesTypeKey) throws Exception {

        if( tf.getFeatureType()== TranscriptFeature.FeatureType.EXON ) {
            logExonsInserted.debug("TR RGD:" + tf.getTranscriptRgdId() + ", exon RGD:" + tf.getRgdId()
                    + " chr" + tf.getChromosome() + ":" + tf.getStartPos() + ".." + tf.getStopPos() + " (" + tf.getStrand() + ")");
        }

        transcriptDAO.createFeature(tf, speciesTypeKey);
    }

    public int unbindFeatureFromTranscript(Transcript tr, TranscriptFeature tf) throws Exception {

        if( tf.getFeatureType()== TranscriptFeature.FeatureType.EXON ) {
            logExonsRemoved.debug(tr.getAccId() + " RGD:" + tr.getRgdId() + ", exon RGD:" + tf.getRgdId()
                + " chr" + tf.getChromosome() + ":" + tf.getStartPos() + ".." + tf.getStopPos() + " (" + tf.getStrand() + ")"
                + " map_key="+tf.getMapKey());
        }
        else if( tf.getFeatureType()== TranscriptFeature.FeatureType.UTR3 || tf.getFeatureType()== TranscriptFeature.FeatureType.UTR5 ) {
            logUtrsRemoved.debug(tr.getAccId() + " RGD:" + tr.getRgdId() + ", "+tf.getCanonicalName()+" RGD:" + tf.getRgdId()
                + " chr" + tf.getChromosome() + ":" + tf.getStartPos() + ".." + tf.getStopPos() + " (" + tf.getStrand() + ")"
                + " map_key="+tf.getMapKey());
        }

        return transcriptDAO.unbindFeatureFromTranscript(tr.getRgdId(), tf.getRgdId());
    }

    ///// STABLE_TRANSCRIPTS - transcript versions

    public String getTranscriptVersionInfo(String acc) throws Exception {
        return transcriptDAO.getTranscriptVersionInfo(acc);
    }

    public void updateTranscriptVersionInfo(String acc, String version) throws Exception {
        transcriptDAO.updateTranscriptVersionInfo(acc, version);
    }

    public void insertTranscriptVersionInfo(String acc, String version, int rgdId) throws Exception {
        transcriptDAO.insertTranscriptVersionInfo(acc, version, rgdId);
    }

    public void insertNomenclatureEvent(NomenclatureEvent event) throws Exception {
        nomenclatureDAO.createNomenEvent(event);
    }

    public List<TranscriptFeature> getFeatures(int transcriptRgdId, int mapKey) throws Exception {
        return transcriptDAO.getFeatures(transcriptRgdId, mapKey);
    }

    public List<String> getRgd_id(String accId, int xdbKey) throws Exception{
        List<String> rgdIds = new ArrayList<>();
        List<Gene> genes = xdbDAO.getActiveGenesByXdbId(xdbKey, accId);
        for( Gene g: genes ) {
            rgdIds.add(Integer.toString(g.getRgdId()));
        }

        if( rgdIds.size()> 0 )
            return rgdIds;
        else return null;
    }

    public List<String> getNcbiRgdId( String accId ) throws Exception{
        String sql = """
            SELECT DISTINCT(rx.rgd_id) FROM rgd_acc_xdb rx, rgd_ids r
            WHERE rx.acc_id=? AND rx.xdb_key = 3 AND r.rgd_id = rx.rgd_id AND r.object_status = 'ACTIVE'
                AND r.object_key = 1 AND rx.src_pipeline = 'ENTREZGENE'
            """;
        return StringListQuery.execute(this, sql, accId);
    }

    public String getEnsemblRgdId( String accId ) throws Exception{
        String sql = """
            SELECT DISTINCT(rx.rgd_id) FROM rgd_acc_xdb rx, rgd_ids r
            WHERE rx.acc_id=? AND rx.xdb_key = 20 AND r.rgd_id = rx.rgd_id and r.object_status = 'ACTIVE'
               AND r.object_key = 1 AND rx.src_pipeline = 'Ensembl'
            """;
        List<String> result =  StringListQuery.execute(this, sql, accId);
        if(result.size() >  0)
            return result.get(0);
        else return null;
    }

    public int insertXdbIds(XdbId id) throws Exception {
        logInsertedXdbs.debug(id.dump("|"));
        return xdbDAO.insertXdb(id);
    }

    public boolean checkXDBRecord(int rgdId, String ensemblId, String srcPipeline) throws Exception{
        String sql = "select count(*) from rgd_acc_xdb where rgd_id = ? and acc_id = ? and src_pipeline = ?";
        return xdbDAO.getCount(sql, rgdId, ensemblId, srcPipeline) != 0;
    }

    public MapData checkrecord_rgdid(String start_pos,String stop_pos,String strand,String chromosome, int mapKey) throws Exception {
     String query="select * from maps_data where map_key=? and start_pos=? and stop_pos=? and strand=? and chromosome=?";
        MapDataQuery q=new MapDataQuery(this.getDataSource(),query);
        List<MapData> mp=execute(q,mapKey,start_pos,stop_pos,strand,chromosome);
        if(!mp.isEmpty())
        {
            return mp.get(0);
        }
        else{
            return null;
        }
    }

    public int insertAlias(Alias alias) throws Exception {

        // ensure alias value does not start with 'LOW QUALITY PROTEIN:' prefix
        String prefix = "LOW QUALITY PROTEIN:";
        if( alias.getValue().startsWith(prefix) ) {
            String newValue = alias.getValue().substring(prefix.length()).trim();
            if( Utils.isStringEmpty(newValue) ) {
                return 0;
            }
            alias.setValue(newValue);
        }
        if( Utils.isStringEmpty(alias.getValue()) ) {
            return 0;
        }
        return aliasDAO.insertAlias(alias);
    }

    public String matchChromosome(String chr, List<Chromosome> chromosomes) {
        // try to match by chromosome directly
        for( Chromosome c: chromosomes ) {
            if(c.getChromosome().equals(chr)) {
                return chr;
            }
        }
        // try to match by genebank id (scaffold assemblies like chinchilla, squirrel)
        for( Chromosome c: chromosomes ) {
            if( c.getGenbankId()!=null && c.getGenbankId().equals(chr)) {
                return c.getChromosome();
            }
        }
        return null;
    }

    public List<Chromosome> getChromosomes(int mapKey) throws Exception {
        return mapDAO.getChromosomes(mapKey);
    }

    public void insertMapData(MapData md, String prefix) throws Exception {
        logInsertedGenePos.info(prefix+md.dump("|"));
        mapDAO.insertMapData(md);
    }

    public List<MapData> getGenesInRegion(int mapKey, String chr, int startPos, int stopPos) throws Exception {
        return mapDAO.getMapDataWithinRange(startPos, stopPos, chr, mapKey, 0);
    }
    public void updateMapData(MapData md, MapData mdOld, String prefix) throws Exception {
        logUpdatedPos.info(prefix+" OLD>"+mdOld.dump("|"));
        logUpdatedPos.info(prefix+" NEW>"+md.dump("|"));
        mapDAO.updateMapData(md);
    }

    public void deleteMapData(List<MapData> posList) throws Exception {
        for( MapData md: posList ) {
            logDeletedGenePos.info(md.dump("|"));
        }
        mapDAO.deleteMapData(posList);
    }

    public Map getAssemblyMap(int mapKey) throws Exception {
        return mapDAO.getMapByKey(mapKey);
    }

    public List<MapData> getGenePositions(int mapKey) throws Exception {
        return mapDAO.getMapDataByMapKeyObject(mapKey, RgdId.OBJECT_KEY_GENES, "Ensembl");
    }
}

