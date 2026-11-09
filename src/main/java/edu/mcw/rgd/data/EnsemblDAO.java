package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.AbstractDAO;
import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.spring.MapDataQuery;
import edu.mcw.rgd.dao.spring.StringListQuery;
import edu.mcw.rgd.dao.spring.StringMapQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.jdbc.object.BatchSqlUpdate;

import javax.sql.DataSource;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
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

    /// preload (once per species) the map of Ensembl gene id -> active gene RGD id, so the transcript
    /// loader does not call getEnsemblRgdId() per transcript. Mirrors getEnsemblRgdId (xdb_key=20, Ensembl).
    public java.util.Map<String,String> getEnsemblGeneRgdIdMap( int speciesTypeKey ) throws Exception {
        String sql = """
            SELECT rx.acc_id, rx.rgd_id FROM rgd_acc_xdb rx, rgd_ids r
            WHERE rx.xdb_key = 20 AND rx.src_pipeline = 'Ensembl'
               AND r.rgd_id = rx.rgd_id AND r.object_status = 'ACTIVE'
               AND r.object_key = 1 AND r.species_type_key = ?
            """;
        java.util.Map<String,String> map = new HashMap<>();
        for( StringMapQuery.MapPair pair: StringMapQuery.execute(this, sql, speciesTypeKey) ) {
            map.putIfAbsent(pair.keyValue, pair.stringValue);
        }
        return map;
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
     String query="select * from maps_data where map_key=? and start_pos=CAST(? AS NUMERIC) and stop_pos=CAST(? AS NUMERIC) and strand=? and chromosome=?";
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

    /// chromosomes for the assembly being loaded: prefer the Ensembl (loading) map's own chromosomes;
    /// fall back to the NCBI map (legacy case: Ensembl and NCBI share the assembly and the chromosomes
    /// are stored only on the NCBI map). An Ensembl-only assembly (no NCBI counterpart, f.e. naked mole-rat
    /// 'Naked mole-rat maternal') keeps its chromosomes on the Ensembl map itself.
    public List<Chromosome> getLoadingChromosomes(int ensemblMapKey, int ncbiMapKey) throws Exception {
        List<Chromosome> chrs = getChromosomes(ensemblMapKey);
        if( chrs!=null && !chrs.isEmpty() ) {
            return chrs;
        }
        return getChromosomes(ncbiMapKey);
    }

    /// true unless the two maps are known to be DIFFERENT assemblies (both carry a GenBank accession and the
    /// accessions differ). Used to decide whether NCBI-gene-by-position matching is valid -- it only is when
    /// the NCBI map shares the Ensembl map's coordinate system. If either accession is missing we cannot tell,
    /// so we assume they are the same (preserves legacy behavior for maps that have no accession populated).
    public boolean mapsShareAssembly(int mapKeyA, int mapKeyB) throws Exception {
        String accA = getAssemblyMap(mapKeyA).getGenBankAssemblyAcc();
        String accB = getAssemblyMap(mapKeyB).getGenBankAssemblyAcc();
        if( Utils.isStringEmpty(accA) || Utils.isStringEmpty(accB) ) {
            return true;
        }
        return Utils.stringsAreEqualIgnoreCase(accA, accB);
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

    ///// POLYPHEN -- variant tables live in the Carpe database

    interface RowHandler {
        void handle(ResultSet rs) throws Exception;
    }

    DataSource getCarpeDataSource() throws Exception {
        return DataSourceFactory.getInstance().getCarpeNovoDataSource();
    }

    /// streams the rows of a query from the Carpe database (Postgres streams only with autocommit off)
    void queryCarpe(String sql, RowHandler handler, Object... params) throws Exception {
        try( Connection conn = getCarpeDataSource().getConnection() ) {
            conn.setAutoCommit(false);
            try( PreparedStatement ps = conn.prepareStatement(sql) ) {
                ps.setFetchSize(10000);
                for( int i=0; i<params.length; i++ ) {
                    if( params[i] instanceof Long[] array ) {
                        ps.setArray(i+1, conn.createArrayOf("bigint", array));
                    } else {
                        ps.setObject(i+1, params[i]);
                    }
                }
                try( ResultSet rs = ps.executeQuery() ) {
                    while( rs.next() ) {
                        handler.handle(rs);
                    }
                }
            } finally {
                conn.rollback();
            }
        }
    }

    /// missense variant transcripts on given assemblies: one reference and one variant amino acid, different;
    /// returns map_key -> rows. All assemblies are read in ONE query: VARIANT_TRANSCRIPT is huge (2 billion rows)
    /// and the planner reads it with a full table scan, so one scan for all assemblies instead of one per assembly
    public java.util.Map<Integer, List<PolyphenLoader.VtRow>> getMissenseVariantTranscripts(Collection<Integer> mapKeys) throws Exception {
        String sql = """
            SELECT map_key, variant_rgd_id, transcript_rgd_id, full_ref_aa_pos, ref_aa, var_aa, full_ref_aa_seq_key
            FROM variant_transcript
            WHERE map_key = ANY(?) AND full_ref_aa_seq_key IS NOT NULL AND full_ref_aa_pos IS NOT NULL
              AND LENGTH(ref_aa)=1 AND LENGTH(var_aa)=1 AND ref_aa<>var_aa
            """;
        java.util.Map<Integer, List<PolyphenLoader.VtRow>> result = new HashMap<>();
        Long[] mapKeyArray = mapKeys.stream().map(Long::valueOf).toArray(Long[]::new);
        queryCarpe(sql, rs -> result.computeIfAbsent(rs.getInt(1), k -> new ArrayList<>())
            .add(new PolyphenLoader.VtRow(rs.getLong(2), rs.getLong(3), rs.getInt(4),
                rs.getString(5).charAt(0), rs.getString(6).charAt(0), rs.getLong(7))), (Object) mapKeyArray);
        return result;
    }

    /// protein sequences by seq_key: seq_key -> md5 and length of the protein sequence without the trailing
    /// stop codon '*' (the same convention as Ensembl's translation_md5); looked up in batches by key,
    /// which uses the seq_key and data_md5 indexes instead of scanning SEQ_DATA
    public java.util.Map<Long, PolyphenLoader.SeqInfo> getProteinSequences(Collection<Long> seqKeys) throws Exception {
        String sql = """
            SELECT s.seq_key, MD5(RTRIM(d.seq_data, '*')), LENGTH(RTRIM(d.seq_data, '*'))
            FROM rgd_sequences s
            JOIN seq_data d ON d.data_md5=s.seq_data_md5
            WHERE s.seq_key = ANY(?)
            """;
        java.util.Map<Long, PolyphenLoader.SeqInfo> result = new HashMap<>();
        for( Long[] batch: batches(seqKeys) ) {
            queryCarpe(sql, rs -> result.put(rs.getLong(1), new PolyphenLoader.SeqInfo(rs.getString(2), rs.getInt(3))), (Object) batch);
        }
        return result;
    }

    /// transcript_rgd_id -> { protein accession, gene symbol }; looked up in batches by key
    public java.util.Map<Long, String[]> getTranscriptProteins(Collection<Long> transcriptRgdIds) throws Exception {
        String sql = """
            SELECT t.transcript_rgd_id, t.protein_acc_id, g.gene_symbol
            FROM transcripts t
            LEFT JOIN genes g ON g.rgd_id=t.gene_rgd_id
            WHERE t.transcript_rgd_id = ANY(?)
            """;
        java.util.Map<Long, String[]> result = new HashMap<>();
        for( Long[] batch: batches(transcriptRgdIds) ) {
            queryCarpe(sql, rs -> result.put(rs.getLong(1), new String[]{rs.getString(2), rs.getString(3)}), (Object) batch);
        }
        return result;
    }

    /// POLYPHEN rows of given variants; looked up in batches by key
    public List<PolyphenRow> getPolyphenRows(Collection<Long> variantRgdIds) throws Exception {
        String sql = """
            SELECT p.polyphen_id, p.variant_rgd_id, p.transcript_rgd_id, p.gene_symbol, p.protein_id, p.position,
              p.aa1, p.aa2, p.prediction, p.pph2_class, p.pph2_prob
            FROM polyphen p
            WHERE p.variant_rgd_id = ANY(?)
            """;
        List<PolyphenRow> rows = new ArrayList<>();
        RowHandler handler = rs -> {
            PolyphenRow r = new PolyphenRow();
            r.polyphenId = rs.getLong(1);
            r.variantRgdId = rs.getLong(2);
            r.transcriptRgdId = rs.getLong(3);
            r.geneSymbol = rs.getString(4);
            r.proteinId = rs.getString(5);
            r.position = Integer.parseInt(rs.getString(6));
            r.aa1 = rs.getString(7);
            r.aa2 = rs.getString(8);
            r.prediction = rs.getString(9);
            r.pph2Class = rs.getString(10);
            r.pph2Prob = rs.getString(11);
            rows.add(r);
        };
        for( Long[] batch: batches(variantRgdIds) ) {
            queryCarpe(sql, handler, (Object) batch);
        }
        return rows;
    }

    /// splits ids into arrays of up to 1000, for '= ANY(?)' lookups by key
    List<Long[]> batches(Collection<Long> ids) {
        List<Long> list = new ArrayList<>(ids);
        List<Long[]> result = new ArrayList<>();
        for( int i=0; i<list.size(); i+=1000 ) {
            result.add(list.subList(i, Math.min(i+1000, list.size())).toArray(new Long[0]));
        }
        return result;
    }

    public void insertPolyphenRows(List<PolyphenRow> rows) throws Exception {
        String sql = """
            INSERT INTO polyphen (polyphen_id, variant_rgd_id, transcript_rgd_id, gene_symbol, protein_id, position,
              aa1, aa2, o_aa1, o_aa2, prediction, pph2_class, pph2_prob, protein_status, creation_date)
            VALUES (NEXTVAL('polyphen_seq'), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, '100 PERC MATCH', LOCALTIMESTAMP(0))
            """;
        BatchSqlUpdate su = new BatchSqlUpdate(getCarpeDataSource(), sql, new int[]{Types.BIGINT, Types.BIGINT,
            Types.VARCHAR, Types.VARCHAR, Types.VARCHAR, Types.VARCHAR, Types.VARCHAR, Types.VARCHAR, Types.VARCHAR,
            Types.VARCHAR, Types.VARCHAR, Types.VARCHAR}, 10000);
        su.compile();
        for( PolyphenRow r: rows ) {
            su.update(r.variantRgdId, r.transcriptRgdId, r.geneSymbol, r.proteinId, Integer.toString(r.position),
                r.aa1, r.aa2, r.aa1, r.aa2, r.prediction, r.pph2Class, r.pph2Prob);
        }
        su.flush();
    }

    public void updatePolyphenRows(List<PolyphenRow> rows) throws Exception {
        String sql = """
            UPDATE polyphen SET gene_symbol=?, prediction=?, pph2_class=?, pph2_prob=?
            WHERE polyphen_id=?
            """;
        BatchSqlUpdate su = new BatchSqlUpdate(getCarpeDataSource(), sql, new int[]{Types.VARCHAR, Types.VARCHAR,
            Types.VARCHAR, Types.VARCHAR, Types.BIGINT}, 10000);
        su.compile();
        for( PolyphenRow r: rows ) {
            su.update(r.geneSymbol, r.prediction, r.pph2Class, r.pph2Prob, r.polyphenId);
        }
        su.flush();
    }

    public void deletePolyphenRows(List<PolyphenRow> rows) throws Exception {
        String sql = """
            DELETE FROM polyphen WHERE polyphen_id=?
            """;
        BatchSqlUpdate su = new BatchSqlUpdate(getCarpeDataSource(), sql, new int[]{Types.BIGINT}, 10000);
        su.compile();
        for( PolyphenRow r: rows ) {
            su.update(r.polyphenId);
        }
        su.flush();
    }
}

