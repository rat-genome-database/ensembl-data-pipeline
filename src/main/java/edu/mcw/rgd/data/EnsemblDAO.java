package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.AbstractDAO;
import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.spring.CountQuery;
import edu.mcw.rgd.dao.spring.MapDataQuery;
import edu.mcw.rgd.dao.spring.StringListQuery;
import edu.mcw.rgd.datamodel.*;
import org.springframework.jdbc.core.SqlParameter;

import java.sql.Types;
import java.util.List;


/**
 * Created by sellanki on 8/6/2019.
 */

public class EnsemblDAO extends AbstractDAO {

    XdbIdDAO xdbDAO = new XdbIdDAO();
    AliasDAO aliasDAO = new AliasDAO();
    GeneDAO geneDAO = new GeneDAO();
    TranscriptDAO transcriptDAO = new TranscriptDAO();
    RGDManagementDAO managementDAO = new RGDManagementDAO();
    NomenclatureDAO nomenclatureDAO = new NomenclatureDAO();
    MapDAO mapDAO = new MapDAO();

    public void insertGeneType(int rgdId,String gene_type) throws Exception
    {
     String sql="update genes set ensembl_gene_type=? where rgd_id=?";
        update(sql,gene_type,rgdId);
    }
    public List<String> getGeneSymbol(int rgdId) throws Exception
    {
        String sql="select gene_symbol from genes where rgd_id=?";
        return StringListQuery.execute(this,sql,rgdId);
    }
    public void insertGeneSymbol(int rgdId,String gene_symbol) throws Exception
    {
        String sql="update genes set ensembl_gene_symbol=? where rgd_id=?";
        update(sql,gene_symbol,rgdId);
    }
    public List<String> getGeneName(int rgdId) throws Exception
    {
        String sql="select full_name from genes where rgd_id=?";
        return StringListQuery.execute(this,sql,rgdId);
    }
    public void insertGeneName(int rgdId,String gene_name) throws Exception
    {
        String sql="update genes set ensembl_full_name=? where rgd_id=?";
        update(sql,gene_name,rgdId);
    }

    public String getGeneStringName(List<String> gene_name) throws Exception
    {
        String result=String.join(" ",gene_name);
        return result;

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
        geneDAO.insertGene(gene);
    }
    public void updateGene(Gene gene) throws Exception {
        geneDAO.updateGene(gene);
    }
    public Gene getGene(int rgdId) throws Exception {
        return  geneDAO.getGene(rgdId);
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
    public void insertTranscriptFeature(TranscriptFeature transcriptFeature,int speciesTypeKey) throws Exception {
        transcriptDAO.createFeature(transcriptFeature,speciesTypeKey);
    }
    public void insertNomenclatureEvent(NomenclatureEvent event) throws Exception {
        nomenclatureDAO.createNomenEvent(event);
    }
    public List<TranscriptFeature> getFeatures(int transcriptRgdId) throws Exception {
        return transcriptDAO.getFeatures(transcriptRgdId);
    }
    public List<String> getRgd_id(String Acc_id, int xdbKey) throws Exception{
        String sql = " select distinct(rx.rgd_id) from rgd_acc_xdb rx, rgd_ids r where rx.acc_id=? and rx.xdb_key = ? and r.rgd_id = rx.rgd_id and r.object_status = 'ACTIVE' and r.object_key = 1";
        List<String> result =  StringListQuery.execute(this,sql,Acc_id, xdbKey);
        if(result.size() >  0)
            return result;
        else return null;
    }
    public List<String> getNcbiRgdId(String Acc_id) throws Exception{
        String sql = " select distinct(rx.rgd_id) from rgd_acc_xdb rx, rgd_ids r where rx.acc_id=? and rx.xdb_key = 3 and r.rgd_id = rx.rgd_id and r.object_status = 'ACTIVE' and r.object_key = 1 " +
                "and rx.src_pipeline = 'ENTREZGENE'";
        List<String> result =  StringListQuery.execute(this,sql,Acc_id);
        return result;
    }
    public String getEnsemblRgdId(String Acc_id) throws Exception{
        String sql = " select distinct(rx.rgd_id) from rgd_acc_xdb rx, rgd_ids r where rx.acc_id=? and rx.xdb_key = 20 and r.rgd_id = rx.rgd_id and r.object_status = 'ACTIVE' and r.object_key = 1 " +
                "and rx.src_pipeline = 'Ensembl'";
        List<String> result =  StringListQuery.execute(this,sql,Acc_id);
        if(result.size() >  0)
            return result.get(0);
        else return null;
    }

    public int insertXdbIds(XdbId xdb) throws Exception {
        return xdbDAO.insertXdb(xdb);
    }

    public MapData checkrecord(int rgd_id,String start_pos,String stop_pos,String strand,String chromosome,int mapKey) throws Exception
    {
        String query="select * from maps_data where map_key=? and rgd_id=? and start_pos=? and stop_pos=? and strand=? and chromosome=?";
        MapDataQuery q=new MapDataQuery(this.getDataSource(),query);
        List<MapData> mp=execute(q,mapKey,rgd_id,start_pos,stop_pos,strand,chromosome);
        if(!mp.isEmpty())
        {
            return mp.get(0);
        }
        else{
            return null;
        }
    }
    public boolean checkXDBRecord(int rgdId,String ensemblId,String srcPipeline) throws Exception{
        String sql = "select count(*) from rgd_acc_xdb where rgd_id = ? and acc_id = ? and src_pipeline = ?";
        CountQuery q = new CountQuery(getDataSource(), sql);
        q.declareParameter(new SqlParameter(Types.INTEGER));
        q.declareParameter(new SqlParameter(Types.VARCHAR));
        q.declareParameter(new SqlParameter(Types.VARCHAR));
        q.compile();
        List<Integer> result = q.execute(rgdId,ensemblId,srcPipeline);
        if(result.get(0) > 0)
            return true;
        else return false;
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

    public int insertAliasType(String aliasType,String notes) throws Exception {
        return aliasDAO.insertAliasType(aliasType,notes);
    }

    public List<String> getAliasTypes() throws Exception {
        return aliasDAO.getAliasTypes();
    }

    public int insertAlias(Alias alias) throws Exception {
       return aliasDAO.insertAlias(alias);
    }

    public List<String> getChromosomes(int mapKey) throws Exception {

        String sql = "SELECT DISTINCT chromosome FROM CHROMOSOMES WHERE map_key=? ";
        return StringListQuery.execute(mapDAO, sql, mapKey);
    }

    public void insertMapData(MapData md) throws Exception {
        mapDAO.insertMapData(md);
    }
}

