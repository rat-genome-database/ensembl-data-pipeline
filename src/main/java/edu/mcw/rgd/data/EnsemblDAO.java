package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.AbstractDAO;
import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.spring.MapDataQuery;
import edu.mcw.rgd.dao.spring.StringListQuery;
import edu.mcw.rgd.datamodel.*;

import java.util.List;


/**
 * Created by sellanki on 8/6/2019.
 */

public class EnsemblDAO extends AbstractDAO {

    XdbIdDAO xdbDAO = new XdbIdDAO();
    AliasDAO aliasDAO = new AliasDAO();
    GeneDAO geneDAO = new GeneDAO();
    RGDManagementDAO managementDAO = new RGDManagementDAO();
    MapDAO mapDAO = new MapDAO();
    AliasDAO aliasdao = new AliasDAO();
    int[] primaryMapKey = new int[4];
    int speciesTypeKey;

    public EnsemblDAO() throws Exception {
        primaryMapKey[1] = mapDAO.getPrimaryRefAssembly(SpeciesType.HUMAN).getKey();
        primaryMapKey[2] = mapDAO.getPrimaryRefAssembly(SpeciesType.MOUSE).getKey();
        primaryMapKey[3] = mapDAO.getPrimaryRefAssembly(SpeciesType.RAT).getKey();
    }

    public List<String> getGeneType(int rgdId) throws Exception
    {
        String sql="select gene_type_lc from genes where rgd_id=?";
        return StringListQuery.execute(this,sql,rgdId);
    }
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
    public void insertGeneSource(int rgdId,String gene_source) throws Exception
    {
        String sql="update genes set gene_source=? where rgd_id=?";
        update(sql,gene_source,rgdId);
    }
    public String getGeneStringName(List<String> gene_name) throws Exception
    {
        String result=String.join(" ",gene_name);
        return result;

    }
    public List<String> checkobjectstatus(int rgdId) throws Exception
    {
        String sql="select object_status from rgd_ids where rgd_id=?";
        return StringListQuery.execute(this,sql,rgdId);
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
public List<String> getXdbIds(int xdbKey, int rgdId) throws Exception
    {
        String sql = "select acc_id from rgd_acc_xdb WHERE xdb_key=? AND rgd_id=?";
        return StringListQuery.execute(this,sql,xdbKey,rgdId);

    }

    public List<String> checkXdbIds(String Acc_id) throws Exception
    {
        String sql = "select acc_id from rgd_acc_xdb WHERE acc_id=?";
        return StringListQuery.execute(this,sql,Acc_id);
    }
    public List<String> checkrgd_id(String Acc_id) throws Exception{
        String sql = " select rgd_id from rgd_acc_xdb where acc_id=?";
        return StringListQuery.execute(this,sql,Acc_id);
    }

    public int insertXdbIds(XdbId xdb) throws Exception {
        return xdbDAO.insertXdb(xdb);
    }

    public MapData checkrecord(int rgd_id,String start_pos,String stop_pos,String strand,String chromosome) throws Exception
    {
        String query="select * from maps_data where map_key=361 and rgd_id=? and start_pos=? and stop_pos=? and strand=? and chromosome=?";
        MapDataQuery q=new MapDataQuery(this.getDataSource(),query);
        List<MapData> mp=execute(q,rgd_id,start_pos,stop_pos,strand,chromosome);
        if(!mp.isEmpty())
        {
            return mp.get(0);
        }
        else{
            return null;
        }
    }
    public MapData checkrecord_rgdid(String start_pos,String stop_pos,String strand,String chromosome) throws Exception {
     String query="select * from maps_data where map_key=361 and start_pos=? and stop_pos=? and strand=? and chromosome=?";
        MapDataQuery q=new MapDataQuery(this.getDataSource(),query);
        List<MapData> mp=execute(q,start_pos,stop_pos,strand,chromosome);
        if(!mp.isEmpty())
        {
            return mp.get(0);
        }
        else{
            return null;
        }
    }

    public int insertAliasType(String aliasType,String notes) throws Exception
    {
        return aliasDAO.insertAliasType(aliasType,notes);
    }
    public List<String> getAliasTypes() throws Exception
    {
        return aliasDAO.getAliasTypes();

    }

    public int insertAlias(Alias alias) throws Exception
    {
       return aliasDAO.insertAlias(alias);
    }

}

