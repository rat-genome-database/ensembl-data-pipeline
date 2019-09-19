package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.datamodel.MapData;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
/**
 * Created by sellanki on 8/29/2019.
 */
public class LogLoader
{
    MapDAO mapDAO = new MapDAO();
    PreparedStatement preparedstatement = null;
    static Logger mismatches_loaded_records_different_rgdid = Logger.getLogger("mismatches_loaded_records_different_rgdid");
    static Logger mismatches_already_loaded_records= Logger.getLogger("mismatches_already_loaded_records");
    static Logger mismatches_ensemblid_doesnotexist_rgdxdb=Logger.getLogger("mismatches_ensemblid_doesnotexist_rgdxdb");
    static Logger test3log=Logger.getLogger("test3");
    static Logger mismatches_newposition_ensemblidexists_rgdxdb=Logger.getLogger("mismatches_newposition_ensemblidexists_rgdxdb");
    EnsemblDAO ensemblDAO = new EnsemblDAO();
    Databaseloading dbloading;
    EnsemblGene gene=null;
    Connection conn;
    public LogLoader() throws Exception {
    }
    public List<EnsemblGene> parseGene(String inputFile) throws Exception
    {
        BufferedReader reader = new BufferedReader(new FileReader(inputFile));
        //conn = dbloading.getConnection();
        List<EnsemblGene> log_genes = new ArrayList<EnsemblGene>();
        String line;
        while ((line = reader.readLine()) != null)
        {
            String[] cols = line.split("\t", -1);
            if (cols.length < 6) {
                System.out.println("\n" + cols);
                throw new Exception("5 columns expected, but found only " + cols.length + " in the file " + inputFile + "\n" +
                        "  offending line: [" + line + "]");
            }
            String ensemblGeneId = cols[0];
            String entrezgene_id = cols[1];
            String rgdid = cols[2];
            String startPos = cols[3];
            String stopPos = cols[4];
            String chromosome = cols[5];
            String strand = cols[6];
            if (gene != null) {
                log_genes.add(gene);
                gene = null;
            }
            if (gene == null) {
                gene = new EnsemblGene();
                gene.setEnsemblGeneId(ensemblGeneId);
            }
            gene.setEnsemblGeneId(ensemblGeneId);
            gene.setEntrezgene_id(entrezgene_id);
            gene.setrgdid(rgdid);
            gene.setStartPos(startPos);
            gene.setStopPos(stopPos);
            gene.setChromosome(chromosome);
            gene.setStrand(strand);
        }
        Map<String,EnsemblGene> load=new HashMap<>();
        for(EnsemblGene gene:log_genes)
        {
            int rgdid_new= Integer.parseInt(gene.getrgdid());
            MapData mapData = new MapData();
            mapData.setSrcPipeline("Ensembl");
            mapData.setChromosome(gene.getChromosome());
            mapData.setMapKey(361);
            mapData.setRgdId(rgdid_new);
            mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
            mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
            mapData.setStrand(gene.getStrand());
            List<String> rgd_ids = ensemblDAO.getrgdid(gene.getStartPos(), gene.getStopPos(), gene.getChromosome(), gene.getStrand());
            if (!rgd_ids.isEmpty() && !rgd_ids.contains(gene.getrgdid()))
            {
                if (ensemblDAO.getXdbIds(20, Integer.parseInt(gene.getrgdid())).contains(gene.getEnsemblGeneId()))
                {
                    load.put(gene.getEnsemblGeneId(), gene);
                    mismatches_loaded_records_different_rgdid.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                } else
                    mismatches_ensemblid_doesnotexist_rgdxdb.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
            }
            else if (rgd_ids.isEmpty())
            {
                if (ensemblDAO.getXdbIds(20, Integer.parseInt(gene.getrgdid())).contains(gene.getEnsemblGeneId()))
                {
                   mapDAO.insertMapData(mapData);
                    mismatches_newposition_ensemblidexists_rgdxdb.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                } else
                    mismatches_ensemblid_doesnotexist_rgdxdb.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
            }
            else {
                mismatches_already_loaded_records.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());

            }
        }
       for(Map.Entry<String,EnsemblGene> entry :load.entrySet())
        {
            int rgdid_new= Integer.parseInt(gene.getrgdid());
            MapData mapData = new MapData();
            mapData.setSrcPipeline("Ensembl");
            mapData.setChromosome(gene.getChromosome());
            mapData.setMapKey(361);
            mapData.setRgdId(rgdid_new);
            mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
            mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
            mapData.setStrand(gene.getStrand());
            test3log.info(entry.getKey()+"\t"+entry.getValue().getEntrezgene_id()+"\t"+entry.getValue().getrgdid()+"\t"+entry.getValue().getStartPos()+"\t"+entry.getValue().getStopPos()+"\t"+entry.getValue().getChromosome()+"\t"+entry.getValue().getStrand());
           mapDAO.insertMapData(mapData);
        }
        reader.close();
        return log_genes;
        }


    }




