package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.impl.GeneDAO;
import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.dao.impl.RGDManagementDAO;
import edu.mcw.rgd.dao.impl.XdbIdDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.PipelineLogger;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */
public class Parser

{
    static Logger statuslog = Logger.getLogger("statuscheck");
    static Logger mismatchlog=Logger.getLogger("mismatchcheck");
    static Logger externalidexsistnorgdid=Logger.getLogger("externalidexsistnorgdid");
    static Logger ncbiexistbutnorgdid=Logger.getLogger("ncbiexistbutnorgdid");
    static Logger test4log=Logger.getLogger("test4");
    static Logger matchedrecords=Logger.getLogger("matchedrecords");
    //matchedrecords.info();
    String strand;
    PipelineLogger dbLogger = PipelineLogger.getInstance();
    EnsemblDAO ensemblDAO = new EnsemblDAO();
    Databaseloading dbloading;
    RGDManagementDAO managementDAO = new RGDManagementDAO();
    GeneDAO geneDAO = new GeneDAO();
    MapDAO mapDAO = new MapDAO();
    XdbIdDAO xdbDAO = new XdbIdDAO();
    EnsemblGene gene;
    Connection conn;
    PreparedStatement preparedstatement = null;
    String ENSEMBL_GENE_ID;
    String CHROMOSE_NAME;
    String START_POSITION;
    String END_POSITION;
    String STRAND;
    String RGD_ID;
    String Map_key = "361";
    int xdbkey;
    String entrezgene_id;
    int speciesTypeKey;

    public Parser() throws Exception {
    }
    public List<EnsemblGene> parseGene(String inputFile) throws Exception {
        dbLogger.log("  parsing gene file ", inputFile, PipelineLogger.INFO);
        BufferedReader reader = new BufferedReader(new FileReader(inputFile));
        conn = dbloading.getConnection();
        List<EnsemblGene> genes = new ArrayList<EnsemblGene>();
        String line;
        while ((line = reader.readLine()) != null) {
            String[] cols = line.split("\t", -1);
            if (cols.length < 6) {
                System.out.println("\n" + cols);
                throw new Exception("5 columns expected, but found only " + cols.length + " in the file " + inputFile + "\n" +
                        "  offending line: [" + line + "]");
            }
            String ensemblGeneId = cols[0];
            String chromosome = cols[1];
            String startPos = cols[2];
            String stopPos = cols[3];
            String strand = cols[4];
            String entrezgene_id = cols[5];
            String gene_name = cols[6];
            String gene_biotype = cols[7];
            String gene_description = cols[8];
            String rgdid = "";
            if (cols.length > 9) {
                rgdid = cols[9];
            }
            if (gene != null) {
                genes.add(gene);
                gene = null;
            }
            if (gene == null) {
                gene = new EnsemblGene();
                gene.setEnsemblGeneId(ensemblGeneId);
            }
            if (rgdid.isEmpty()) {
                gene.setrgdid("0");
            } else {
                gene.setrgdid(rgdid);
            }
            if (entrezgene_id.isEmpty()) {
                gene.setEntrezgene_id("NIL");
            } else {
                gene.setEntrezgene_id(entrezgene_id);
            }
            gene.setChromosome(chromosome);
            gene.setgene_biotype(gene_biotype);
            gene.setgene_name(gene_name);
            gene.setStartPos(startPos);
            gene.setStopPos(stopPos);
            gene.setgene_description(gene_description);
            if (strand.equals("1"))
                gene.setStrand("+");
            else if (strand.equals("-1"))
                gene.setStrand("-");
            else
                gene.setStrand(strand);

        }
        //mismatchlog.info("EnsemblGene_Id"+ "\t" + "Entrezgene_id" + "\t" + "rgd_id" + "\t" + "Start_Pos" + "\t" + "Stop_Pos" + "\t" + "Chromosome" + "\t" + "Strand");
        matchedrecords.info("EnsemblGene_Id"+ "\t" + "Entrezgene_id" + "\t" + "rgd_id" + "\t" + "Start_Pos" + "\t" + "Stop_Pos" + "\t" + "Chromosome" + "\t" + "Strand");
        for (EnsemblGene gene : genes)
        {
            int rgdid_new = Integer.parseInt(gene.getrgdid());
            MapData mapData = new MapData();
            mapData.setSrcPipeline("Ensembl");
            mapData.setChromosome(gene.getChromosome());
            mapData.setMapKey(361);
            mapData.setRgdId(rgdid_new);
            mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
            mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
            mapData.setStrand(gene.getStrand());
            if (gene.getrgdid() != "0")
                        try {
                            if (managementDAO.getRgdId2(rgdid_new) != null)
                            {
                                if (!ensemblDAO.checkobjectstatus(rgdid_new).contains("ACTIVE"))
                                {
                                    statuslog.info(gene.getrgdid() + "\t" + gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + ensemblDAO.checkobjectstatus(rgdid_new));
                                }
                                if (!gene.getEntrezgene_id().contains("NIL")) {
                                    if (!(ensemblDAO.getXdbIds(3, rgdid_new).isEmpty()) && !(ensemblDAO.getXdbIds(20, rgdid_new).isEmpty())) {
                                        if (ensemblDAO.getXdbIds(3, rgdid_new).contains(gene.getEntrezgene_id()) && ensemblDAO.getXdbIds(20, rgdid_new).contains(gene.getEnsemblGeneId()))
                                        {
                                          if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                            {
                                                matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                                mapDAO.insertMapData(mapData);
                                            }
                                        } else
                                            mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());

                                    }
                                    if (!(ensemblDAO.getXdbIds(3, rgdid_new).isEmpty()) && (ensemblDAO.getXdbIds(20, rgdid_new).isEmpty())) {
                                        if (ensemblDAO.getXdbIds(3, rgdid_new).contains(gene.getEntrezgene_id()))
                                        {
                                         if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                            {
                                                matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                                mapDAO.insertMapData(mapData);
                                          }
                                        } else
                                            mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                    }
                                    if ((ensemblDAO.getXdbIds(3, rgdid_new).isEmpty()) && !(ensemblDAO.getXdbIds(20, rgdid_new).isEmpty())) {
                                        if (ensemblDAO.getXdbIds(20, rgdid_new).contains(gene.getEnsemblGeneId()))
                                        {
                                            if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                            {
                                                matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                                mapDAO.insertMapData(mapData);
                                          }
                                        } else
                                            mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                    }

                                }
                                if (gene.getEntrezgene_id().contains("NIL")) {
                                    if (ensemblDAO.getXdbIds(20, rgdid_new).contains(gene.getEnsemblGeneId()))
                                    {
                                      if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                        {
                                            matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                            mapDAO.insertMapData(mapData);
                                        }
                                    } else
                                        mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());

                                }
                            }
                        } catch (Exception e) {
                            System.out.println("RGD id is" + gene.getrgdid());
                            e.printStackTrace();
                        }


                if (gene.getrgdid() == "0")
                {
                    //dbloading.addgene(conn, preparedstatement, gene.getEnsemblGeneId(), gene.getChromosome(), gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getrgdid(), gene.getMap_key(), gene.getEntrezgene_id());
                    if ((ensemblDAO.checkXdbIds(gene.getEntrezgene_id()).contains(gene.getEntrezgene_id())) || (ensemblDAO.checkXdbIds(gene.getEnsemblGeneId()).contains(gene.getEnsemblGeneId()))) {
                 if ((ensemblDAO.checkXdbIds(gene.getEntrezgene_id()).contains(gene.getEntrezgene_id())) && (ensemblDAO.checkXdbIds(gene.getEnsemblGeneId()).contains(gene.getEnsemblGeneId()))) {
                     externalidexsistnorgdid.info("Record with ncbi id and ensembl id exsists in rgd_acc_xdb table with, Ncbi id:" + gene.getEntrezgene_id() + ", Ensembl id:" + gene.getEnsemblGeneId() + " but rgd_id is not found in the file");
                 } else if ((ensemblDAO.checkXdbIds(gene.getEntrezgene_id()).contains(gene.getEntrezgene_id()))) {
                     externalidexsistnorgdid.info("Record with ncbi id exsists in rgd_acc_xdb table with, Ncbi id:" + gene.getEntrezgene_id() + ", Ensembl id:" + gene.getEnsemblGeneId() + " but rgd_id is not found in the file");
                 } else if ((ensemblDAO.checkXdbIds(gene.getEnsemblGeneId()).contains(gene.getEnsemblGeneId()))) {
                     externalidexsistnorgdid.info("Record with ensemblgene id exsists in rgd_acc_xdb table with, Ncbi id:" + gene.getEntrezgene_id() + ", Ensembl id:" + gene.getEnsemblGeneId() + " but rgd_id is not found in the file");
                 }
                    }
                    // else
                    //{dbloading.addgene(conn, preparedstatement, gene.getEnsemblGeneId(), gene.getChromosome(), gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getrgdid(), gene.getMap_key(), gene.getEntrezgene_id(),gene.getgene_name(),gene.getgene_biotype(),gene.getgene_description());}
                    else
                        {
                        if (gene.getEntrezgene_id().equals("NIL"))
                        {
                            if (ensemblDAO.checkrecord_rgdid(gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                            {
                                test4log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                /*String geneTypeLc = gene.getgene_biotype().toLowerCase();
                                 if (!ensemblDAO.existsGeneType(geneTypeLc))
                                 ensemblDAO.createGeneType(geneTypeLc);
                                RgdId newRgdId = ensemblDAO.createRgdId(RgdId.OBJECT_KEY_GENES, SpeciesType.RAT);
                                rgdid_new = newRgdId.getRgdId();
                                Gene newGene = new Gene();
                                newGene.setSymbol(gene.getgene_name());
                                newGene.setRgdId(rgdid_new);
                                newGene.setType(geneTypeLc);
                                newGene.setName(gene.getgene_description());
                                ensemblDAO.insertGene(newGene);
                                mapData.setRgdId(rgdid_new);
                                mapDAO.insertMapData(mapData);
                               //matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                //List<XdbId> ensemblXdbIds = new ArrayList<XdbId>();
                                XdbId xdbId = new XdbId();
                                xdbId.setRgdId(rgdid_new);
                                xdbId.setSrcPipeline("Ensembl");
                                xdbId.setAccId(gene.getEnsemblGeneId());
                                xdbId.setXdbKey(XdbId.XDB_KEY_ENSEMBL_GENES);
                                ensemblDAO.insertXdbIds(xdbId);*/
                            }

                        }
                        else
                        {
                            ncbiexistbutnorgdid.info("Record has an ncbi_id  but no rgdid:" + gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                        }
                    }
                }

                //dbloading.addgene(conn,preparedstatement,genes.get(s).getEnsemblGeneId(),genes.get(s).getChromosome(),genes.get(s).getStartPos(),genes.get(s).getStopPos(),genes.get(s).getStrand(),genes.get(s).getrgdid(),genes.get(s).getMap_key(),genes.get(s).getEntrezgene_id());


        }
            //dbloading.addgene(conn,preparedstatement,genes.get(s).getEnsemblGeneId(),genes.get(s).getChromosome(),genes.get(s).getStartPos(),genes.get(s).getStopPos(),genes.get(s).getStrand(),genes.get(s).getrgdid(),genes.get(s).getMap_key());
            reader.close();
            return genes;
        }



    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }



    //dbloading.addgene(conn, preparedstatement, gene.getEnsemblGeneId(), gene.getChromosome(), gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getrgdid(), gene.getMap_key());

    }








