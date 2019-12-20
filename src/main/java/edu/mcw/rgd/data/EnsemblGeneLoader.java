package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.dao.impl.RGDManagementDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Created by sellanki on 10/18/2019.
 */
public class EnsemblGeneLoader {
    EnsemblDAO ensemblDAO;
    MapDAO mapDAO = new MapDAO();
    static Logger statuslog = Logger.getLogger("statuscheck");

    Map<String,String> matches=new HashMap<String,String>();
    List mismatches=new ArrayList();
    List newGenes = new ArrayList<>();
    public EnsemblGeneLoader() throws Exception {
        ensemblDAO = new EnsemblDAO();
    }
    public void run(Collection<EnsemblGene> genes,int speciesTypeKey) throws Exception
    {
        System.out.println("Loading the file");
        int mapKey = 0;
        List<String> chromosomes = ensemblDAO.getChromosomes(MapManager.getInstance().getReferenceAssembly(speciesTypeKey).getKey());

        if(speciesTypeKey == SpeciesType.RAT) {
            mapKey = 361;
        }
        else if(speciesTypeKey == SpeciesType.MOUSE) {
            mapKey =  39;
        }
        else if(speciesTypeKey == SpeciesType.HUMAN){
            mapKey = 40;
        }


        for (EnsemblGene gene : genes) {
            if(chromosomes.contains(gene.getChromosome())) {
                String ncbiRgdId = null;
                String rgdId = null;

                List<String> ensembleRgdIds = ensemblDAO.getRgd_id(gene.getEnsemblGeneId(), 20);
                if(!gene.getEntrezgene_id().equals("NIL")) {
                    List<String> ncbiIds = ensemblDAO.getRgd_id(gene.getEntrezgene_id(), 3);
                    if (ncbiIds != null && ncbiIds.size() == 1)
                        ncbiRgdId = ncbiIds.get(0);
                    else {
                        //This indicates multiple rgdIds for a ncbi
                        if (ncbiIds != null) {
                            statuslog.info("Check the ncbi Id: " + gene.getEntrezgene_id());
                        } else ncbiRgdId = null; //Ncbi Id doesnt exist in rgd database
                    }
                }
                // Get rgdId based on the id from the file. RgdiD for Rat, MGI Id for Mouse and HGNC id for human are the ids from file.
                if (!gene.getrgdid().equals("0")) {
                    if (speciesTypeKey == SpeciesType.RAT)
                        rgdId = gene.getrgdid();
                    else if (speciesTypeKey == SpeciesType.MOUSE) {
                        List<String> rgdIds = ensemblDAO.getRgd_id(gene.getrgdid(), 5);
                        if(rgdIds != null && rgdIds.size() == 1)
                            rgdId = rgdIds.get(0);
                        else {
                            if(rgdId != null)
                                statuslog.info("Check these ids: Multiple rgdids for MGI Id -" +gene.getrgdid());
                            rgdId = null;
                        }
                    }
                    else {
                        List<String> rgdIds = ensemblDAO.getRgd_id(gene.getrgdid(), 21);
                        if(rgdIds != null && rgdIds.size() == 1)
                            rgdId = rgdIds.get(0);
                        else {
                            if(rgdId != null)
                                statuslog.info("Check these ids: Multiple rgdids for HGNC Id -" +gene.getrgdid());
                            rgdId = null;
                        }
                    }
                }

                // Ncbi id missing in the db
                // check for rgdId in file and rgdId in database for ensemble gene id
                if (ncbiRgdId == null) {
                    if (rgdId == null) {
                        // Case 1: No rgdid and no ncbi id in the file.
                        if (ensembleRgdIds != null && ensembleRgdIds.size() > 1)
                            statuslog.info("Check this out: Multiple RgdIds for EnsembleGene Id: " + gene.getEnsemblGeneId());
                        else {
                            if (ensembleRgdIds == null || ensembleRgdIds.isEmpty())
                                createNewEnsemblGene(gene, mapKey,null);
                            else
                                updateData(gene, ensembleRgdIds.get(0), mapKey);

                        }
                    } else {

                        // Case 2: No ncbi id in the file
                        if (ensembleRgdIds != null && ensembleRgdIds.contains(rgdId)) {
                            updateData(gene, rgdId, mapKey);
                        } else {
                            // Ignore the duplicate entries which ensemble sends in the file with wrong ncbi ids and rgd ids
                            if (matches.containsKey(gene.getEnsemblGeneId()) )
                                continue;
                            else {
                                    if(ensembleRgdIds != null) {
                                        mismatches.add(gene.getEnsemblGeneId());
                                        statuslog.info(" Ensemble Rgd ID and RgdId in file mismatch: " + gene.getEnsemblGeneId());
                                    } else {
                                        createNewEnsemblGene(gene, mapKey,rgdId);
                                    }
                                }
                            }
                        }
                    }
                 else {

                    // Ncbi Rgd Id and ensemble RgdId matches
                    if (ensembleRgdIds != null && ensembleRgdIds.contains(ncbiRgdId)) {
                        updateData(gene, ncbiRgdId, mapKey);
                        matches.put(gene.getEnsemblGeneId(), ncbiRgdId);
                    } else {
                        // Check if ncbi rgdId and rgdId from file matches
                        if (ensembleRgdIds == null && rgdId == null) {
                            createNewEnsemblGene(gene, mapKey,ncbiRgdId);
                        } else {
                            if (rgdId == null) {
                                if (matches.containsKey(gene.getEnsemblGeneId()))
                                    continue;
                                else {
                                    mismatches.add(gene.getEnsemblGeneId());
                                    statuslog.info(" Ensemble Rgd ID and Ncbi RgdId in db mismatch: " + gene.getEnsemblGeneId());
                                }
                            } else if (rgdId.equals(ncbiRgdId)) {
                                updateData(gene, ncbiRgdId, mapKey);
                                matches.put(gene.getEnsemblGeneId(), rgdId);
                            }
                        }
                    }

                }
            }
        }

        statuslog.info("Total mismatches: "+mismatches.size());
        statuslog.info("Total matches: "+ matches.size());
        statuslog.info("Total new Genes: "+ newGenes.size());
        statuslog.info("Total genes in file: "+genes.size());
    }
   public void aliasesinsert(int rgdid_new, EnsemblGene gene) throws Exception {
        List<String> gene_name = ensemblDAO.getGeneName(rgdid_new);
        String gene_name_str = ensemblDAO.getGeneStringName(gene_name);
        Alias aliasData = new Alias();
        aliasData.setNotes("Added by Ensembl pipeline");


        if (!ensemblDAO.getGeneSymbol(rgdid_new).contains(gene.getgene_name())) {
            aliasData.setRgdId(rgdid_new);
            aliasData.setValue(gene.getgene_name());
            aliasData.setTypeName("ensembl_gene_symbol");
            ensemblDAO.insertAlias(aliasData);

        }
        if (gene.getgene_description() != null && !gene.getgene_description().isEmpty() && !gene.getgene_description().contains(gene_name_str)) {
            aliasData.setRgdId(rgdid_new);
            aliasData.setValue(gene.getgene_description());
            aliasData.setTypeName("ensembl_full_name");
            ensemblDAO.insertAlias(aliasData);

        }

        ensemblDAO.insertGeneType(rgdid_new, gene.getgene_biotype());
        ensemblDAO.insertGeneSymbol(rgdid_new, gene.getgene_name());
        ensemblDAO.insertGeneName(rgdid_new, gene.getgene_description());
    }

   public void updateData(EnsemblGene gene,String rgdId, int mapKey) throws Exception{

       MapData mapData = new MapData();
       mapData.setSrcPipeline("Ensembl");
       mapData.setChromosome(gene.getChromosome());
       mapData.setMapKey(mapKey);
       mapData.setRgdId(Integer.parseInt(rgdId));
       mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
       mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
       mapData.setStrand(gene.getStrand());


       if (ensemblDAO.checkrecord(Integer.parseInt(rgdId), gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome(),mapKey) == null) {
           mapDAO.insertMapData(mapData);
           aliasesinsert(Integer.parseInt(rgdId), gene);
       }

       if(!ensemblDAO.checkXDBRecord(Integer.parseInt(rgdId),gene.getEnsemblGeneId(),"Ensembl")) {
           XdbId xdbId = new XdbId();
           xdbId.setRgdId(Integer.parseInt(rgdId));
           xdbId.setSrcPipeline("Ensembl");
           xdbId.setAccId(gene.getEnsemblGeneId());
           xdbId.setXdbKey(XdbId.XDB_KEY_ENSEMBL_GENES);
           ensemblDAO.insertXdbIds(xdbId);
       }
   }

    public void createNewEnsemblGene(EnsemblGene gene, int mapKey,String rgdId) throws Exception {


        int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();
        if (ensemblDAO.checkrecord_rgdid(gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome(),mapKey) == null) {
            String geneTypeLc = gene.getgene_biotype().toLowerCase();
            if (!ensemblDAO.existsGeneType(geneTypeLc))
                ensemblDAO.createGeneType(geneTypeLc);
            if(rgdId == null) {
                newGenes.add(gene);
                RgdId newRgdId = ensemblDAO.createRgdId(RgdId.OBJECT_KEY_GENES, speciesTypeKey);
                rgdId = String.valueOf(newRgdId.getRgdId());
                Gene newGene = new Gene();
                newGene.setSymbol(gene.getgene_name());
                newGene.setRgdId(Integer.parseInt(rgdId));
                newGene.setType(geneTypeLc);
                newGene.setName(gene.getgene_description());
                newGene.setGeneSource("Ensembl");
                newGene.setEnsemblFullName(gene.getgene_description());
                newGene.setEnsemblGeneSymbol(gene.getgene_name());
                newGene.setEnsemblGeneType(gene.getgene_biotype());
                ensemblDAO.insertGene(newGene);
            } else aliasesinsert(Integer.parseInt(rgdId), gene);


            MapData mapData = new MapData();
            mapData.setSrcPipeline("Ensembl");
            mapData.setChromosome(gene.getChromosome());
            mapData.setMapKey(mapKey);
            mapData.setRgdId(Integer.parseInt(rgdId));
            mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
            mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
            mapData.setStrand(gene.getStrand());
            mapDAO.insertMapData(mapData);


            XdbId xdbId = new XdbId();
            xdbId.setRgdId(Integer.parseInt(rgdId));
            xdbId.setSrcPipeline("Ensembl");
            xdbId.setAccId(gene.getEnsemblGeneId());
            xdbId.setXdbKey(XdbId.XDB_KEY_ENSEMBL_GENES);
            ensemblDAO.insertXdbIds(xdbId);
        }
    }
}
