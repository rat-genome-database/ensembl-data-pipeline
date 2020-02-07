package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.impl.GeneDAO;
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
    List nomenEvents = new ArrayList<>();
    public EnsemblGeneLoader() throws Exception {
        ensemblDAO = new EnsemblDAO();
    }
    public void run(Collection<EnsemblGene> genes,int speciesTypeKey) throws Exception
    {
        System.out.println("Loading the file");
        int mapKey = 0;
        edu.mcw.rgd.datamodel.Map reference = MapManager.getInstance().getReferenceAssembly(speciesTypeKey);
        List<String> chromosomes = ensemblDAO.getChromosomes(reference.getKey());

        if(speciesTypeKey == SpeciesType.MOUSE) {
            mapKey =  reference.getKey() + 4;
        }
        else if(speciesTypeKey == SpeciesType.HUMAN){
            mapKey = reference.getKey() + 2;
        }else mapKey = reference.getKey() + 1;

        for (EnsemblGene gene : genes) {
            if(chromosomes.contains(gene.getChromosome())) {
                String ncbiRgdId = null;
                String rgdId = null;

                List<String> ensembleRgdIds = ensemblDAO.getRgd_id(gene.getEnsemblGeneId(), XdbId.XDB_KEY_ENSEMBL_GENES);
                if(!gene.getEntrezgene_id().equals("NIL")) {
                    List<String> ncbiIds = ensemblDAO.getNcbiRgdId(gene.getEntrezgene_id());
                    if (ncbiIds != null && ncbiIds.size() == 1)
                        ncbiRgdId = ncbiIds.get(0);
                    else {
                        //This indicates multiple rgdIds for a ncbi
                        if (ncbiIds != null) {
                            statuslog.info("Check the ncbi Id: " + gene.getEntrezgene_id() +"\n");
                        } else ncbiRgdId = null; //Ncbi Id doesnt exist in rgd database
                    }
                }
                // Get rgdId based on the id from the file. RgdiD for Rat, MGI Id for Mouse and HGNC id for human are the ids from file.
                if (!gene.getrgdid().equals("0")) {
                    if (speciesTypeKey == SpeciesType.RAT)
                        rgdId = gene.getrgdid();
                    else if (speciesTypeKey == SpeciesType.MOUSE) {
                        List<String> rgdIds = ensemblDAO.getRgd_id(gene.getrgdid(), XdbId.XDB_KEY_MGD);
                        if(rgdIds != null && rgdIds.size() == 1)
                            rgdId = rgdIds.get(0);
                        else {
                            if(rgdId != null)
                                statuslog.info("Check these ids: Multiple rgdids for MGI Id -" +gene.getrgdid());
                            rgdId = null;
                        }
                    }
                    else {
                        List<String> rgdIds = ensemblDAO.getRgd_id(gene.getrgdid(), XdbId.XDB_KEY_HGNC);
                        if(rgdIds != null && rgdIds.size() == 1)
                            rgdId = rgdIds.get(0);
                        else {
                            if(rgdId != null)
                                statuslog.info("Check these ids: Multiple rgdids for HGNC Id -" +gene.getrgdid()+"\n");
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
                            statuslog.info("Check this out: Multiple RgdIds for EnsembleGene Id: " + gene.getEnsemblGeneId()+"\n");
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
                                        statuslog.info(" Ensemble Rgd ID and RgdId in file mismatch: " + gene.getEnsemblGeneId()+"\n");
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
                                    statuslog.info(" Ensemble Rgd ID and Ncbi RgdId in db mismatch: " + gene.getEnsemblGeneId()+"\n");
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

        statuslog.info("Total mismatches: "+mismatches.size()+"\n");
        statuslog.info("Total matches: "+ matches.size()+"\n");
        statuslog.info("Total new Genes: "+ newGenes.size()+"\n");
        statuslog.info("Total genes in file: "+genes.size()+"\n");
        statuslog.info("Total nomenEvents in file: "+nomenEvents.size()+"\n");
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

       if (ensemblDAO.checkrecord(Integer.parseInt(rgdId), gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome(),mapKey) == null) {
           MapData mapData = new MapData();
           mapData.setSrcPipeline("Ensembl");
           mapData.setChromosome(gene.getChromosome());
           mapData.setMapKey(mapKey);
           mapData.setRgdId(Integer.parseInt(rgdId));
           mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
           mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
           mapData.setStrand(gene.getStrand());
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
       } else {
            updateNomenEvents(gene,Integer.parseInt(rgdId));
       }
   }

    public void updateNomenEvents(EnsemblGene gene,int rgdId) throws  Exception{

        Gene existing = ensemblDAO.getGene(rgdId);
        if(existing.getGeneSource().equals("Ensembl") && (existing.getNomenSource() == null || !existing.getNomenSource().equals("HGNC"))){

            if(!existing.getSymbol().equalsIgnoreCase(gene.getgene_name()) || (existing.getName() != null && !existing.getName().equalsIgnoreCase(gene.getgene_description()) )){

                NomenclatureEvent event = new NomenclatureEvent();
                event.setRgdId(rgdId);
                event.setSymbol(gene.getgene_name());
                event.setName(gene.getgene_description());
                event.setRefKey("133850");
                event.setNomenStatusType("PROVISIONAL");
                event.setDesc("Symbol and/or name change");
                event.setEventDate(new Date());
                event.setOriginalRGDId(rgdId);
                event.setPreviousName(existing.getName());
                event.setPreviousSymbol(existing.getSymbol());
                ensemblDAO.insertNomenclatureEvent(event);

                Alias aliasData = new Alias();
                aliasData.setNotes("Added by Ensembl pipeline");
                aliasData.setRgdId(rgdId);
                if(!existing.getSymbol().equalsIgnoreCase(gene.getgene_name())){
                    aliasData.setValue(existing.getSymbol());
                    aliasData.setTypeName("old_gene_symbol");
                    ensemblDAO.insertAlias(aliasData);
                }
                if(existing.getName() != null && !existing.getName().equalsIgnoreCase(gene.getgene_description())) {
                    aliasData.setValue(existing.getName());
                    aliasData.setTypeName("old_gene_name");
                    ensemblDAO.insertAlias(aliasData);
                }

                existing.setSymbol(gene.getgene_name());
                existing.setName(gene.getgene_description());
                existing.setEnsemblGeneSymbol(gene.getgene_name());
                existing.setEnsemblFullName(gene.getgene_description());
                existing.setNomenSource("Ensembl");
                ensemblDAO.updateGene(existing);
                nomenEvents.add(rgdId);
            }

        } else {
            if(!existing.getEnsemblGeneSymbol().equalsIgnoreCase(gene.getgene_name()) || (existing.getEnsemblFullName() != null && !existing.getEnsemblFullName().equalsIgnoreCase(gene.getgene_description()) )){
                existing.setEnsemblGeneSymbol(gene.getgene_name());
                existing.setEnsemblFullName(gene.getgene_description());
                ensemblDAO.updateGene(existing);
            }
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
                newGene.setNomenSource("Ensembl");
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
