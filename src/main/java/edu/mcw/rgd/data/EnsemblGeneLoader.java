package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.Map;

/**
 * Created by sellanki on 10/18/2019.
 */
public class EnsemblGeneLoader {
    EnsemblDAO ensemblDAO = new EnsemblDAO();
    Logger statuslog = LogManager.getLogger("status");
    Logger conflictLog = LogManager.getLogger("conflicts");

    Map<String,String> matches;
    List mismatches;
    List<EnsemblGene> newGenes;
    List nomenEvents;
    AssemblyPositions genePositions;

    public void run(Collection<EnsemblGene> genes, int speciesTypeKey, int ensemblMapKey, int ncbiAssemblyMapKey) throws Exception {

        matches=new HashMap<>();
        mismatches=new ArrayList();
        newGenes = new ArrayList<>();
        nomenEvents = new ArrayList<>();
        genePositions = new AssemblyPositions();

        genePositions.loadPositionsInRgd(ensemblMapKey, ensemblDAO);

        // we have chromosome data only for NCBI assemblies
        List<Chromosome> chromosomes = ensemblDAO.getChromosomes(ncbiAssemblyMapKey);

        int skippedGenes = 0;
        int genesNoSymbolSkipped = 0;

        for (EnsemblGene gene : genes) {
            String chr = ensemblDAO.matchChromosome(gene.getChromosome(), chromosomes);
            if (chr == null) {
                skippedGenes++;
                statuslog.debug("gene_skipped: unexpected chromosome " + gene.getChromosome());
                continue;
            }
            gene.setChromosome(chr); // replace chr GeneBank id with NCBI scaffold acc

            // make sure gene name does not end with whitespace
            if( gene.getgene_description()!=null ) {
                gene.setgene_description(gene.getgene_description().trim());
            }

            String ncbiRgdId = null;

            List<String> ensembleRgdIds = ensemblDAO.getRgd_id(gene.getEnsemblGeneId(), XdbId.XDB_KEY_ENSEMBL_GENES);

            // some new genes loaded from rapid Release files do not have gene symbols;
            // we will use gene_ensembl_id as the gene symbol
            if( Utils.isStringEmpty(gene.getGeneSymbol()) ) {
                if( gene.getRefseqAccIds()==null ) {
                    //System.out.println("problem");
                }
            }

            String ncbiId = Utils.NVL(gene.getEntrezgene_id(), "NIL");
            if( !ncbiId.equals("NIL") ) {
                List<String> ncbiIds = ensemblDAO.getNcbiRgdId(ncbiId);
                if (ncbiIds != null && ncbiIds.size() == 1)
                    ncbiRgdId = ncbiIds.get(0);
                else {
                    //This indicates multiple rgdIds for a ncbi
                    if (ncbiIds != null) {
                        conflictLog.info("Check the ncbi Id: " + ncbiId);
                    } else ncbiRgdId = null; //Ncbi Id doesnt exist in rgd database
                }
            }

            // RgdId for Rat, MGI Id for Mouse and HGNC id for human are the ids from file.
            String accId = null;
            String incomingAcc = Utils.NVL(gene.getrgdid(), "0");
            if( !incomingAcc.equals("0")) {
                if (speciesTypeKey == SpeciesType.RAT)
                    accId = incomingAcc;
                else if (speciesTypeKey == SpeciesType.MOUSE) {
                    List<String> rgdIds = ensemblDAO.getRgd_id(incomingAcc, XdbId.XDB_KEY_MGD);
                    if(rgdIds != null && rgdIds.size() == 1)
                        accId = rgdIds.get(0);
                    else {
                        if(accId != null)
                            conflictLog.info("Check these ids: Multiple rgdids for MGI Id -" +incomingAcc);
                        accId = null;
                    }
                }
                else {
                    List<String> rgdIds = ensemblDAO.getRgd_id(incomingAcc, XdbId.XDB_KEY_HGNC);
                    if(rgdIds != null && rgdIds.size() == 1)
                        accId = rgdIds.get(0);
                    else {
                        if(accId != null)
                            conflictLog.info("Check these ids: Multiple rgdids for HGNC Id -" +incomingAcc);
                        accId = null;
                    }
                }
            }

            // Ncbi id missing in the db
            // check for rgdId in file and rgdId in database for ensemble gene id
            if (ncbiRgdId == null) {
                if (accId == null) {
                    // Case 1: No rgdid and no ncbi id in the file.
                    if (ensembleRgdIds != null && ensembleRgdIds.size() > 1)
                        conflictLog.info("  check this out: multiple RgdIds for EnsembleGene Id: " + gene.getEnsemblGeneId());
                    else {
                        if (ensembleRgdIds == null || ensembleRgdIds.isEmpty()) {
                            if (!createNewEnsemblGene(gene, ensemblMapKey, null, speciesTypeKey)) {
                                genesNoSymbolSkipped++;
                            }
                        } else {
                            updateData(gene, ensembleRgdIds.get(0), ensemblMapKey);
                        }

                    }
                } else {

                    // Case 2: No ncbi id in the file
                    if (ensembleRgdIds != null && ensembleRgdIds.contains(accId)) {
                        updateData(gene, accId, ensemblMapKey);
                    } else {
                        // Ignore the duplicate entries which ensemble sends in the file with wrong ncbi ids and rgd ids
                        if (matches.containsKey(gene.getEnsemblGeneId()) )
                            continue;
                        else {
                                if(ensembleRgdIds != null) {
                                    mismatches.add(gene.getEnsemblGeneId());
                                    conflictLog.info("Ensemble Rgd ID and RgdId in file mismatch: " + gene.getEnsemblGeneId());
                                } else {
                                    if( !createNewEnsemblGene(gene, ensemblMapKey, accId, speciesTypeKey) ) {
                                        genesNoSymbolSkipped++;
                                    }
                                }
                            }
                        }
                    }
                }
            else {

                // Ncbi Rgd Id and ensemble RgdId matches
                if (ensembleRgdIds != null && ensembleRgdIds.contains(ncbiRgdId)) {
                    updateData(gene, ncbiRgdId, ensemblMapKey);
                    matches.put(gene.getEnsemblGeneId(), ncbiRgdId);
                } else {
                    // Check if ncbi rgdId and rgdId from file matches
                    if (ensembleRgdIds == null && accId == null) {
                        if( !createNewEnsemblGene(gene, ensemblMapKey, ncbiRgdId, speciesTypeKey) ) {
                            genesNoSymbolSkipped++;
                        }
                    } else {
                        if (accId == null) {
                            if (matches.containsKey(gene.getEnsemblGeneId()))
                                continue;
                            else {
                                mismatches.add(gene.getEnsemblGeneId());
                                conflictLog.info("Ensemble Rgd ID and Ncbi RgdId in db mismatch: " + gene.getEnsemblGeneId());
                            }
                        } else if (accId.equals(ncbiRgdId)) {
                            updateData(gene, ncbiRgdId, ensemblMapKey);
                            matches.put(gene.getEnsemblGeneId(), accId);
                        }
                    }
                }

            }
        }

        statuslog.info("Total mismatches: "+mismatches.size());
        statuslog.info("Total matches: "+ matches.size());
        if( newGenes.size()>0 ) {
            statuslog.info("Total new Genes: " + newGenes.size());
        }
        statuslog.info("Total genes in file: "+genes.size());
        if( nomenEvents.size()>0 ) {
            statuslog.info("Total nomenEvents in file: " + nomenEvents.size());
        }
        if( skippedGenes>0 ) {
            statuslog.info("Genes skipped, unexpected chromosome: " + skippedGenes);
        }
        if( genesNoSymbolSkipped>0 ) {
            statuslog.info("Genes skipped, no symbol, cannot match with NCBI gene: " + genesNoSymbolSkipped);
        }

        genePositions.qcAndLoad(statuslog, ensemblDAO);
    }

    public void aliasesinsert(int newRgdId, EnsemblGene gene) throws Exception {

        Gene geneInRgd = ensemblDAO.getGene(newRgdId);
        String geneName = geneInRgd.getName();

        String geneSymbolIncoming = Utils.NVL(gene.getGeneSymbol(), gene.getEnsemblGeneId());

        Alias aliasData = new Alias();
        aliasData.setNotes("Added by Ensembl pipeline");

        if( !geneInRgd.getSymbol().equals(geneSymbolIncoming) ) {
            aliasData.setRgdId(newRgdId);
            aliasData.setValue(geneSymbolIncoming);
            aliasData.setTypeName("ensembl_gene_symbol");
            ensemblDAO.insertAlias(aliasData);
        }

        if (gene.getgene_description() != null && !gene.getgene_description().isEmpty() && !gene.getgene_description().contains(geneName)) {
            aliasData.setRgdId(newRgdId);
            aliasData.setValue(gene.getgene_description());
            aliasData.setTypeName("ensembl_full_name");
            ensemblDAO.insertAlias(aliasData);
        }

        geneInRgd.setEnsemblFullName(gene.getgene_description());
        geneInRgd.setEnsemblGeneSymbol(geneSymbolIncoming);
        geneInRgd.setEnsemblGeneType(gene.getgene_biotype());
        ensemblDAO.updateGene(geneInRgd);
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
       genePositions.addIncomingPosition(mapData);

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

            String geneSymbolIncoming = Utils.NVL(gene.getGeneSymbol(), gene.getEnsemblGeneId());

            if( !existing.getSymbol().equalsIgnoreCase(geneSymbolIncoming) || (existing.getName() != null && !existing.getName().equalsIgnoreCase(gene.getgene_description()) )){

                NomenclatureEvent event = new NomenclatureEvent();
                event.setRgdId(rgdId);
                event.setSymbol(geneSymbolIncoming);
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
                if(!existing.getSymbol().equalsIgnoreCase(geneSymbolIncoming)){
                    aliasData.setValue(existing.getSymbol());
                    aliasData.setTypeName("old_gene_symbol");
                    ensemblDAO.insertAlias(aliasData);
                }
                if(existing.getName() != null && !existing.getName().equalsIgnoreCase(gene.getgene_description())) {
                    aliasData.setValue(existing.getName());
                    aliasData.setTypeName("old_gene_name");
                    ensemblDAO.insertAlias(aliasData);
                }

                existing.setSymbol(geneSymbolIncoming);
                existing.setName(gene.getgene_description());
                existing.setEnsemblGeneSymbol(geneSymbolIncoming);
                existing.setEnsemblFullName(gene.getgene_description());
                existing.setNomenSource("Ensembl");
                ensemblDAO.updateGene(existing);
                nomenEvents.add(rgdId);
            }

        } else {
            if( !Utils.stringsAreEqualIgnoreCase(existing.getEnsemblGeneSymbol(), gene.getGeneSymbol())
             || !Utils.stringsAreEqualIgnoreCase(existing.getEnsemblFullName(), gene.getgene_description()) ){

                existing.setEnsemblGeneSymbol(gene.getGeneSymbol());
                existing.setEnsemblFullName(gene.getgene_description());
                ensemblDAO.updateGene(existing);
            }
        }
    }

    public boolean createNewEnsemblGene(EnsemblGene gene, int mapKey, String rgdId, int speciesTypeKey) throws Exception {

        if (ensemblDAO.checkrecord_rgdid(gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome(),mapKey) == null) {

            String geneSymbol = gene.getGeneSymbol();
            if( false ) {
                // gene symbol must be non-null
                if (Utils.isStringEmpty(gene.getGeneSymbol())) {
                    conflictLog.info("gene " + gene.getEnsemblGeneId() + " does not have a symbol! gene skipped -- not inserted");
                    return false;
                }
            } else {
                if( Utils.isStringEmpty(geneSymbol) ) {
                    geneSymbol = gene.getEnsemblGeneId();
                }
            }

            String geneTypeLc = gene.getgene_biotype().toLowerCase();
            if (!ensemblDAO.existsGeneType(geneTypeLc)) {
                ensemblDAO.createGeneType(geneTypeLc);
            }

            if(rgdId == null) {

                newGenes.add(gene);
                RgdId newRgdId = ensemblDAO.createRgdId(RgdId.OBJECT_KEY_GENES, speciesTypeKey);
                rgdId = String.valueOf(newRgdId.getRgdId());
                Gene newGene = new Gene();
                newGene.setSymbol(geneSymbol);
                newGene.setRgdId(Integer.parseInt(rgdId));
                newGene.setType(geneTypeLc);
                newGene.setName(gene.getgene_description());
                newGene.setGeneSource("Ensembl");
                newGene.setNomenSource("Ensembl");
                newGene.setEnsemblFullName(gene.getgene_description());
                newGene.setEnsemblGeneSymbol(geneSymbol);
                newGene.setEnsemblGeneType(gene.getgene_biotype());
                ensemblDAO.insertGene(newGene);

                // always create PROVISIONAL nomenclature event for newly created rat gene
                if( speciesTypeKey==SpeciesType.RAT ) {
                    NomenclatureEvent event = new NomenclatureEvent();
                    event.setRgdId(newRgdId.getRgdId());
                    event.setSymbol(newGene.getSymbol());
                    event.setName(newGene.getName());
                    event.setRefKey("20683");
                    event.setNomenStatusType("PROVISIONAL");
                    event.setDesc("Symbol and Name status set to provisional");
                    event.setEventDate(new Date());
                    event.setOriginalRGDId(newRgdId.getRgdId());
                    ensemblDAO.insertNomenclatureEvent(event);
                }
            } else {
                aliasesinsert(Integer.parseInt(rgdId), gene);
            }


            MapData mapData = new MapData();
            mapData.setSrcPipeline("Ensembl");
            mapData.setChromosome(gene.getChromosome());
            mapData.setMapKey(mapKey);
            mapData.setRgdId(Integer.parseInt(rgdId));
            mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
            mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
            mapData.setStrand(gene.getStrand());
            ensemblDAO.insertMapData(mapData, "GENE ");


            XdbId xdbId = new XdbId();
            xdbId.setRgdId(Integer.parseInt(rgdId));
            xdbId.setSrcPipeline("Ensembl");
            xdbId.setAccId(gene.getEnsemblGeneId());
            xdbId.setXdbKey(XdbId.XDB_KEY_ENSEMBL_GENES);
            ensemblDAO.insertXdbIds(xdbId);
        }
        return true;
    }
}
