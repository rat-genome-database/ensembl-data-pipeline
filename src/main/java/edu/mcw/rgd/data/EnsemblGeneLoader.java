package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
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
    Logger genesUpdatedLog = LogManager.getLogger("genes_updated");

    Map<String,String> matches;
    List mismatches;
    List<EnsemblGene> newGenes;
    List nomenEvents;
    AssemblyPositions genePositions;
    CounterPool counters;

    public void run(Collection<EnsemblGene> genes, int speciesTypeKey, int ensemblMapKey, int ncbiAssemblyMapKey) throws Exception {

        matches=new HashMap<>();
        mismatches=new ArrayList();
        newGenes = new ArrayList<>();
        nomenEvents = new ArrayList<>();
        genePositions = new AssemblyPositions();
        counters = new CounterPool();

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
            if( gene.getGeneName()!=null ) {
                gene.setGeneName(gene.getGeneName().trim());
            }

            String ncbiRgdId = qcNcbiId(gene.getEntrezgene_id());

            List<String> ensembleRgdIds = ensemblDAO.getRgd_id(gene.getEnsemblGeneId(), XdbId.XDB_KEY_ENSEMBL_GENES);

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
                                    conflictLog.info("NO NCBI rgd ids; incoming " + gene.getEnsemblGeneId()+" "+gene.getGeneSymbol()+"  has  RGD:"+accId+" and Ensembl RGD IDS: "+Utils.concatenate(ensembleRgdIds, ","));
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
                                conflictLog.info("NCBI mismatch for " + gene.getEnsemblGeneId()+": NCBI RGD:"+ncbiRgdId+", acc="+accId+", ensembl rgd ids:"+Utils.concatenate(ensembleRgdIds,","));
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
            statuslog.info("  Total nomenEvents in file: " + nomenEvents.size());
        }
        int count = counters.get("ENSEMBL_GENE_SYMBOL_CHANGED");
        if( count>0 ){
            statuslog.info("  genes with Ensembl symbol change: " + count);
        }
        count = counters.get("ENSEMBL_GENE_NAME_CHANGED");
        if( count>0 ){
            statuslog.info("  genes with Ensembl name change: " + count);
        }
        count = counters.get("ENSEMBL_GENE_TYPE_CHANGED");
        if( count>0 ){
            statuslog.info("  genes with Ensembl biotype change: " + count);
        }

        if( genesNoSymbolSkipped>0 ) {
            statuslog.info("Genes skipped, no symbol, cannot match with NCBI gene: " + genesNoSymbolSkipped);
        }

        genePositions.qcAndLoad(statuslog, ensemblDAO);
    }

    String qcNcbiId(String ncbiGeneId) throws Exception {

        String ncbiRgdId = null;

        String ncbiId = Utils.NVL(ncbiGeneId, "NIL");
        if( !ncbiId.equals("NIL") ) {
            List<String> ncbiRgdIds = ensemblDAO.getNcbiRgdId(ncbiId);
            if( ncbiRgdIds.size() == 1 ) {
                ncbiRgdId = ncbiRgdIds.get(0);
            } else if( ncbiRgdIds.size()>1 ) {
                //This indicates multiple rgdIds for a ncbi
                conflictLog.info("MULTIs: NCBI Id:" + ncbiId+" resolves to RGD ids: "+Utils.concatenate(ncbiRgdIds,","));
            } else {
                //Ncbi Id doesnt exist in rgd database
                conflictLog.info("NCBI Id:" + ncbiId+" is inactive in RGD or is not present in RGD");
            }
        }
        return ncbiRgdId;
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

        if( gene.getGeneName() != null && !gene.getGeneName().contains(geneName) ) {
            aliasData.setRgdId(newRgdId);
            aliasData.setValue(gene.getGeneName());
            aliasData.setTypeName("ensembl_full_name");
            ensemblDAO.insertAlias(aliasData);
        }

        geneInRgd.setEnsemblFullName(gene.getGeneName());
        geneInRgd.setEnsemblGeneSymbol(geneSymbolIncoming);
        geneInRgd.setEnsemblGeneType(gene.getGeneBioType());
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

            if( !existing.getSymbol().equalsIgnoreCase(geneSymbolIncoming) || (existing.getName() != null && !existing.getName().equalsIgnoreCase(gene.getGeneName()) )){

                NomenclatureEvent event = new NomenclatureEvent();
                event.setRgdId(rgdId);
                event.setSymbol(geneSymbolIncoming);
                event.setName(gene.getGeneName());
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
                if(existing.getName() != null && !existing.getName().equalsIgnoreCase(gene.getGeneName())) {
                    aliasData.setValue(existing.getName());
                    aliasData.setTypeName("old_gene_name");
                    ensemblDAO.insertAlias(aliasData);
                }

                existing.setSymbol(geneSymbolIncoming);
                existing.setName(gene.getGeneName());
                existing.setEnsemblGeneSymbol(geneSymbolIncoming);
                existing.setEnsemblFullName(gene.getGeneName());
                existing.setEnsemblGeneType(gene.getGeneBioType());
                existing.setNomenSource("Ensembl");
                ensemblDAO.updateGene(existing);
                nomenEvents.add(rgdId);
            }

        } else {
            // non-Ensembl genes or Ensembl genes with HGNC authority (reviewed May 20, 2022)
            boolean symbolChanged = !Utils.stringsAreEqualIgnoreCase(existing.getEnsemblGeneSymbol(), gene.getGeneSymbol());
            boolean nameChanged = !Utils.stringsAreEqualIgnoreCase(existing.getEnsemblFullName(), gene.getGeneName());
            boolean typeChanged = !Utils.stringsAreEqualIgnoreCase(existing.getEnsemblGeneType(), gene.getGeneBioType());

            if( symbolChanged || nameChanged || typeChanged ) {

                if( symbolChanged ) {
                    counters.increment("ENSEMBL_GENE_SYMBOL_CHANGED");
                }
                if( nameChanged ) {
                    counters.increment("ENSEMBL_GENE_NAME_CHANGED");
                }
                if( typeChanged ) {
                    counters.increment("ENSEMBL_GENE_TYPE_CHANGED");
                }

                String msg = "RGD:" + existing.getRgdId() + "\n" +
                        "OLD: SYMBOL [" + existing.getEnsemblGeneSymbol() + "] NAME [" + existing.getEnsemblFullName() + "] TYPE [" + existing.getEnsemblGeneType() + "]\n" +
                        "NEW: SYMBOL [" + gene.getGeneSymbol() + "] NAME [" + gene.getGeneName() + "] TYPE [" + gene.getGeneBioType() + "]";
                genesUpdatedLog.debug(msg);

                existing.setEnsemblGeneSymbol(gene.getGeneSymbol());
                existing.setEnsemblFullName(gene.getGeneName());
                existing.setEnsemblGeneType(gene.getGeneBioType());
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

            String geneTypeLc = gene.getGeneBioType().toLowerCase();
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
                newGene.setName(gene.getGeneName());
                newGene.setGeneSource("Ensembl");
                newGene.setNomenSource("Ensembl");
                newGene.setEnsemblFullName(gene.getGeneName());
                newGene.setEnsemblGeneSymbol(geneSymbol);
                newGene.setEnsemblGeneType(gene.getGeneBioType());
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
