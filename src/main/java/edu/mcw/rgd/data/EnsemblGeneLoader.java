package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.dao.impl.RGDManagementDAO;
import edu.mcw.rgd.datamodel.Alias;
import edu.mcw.rgd.datamodel.MapData;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by sellanki on 10/18/2019.
 */
public class EnsemblGeneLoader {
    EnsemblDAO ensemblDAO;
    static Logger statuslog = Logger.getLogger("statuscheck");
    static Logger mismatchlog = Logger.getLogger("mismatchcheck");
    static Logger externalidexsistnorgdid = Logger.getLogger("externalidexsistnorgdid");
    static Logger ncbiexistbutnorgdid = Logger.getLogger("ncbiexistbutnorgdid");
    static Logger test4log = Logger.getLogger("test4");
    static Logger test1log = Logger.getLogger("test1");
    static Logger test2log = Logger.getLogger("test2");
    static Logger test3log = Logger.getLogger("test3");
    static Logger test5log = Logger.getLogger("test5");
    static Logger matchedrecords = Logger.getLogger("matchedrecords");
    RGDManagementDAO managementDAO = new RGDManagementDAO();
    MapDAO mapDAO = new MapDAO();
    List totalgenesloaded=new ArrayList();
    List totalmismatches=new ArrayList();
    List external_idexistnorgd_id=new ArrayList();
    Set toalaliasesinserted = new HashSet();
    List newgenerecords = new ArrayList();
    List ncbi_idexistnorgd_id=new ArrayList();

    public EnsemblGeneLoader() throws Exception {
        ensemblDAO = new EnsemblDAO();
    }

    public void run(Collection<EnsemblGene> genes) throws Exception
    {
        for (EnsemblGene gene : genes) {
            int rgdid_new = Integer.parseInt(gene.getrgdid());
            MapData mapData = new MapData();
            mapData.setSrcPipeline("Ensembl");
            mapData.setChromosome(gene.getChromosome());
            mapData.setMapKey(361);
            mapData.setRgdId(rgdid_new);
            mapData.setStartPos(Integer.parseInt(gene.getStartPos()));
            mapData.setStopPos(Integer.parseInt(gene.getStopPos()));
            mapData.setStrand(gene.getStrand());
            Alias aliasData = new Alias();
            aliasData.setNotes("Added by Ensembl pipeline");
            List<String> gene_name = ensemblDAO.getGeneName(rgdid_new);
            String gene_name_str = ensemblDAO.getGeneStringName(gene_name);

            if (gene.getrgdid() != "0")
                try {

                    if (managementDAO.getRgdId2(rgdid_new) != null) {
                        if (!ensemblDAO.checkobjectstatus(rgdid_new).contains("ACTIVE")) {
                            statuslog.info(gene.getrgdid() + "\t" + gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + ensemblDAO.checkobjectstatus(rgdid_new));
                        }
                        if (!gene.getEntrezgene_id().contains("NIL"))
                        {
                            if (!(ensemblDAO.getXdbIds(3, rgdid_new).isEmpty()) && !(ensemblDAO.getXdbIds(20, rgdid_new).isEmpty()))
                            {
                                if (ensemblDAO.getXdbIds(3, rgdid_new).contains(gene.getEntrezgene_id()) && ensemblDAO.getXdbIds(20, rgdid_new).contains(gene.getEnsemblGeneId()))
                                {
                                if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                    {
                                        matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                        mapDAO.insertMapData(mapData);
                                        totalgenesloaded.add(gene.getEnsemblGeneId());
                                        aliasesinsert(rgdid_new, gene);
                                    }
                                }
                                else {
                                    mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                    totalmismatches.add(gene.getEnsemblGeneId());
                                }
                            }
                            if (!(ensemblDAO.getXdbIds(3, rgdid_new).isEmpty()) && (ensemblDAO.getXdbIds(20, rgdid_new).isEmpty()))
                            {
                                if (ensemblDAO.getXdbIds(3, rgdid_new).contains(gene.getEntrezgene_id()))
                                    {
                                        if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                        {
                                        matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                        mapDAO.insertMapData(mapData);
                                        totalgenesloaded.add(gene.getEnsemblGeneId());
                                        aliasesinsert(rgdid_new, gene);
                                        }
                                }
                                else {
                                    mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                    totalmismatches.add(gene.getEnsemblGeneId());
                                }
                            }
                            if ((ensemblDAO.getXdbIds(3, rgdid_new).isEmpty()) && !(ensemblDAO.getXdbIds(20, rgdid_new).isEmpty()))
                            {

                                if (ensemblDAO.getXdbIds(20, rgdid_new).contains(gene.getEnsemblGeneId()))
                                    {
                                        if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                        {
                                        matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                        mapDAO.insertMapData(mapData);
                                        totalgenesloaded.add(gene.getEnsemblGeneId());
                                        aliasesinsert(rgdid_new, gene);
                                    }
                                }
                                else {
                                    mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                    totalmismatches.add(gene.getEnsemblGeneId());
                                }
                            }
                        }
                        if (gene.getEntrezgene_id().contains("NIL"))
                        {

                                if (ensemblDAO.getXdbIds(20, rgdid_new).contains(gene.getEnsemblGeneId()))
                                {
                                    if (ensemblDAO.checkrecord(rgdid_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                    {

                                    matchedrecords.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                    mapDAO.insertMapData(mapData);
                                    totalgenesloaded.add(gene.getEnsemblGeneId());
                                    aliasesinsert(rgdid_new, gene);
                                }

                            }
                                else {
                                    mismatchlog.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                                    totalmismatches.add(gene.getEnsemblGeneId());
                                }
                        }
                    }
                } catch (Exception e) {
                    System.out.println("RGD id is" + gene.getrgdid());
                    e.printStackTrace();
                }
            if (gene.getrgdid() == "0")
            {
                if ((ensemblDAO.checkXdbIds(gene.getEntrezgene_id()).contains(gene.getEntrezgene_id())) || (ensemblDAO.checkXdbIds(gene.getEnsemblGeneId()).contains(gene.getEnsemblGeneId())))
                {
                    if ((ensemblDAO.checkXdbIds(gene.getEntrezgene_id()).contains(gene.getEntrezgene_id())) && (ensemblDAO.checkXdbIds(gene.getEnsemblGeneId()).contains(gene.getEnsemblGeneId())))
                    {
                        MapData[] array = new MapData[ensemblDAO.checkrgd_id(gene.getEnsemblGeneId()).size()];
                        List<String> list3 = new ArrayList<String>();
                        list3 = ensemblDAO.checkrgd_id(gene.getEntrezgene_id()).stream()
                                .filter(ensemblDAO.checkrgd_id(gene.getEnsemblGeneId())::contains)
                                .collect(Collectors.toList());
                        if (!list3.isEmpty()) {
                            String rgd1_new = "";
                            for (String s : list3) {
                                rgd1_new += s;
                                int rgd_id_new = Integer.parseInt(rgd1_new);
                                mapData.setRgdId(rgd_id_new);
                                if (ensemblDAO.checkrecord(rgd_id_new, gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null)
                                {
                                    mapDAO.insertMapData(mapData);
                                    totalgenesloaded.add(gene.getEnsemblGeneId());
                                    if (!ensemblDAO.getGeneType(rgd_id_new).contains(gene.getgene_biotype()))
                                    {
                                        aliasData.setRgdId(rgd_id_new);
                                        aliasData.setValue(gene.getgene_biotype());
                                        aliasData.setTypeName("ensembl_gene_type");
                                        //ensemblDAO.insertAlias(aliasData);
                                        toalaliasesinserted.add(aliasData);
                                        test1log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand() + "\t" + gene.getgene_biotype() + "\t" + ensemblDAO.getGeneType(rgdid_new) + "\t" + gene.getrgdid());
                                    }
                                    if (!ensemblDAO.getGeneSymbol(rgd_id_new).contains(gene.getgene_name()))
                                    {
                                        aliasData.setRgdId(rgd_id_new);
                                        aliasData.setValue(gene.getgene_name());
                                        aliasData.setTypeName("ensembl_gene_symbol");
                                        ensemblDAO.insertAlias(aliasData);
                                        toalaliasesinserted.add(aliasData);
                                        test2log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand() + "\t" + gene.getgene_name() + "\t" + ensemblDAO.getGeneSymbol(rgdid_new) + "\t" + gene.getrgdid());
                                    }
                                    if (!gene.getgene_description().contains(gene_name_str)) {
                                        aliasData.setRgdId(rgd_id_new);
                                        aliasData.setValue(gene.getgene_description());
                                        aliasData.setTypeName("ensembl_full_name");
                                        ensemblDAO.insertAlias(aliasData);
                                        toalaliasesinserted.add(aliasData);
                                        test3log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand() + "\t" + gene.getgene_description() + "\t" + ensemblDAO.getGeneName(rgdid_new) + "\t" + gene.getrgdid());
                                    }
                                    ensemblDAO.insertGeneType(rgd_id_new, gene.getgene_biotype());
                                    ensemblDAO.insertGeneSymbol(rgd_id_new, gene.getgene_name());
                                    ensemblDAO.insertGeneName(rgd_id_new, gene.getgene_description());
                                    test5log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand() + "\t" + gene.getgene_description() + "\t" + rgd1_new);
                                }
                            }
                        } else {
                            externalidexsistnorgdid.info("Record with ncbi id and ensembl id exsists in rgd_acc_xdb table with, Ncbi id:" + gene.getEntrezgene_id() + ", Ensembl id:" + gene.getEnsemblGeneId() + " but rgd_id is not found in the file");
                            external_idexistnorgd_id.add(gene.getEnsemblGeneId());
                        }

                    } else if ((ensemblDAO.checkXdbIds(gene.getEntrezgene_id()).contains(gene.getEntrezgene_id()))) {
                        externalidexsistnorgdid.info("Record with ncbi id exsists in rgd_acc_xdb table with, Ncbi id:" + gene.getEntrezgene_id() + ", Ensembl id:" + gene.getEnsemblGeneId() + " but rgd_id is not found in the file");
                        external_idexistnorgd_id.add(gene.getEnsemblGeneId());
                    } else if ((ensemblDAO.checkXdbIds(gene.getEnsemblGeneId()).contains(gene.getEnsemblGeneId()))) {
                        externalidexsistnorgdid.info("Record with ensemblgene id exsists in rgd_acc_xdb table with, Ncbi id:" + gene.getEntrezgene_id() + ", Ensembl id:" + gene.getEnsemblGeneId() + " but rgd_id is not found in the file");
                        external_idexistnorgd_id.add(gene.getEnsemblGeneId());
                    }
                } else {
                    if (gene.getEntrezgene_id().equals("NIL")) {
                        if (ensemblDAO.checkrecord_rgdid(gene.getStartPos(), gene.getStopPos(), gene.getStrand(), gene.getChromosome()) == null) {
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
                                newgenerecords.add(gene.getEnsemblGeneId());
                                ensemblDAO.insertGeneSource(rgdid_new,"Ensembl");
                                ensemblDAO.insertGeneType(rgdid_new,gene.getgene_biotype());
                                ensemblDAO.insertGeneSymbol(rgdid_new,gene.getgene_name());
                                ensemblDAO.insertGeneName(rgdid_new,gene.getgene_description());
                                mapData.setRgdId(rgdid_new);
                                mapDAO.insertMapData(mapData);
                                totalgenesloaded.add(gene.getEnsemblGeneId());
                                XdbId xdbId = new XdbId();
                                xdbId.setRgdId(rgdid_new);
                                xdbId.setSrcPipeline("Ensembl");
                                xdbId.setAccId(gene.getEnsemblGeneId());
                                xdbId.setXdbKey(XdbId.XDB_KEY_ENSEMBL_GENES);
                                ensemblDAO.insertXdbIds(xdbId);*/
                        }
                    } else {
                        ncbiexistbutnorgdid.info("Record has an ncbi_id  but no rgdid:" + gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getrgdid() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand());
                        ncbi_idexistnorgd_id.add(gene.getEnsemblGeneId());
                    }
                }
            }
        }
        statuslog.info("Total loaded genes are \t"+totalgenesloaded.size());
        statuslog.info("\nTotal no of mismatched genes\t"+totalmismatches.size());
        statuslog.info("\nTotal no of genes external id exist and no rgdid\t"+external_idexistnorgd_id.size());
        statuslog.info("\nTotal number of Aliases inserted\t"+toalaliasesinserted.size());
        statuslog.info("\nTotal number of new gene records created\t"+newgenerecords.size());
        statuslog.info("\nTotal number of records in which ncbi exist and no rgdid\t"+ncbi_idexistnorgd_id.size());
    }
    public void aliasesinsert(int rgdid_new, EnsemblGene gene) throws Exception {
        List<String> gene_name = ensemblDAO.getGeneName(rgdid_new);
        String gene_name_str = ensemblDAO.getGeneStringName(gene_name);
        Alias aliasData = new Alias();
        aliasData.setNotes("Added by Ensembl pipeline");
        //System.out.println(gene.getgene_biotype());
        if (!ensemblDAO.getGeneType(rgdid_new).contains(gene.getgene_biotype())) {
            aliasData.setRgdId(rgdid_new);
            aliasData.setValue(gene.getgene_biotype());
            aliasData.setTypeName("ensembl_gene_type");
            ensemblDAO.insertAlias(aliasData);
            toalaliasesinserted.add(aliasData);
           // test1log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand() + "\t" + gene.getgene_biotype() + "\t" + ensemblDAO.getGeneType(rgdid_new) + "\t" + gene.getrgdid());
        }
        if (!ensemblDAO.getGeneSymbol(rgdid_new).contains(gene.getgene_name())) {
            aliasData.setRgdId(rgdid_new);
            aliasData.setValue(gene.getgene_name());
            aliasData.setTypeName("ensembl_gene_symbol");
            ensemblDAO.insertAlias(aliasData);
            toalaliasesinserted.add(aliasData);
           // test2log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand() + "\t" + gene.getgene_name() + "\t" + ensemblDAO.getGeneSymbol(rgdid_new) + "\t" + gene.getrgdid());
        }
        if (!gene.getgene_description().contains(gene_name_str)) {
            aliasData.setRgdId(rgdid_new);
            aliasData.setValue(gene.getgene_description());
            aliasData.setTypeName("ensembl_full_name");
            ensemblDAO.insertAlias(aliasData);
            toalaliasesinserted.add(aliasData);
           // test3log.info(gene.getEnsemblGeneId() + "\t" + gene.getEntrezgene_id() + "\t" + gene.getStartPos() + "\t" + gene.getStopPos() + "\t" + gene.getChromosome() + "\t" + gene.getStrand() + "\t" + gene.getgene_description() + "\t" + ensemblDAO.getGeneName(rgdid_new) + "\t" + gene.getrgdid());
        }
        ensemblDAO.insertGeneType(rgdid_new, gene.getgene_biotype());
        ensemblDAO.insertGeneSymbol(rgdid_new, gene.getgene_name());
        ensemblDAO.insertGeneName(rgdid_new, gene.getgene_description());
    }
}
