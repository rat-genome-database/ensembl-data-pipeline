package edu.mcw.rgd.data;


import edu.mcw.rgd.process.PipelineLogger;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.util.Collection;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */

public class EnsemblPipelinePreprocessor {

    EnsemblDataPuller dataPuller;
    EnsemblDAO ensemblDAO;
    //Parser parser;
    Databaseloading dbloading;
     Parser dataParser;
    LogLoader loader;
    int speciesTypeKey;
    PipelineLogger dbLogger = PipelineLogger.getInstance();
     EnsemblGene gene;
    private Connection conn;
    PreparedStatement preparedstatement = null;
    String ENSEMBL_GENE_ID;
    String CHROMOSE_NAME;
    String START_POSITION;
    String END_POSITION;
    String STRAND;
    String RGD_ID;
    String entrezgene_id;

    public Collection<EnsemblGene> run() throws Exception {
//public void run() throws Exception{
        // download genes data from Ensembl biomart and store it locally in data folder
       dataPuller.setSpeciesTypeKey(speciesTypeKey);
        String dataFile = dataPuller.downloadGenesFile();
       // String dataFile1="C:/git/ensembl-data/src/main/dist/logs/mismatchofexternalidwithrgdid.log";
        // parse the file; map key is 'ensembl gene id'
        //System.out.println(dataFile);
        //Map<String, EnsemblGene> genes = dataParser.parseGene(dataFile);
       List<EnsemblGene> genes=dataParser.parseGene(dataFile);
        //List<EnsemblGene> log_genes = loader.parseGene(dataFile1);
        //List<EnsemblGene> genes =parser.parseGene(dataFile);
        //traverseallgenes(genes);
        /*
        */
        //dbLogger.log(null, Integer.toString(log_genes.size()), PipelineLog.LOGPROP_RECCOUNT);
    //process all ensembl genes*/
         //return genes.values();
        return genes;
    }
    void traverseallgenes(Collection<EnsemblGene> genes) throws Exception
    {
        for(EnsemblGene gene:genes)
        {
            System.out.println(gene.getChromosome());
        }
    }
    public LogLoader getLoader()
    {
        return loader;
    }
    public void setLoader(LogLoader loader)
    {
        this.loader=loader;
    }

    public EnsemblDataPuller getDataPuller() {
        return dataPuller;
    }

    public void setDataPuller(EnsemblDataPuller dataPuller) {
        this.dataPuller = dataPuller;
        System.out.println(dataPuller.getVersion());
    }

    public Parser getDataParser() {
        return dataParser;
    }

    public void setDataParser(Parser dataParser) {
        this.dataParser = dataParser;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }
    /*public EnsemblParser getPaser()
    {
        return parser;
    }
    public void setParser(EnsemblParser parser)
    {
        this.parser=parser;
    }*/
}
