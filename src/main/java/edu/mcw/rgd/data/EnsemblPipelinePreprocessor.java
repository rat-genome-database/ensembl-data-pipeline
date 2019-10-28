package edu.mcw.rgd.data;


import edu.mcw.rgd.process.PipelineLogger;

import java.util.Collection;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */

public class EnsemblPipelinePreprocessor {

    EnsemblDataPuller dataPuller;
    EnsemblDAO ensemblDAO;
     Parser dataParser;

    int speciesTypeKey;
    PipelineLogger dbLogger = PipelineLogger.getInstance();
     EnsemblGene gene;
    private EnsemblGeneLoader geneLoader;

    public Collection<EnsemblGene> run() throws Exception {
        // download genes data from Ensembl biomart and store it locally in data folder
       dataPuller.setSpeciesTypeKey(speciesTypeKey);

        String dataFile = dataPuller.downloadGenesFile();
       List<EnsemblGene> genes=dataParser.parseGene(dataFile);
        geneLoader.run(genes);
        return genes;
    }
    void traverseallgenes(Collection<EnsemblGene> genes) throws Exception
    {
        for(EnsemblGene gene:genes)
        {
            System.out.println(gene.getChromosome());
        }
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

    public void setGeneLoader(EnsemblGeneLoader geneLoader) {
        this.geneLoader = geneLoader;
    }

    public EnsemblGeneLoader getGeneLoader() {
        return geneLoader;
    }

}
