package edu.mcw.rgd.data;

import java.util.Collection;
import java.util.List;

/**
 * Created by sellanki on 9/9/2019.
 */
public class PipelinePreprocessor2
{
    LogLoader loader;
    int speciesTypeKey;
    EnsemblGene gene;
    public Collection<EnsemblGene> run() throws Exception
    {
        String dataFile1="C:/git/ensembl-data/src/main/dist/logs/mismatchofexternalidwithrgdid.log";
        List<EnsemblGene> log_genes = loader.parseGene(dataFile1);
        //Files.delete(Paths.get("C:\\git\\ensembl-data\\src\\main\\dist\\logs\\ncbiexistbutnorgdid.log"));
        return log_genes;

    }
    public LogLoader getLoader()
    {
        return loader;
    }
    public void setLoader(LogLoader loader)
    {
        this.loader=loader;
    }
}
