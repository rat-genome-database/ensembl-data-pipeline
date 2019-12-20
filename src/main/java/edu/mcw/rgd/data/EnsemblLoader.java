package edu.mcw.rgd.data;
import edu.mcw.rgd.dao.impl.RGDManagementDAO;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.PipelineLogFlagManager;
import edu.mcw.rgd.process.PipelineLogger;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;
import java.util.Collection;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */
public class EnsemblLoader {


    private String version;
    static Logger log = Logger.getLogger("statuscheck");
    EnsemblDataPuller dataPuller;
    Parser dataParser;
    EnsemblGeneLoader geneLoader;
    EnsemblTranscriptLoader transcriptLoader;

    /**
     * starts the pipeline; properties are read from properties/AppConfigure.xml file
     * @param args cmd line arguments, like species
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        DefaultListableBeanFactory bf= new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        edu.mcw.rgd.data.EnsemblLoader loader=(edu.mcw.rgd.data.EnsemblLoader) (bf.getBean("loader"));

        // parse cmd line params
        if( args.length<2 ) {
            usage();
            return;
        }
        int speciesTypeKey = -1;
        if( args[0].equals("-species") ) {
            speciesTypeKey = SpeciesType.parse(args[1]);
        }
        // if species type key is all, run for all species
        if( speciesTypeKey<0 ) {
            throw new Exception("Aborted: please specify the species in cmd line");
        }
        if( speciesTypeKey== SpeciesType.ALL ) {
            loader.run(SpeciesType.RAT);
        }
        else {
            loader.run(speciesTypeKey);
        }
    }
    /**
     * run the Ensembl pipeline in download+process mode;
     * <ol>
     *     <li>download genes data from Ensembl biomart and store it locally in data folder</li>
     *     <li>download file with NcbiGene ids mapped to Ensembl ids</li>
     * </ol>
     * @param speciesTypeKey species type key
     * @throws Exception
     */
    public void run(int speciesTypeKey) throws Exception {
        log.info(SpeciesType.getCommonName(speciesTypeKey)+" " +getVersion());
        try {
          //  Collection<EnsemblGene> genes = pipelinePreprocessor.run();

           dataPuller.setSpeciesTypeKey(speciesTypeKey);
            String dataFile = dataPuller.downloadGenesFile();
            List<EnsemblGene> genes=dataParser.parseGene(dataFile);
            System.out.println("Total genes parsed from file: "+genes.size());
            geneLoader.run(genes,speciesTypeKey);

           String transcriptsFile = dataPuller.downloadTranscriptsFile();
           List<EnsemblTranscript> transcripts = dataParser.parseTranscript(transcriptsFile);
            transcriptLoader.run(transcripts,speciesTypeKey);
           System.out.println("Total transcripts parsed from file: "+transcripts.size());

       }
        catch(Exception e)
        {
            e.printStackTrace();
            throw e;
        }

    }
    /**
     * print to stdout the information about command line parameters
     */
    static public void usage() {
        log.info("Command line parameters required:");
        log.info(" -species 0|1|2|3|Rat|Mouse|Human|All");
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public EnsemblDataPuller getDataPuller() {
        return dataPuller;
    }

    public void setDataPuller(EnsemblDataPuller dataPuller) {
        this.dataPuller = dataPuller;
    }

    public EnsemblGeneLoader getGeneLoader() {
        return geneLoader;
    }

    public void setGeneLoader(EnsemblGeneLoader geneLoader) {
        this.geneLoader = geneLoader;
    }

    public Parser getDataParser() {
        return dataParser;
    }

    public void setDataParser(Parser dataParser) {
        this.dataParser = dataParser;
    }

    public EnsemblTranscriptLoader getTranscriptLoader() {
        return transcriptLoader;
    }

    public void setTranscriptLoader(EnsemblTranscriptLoader transcriptLoader) {
        this.transcriptLoader = transcriptLoader;
    }
}
