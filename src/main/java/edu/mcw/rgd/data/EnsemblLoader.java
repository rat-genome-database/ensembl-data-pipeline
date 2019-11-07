package edu.mcw.rgd.data;
import edu.mcw.rgd.dao.impl.RGDManagementDAO;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.PipelineLogFlagManager;
import edu.mcw.rgd.process.PipelineLogger;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;
import java.util.Collection;

/**
 * Created by sellanki on 8/6/2019.
 */
public class EnsemblLoader {
    RGDManagementDAO managementDAO = new RGDManagementDAO();
    EnsemblDAO ensemblDAO;
    PipelineLogger dbLogger = PipelineLogger.getInstance();
    PipelineLogFlagManager dbFlagManager = new PipelineLogFlagManager(dbLogger);
    EnsemblPipelinePreprocessor pipelinePreprocessor;
    private String version;
    private Log log = LogFactory.getLog("status");
    EnsemblGene gene;

    /**
     * starts the pipeline; properties are read from properties/AppConfigure.xml file
     * @param args cmd line arguments, like species
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        DefaultListableBeanFactory bf= new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        edu.mcw.rgd.data.EnsemblLoader loader=(edu.mcw.rgd.data.EnsemblLoader) (bf.getBean("loader"));
        loader.ensemblDAO = new EnsemblDAO();
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
        CounterPool counters = new CounterPool();
        log.info(SpeciesType.getCommonName(speciesTypeKey)+" " +getVersion());
        dbLogger.init(speciesTypeKey, "download+process", "Ensembl");
        pipelinePreprocessor.setSpeciesTypeKey(speciesTypeKey);
        try {
            Collection<EnsemblGene> genes = pipelinePreprocessor.run();
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
        System.out.println("Command line parameters required:");
        System.out.println(" -species 0|1|2|3|Rat|Mouse|Human|All");
    }
    public EnsemblPipelinePreprocessor getPipelinePreprocessor() {
        return pipelinePreprocessor;
    }

    public void setPipelinePreprocessor(EnsemblPipelinePreprocessor pipelinePreprocessor) {
        this.pipelinePreprocessor = pipelinePreprocessor;
    }
    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }
}
