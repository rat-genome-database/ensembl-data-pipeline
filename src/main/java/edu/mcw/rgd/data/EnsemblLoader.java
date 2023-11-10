package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.util.*;

/**
 * Created by sellanki on 8/6/2019.
 */
public class EnsemblLoader {

    private String version;
    Logger log = LogManager.getLogger("status");
    EnsemblDataPuller dataPuller;
    Parser dataParser;
    EnsemblGff3Parser dataGff3Parser;
    EnsemblGeneLoader geneLoader;
    EnsemblTranscriptLoader transcriptLoader;
    private Map<Integer, Integer> ensemblAssemblyMap;
    private Map<Integer, Integer> ncbiAssemblyMap;
    private Map<Integer, String> assemblyMapNames;

    private boolean skipGeneLoader = false;
    private boolean skipTranscriptLoader = false;
    private boolean useGff3Loader = false; // by default load data from Ensembl BioMart

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
        for( int argc=0; argc<args.length; argc++ ) {
            String arg = args[argc];
            if (arg.equals("-species")) {
                speciesTypeKey = SpeciesType.parse(args[++argc]);
            }
            else if( arg.equals("-skipGenes") ) {
                loader.skipGeneLoader = true;
            }
            else if( arg.equals("-skipTranscripts") ) {
                loader.skipTranscriptLoader = true;
            }
            else if( arg.equals("-useGff3Loader") ) {
                loader.useGff3Loader = true;
            }
        }

        // if species type key is all, run for all species
        if( speciesTypeKey<0 ) {
            throw new Exception("Aborted: please specify the species in cmd line");
        }

        if( speciesTypeKey== SpeciesType.ALL ) {
            List<Integer> speciesList = new ArrayList<>(loader.getEnsemblAssemblyMap().keySet());
            Collections.shuffle(speciesList);
            for( Integer spKey: speciesList){
                loader.run(spKey);
            }
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
        long time0 = System.currentTimeMillis();
        String speciesName = SpeciesType.getCommonName(speciesTypeKey);
        int ensemblMapKey = getEnsemblAssemblyMap().get(speciesTypeKey);
        int ncbiAssemblyMapKey = getNcbiAssemblyMap().get(speciesTypeKey);
        log.info(speciesName+" " +getVersion());
        try {

            dataPuller.setSpeciesTypeKey(speciesTypeKey);
            validateAssemblyName(ensemblMapKey);
            if( skipGeneLoader ) {
                log.warn("WARNING: gene processing skipped!");
            } else {
                List<EnsemblGene> genes;
                if( useGff3Loader ) {
                    genes = dataGff3Parser.parseGenes();
                    ensemblMapKey = dataGff3Parser.getEnsemblAssemblyMapKey();
                    ncbiAssemblyMapKey = dataGff3Parser.getNcbiAssemblyMapKey();
                } else {
                    String dataFile = dataPuller.downloadGenesFile();
                    genes = dataParser.parseGene(dataFile);
                }
                log.info("Total genes parsed from file: " + genes.size());

                Collections.shuffle(genes);
                geneLoader.run(genes, speciesTypeKey, ensemblMapKey, ncbiAssemblyMapKey);
            }

            if( skipTranscriptLoader ) {
                log.warn("WARNING: transcript processing skipped!");
            } else {
                String transcriptsFile;
                if( useGff3Loader ) {
                    ensemblMapKey = dataGff3Parser.getEnsemblAssemblyMapKey();
                    ncbiAssemblyMapKey = dataGff3Parser.getNcbiAssemblyMapKey();
                    transcriptsFile = dataGff3Parser.generateTranscriptFile();
                } else {
                    transcriptsFile = dataPuller.downloadTranscriptsFile();
                }
                TranscriptVersionManager.getInstance().init();
                Collection<EnsemblTranscript> transcripts = dataParser.parseTranscript(transcriptsFile);
                log.info("Total transcripts parsed from file: " + transcripts.size());
                transcriptLoader.run(transcripts, speciesTypeKey, ensemblMapKey, ncbiAssemblyMapKey);
            }

            log.info(speciesName.toUpperCase()+" DONE -- TIME ELAPSED "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
            log.info("===");
        }
        catch(Exception e) {
            Utils.printStackTrace(e, log);
            throw e;
        }
    }

    void validateAssemblyName(int ensemblMapKey) throws Exception {

        String assemblyName = dataPuller.getAssemblyNameFromEnsemblRest();
        log.info("Ensembl Rest service: assembly name: "+assemblyName);
        String expectedAssemblyNameInRgd = assemblyName+" Ensembl";

        edu.mcw.rgd.datamodel.Map ensemblMap = new EnsemblDAO().getAssemblyMap(ensemblMapKey);

        if( !expectedAssemblyNameInRgd.equalsIgnoreCase(ensemblMap.getName()) ) {

            // try to use alternate map names, if available
            String altAssemblyName = getAssemblyMapNames().get(ensemblMapKey);
            if( altAssemblyName!=null ) {
                if (expectedAssemblyNameInRgd.equalsIgnoreCase(altAssemblyName)) {
                    return; // validation succeeded on alt assembly name
                } else {
                    throw new Exception("Assembly map mismatch: expected [" + altAssemblyName + "] found [" + expectedAssemblyNameInRgd + "]");
                }
            }

            throw new Exception("Assembly map mismatch: expected ["+ensemblMap.getName()+"] found ["+expectedAssemblyNameInRgd+"]");
        }
    }

    /**
     * print to stdout the information about command line parameters
     */
    static public void usage() {
        System.out.println("Command line parameters required:");
        System.out.println(" -species 0|1|2|3|...|Rat|Mouse|Human|...|All");
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

    public void setEnsemblAssemblyMap(Map<Integer, Integer> ensemblAssemblyMap) {
        this.ensemblAssemblyMap = ensemblAssemblyMap;
    }

    public Map<Integer, Integer> getEnsemblAssemblyMap() {
        return ensemblAssemblyMap;
    }

    public Map<Integer, Integer> getNcbiAssemblyMap() {
        return ncbiAssemblyMap;
    }

    public void setNcbiAssemblyMap(Map<Integer, Integer> ncbiAssemblyMap) {
        this.ncbiAssemblyMap = ncbiAssemblyMap;
    }

    public Map<Integer, String> getAssemblyMapNames() {
        return assemblyMapNames;
    }

    public void setAssemblyMapNames(Map<Integer, String> assemblyMapNames) {
        this.assemblyMapNames = assemblyMapNames;
    }

    public EnsemblGff3Parser getDataGff3Parser() {
        return dataGff3Parser;
    }

    public void setDataGff3Parser(EnsemblGff3Parser dataGff3Parser) {
        this.dataGff3Parser = dataGff3Parser;
    }
}
