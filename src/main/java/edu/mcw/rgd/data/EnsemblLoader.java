package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.MemoryMonitor;
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
    private List<AssemblyConfig> assemblies; // all assemblies processed by the pipeline, grouped per species in AppConfigure.xml

    private boolean skipGeneLoader = false;
    private boolean skipTranscriptLoader = false;
    // data source -- exactly one must be selected on the command line (no default)
    private boolean useGff3Loader = false;
    private boolean useBioMart = false;

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
        if( args.length<1 ) {
            usage();
            return;
        }

        int speciesTypeKey = SpeciesType.ALL; // default: all configured species
        List<Integer> mapKeys = new ArrayList<>();
        for( int argc=0; argc<args.length; argc++ ) {
            String arg = args[argc];
            if (arg.equals("-species")) {
                speciesTypeKey = SpeciesType.parse(args[++argc]);
            }
            else if( arg.equals("-mapKey") ) {
                for( String mapKeyStr: args[++argc].split(",") ) {
                    mapKeys.add(Integer.parseInt(mapKeyStr.trim()));
                }
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
            else if( arg.equals("-useBioMart") ) {
                loader.useBioMart = true;
            }
        }

        // the data source must be chosen explicitly -- exactly one of -useBioMart / -useGff3Loader
        if( loader.useBioMart == loader.useGff3Loader ) {
            if( loader.useBioMart ) {
                System.out.println("Aborted: specify only ONE data source: -useBioMart OR -useGff3Loader");
            } else {
                System.out.println("Aborted: please specify the data source: -useBioMart OR -useGff3Loader");
            }
            return;
        }

        List<AssemblyConfig> assembliesToRun = loader.selectAssemblies(speciesTypeKey, mapKeys);
        Collections.shuffle(assembliesToRun);
        for( AssemblyConfig assembly: assembliesToRun ) {
            loader.run(assembly);
        }
    }

    /**
     * select the configured assemblies to be processed in this run:
     * all of them by default, restricted by the optional '-species' and '-mapKey' cmdline filters
     */
    List<AssemblyConfig> selectAssemblies(int speciesTypeKey, List<Integer> mapKeys) throws Exception {

        List<AssemblyConfig> result = new ArrayList<>(getAssemblies());

        if( speciesTypeKey!=SpeciesType.ALL ) {
            result.removeIf(a -> a.getSpeciesTypeKey()!=speciesTypeKey);
        }

        if( !mapKeys.isEmpty() ) {
            result.removeIf(a -> !mapKeys.contains(a.getEnsemblMapKey()));

            for( int mapKey: mapKeys ) {
                boolean found = false;
                for( AssemblyConfig a: result ) {
                    if( a.getEnsemblMapKey()==mapKey ) {
                        found = true;
                        break;
                    }
                }
                if( !found ) {
                    throw new Exception("Aborted: map_key "+mapKey+" is not configured"
                        + (speciesTypeKey!=SpeciesType.ALL ? " for species "+SpeciesType.getCommonName(speciesTypeKey) : "")
                        + " -- see the 'assemblies' property in AppConfigure.xml");
                }
            }
        }

        // strain assemblies (Ensembl and NCBI positions share one map key) are loadable only from GFF3 files
        if( !useGff3Loader ) {
            for( AssemblyConfig a: new ArrayList<>(result) ) {
                if( a.isSharedMapKey() ) {
                    log.info("  "+a.describe()+" skipped -- available only with -useGff3Loader");
                    result.remove(a);
                }
            }
        }

        if( result.isEmpty() ) {
            throw new Exception("Aborted: no assemblies selected to run");
        }
        return result;
    }

    /**
     * run the Ensembl pipeline for one configured assembly in download+process mode;
     * <ol>
     *     <li>download genes data from Ensembl (gff3 file, or biomart) and store it locally in data folder</li>
     *     <li>download file with NcbiGene ids mapped to Ensembl ids</li>
     * </ol>
     * @param assembly assembly to be processed
     * @throws Exception
     */
    public void run(AssemblyConfig assembly) throws Exception {
        long time0 = System.currentTimeMillis();
        int speciesTypeKey = assembly.getSpeciesTypeKey();
        String speciesName = SpeciesType.getCommonName(speciesTypeKey);
        int ensemblMapKey = assembly.getEnsemblMapKey();
        int ncbiAssemblyMapKey = assembly.getNcbiMapKey(); // 0 = no NCBI assembly (f.e. naked mole-rat)
        log.info(speciesName+" " +getVersion());
        log.info("  assembly: "+assembly.describe());

        // QC pre-check: the assembly we load onto must be an Ensembl-source assembly in RGD;
        // a strain assembly shares one NCBI-source map with the NCBI positions -- allowed for the GFF3
        // path, because the assembly identity is verified by GenBank accession right after the download
        checkLoadingAssemblyIsEnsembl(ensemblMapKey, speciesName, useGff3Loader && assembly.isSharedMapKey());

        // GFF3 path: download the assembly's gff3 + entrez files and configure the parser
        if( useGff3Loader ) {
            prepareGff3Parser(assembly);
        }

        MemoryMonitor memoryMonitor = new MemoryMonitor();
        memoryMonitor.start();
        try {

            dataPuller.setSpeciesTypeKey(speciesTypeKey);
            if( skipGeneLoader ) {
                log.warn("WARNING: gene processing skipped!");
            } else {
                List<EnsemblGene> genes;
                if( useGff3Loader ) {
                    genes = dataGff3Parser.parseGenes();
                    ensemblMapKey = dataGff3Parser.getEnsemblAssemblyMapKey();
                    ncbiAssemblyMapKey = dataGff3Parser.getNcbiAssemblyMapKey();
                } else {
                    validateAssemblyName(assembly);

                    String dataFile = dataPuller.downloadGenesFile();
                    genes = dataParser.parseGene(dataFile);

                    if( genes==null ) {
                        log.info(speciesName.toUpperCase() + " SKIPPED -- TIME ELAPSED " + Utils.formatElapsedTime(time0, System.currentTimeMillis()));
                        log.info("===");
                        return;
                    }
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
        finally {
            memoryMonitor.stop();
            log.info(memoryMonitor.getSummary());
        }
    }

    /**
     * QC pre-check run before processing an assembly: the assembly map we load Ensembl positions onto
     * must itself be an Ensembl-source assembly in RGD. Guards against a misconfigured map key
     * that would load Ensembl data onto an NCBI (or other) assembly. For a strain assembly the map is
     * typically NCBI-source: allowed (GFF3 path only), because the Ensembl positions coexist on the same
     * map key, distinguished by SRC_PIPELINE='Ensembl', and the assembly identity is verified by
     * GenBank accession.
     */
    void checkLoadingAssemblyIsEnsembl(int ensemblMapKey, String speciesName, boolean allowNonEnsemblSource) throws Exception {

        edu.mcw.rgd.datamodel.Map map = new EnsemblDAO().getAssemblyMap(ensemblMapKey);
        if( map==null ) {
            throw new Exception("QC pre-check failed for "+speciesName+": assembly map_key "+ensemblMapKey+" not found in RGD");
        }
        if( !Utils.stringsAreEqualIgnoreCase(map.getSource(), "Ensembl") ) {
            if( allowNonEnsemblSource ) {
                log.warn("  QC: loading assembly map_key "+ensemblMapKey+" ["+map.getName()+"] has source '"+map.getSource()
                        +"'; Ensembl positions will coexist on this map, distinguished by SRC_PIPELINE='Ensembl';"
                        +" assembly identity is verified by GenBank accession");
                return;
            }
            throw new Exception("QC pre-check failed for "+speciesName+": loading assembly map_key "+ensemblMapKey
                    +" ["+map.getName()+"] has source '"+map.getSource()+"', expected 'Ensembl'");
        }
        log.info("  QC: loading assembly map_key "+ensemblMapKey+" ["+map.getName()+"] source=Ensembl -- OK");
    }

    /**
     * QC pre-check (GFF3 path): make certain the map_key we load onto is the SAME assembly as the downloaded
     * GFF3 -- so there is no doubt which assembly the data goes into. The check is by GenBank assembly
     * accession only (the definitive identity of an assembly); a missing or mismatched accession aborts the
     * assembly rather than guessing by assembly name.
     * @param ensemblMapKey the map_key the pipeline loads Ensembl positions onto
     * @param gff3AssemblyName  GFF3 '#!genome-build' value, f.e. 'Naked_mole-rat_maternal' (for messages)
     * @param gff3Accession     GFF3 '#!genome-build-accession' value, f.e. 'GCA_944319715.1'
     */
    void verifyLoadingAssembly(int ensemblMapKey, String gff3AssemblyName, String gff3Accession) throws Exception {

        edu.mcw.rgd.datamodel.Map map = new EnsemblDAO().getAssemblyMap(ensemblMapKey);
        if( map==null ) {
            throw new Exception("QC pre-check failed: assembly map_key "+ensemblMapKey+" not found in RGD");
        }
        String mapAcc = map.getGenBankAssemblyAcc();

        if( Utils.isStringEmpty(gff3Accession) ) {
            throw new Exception("QC pre-check failed: GFF3 for "+gff3AssemblyName+" has no '#!genome-build-accession' line");
        }
        if( Utils.isStringEmpty(mapAcc) ) {
            throw new Exception("QC pre-check failed: map_key "+ensemblMapKey+" ["+map.getName()
                    +"] has no GenBank accession in RGD -- set maps.genbank_assembly_acc to "+gff3Accession);
        }
        if( !Utils.stringsAreEqualIgnoreCase(gff3Accession, mapAcc) ) {
            throw new Exception("QC pre-check failed: WRONG ASSEMBLY -- map_key "+ensemblMapKey+" ["+map.getName()
                    +"] GenBank acc "+mapAcc+" != GFF3 genome-build-accession "+gff3Accession+" ("+gff3AssemblyName+")");
        }
        log.info("  QC: assembly map_key "+ensemblMapKey+" ["+map.getName()+"] GenBank acc "+mapAcc
                +" == GFF3 "+gff3Accession+" -- OK");
    }

    /// download the assembly's gff3 + entrez source files (full URLs configured per assembly in
    /// AppConfigure.xml) and configure the gff3 parser
    void prepareGff3Parser(AssemblyConfig assembly) throws Exception {

        dataPuller.setSpeciesTypeKey(assembly.getSpeciesTypeKey());
        if( Utils.isStringEmpty(assembly.getGff3Url()) || Utils.isStringEmpty(assembly.getEntrezUrl()) ) {
            throw new Exception("no GFF3/entrez source file configured for "+assembly.describe()
                    +" -- see the 'assemblies' property in AppConfigure.xml");
        }

        String gff3File = dataPuller.downloadEnsemblFile(assembly.getGff3Url());
        String entrezFile = dataPuller.downloadEnsemblFile(assembly.getEntrezUrl());

        // pre-check: the map_key we load onto must be the SAME assembly as the GFF3 we just downloaded
        String[] gff3Assembly = EnsemblGff3Parser.readAssemblyHeader(gff3File); // { name, GenBank accession }
        verifyLoadingAssembly(assembly.getEnsemblMapKey(), gff3Assembly[0], gff3Assembly[1]);

        dataGff3Parser.setGff3File(gff3File);
        dataGff3Parser.setEntrezFile(entrezFile);
        dataGff3Parser.setGenomeBuild(resolveEnsemblAssembly(assembly));
        dataGff3Parser.setXrefAuthority(assembly.getXrefAuthority());
        dataGff3Parser.setEnsemblAssemblyMapKey(assembly.getEnsemblMapKey());
        dataGff3Parser.setNcbiAssemblyMapKey(assembly.getNcbiMapKey());

        log.info("  GFF3: "+assembly.getGff3Url());
    }

    /// Ensembl assembly name for a configured assembly: the assemblyName override if present, else the RGD
    /// map name, with the trailing ' Ensembl' stripped. f.e. map 381 -> 'GRCr8';
    /// naked mole-rat (override) -> 'Naked_mole-rat_maternal'
    String resolveEnsemblAssembly(AssemblyConfig assembly) throws Exception {
        String name = assembly.getAssemblyName();
        if( name==null ) {
            name = new EnsemblDAO().getAssemblyMap(assembly.getEnsemblMapKey()).getName();
        }
        return name.replace(" Ensembl", "").trim();
    }

    void validateAssemblyName(AssemblyConfig assembly) throws Exception {

        String assemblyName = dataPuller.getAssemblyNameFromEnsemblRest();
        log.info("Ensembl Rest service: assembly name: "+assemblyName);
        String expectedAssemblyNameInRgd = assemblyName+" Ensembl";

        edu.mcw.rgd.datamodel.Map ensemblMap = new EnsemblDAO().getAssemblyMap(assembly.getEnsemblMapKey());

        if( !expectedAssemblyNameInRgd.equalsIgnoreCase(ensemblMap.getName()) ) {

            // try to use alternate assembly name, if available
            String altAssemblyName = assembly.getAssemblyName();
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
        System.out.println("Command line parameters:");
        System.out.println(" data source (pick one, required): -useBioMart | -useGff3Loader");
        System.out.println(" optional: -species 1|2|3|...|Rat|Mouse|Human|... -- process only the given species");
        System.out.println("           (default: all species and assemblies configured in AppConfigure.xml)");
        System.out.println(" optional: -mapKey <k1[,k2,...]> -- process only the given assembly map(s), f.e. -mapKey 301,302,303");
        System.out.println(" optional: -skipGenes -skipTranscripts");
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

    public List<AssemblyConfig> getAssemblies() {
        return assemblies;
    }

    public void setAssemblies(List<AssemblyConfig> assemblies) {
        this.assemblies = assemblies;
    }

    public EnsemblGff3Parser getDataGff3Parser() {
        return dataGff3Parser;
    }

    public void setDataGff3Parser(EnsemblGff3Parser dataGff3Parser) {
        this.dataGff3Parser = dataGff3Parser;
    }
}
