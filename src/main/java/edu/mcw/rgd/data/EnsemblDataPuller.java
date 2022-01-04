package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.process.FileDownloader;
import edu.mcw.rgd.process.FileExternalSort;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.net.URLEncoder;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */
public class EnsemblDataPuller {
    String websiteUrl; // ensembl biomart url
    String restGenomeUrl;
    String biomartQueryTemplate; // xml config file containing the query template to be POST-ed to biomart
    List<String> biomartQueryAttrsGenes; // list of attributes to be loaded for main genes query
    List<String> biomartQueryAttrsTranscripts; // list of attributes to be loaded for main transcripts query
    List<String> biomartQueryAttrsRatGenes; // extra attributes to be loaded for rat genes
    List<String> biomartQueryAttrsMouseGenes;
    List<String> biomartQueryAttrsHumanGenes;
    Logger statuslog = LogManager.getLogger("status");
    int speciesTypeKey;
    SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMdd");

    /**
     * download data from Ensembl biomart
     * @return name of the local copy of sorted data file
     * @throws Exception
     */
    public String downloadGenesFile() throws Exception {

        // create a list with all attributes needed by biomart query
        // it consist of attributes shared for any species, and attrs specific to species
        List<String> attrs = new ArrayList<String>(this.getBiomartQueryAttrsGenes());
        if( speciesTypeKey== SpeciesType.RAT )
            attrs.addAll(this.getBiomartQueryAttrsRatGenes());
        else if(speciesTypeKey == SpeciesType.MOUSE)
            attrs.addAll(this.getBiomartQueryAttrsMouseGenes());
        else if(speciesTypeKey == SpeciesType.HUMAN)
            attrs.addAll(this.getBiomartQueryAttrsHumanGenes());

        // build biomart url
        // and download the file from biomart given the url
        String speciesName = SpeciesType.getCommonName(speciesTypeKey).toLowerCase();
        statuslog.info("Downloading the genes file for "+ speciesName);

        String finalOutputFile = "data/" + sdf.format(new Date()) +"_"+speciesName + "_genes.txt.sorted.gz";
        if( !(new File(finalOutputFile).exists()) ) {

            String fileName = "data/" + speciesName + "_genes.txt";
            String downloadedFile = downloadFile(attrs, fileName);
            sortFile(downloadedFile, finalOutputFile);
        }
        return finalOutputFile;
    }

    /**
     * download data from Ensembl biomart
     * @return name of the local copy of sorted data file
     * @throws Exception
     */
    public String downloadTranscriptsFile() throws Exception {

        // create a list with all attributes needed by biomart query
        // it consist of attributes shared for any species, and attrs specific to species
        List<String> attrs = new ArrayList<String>(this.getBiomartQueryAttrsTranscripts());

        // build biomart url
        // and download the file from biomart given the url
        String speciesName = SpeciesType.getCommonName(speciesTypeKey).toLowerCase();
        statuslog.info("Downloading the transcripts file for "+ speciesName);

        String finalOutputFile = "data/" + sdf.format(new Date()) +"_"+speciesName + "_transcripts.txt.sorted.gz";
        if( !(new File(finalOutputFile).exists()) ) {

            String fileName = "data/" + speciesName + "_transcripts.txt";
            String downloadedFile = downloadFile(attrs, fileName);
            sortFile(downloadedFile, finalOutputFile);
        }
        return finalOutputFile;
    }

    // return the name of output file
    String downloadFile(List<String> attributes, String outFile) throws Exception {

        String data = buildBiomartQuery(attributes);

        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile(websiteUrl + "?" + data);
        downloader.setLocalFile(outFile);
        String outPath = downloader.download();

        return outPath;
    }

    // return url encoded biomart xml query
    String buildBiomartQuery(List<String> attributes) throws Exception {

        // load the biomart config file and remove all new lines
        BufferedReader fileReader = new BufferedReader(new FileReader(this.getBiomartQueryTemplate()));
        String line;
        StringBuilder buf = new StringBuilder();
        while( (line=fileReader.readLine())!=null ) {
            line = line.trim();
            // replace any line with '#SPECIES#' tag
            if( line.contains("#SPECIES#") ) {
                // from species taxonomic name take first letter of 1st word and entire second word
                // f.e. 'Homo sapiens' --> 'hsapiens'


                String[] taxoNameWords = SpeciesType.getTaxonomicName(speciesTypeKey).split(" ");
                String taxoName = taxoNameWords[0].substring(0, 1).toLowerCase() + taxoNameWords[1];
                if(speciesTypeKey == SpeciesType.DOG)
                    taxoName = "clfamiliaris";
                line = line.replaceAll("\\#SPECIES\\#", taxoName);
            }
            buf.append(line);

            // if we are processing <Attributes> line, add all custom attributes as well
            if( line.contains("Attribute") ) {
                int quoteStart = line.indexOf('"');
                int quoteStop = line.lastIndexOf('"');
                String prefix = line.substring(0, quoteStart+1);
                String suffix = line.substring(quoteStop);
                for( String attr: attributes ) {
                    buf.append(prefix).append(attr).append(suffix);
                }
            }
        }
        fileReader.close();

        return "query="+URLEncoder.encode(buf.toString(), "UTF-8");
    }

    void sortFile(String inputFile, String finalOutputFile) throws IOException {

        // now sort the file using external merge sort
        String outputFile = inputFile + ".tmp";
        Comparator<String> comparator = new Comparator<String>() {
            public int compare(String r1, String r2){
                return r1.compareTo(r2);
            }
        };
        File inFile = new File(inputFile);
        List<File> l = FileExternalSort.sortInBatch(inFile, comparator) ;
        FileExternalSort.mergeSortedFiles(l, new File(outputFile), comparator, false);

        // to save disk space, compress the file
        compressFile(outputFile, finalOutputFile);
    }

    void compressFile(String inFileName, String outFileName) throws IOException {
        BufferedReader in = Utils.openReader(inFileName);
        BufferedWriter out = Utils.openWriter(outFileName);
        String line;
        while( (line=in.readLine())!=null ) {
            out.write(line);
            out.write("\n");
        }
        in.close();
        out.close();
    }

    String getAssemblyNameFromEnsemblRest() throws Exception {

        FileDownloader fd = new FileDownloader();
        String url = getRestGenomeUrl() + SpeciesType.getTaxonomicName(speciesTypeKey).replace(" ", "_");
        fd.setExternalFile(url);
        fd.setLocalFile(null);
        String content = fd.download();

        // sexample line that interests us: "assembly_name: Sscrofa11.1"
        String[] lines = content.split("[\\n\\r]");
        String pattern = "assembly_name:";
        for( String line: lines ) {
            int pos = line.indexOf(pattern);
            if( pos>=0 ) {
                String assemblyName = line.substring(pos+pattern.length()).trim();
                return assemblyName;
            }
        }

        throw new Exception("Cannot retrieve an assembly name from "+url);
    }

    public String getWebsiteUrl() {
        return websiteUrl;
    }

    public void setWebsiteUrl(String websiteUrl) {
        this.websiteUrl = websiteUrl;
    }

    public String getBiomartQueryTemplate() {
        return biomartQueryTemplate;
    }

    public void setBiomartQueryTemplate(String biomartQueryTemplate) {
        this.biomartQueryTemplate = biomartQueryTemplate;
    }

    public List<String> getBiomartQueryAttrsGenes() {
        return biomartQueryAttrsGenes;
    }

    public void setBiomartQueryAttrsGenes(List<String> biomartQueryAttrsGenes) {
        this.biomartQueryAttrsGenes = biomartQueryAttrsGenes;
    }

    public List<String> getBiomartQueryAttrsTranscripts() {
        return biomartQueryAttrsTranscripts;
    }

    public void setBiomartQueryAttrsTranscripts(List<String> biomartQueryAttrsTranscripts) {
        this.biomartQueryAttrsTranscripts = biomartQueryAttrsTranscripts;
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public List<String> getBiomartQueryAttrsRatGenes() {
        return biomartQueryAttrsRatGenes;
    }

    public void setBiomartQueryAttrsRatGenes(List<String> biomartQueryAttrsRatGenes) {
        this.biomartQueryAttrsRatGenes = biomartQueryAttrsRatGenes;
    }

    public List<String> getBiomartQueryAttrsMouseGenes() {
        return biomartQueryAttrsMouseGenes;
    }

    public void setBiomartQueryAttrsMouseGenes(List<String> biomartQueryAttrsMouseGenes) {
        this.biomartQueryAttrsMouseGenes = biomartQueryAttrsMouseGenes;
    }

    public List<String> getBiomartQueryAttrsHumanGenes() {
        return biomartQueryAttrsHumanGenes;
    }

    public void setBiomartQueryAttrsHumanGenes(List<String> biomartQueryAttrsHumanGenes) {
        this.biomartQueryAttrsHumanGenes = biomartQueryAttrsHumanGenes;
    }

    public String getRestGenomeUrl() {
        return restGenomeUrl;
    }

    public void setRestGenomeUrl(String restGenomeUrl) {
        this.restGenomeUrl = restGenomeUrl;
    }
}
