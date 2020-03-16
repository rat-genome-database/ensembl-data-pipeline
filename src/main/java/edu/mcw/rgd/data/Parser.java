package edu.mcw.rgd.data;
import edu.mcw.rgd.datamodel.SpeciesType;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */
public class Parser {

    private int speciesTypeKey;
    Logger log = Logger.getLogger("status");

    public List<EnsemblGene> parseGene(String inputFile) throws Exception {

        log.info(SpeciesType.getCommonName(speciesTypeKey)+": parsing gene file");
        BufferedReader reader = new BufferedReader(new FileReader(inputFile));

        List<EnsemblGene> genes = new ArrayList<>();
        String line;
        while( (line = reader.readLine()) != null ) {

            EnsemblGene gene = new EnsemblGene();
            String[] cols = line.split("\t", -1);
            if (cols.length < 6) {
                String msg = "5 columns expected, but found only " + cols.length + " in the file " + inputFile + "\n" +
                        "  offending line: [" + line + "]";
                log.warn(msg);
                throw new Exception(msg);
            }
            String ensemblGeneId = cols[0];
            String chromosome = cols[1];
            String startPos = cols[2];
            String stopPos = cols[3];
            String strand = cols[4];
            String entrezgene_id = cols[5];
            String gene_name = cols[6];
            String gene_biotype = cols[7];
            String gene_description = cols[8];
            String rgdid = "";
            if (cols.length > 9) {
                    rgdid = cols[9];
            }
            gene.setEnsemblGeneId(ensemblGeneId);

            if (rgdid.isEmpty()) {
                gene.setrgdid("0");
            } else {
                gene.setrgdid(rgdid);
            }
            if (entrezgene_id.isEmpty()) {
                gene.setEntrezgene_id("NIL");
            } else {
                gene.setEntrezgene_id(entrezgene_id);
            }
            gene.setChromosome(chromosome);
            if (gene_biotype.equals("protein_coding"))
                gene.setgene_biotype("protein-coding");
            else
                gene.setgene_biotype(gene_biotype);
            gene.setGeneSymbol(gene_name);
            gene.setStartPos(startPos);
            gene.setStopPos(stopPos);
            if (!gene_description.isEmpty()) {
                String gene_desc = "";
                if(gene_description.contains("[")) {
                    gene_desc = gene_description.substring(0, gene_description.indexOf('[')).trim();
                }
                else gene_desc = gene_description;
                gene.setgene_description(gene_desc);
            } else {
                gene.setgene_description(gene_description);
            }
            if (strand.equals("1"))
                gene.setStrand("+");
            else if (strand.equals("-1"))
                gene.setStrand("-");
            else
                gene.setStrand(strand);

            if(!gene.getGeneSymbol().isEmpty())
                genes.add(gene);
        }

        reader.close();
        return genes;
    }

    public Collection<EnsemblTranscript> parseTranscript(String inputFile) throws Exception {

        log.info(SpeciesType.getCommonName(speciesTypeKey)+": parsing transcript file");
        BufferedReader reader = new BufferedReader(new FileReader(inputFile));

        HashMap<String,EnsemblTranscript> transcripts = new HashMap<String,EnsemblTranscript>();
        String line;
        while ((line = reader.readLine()) != null)
        {
            String[] cols = line.split("\t", -1);
            if (cols.length < 10) {
                String msg = "10 columns expected, but found only " + cols.length + " in the file " + inputFile + "\n" +
                        "  offending line: [" + cols + "]";
                log.warn(msg);
                throw new Exception(msg);
            }

            EnsemblTranscript transcript = new EnsemblTranscript();
            transcript.setEnsGeneId(cols[0]);
            transcript.setEnsTranscriptId(cols[1]);
            transcript.setChromosome(cols[3]);
            transcript.setStart(Integer.valueOf(cols[4]));
            transcript.setStop(Integer.valueOf(cols[5]));
            if(cols[6].equals("1"))
                transcript.setStrand("+");
            else if(cols[6].equals("-1"))
                transcript.setStrand("-");


            EnsemblExon exon = new EnsemblExon();
            exon.setExonChromosome(cols[3]);
            exon.setExonStart(Integer.valueOf(cols[7]));
            exon.setExonStop(Integer.valueOf(cols[8]));
            if(cols[6].equals("1"))
                exon.setStrand("+");
            else if(cols[6].equals("-1"))
                exon.setStrand("-");
            exon.setExonTranscriptAccId(cols[1]);
            exon.setExonNumber(cols[9]);

            if(cols[10] != null && cols[11] != null)
                transcript.setNonCodingInd(false);

            transcript.setProteinId(cols[12]);
            transcript.setType(cols[13]);
            ArrayList<EnsemblExon> exons = new ArrayList<>();

            if(!transcripts.isEmpty() && transcripts.containsKey(transcript.getEnsTranscriptId())) {
                transcript = transcripts.get(transcript.getEnsTranscriptId());
                exons = transcript.getExonList();
            }

            exons.add(exon);
            transcript.setExonList(exons);
            transcripts.put(transcript.getEnsTranscriptId(),transcript);
        }

        reader.close();

        return transcripts.values();
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }
}
