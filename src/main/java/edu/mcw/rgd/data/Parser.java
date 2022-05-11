package edu.mcw.rgd.data;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */
public class Parser {

    Logger log = LogManager.getLogger("status");

    public List<EnsemblGene> parseGene(String inputFile) throws Exception {

        BufferedReader reader = Utils.openReader(inputFile);

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
            String geneName = cols[8];
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
            if (!geneName.isEmpty()) {
                String gene_desc = "";
                if(geneName.contains("[")) {
                    gene_desc = geneName.substring(0, geneName.indexOf('[')).trim();
                }
                else gene_desc = geneName;
                gene.setGeneName(gene_desc);
            } else {
                gene.setGeneName(geneName);
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

        BufferedReader reader = Utils.openReader(inputFile);

        HashMap<String,EnsemblTranscript> transcripts = new HashMap<>();
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

            String trAcc = cols[1];
            EnsemblTranscript transcript = transcripts.get(trAcc);
            if( transcript==null ) {
                transcript = new EnsemblTranscript();
                transcripts.put(trAcc, transcript);

                transcript.setEnsGeneId(cols[0]);
                transcript.setEnsTranscriptId(cols[1]);
                transcript.setEnsTranscriptVer(cols[2]);
                transcript.setChromosome(cols[3]);
                transcript.setStart(Integer.parseInt(cols[4]));
                transcript.setStop(Integer.parseInt(cols[5]));

                if(cols[6].equals("1"))
                    transcript.setStrand("+");
                else if(cols[6].equals("-1"))
                    transcript.setStrand("-");
                if(cols[10] != null && cols[11] != null)
                    transcript.setNonCodingInd(false);

                transcript.setProteinId(cols[12]);
                transcript.setType(cols[13]);
            }

            EnsemblExon exon = new EnsemblExon();
            exon.setExonChromosome(cols[3]);
            exon.setExonStart(Integer.parseInt(cols[7]));
            exon.setExonStop(Integer.parseInt(cols[8]));
            if(cols[6].equals("1"))
                exon.setStrand("+");
            else if(cols[6].equals("-1"))
                exon.setStrand("-");
            exon.setExonTranscriptAccId(cols[1]);
            exon.setExonNumber(cols[9]);

            if( !Utils.isStringEmpty(cols[10]) ) {
                exon.setCdsStart(Integer.parseInt(cols[10]));
            }
            if( !Utils.isStringEmpty(cols[11]) ) {
                exon.setCdsStop(Integer.parseInt(cols[11]));
            }

            transcript.getExonList().add(exon);
        }

        reader.close();

        return transcripts.values();
    }
}
