package edu.mcw.rgd.data;
import edu.mcw.rgd.process.PipelineLogger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.sql.Connection;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by sellanki on 8/6/2019.
 */
public class Parser
{
    PipelineLogger dbLogger = PipelineLogger.getInstance();
    EnsemblDAO ensemblDAO = new EnsemblDAO();

    EnsemblGene gene;
    Connection conn;
    int speciesTypeKey;
    public Parser() throws Exception
    {
    }
    public List<EnsemblGene> parseGene(String inputFile) throws Exception
    {
        dbLogger.log("  parsing gene file ", inputFile, PipelineLogger.INFO);
        BufferedReader reader = new BufferedReader(new FileReader(inputFile));

        List<EnsemblGene> genes = new ArrayList<EnsemblGene>();
        String line;
        while ((line = reader.readLine()) != null)
        {
            String[] cols = line.split("\t", -1);
            if (cols.length < 6) {
                System.out.println("\n" + cols);
                throw new Exception("5 columns expected, but found only " + cols.length + " in the file " + inputFile + "\n" +
                        "  offending line: [" + line + "]");
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
            if (gene != null) {
                genes.add(gene);
                gene = null;
            }
            if (gene == null) {
                gene = new EnsemblGene();
                gene.setEnsemblGeneId(ensemblGeneId);
            }
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
            gene.setgene_name(gene_name);
            gene.setStartPos(startPos);
            gene.setStopPos(stopPos);
            if (!gene_description.isEmpty()) {
                String gene_desc = gene_description.substring(0, gene_description.indexOf('[')).trim();
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
        }
        if(!ensemblDAO.getAliasTypes().contains("ensembl_gene_symbol"))
        {
            ensemblDAO.insertAliasType("ensembl_gene_symbol", "Created by Ensembl pipeline");
        }
        if(!ensemblDAO.getAliasTypes().contains("ensembl_full_name"))
        {
            ensemblDAO.insertAliasType("ensembl_full_name", "created by Ensembl pipeline");
        }
        if(!ensemblDAO.getAliasTypes().contains("ensembl_gene_type"))
        {
            ensemblDAO.insertAliasType("ensembl_gene_type", "created by ensembl pipeline");
        }
        reader.close();
            return genes;
        }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey)
    {
        this.speciesTypeKey = speciesTypeKey;
    }

    }








