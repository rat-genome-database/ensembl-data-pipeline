package edu.mcw.rgd.data;

import java.util.List;
import java.util.Set;

/**
 * Created by sellanki on 8/6/2019.
 */
public class EnsemblGene {

    private String rgdid;
    private String entrezgene_id;
    // incoming data from Ensembl
    private String ensemblGeneId;
    // genomic position
    private String chromosome;
    private String startPos, stopPos;
    private String strand;
    private String geneSymbol;
    private String gene_biotype;
    private String gene_description;

    // transcript RefSeq acc ids, used by EnsemblGff3Parser
    private Set<String> refseqAccIds;

    public String getgene_description(){return gene_description;}
    public void setgene_description(String gene_description){this.gene_description=gene_description;}

    public String getGeneSymbol()
    {
        return geneSymbol;
    }
    public void setGeneSymbol(String geneSymbol)
    {
        this.geneSymbol=geneSymbol;
    }

    public String getgene_biotype()
    {
        return gene_biotype;
    }
    public void setgene_biotype(String gene_biotype)
    {
        this.gene_biotype=gene_biotype;
    }

    public String getrgdid() {
        return rgdid;
    }
    public void setrgdid(String rgdid) {
        this.rgdid = rgdid;
    }

    public String getEntrezgene_id()
    {
        return entrezgene_id;
    }
    public void setEntrezgene_id(String entrezgene_id)
    {
        this.entrezgene_id=entrezgene_id;
    }

    public String getChromosome() {
        return chromosome;
    }
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getStartPos() {
        return startPos;
    }
    public void setStartPos(String startPos) {
        this.startPos = startPos;
    }

    public String getStopPos() {
        return stopPos;
    }
    public void setStopPos(String stopPos) {
        this.stopPos = stopPos;
    }

    public String getStrand() {
        return strand;
    }
    public void setStrand(String strand) {
        this.strand = strand;
    }

    public String getEnsemblGeneId() {
        return ensemblGeneId;
    }

    public void setEnsemblGeneId(String ensemblGeneId) {
        this.ensemblGeneId = ensemblGeneId;
    }

    public Set<String> getRefseqAccIds() {
        return refseqAccIds;
    }

    public void setRefseqAccIds(Set<String> refseqAccIds) {
        this.refseqAccIds = refseqAccIds;
    }
}