package edu.mcw.rgd.data;

import java.util.Set;

/**
 * Created by sellanki on 8/6/2019.
 */
public class EnsemblGene {

    private String rgdid;
    private String entrezGeneId;
    // incoming data from Ensembl
    private String ensemblGeneId;
    // genomic position
    private String chromosome;
    private String startPos, stopPos;
    private String strand;
    private String geneSymbol;
    private String geneBioType;
    private String geneName;

    // transcript RefSeq acc ids, used by EnsemblGff3Parser
    private Set<String> refseqAccIds;

    public String getGeneName() {
        return geneName;
    }
    public void setGeneName( String geneName ) {
        this.geneName = geneName;
    }

    public String getGeneSymbol()
    {
        return geneSymbol;
    }
    public void setGeneSymbol(String geneSymbol)
    {
        this.geneSymbol=geneSymbol;
    }

    public String getGeneBioType() {
        return geneBioType;
    }
    public void setGeneBioType(String geneBioType) {
        this.geneBioType = geneBioType;
    }

    public String getrgdid() {
        return rgdid;
    }
    public void setrgdid(String rgdid) {
        this.rgdid = rgdid;
    }

    public String getEntrezGeneId()
    {
        return entrezGeneId;
    }
    public void setEntrezGeneId(String entrezGeneId)
    {
        this.entrezGeneId=entrezGeneId;
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