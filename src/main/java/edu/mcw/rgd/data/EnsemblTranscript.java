package edu.mcw.rgd.data;



import java.util.ArrayList;

/**
 * Created by hsnalabolu on 12/4/2019.
 */
public class EnsemblTranscript {

    int rgdId;
    String ensTranscriptId;
    String ensGeneId;
    String proteinId;
    String chromosome;
    int start;
    int stop;
    String strand;
    boolean nonCodingInd;
    String type;

    ArrayList<EnsemblExon> exonList = new ArrayList<EnsemblExon>();

    public int getRgdId() {
        return rgdId;
    }

    public void setRgdId(int rgdId) {
        this.rgdId = rgdId;
    }

    public String getEnsTranscriptId() {
        return ensTranscriptId;
    }

    public void setEnsTranscriptId(String ensTranscriptId) {
        this.ensTranscriptId = ensTranscriptId;
    }

    public String getEnsGeneId() {
        return ensGeneId;
    }

    public void setEnsGeneId(String ensGeneId) {
        this.ensGeneId = ensGeneId;
    }

    public String getProteinId() {
        return proteinId;
    }

    public void setProteinId(String proteinId) {
        this.proteinId = proteinId;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public int getStart() {
        return start;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getStop() {
        return stop;
    }

    public void setStop(int stop) {
        this.stop = stop;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public boolean isNonCodingInd() {
        return nonCodingInd;
    }

    public void setNonCodingInd(boolean nonCodingInd) {
        this.nonCodingInd = nonCodingInd;
    }
    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }


    public ArrayList<EnsemblExon> getExonList() {
        return exonList;
    }

    public void setExonList(ArrayList<EnsemblExon> exonList) {
        this.exonList = exonList;
    }
}
