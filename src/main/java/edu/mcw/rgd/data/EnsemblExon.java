package edu.mcw.rgd.data;


/**
 * Created by hsnalabolu on 12/4/2019.
 */
public class EnsemblExon {

    int rgdId;
    int exonStart;
    int exonStop;
    String strand;
    String exonChromosome;
    String exonNumber;
    String exonTranscriptAccId;

    public int getRgdId() {
        return rgdId;
    }

    public void setRgdId(int rgdId) {
        this.rgdId = rgdId;
    }

    public int getExonStop() {
        return exonStop;
    }

    public void setExonStop(int exonStop) {
        this.exonStop = exonStop;
    }

    public int getExonStart() {
        return exonStart;
    }

    public void setExonStart(int exonStart) {
        this.exonStart = exonStart;
    }

    public String getExonChromosome() {
        return exonChromosome;
    }

    public void setExonChromosome(String exonChromosome) {
        this.exonChromosome = exonChromosome;
    }

    public String getExonNumber() {
        return exonNumber;
    }

    public void setExonNumber(String exonNumber) {
        this.exonNumber = exonNumber;
    }

    public String getExonTranscriptAccId() {
        return exonTranscriptAccId;
    }

    public void setExonTranscriptAccId(String exonTranscriptAccId) {
        this.exonTranscriptAccId = exonTranscriptAccId;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }
}