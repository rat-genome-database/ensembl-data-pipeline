package edu.mcw.rgd.data;


import edu.mcw.rgd.datamodel.TranscriptFeature;
import edu.mcw.rgd.process.Utils;

/**
 * Created by hsnalabolu on 12/4/2019.
 */
public class EnsemblExon {

    int exonStart;
    int exonStop;
    int cdsStart;
    int cdsStop;
    String strand;
    String exonChromosome;
    String exonNumber;
    String exonTranscriptAccId;

    public boolean matchesByPos(TranscriptFeature tf) {
        return Utils.intsAreEqual(this.exonStart, tf.getStartPos())
        && Utils.intsAreEqual(this.exonStop, tf.getStopPos())
        && Utils.stringsAreEqual(this.strand, tf.getStrand())
        && Utils.stringsAreEqual(this.exonChromosome, tf.getChromosome());
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

    public int getCdsStart() {
        return cdsStart;
    }

    public void setCdsStart(int cdsStart) {
        this.cdsStart = cdsStart;
    }

    public int getCdsStop() {
        return cdsStop;
    }

    public void setCdsStop(int cdsStop) {
        this.cdsStop = cdsStop;
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