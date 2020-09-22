package edu.mcw.rgd.data;



import edu.mcw.rgd.datamodel.TranscriptFeature;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by hsnalabolu on 12/4/2019.
 */
public class EnsemblTranscript {

    int rgdId;
    String ensTranscriptId;
    String ensTranscriptVer;
    String ensGeneId;
    String proteinId;
    String chromosome;
    int start;
    int stop;
    String strand;
    boolean nonCodingInd;
    String type;

    ArrayList<EnsemblExon> exonList = new ArrayList<EnsemblExon>();

    public List<TranscriptFeature> getUtrs() {
        int cdsStart = 0;
        int cdsStop = 0;
        for( EnsemblExon exon: exonList ) {
            if( exon.getCdsStart()>0 ) {
                if( cdsStart==0 || exon.getCdsStart()<cdsStart ) {
                    cdsStart = exon.getCdsStart();
                }
            }
            if( exon.getCdsStop()>0 ) {
                if( cdsStop==0 || exon.getCdsStop()>cdsStop ) {
                    cdsStop = exon.getCdsStop();
                }
            }
        }

        if( cdsStart==0 || cdsStop==0 ) {
            return null; // no CDS
        }
        if( cdsStart==start && cdsStop==stop ) {
            return null; // no CDS
        }

        // there are some UTRs
        List<TranscriptFeature> utrs = new ArrayList<>();
        if( getStrand().equals("+") ) {
            if( cdsStart>start ) {
                TranscriptFeature tf = new TranscriptFeature();
                tf.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                tf.setChromosome(chromosome);
                tf.setStrand(strand);
                tf.setStartPos(start);
                tf.setStopPos(cdsStart-1);
                utrs.add(tf);
            }

            if( cdsStop<stop ) {
                TranscriptFeature tf = new TranscriptFeature();
                tf.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                tf.setChromosome(chromosome);
                tf.setStrand(strand);
                tf.setStartPos(cdsStop+1);
                tf.setStopPos(stop);
                utrs.add(tf);
            }
        }
        else if( getStrand().equals("-") ) {
            if( cdsStart>start ) {
                TranscriptFeature tf = new TranscriptFeature();
                tf.setFeatureType(TranscriptFeature.FeatureType.UTR3);
                tf.setChromosome(chromosome);
                tf.setStrand(strand);
                tf.setStartPos(start);
                tf.setStopPos(cdsStart-1);
                utrs.add(tf);
            }

            if( cdsStop<stop ) {
                TranscriptFeature tf = new TranscriptFeature();
                tf.setFeatureType(TranscriptFeature.FeatureType.UTR5);
                tf.setChromosome(chromosome);
                tf.setStrand(strand);
                tf.setStartPos(cdsStop+1);
                tf.setStopPos(stop);
                utrs.add(tf);
            }
        } else {
            System.out.println("unknown strand");
        }
        return utrs;
    }

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

    public String getEnsTranscriptVer() {
        return ensTranscriptVer;
    }

    public void setEnsTranscriptVer(String ensTranscriptVer) {
        this.ensTranscriptVer = ensTranscriptVer;
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
