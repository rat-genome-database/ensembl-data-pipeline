package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.CounterPool;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Created by hsnalabolu on 12/4/2019.
 */
public class EnsemblTranscriptLoader {

    EnsemblDAO ensemblDAO = new EnsemblDAO();

    Logger log = Logger.getLogger("status");

    List genesNotFound=new ArrayList();

    public void run(Collection<EnsemblTranscript> transcriptCollection, int speciesTypeKey, int ensemblMapKey) throws Exception {

        log.debug("Loading the transcripts file");
        CounterPool counters = new CounterPool();

        // we have chromosome data only for NCBI assemblies
        edu.mcw.rgd.datamodel.Map referenceAssembly = MapManager.getInstance().getReferenceAssembly(speciesTypeKey);
        List<Chromosome> chromosomes = ensemblDAO.getChromosomes(referenceAssembly.getKey());

        List<EnsemblTranscript> transcripts = new ArrayList<>(transcriptCollection);
        Collections.shuffle(transcripts);

        int t=0;
        for (EnsemblTranscript transcript : transcripts) {
            log.debug((++t)+". "+transcript.getEnsTranscriptId());
            counters.increment("TRANSCRIPTS_INCOMING");

            String chr = ensemblDAO.matchChromosome(transcript.getChromosome(), chromosomes);
            if( chr==null ) {
                counters.increment("TRANSCRIPTS_SKIPPED_UNEXPECTED_CHROMOSOME");
                log.debug("transcript_skipped: unexpected chromosome "+transcript.getChromosome());
                continue;
            }
            List<EnsemblExon> exons = transcript.getExonList();

            if( !transcript.getChromosome().equals(chr) ) {
                transcript.setChromosome(chr); // replace chr GeneBank id with NCBI scaffold acc

                // ensure exons have the same chr as transcript
                for (EnsemblExon exon : exons) {
                    exon.setExonChromosome(chr);
                }
            }

            boolean transcriptMatch = false;
            String rgdId = ensemblDAO.getEnsemblRgdId(transcript.getEnsGeneId());
            if (rgdId != null) {

                TranscriptVersionManager.getInstance().addVersion(transcript.getEnsTranscriptId(), transcript.getEnsTranscriptVer());

                List<TranscriptFeature> utrs = transcript.getUtrs();

                int geneRgdId = Integer.parseInt(ensemblDAO.getEnsemblRgdId(transcript.getEnsGeneId()));
                List<Transcript> transcriptsForGene = ensemblDAO.getTranscriptsForGene(geneRgdId);
                for (Transcript oldTranscript : transcriptsForGene) {
                    if (oldTranscript.getAccId().equalsIgnoreCase(transcript.getEnsTranscriptId()) &&
                            oldTranscript.isNonCoding() == transcript.isNonCodingInd()) {

                        transcriptMatch = true;
                        transcript.setRgdId(oldTranscript.getRgdId());
                        List<TranscriptFeature> tfsInRgd = ensemblDAO.getFeatures(oldTranscript.getRgdId());
                        List<TranscriptFeature> obsoleteExons = new ArrayList<>();

                        List<EnsemblExon> matchedPositions = new ArrayList<>();
                        for (TranscriptFeature oldPos : tfsInRgd) {
                            boolean matchFound = false;
                            for (EnsemblExon exon: exons) {
                                if( exon.matchesByPos(oldPos) && oldPos.getFeatureType()==TranscriptFeature.FeatureType.EXON) {
                                    matchedPositions.add(exon);
                                    matchFound = true;
                                    break;
                                }
                            }
                            if( !matchFound && oldPos.getFeatureType()==TranscriptFeature.FeatureType.EXON ) {
                                obsoleteExons.add(oldPos);
                            }
                        }
                        counters.add("TRANSCRIPT_EXONS_MATCHED", matchedPositions.size());
                        exons.removeAll(matchedPositions);

                        deleteExons(oldTranscript, obsoleteExons, counters);

                        qcUtrs(utrs, tfsInRgd, counters, oldTranscript, ensemblMapKey, speciesTypeKey);
                    }
                }


                if (transcriptMatch) {
                    insertExons(transcript.getRgdId(), exons, ensemblMapKey, speciesTypeKey, counters);
                    updateTranscriptType(transcript);
                } else {
                    createNewEnsemblTranscript(transcript, ensemblMapKey, speciesTypeKey, counters);
                }

            } else {
                genesNotFound.add(transcript.getEnsGeneId());
                counters.increment("TRANSCRIPTS_WITHOUT_MATCHING_GENE");
            }
        }

        TranscriptVersionManager.getInstance().qcAndLoad(counters);
        log.info(counters.dumpAlphabetically());
    }

    void qcUtrs(List<TranscriptFeature> utrs, List<TranscriptFeature> tfsInRgd, CounterPool counters, Transcript tr, int ensemblMapKey, int speciesTypeKey) throws Exception {

        List<TranscriptFeature> utrsInRgd = new ArrayList<>();
        for( TranscriptFeature tf: tfsInRgd ) {
            if( tf.getFeatureType()==TranscriptFeature.FeatureType.UTR3 || tf.getFeatureType()==TranscriptFeature.FeatureType.UTR5 ) {
                utrsInRgd.add(tf);
            }
        }

        // qc utrs
        if( utrs!=null ) {

            // remove matching utrs from both incoming utrs and utrs in rgd
            boolean matchingUtrFound;
            do {
                matchingUtrFound = false;
                for (int u = 0; u < utrsInRgd.size(); u++) {
                    TranscriptFeature utrInRgd = utrsInRgd.get(u);
                    for (int i = 0; i < utrs.size(); i++) {
                        TranscriptFeature utr = utrs.get(i);
                        if (utr.getFeatureType() == utrInRgd.getFeatureType()
                                && utr.getStartPos().equals(utrInRgd.getStartPos())
                                && utr.getStopPos().equals(utrInRgd.getStopPos())) {

                            utrs.remove(i);
                            utrsInRgd.remove(u);
                            counters.increment("TRANSCRIPT_UTRS_MATCHED");
                            matchingUtrFound = true;
                            break;
                        }
                    }
                    if( matchingUtrFound ) {
                        break;
                    }
                }
            } while(matchingUtrFound);

            // remaining entries in 'utrs' must be inserted
            insertUtrs(utrs, tr.getRgdId(), ensemblMapKey, speciesTypeKey, counters);
        }

        // remaining entries in 'utrsInRgd' table must be deleted
        removeUtrs(tr, utrsInRgd, counters);
    }

    public void updateTranscriptType(EnsemblTranscript transcript) throws Exception{
        Transcript t = ensemblDAO.getTranscript(transcript.getRgdId());
        t.setType(transcript.getType());
        ensemblDAO.updateTranscript(t);

        TranscriptVersionManager.getInstance().addRgdId(t.getAccId(), t.getRgdId());
    }

    public void insertExons(int transcriptRgdId, List<EnsemblExon> exons, int mapKey, int speciesTypeKey, CounterPool counters) throws Exception{
        for(EnsemblExon exon: exons) {
            TranscriptFeature newFeature = new TranscriptFeature();
            newFeature.setTranscriptRgdId(transcriptRgdId);
            newFeature.setChromosome(exon.getExonChromosome());
            newFeature.setStartPos(exon.getExonStart());
            newFeature.setStopPos(exon.getExonStop());
            newFeature.setStrand(exon.getStrand());
            newFeature.setSrcPipeline("Ensembl");
            newFeature.setMapKey(mapKey);
            newFeature.setFeatureType(TranscriptFeature.FeatureType.EXON);
            ensemblDAO.insertTranscriptFeature(newFeature,speciesTypeKey);
            counters.increment("TRANSCRIPT_EXONS_INSERTED");
        }
    }

    public void deleteExons(Transcript tr, List<TranscriptFeature> exons, CounterPool counters) throws Exception{

        for(TranscriptFeature exon: exons) {

            // sanity check
            if( exon.getFeatureType()!=TranscriptFeature.FeatureType.EXON ) {
                throw new Exception("deleteExons() can process only EXONs!");
            }

            ensemblDAO.unbindFeatureFromTranscript(tr, exon);
            counters.increment("TRANSCRIPT_EXONS_REMOVED");
        }
    }
    public void removeUtrs(Transcript tr, List<TranscriptFeature> utrsInRgd, CounterPool counters) throws Exception{

        for(TranscriptFeature utrInRgd: utrsInRgd) {

            // sanity check
            if( utrInRgd.getFeatureType()==TranscriptFeature.FeatureType.UTR3 ) {
                ensemblDAO.unbindFeatureFromTranscript(tr, utrInRgd);
                counters.increment("TRANSCRIPT_UTR3_REMOVED");
            }
            else if( utrInRgd.getFeatureType()==TranscriptFeature.FeatureType.UTR5 ) {
                ensemblDAO.unbindFeatureFromTranscript(tr, utrInRgd);
                counters.increment("TRANSCRIPT_UTR5_REMOVED");
            }
            else {
                throw new Exception("removeUtrs() can process only UTRs!");
            }
        }
    }

    public void insertUtrs(List<TranscriptFeature> utrs, int transcriptRgdId, int mapKey, int speciesTypeKey, CounterPool counters) throws Exception {

        if( utrs!=null)
        for(TranscriptFeature utr: utrs) {
            TranscriptFeature newFeature = new TranscriptFeature();
            newFeature.setTranscriptRgdId(transcriptRgdId);
            newFeature.setChromosome(utr.getChromosome());
            newFeature.setStartPos(utr.getStartPos());
            newFeature.setStopPos(utr.getStopPos());
            newFeature.setStrand(utr.getStrand());
            newFeature.setSrcPipeline("Ensembl");
            newFeature.setMapKey(mapKey);
            newFeature.setFeatureType(utr.getFeatureType());
            ensemblDAO.insertTranscriptFeature(newFeature,speciesTypeKey);
            counters.increment("TRANSCRIPT_UTRS_INSERTED");
        }
    }

    public void createNewEnsemblTranscript(EnsemblTranscript transcript, int mapKey, int speciesTypeKey, CounterPool counters) throws Exception {

        Transcript newTranscript = new Transcript();
        newTranscript.setAccId(transcript.getEnsTranscriptId());
        int geneRgdId = Integer.parseInt(ensemblDAO.getEnsemblRgdId(transcript.getEnsGeneId()));
        newTranscript.setGeneRgdId(geneRgdId);
        newTranscript.setNonCoding(transcript.isNonCodingInd());
        newTranscript.setProteinAccId(transcript.getProteinId());
        newTranscript.setType(transcript.getType());

        MapData mapData = new MapData();
        mapData.setSrcPipeline("Ensembl");
        mapData.setChromosome(transcript.getChromosome());
        mapData.setMapKey(mapKey);
        mapData.setStartPos(transcript.getStart());
        mapData.setStopPos(transcript.getStop());
        mapData.setStrand(transcript.getStrand());

        List<MapData> genomicPos = new ArrayList<>();
        genomicPos.add(mapData);

        newTranscript.setGenomicPositions(genomicPos);

        ensemblDAO.insertTranscript(newTranscript, speciesTypeKey);
        counters.increment("TRANSCRIPTS_INSERTED");

        insertExons(newTranscript.getRgdId(), transcript.getExonList(), mapKey, speciesTypeKey, counters);

        List<TranscriptFeature> utrs = transcript.getUtrs();
        insertUtrs(utrs, newTranscript.getRgdId(), mapKey, speciesTypeKey, counters);
    }
}
