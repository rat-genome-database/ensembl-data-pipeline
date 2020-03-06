package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * Created by hsnalabolu on 12/4/2019.
 */
public class EnsemblTranscriptLoader {

    EnsemblDAO ensemblDAO;

    static Logger statuslog = Logger.getLogger("statuscheck");


    List loaded=new ArrayList();
    List genesNotFound=new ArrayList();

    public EnsemblTranscriptLoader() throws Exception {
        ensemblDAO = new EnsemblDAO();
    }

    public void run(Collection<EnsemblTranscript> transcripts, int speciesTypeKey, int ensemblMapKey) throws Exception {

        System.out.println("Loading the transcripts file");

        // we have chromosome data only for NCBI assemblies
        edu.mcw.rgd.datamodel.Map referenceAssembly = MapManager.getInstance().getReferenceAssembly(speciesTypeKey);
        List<String> chromosomes = ensemblDAO.getChromosomes(referenceAssembly.getKey());

        for (EnsemblTranscript transcript : transcripts) {
            if(chromosomes.contains(transcript.getChromosome())) {
                boolean transcriptMatch = false;
                String rgdId = ensemblDAO.getEnsemblRgdId(transcript.getEnsGeneId());
                if (rgdId != null) {
                    int geneRgdId = Integer.parseInt(ensemblDAO.getEnsemblRgdId(transcript.getEnsGeneId()));
                    List<EnsemblExon> features = transcript.getExonList();
                    List<Transcript> transcriptsForGene = ensemblDAO.getTranscriptsForGene(geneRgdId);
                    for (Transcript oldTranscript : transcriptsForGene) {
                        if (oldTranscript.getAccId().equalsIgnoreCase(transcript.getEnsTranscriptId()) &&
                                oldTranscript.isNonCoding() == transcript.isNonCodingInd()) {
                            transcriptMatch = true;
                            transcript.setRgdId(oldTranscript.getRgdId());
                            List<TranscriptFeature> positions = ensemblDAO.getFeatures(oldTranscript.getRgdId());
                            List<EnsemblExon> matchedPositions = new ArrayList<>();
                            for (TranscriptFeature oldPos : positions) {
                                for (EnsemblExon feature : features) {
                                    if (oldPos.getChromosome().equalsIgnoreCase(feature.getExonChromosome()) && oldPos.getStartPos() == feature.getExonStart() && oldPos.getStopPos() == feature.getExonStop()
                                            && oldPos.getStrand().equalsIgnoreCase(feature.getStrand())) {
                                     
                                        matchedPositions.add(feature);
                                    }
                                }
                            }

                            features.removeAll(matchedPositions);
                        }
                    }
                    if (transcriptMatch) {
                        updateTranscriptData(transcript.getRgdId(), features, ensemblMapKey, speciesTypeKey);
                        updateTranscriptType(transcript);
                    } else {
                        createNewEnsemblTranscript(transcript, ensemblMapKey, speciesTypeKey);
                    }

                } else {
                    genesNotFound.add(transcript.getEnsGeneId());
                }
            }
        }

        statuslog.info("Total loaded: "+loaded.size()+"\n");
        statuslog.info("Total genes not found: "+genesNotFound.size()+"\n");
        statuslog.info("Total in file: "+ transcripts.size()+"\n");

    }

    public void updateTranscriptType(EnsemblTranscript transcript) throws Exception{
        Transcript t = ensemblDAO.getTranscript(transcript.getRgdId());
        t.setType(transcript.getType());
        ensemblDAO.updateTranscript(t);
    }

    public void updateTranscriptData(int transcriptRgdId, List<EnsemblExon> features, int mapKey, int speciesTypeKey) throws Exception{
        for(EnsemblExon feature: features) {
            TranscriptFeature newFeature = new TranscriptFeature();
            newFeature.setTranscriptRgdId(transcriptRgdId);
            newFeature.setChromosome(feature.getExonChromosome());
            newFeature.setStartPos(feature.getExonStart());
            newFeature.setStopPos(feature.getExonStop());
            newFeature.setStrand(feature.getStrand());
            newFeature.setSrcPipeline("Ensembl");
            newFeature.setMapKey(mapKey);
            newFeature.setFeatureType(TranscriptFeature.FeatureType.EXON);
            ensemblDAO.insertTranscriptFeature(newFeature,speciesTypeKey);
        }
    }

    public void createNewEnsemblTranscript(EnsemblTranscript transcript, int mapKey, int speciesTypeKey) throws Exception {
        loaded.add(transcript);
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

        updateTranscriptData(newTranscript.getRgdId(), transcript.getExonList(), mapKey, speciesTypeKey);
    }
}
