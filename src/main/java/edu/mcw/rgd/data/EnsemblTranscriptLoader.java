package edu.mcw.rgd.data;

import edu.mcw.rgd.dao.impl.MapDAO;
import edu.mcw.rgd.dao.impl.TranscriptDAO;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.process.mapping.MapManager;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.Map;

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
    public void run(Collection<EnsemblTranscript> transcripts,int speciesTypeKey) throws Exception
    {
        System.out.println("Loading the transcripts file");
        int mapKey = 0;
        edu.mcw.rgd.datamodel.Map referenceAssembly = MapManager.getInstance().getReferenceAssembly(speciesTypeKey);
        List<String> chromosomes = ensemblDAO.getChromosomes(referenceAssembly.getKey());

        if(speciesTypeKey == SpeciesType.RAT || speciesTypeKey == SpeciesType.DOG || speciesTypeKey == SpeciesType.PIG || speciesTypeKey == SpeciesType.BONOBO
                || speciesTypeKey == SpeciesType.CHINCHILLA || speciesTypeKey == SpeciesType.SQUIRREL) {
            mapKey = referenceAssembly.getKey() + 1;
        }
        else if(speciesTypeKey == SpeciesType.MOUSE) {
            mapKey =  referenceAssembly.getKey() + 4;
        }
        else if(speciesTypeKey == SpeciesType.HUMAN){
            mapKey = referenceAssembly.getKey() + 2;
        }
        TranscriptDAO transcriptDAO = new TranscriptDAO();

        for (EnsemblTranscript transcript : transcripts) {
            if(chromosomes.contains(transcript.getChromosome())) {
                boolean transcriptMatch = false;
                String rgdId = ensemblDAO.getEnsemblRgdId(transcript.getEnsGeneId());
                if (rgdId != null) {
                    int geneRgdId = Integer.parseInt(ensemblDAO.getEnsemblRgdId(transcript.getEnsGeneId()));
                    List<EnsemblExon> features = transcript.getExonList();
                    List<Transcript> transcriptsForGene = transcriptDAO.getTranscriptsForGene(geneRgdId);
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
                        updateTranscriptData(transcript.getRgdId(), features, mapKey);
                        updateTranscriptType(transcript);
                    }else createNewEnsemblTranscript(transcript, mapKey);

                }else genesNotFound.add(transcript.getEnsGeneId());
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
   public void updateTranscriptData(int transcriptRgdId,List<EnsemblExon> features,int mapKey) throws Exception{
       int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();
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


    public void createNewEnsemblTranscript(EnsemblTranscript transcript, int mapKey) throws Exception {
        loaded.add(transcript);
        int speciesTypeKey = MapManager.getInstance().getMap(mapKey).getSpeciesTypeKey();
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

        ensemblDAO.insertTranscript(newTranscript,speciesTypeKey);

        updateTranscriptData(newTranscript.getRgdId(),transcript.getExonList(),mapKey);
    }

}
