package edu.mcw.rgd.data;

/// Ensembl variation database dump files providing PolyPhen-2 predictions for one RGD assembly;
/// see the 'sources' property of the 'polyphenLoader' bean in AppConfigure.xml
public class PolyphenSource {

    private int mapKey;               // RGD map key of the variants, f.e. 38 (GRCh38)
    private String attribUrl;         // attrib.txt[.gz] -- resolves the analysis attrib id of the PolyPhen model
    private String translationMd5Url; // translation_md5.txt[.gz] -- md5 of protein sequence -> translation_md5_id
    private String predictionsUrl;    // protein_function_predictions.txt[.gz] -- gzipped prediction matrix per protein and analysis

    public int getMapKey() {
        return mapKey;
    }

    public void setMapKey(int mapKey) {
        this.mapKey = mapKey;
    }

    public String getAttribUrl() {
        return attribUrl;
    }

    public void setAttribUrl(String attribUrl) {
        this.attribUrl = attribUrl;
    }

    public String getTranslationMd5Url() {
        return translationMd5Url;
    }

    public void setTranslationMd5Url(String translationMd5Url) {
        this.translationMd5Url = translationMd5Url;
    }

    public String getPredictionsUrl() {
        return predictionsUrl;
    }

    public void setPredictionsUrl(String predictionsUrl) {
        this.predictionsUrl = predictionsUrl;
    }
}
