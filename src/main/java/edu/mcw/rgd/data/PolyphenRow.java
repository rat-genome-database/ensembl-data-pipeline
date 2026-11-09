package edu.mcw.rgd.data;

import java.util.Objects;

/// one row of the POLYPHEN table: PolyPhen-2 prediction for a variant on a transcript (only the columns
/// the Ensembl-sourced predictions provide)
public class PolyphenRow {

    long polyphenId;       // 0 for rows not in the database yet
    long variantRgdId;
    long transcriptRgdId;
    String geneSymbol;
    String proteinId;
    int position;          // 1-based position in the reference protein
    String aa1;            // reference amino acid
    String aa2;            // variant amino acid
    String prediction;     // 'benign', 'possibly damaging', 'probably damaging', 'unknown'
    String pph2Class;      // 'neutral', 'deleterious', 'none'
    String pph2Prob;       // score 0..1, null for 'unknown'

    /// unique key of a prediction: the same variant/transcript pair may have different reference proteins
    /// on GRCh37 and GRCh38
    public String getKey() {
        return variantRgdId+"|"+transcriptRgdId+"|"+proteinId+"|"+position+"|"+aa1+"|"+aa2;
    }

    public boolean hasSameData(PolyphenRow r) {
        return Objects.equals(geneSymbol, r.geneSymbol)
            && Objects.equals(prediction, r.prediction)
            && Objects.equals(pph2Class, r.pph2Class)
            && Objects.equals(pph2Prob, r.pph2Prob);
    }

    public String dump() {
        return "POLYPHEN_ID="+polyphenId+" VARIANT_RGD_ID="+variantRgdId+" TRANSCRIPT_RGD_ID="+transcriptRgdId
            +" GENE="+geneSymbol+" PROTEIN="+proteinId+" "+aa1+position+aa2
            +" PREDICTION="+prediction+" PPH2_PROB="+pph2Prob;
    }
}
