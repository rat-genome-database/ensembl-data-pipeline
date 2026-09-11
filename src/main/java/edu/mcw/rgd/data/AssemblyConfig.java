package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.SpeciesType;

/// configuration of a single assembly processed by the pipeline: the species, the map keys the data
/// is loaded onto, and the Ensembl source files; see the 'assemblies' property in AppConfigure.xml
public class AssemblyConfig {

    private int speciesTypeKey;
    private int ensemblMapKey;    // map key the Ensembl positions are loaded onto
    private int ncbiMapKey;       // map key of the NCBI assembly used for gene-by-position matching; 0 = none;
                                  //   for strain assemblies it is the SAME map key as ensemblMapKey
    private String gff3Url;       // full Ensembl GFF3 source URL (gff3 loader)
    private String entrezUrl;     // full Ensembl entrez tsv source URL (gff3 loader)
    private String xrefAuthority; // 'RGD' | 'MGI' | 'HGNC' | null -- xref authority in the gff3 gene description
    private String assemblyName;  // optional override of the RGD map name, f.e. 'Naked_mole-rat_maternal Ensembl'

    /// true for strain assemblies: Ensembl and NCBI positions share one map key (SRC_PIPELINE tells them apart)
    public boolean isSharedMapKey() {
        return ensemblMapKey==ncbiMapKey;
    }

    public String describe() {
        return SpeciesType.getCommonName(speciesTypeKey)+" map_key "+ensemblMapKey
            + (isSharedMapKey() ? " (shared with NCBI)" : ncbiMapKey>0 ? " (NCBI map_key "+ncbiMapKey+")" : " (no NCBI map)");
    }

    public int getSpeciesTypeKey() {
        return speciesTypeKey;
    }

    public void setSpeciesTypeKey(int speciesTypeKey) {
        this.speciesTypeKey = speciesTypeKey;
    }

    public int getEnsemblMapKey() {
        return ensemblMapKey;
    }

    public void setEnsemblMapKey(int ensemblMapKey) {
        this.ensemblMapKey = ensemblMapKey;
    }

    public int getNcbiMapKey() {
        return ncbiMapKey;
    }

    public void setNcbiMapKey(int ncbiMapKey) {
        this.ncbiMapKey = ncbiMapKey;
    }

    public String getGff3Url() {
        return gff3Url;
    }

    public void setGff3Url(String gff3Url) {
        this.gff3Url = gff3Url;
    }

    public String getEntrezUrl() {
        return entrezUrl;
    }

    public void setEntrezUrl(String entrezUrl) {
        this.entrezUrl = entrezUrl;
    }

    public String getXrefAuthority() {
        return xrefAuthority;
    }

    public void setXrefAuthority(String xrefAuthority) {
        this.xrefAuthority = xrefAuthority;
    }

    public String getAssemblyName() {
        return assemblyName;
    }

    public void setAssemblyName(String assemblyName) {
        this.assemblyName = assemblyName;
    }
}
