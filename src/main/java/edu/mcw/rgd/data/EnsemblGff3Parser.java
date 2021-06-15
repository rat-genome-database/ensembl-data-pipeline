package edu.mcw.rgd.data;

import edu.mcw.rgd.process.Utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Ensembl Gff3 file has genes separated by '###' lines
 * ###
 * Y	ensembl	gene	7732918	7746199	.	+	.	ID=gene:ENSRNOG00000065605;biotype=protein_coding;gene_id=ENSRNOG00000065605;version=1
 * Y	ensembl	mRNA	7732918	7746199	.	+	.	ID=transcript:ENSRNOT00000096124;Parent=gene:ENSRNOG00000065605;biotype=protein_coding;transcript_id=ENSRNOT00000096124;version=1
 * Y	ensembl	exon	7732918	7732964	.	+	.	Parent=transcript:ENSRNOT00000096124;constitutive=1;exon_id=ENSRNOE00000621907;rank=1;version=1
 * Y	ensembl	five_prime_UTR	7732918	7732964	.	+	.	Parent=transcript:ENSRNOT00000096124
 * Y	ensembl	five_prime_UTR	7745762	7745823	.	+	.	Parent=transcript:ENSRNOT00000096124
 * Y	ensembl	exon	7745762	7746199	.	+	.	Parent=transcript:ENSRNOT00000096124;constitutive=1;exon_id=ENSRNOE00000651328;rank=2;version=1
 * Y	ensembl	CDS	7745824	7746039	.	+	0	ID=CDS:ENSRNOP00000080992;Parent=transcript:ENSRNOT00000096124;protein_id=ENSRNOP00000080992;version=1
 * Y	ensembl	three_prime_UTR	7746040	7746199	.	+	.	Parent=transcript:ENSRNOT00000096124
 * ###
 */
public class EnsemblGff3Parser {

    private String gff3File;
    private String xrefFile;

    public List<EnsemblGene> parseGenes() throws IOException {

        final int mapKey = 372; // rn7_2
/*
        // load primary information from gff3 genes file
        List<EnsemblGene> genes = parseGff3File();

        // load supplemental information from xref file
        parseXrefFile(genes);

        BufferedReader in = Utils.openReader(getSrcFile());
        String line;
        while( (line=in.readLine())!=null ) {

        }

        in.close();
        return genes;
*/
        return null;
    }


    public String getGff3File() {
        return gff3File;
    }

    public void setGff3File(String gff3File) {
        this.gff3File = gff3File;
    }

    public String getXrefFile() {
        return xrefFile;
    }

    public void setXrefFile(String xrefFile) {
        this.xrefFile = xrefFile;
    }
}
