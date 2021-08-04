package edu.mcw.rgd.data;

import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.Transcript;
import edu.mcw.rgd.process.Utils;
import oracle.ucp.util.MapChain;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

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

    private int ncbiAssemblyMapKey;
    private int ensemblAssemblyMapKey;

    private String gff3File;
    private String xrefFile;

    public List<EnsemblGene> parseGenes() throws Exception {

        // rn7
        setNcbiAssemblyMapKey(372);
        setEnsemblAssemblyMapKey(373);

        // load primary information from gff3 genes file
        Map<String, EnsemblGene> genes = parseGff3File();
        System.out.println("genes loaded from gff3 file");

        // load supplemental information from xref file
        parseXrefFile(genes);

        /*

        BufferedReader in = Utils.openReader(getSrcFile());
        String line;
        while( (line=in.readLine())!=null ) {

        }

        in.close();

         */
        return new ArrayList<>(genes.values());
    }

    Map<String, EnsemblGene> parseGff3File() throws IOException {

        Map<String, EnsemblGene> results = new HashMap<>();

        String genomeBuild = null;
        if( getNcbiAssemblyMapKey()==372 || getEnsemblAssemblyMapKey()==373 ) {
            genomeBuild = "mRatBN7.2";
        }
        boolean genomeBuildVerified = false;

        BufferedReader in = Utils.openReader(getGff3File());
        String line;
        while( (line=in.readLine())!=null ) {
            // header parsing
            if( line.startsWith("#") ) {
                if( line.contains("genome-build") ) {
                    if( line.contains(genomeBuild) ) {
                        genomeBuildVerified = true;
                    }
                }
                continue;
            }

            // data lines
            if( !genomeBuildVerified ) {
                // assembly mismatch
                System.out.println("assembly "+genomeBuild+" expected! not found");
                break;
            }

            // skip non-gene lines
            // Y	ensembl	gene	7732918	7746199	.	+	.	ID=gene:ENSRNOG00000065605;biotype=protein_coding;gene_id=ENSRNOG00000065605;version=1
            String[] cols = line.split("[\\t]", -1);
            String info = cols[8];
            if( !info.startsWith("ID=gene:ENSRNOG") ) {
                continue;
            }
            String chr = cols[0];
            String startPosStr = cols[3];
            String stopPosStr = cols[4];
            String strand = cols[6];
            String geneId = null, bioType = null;

            // parse gene data
            String[] infos = info.split("[\\;]");
            for( String inf: infos ) {
                if( inf.startsWith("biotype=") ) {
                    bioType = inf.substring(8);
                } else if( inf.startsWith("gene_id=") ) {
                    geneId = inf.substring(8);
                }
            }
            if( geneId==null || bioType==null ) {
                System.out.println("unexpected gene line: "+line);
                continue;
            }

            EnsemblGene g = new EnsemblGene();
            g.setChromosome(chr);
            g.setStartPos(startPosStr);
            g.setStopPos(stopPosStr);
            g.setStrand(strand);
            g.setgene_biotype(bioType);
            g.setEnsemblGeneId(geneId);
            results.put(geneId, g);
        }
        in.close();

        return results;
    }

    void parseXrefFile( Map<String, EnsemblGene> genes ) throws Exception {

        // gene_stable_id	transcript_stable_id	protein_stable_id	xref_id	xref_label	description	db_name	info_type	source	ensembl_identity	xref_identity
        //
        // we will parse only these two lines
        //ENSRNOG00000000001	ENSRNOT00000055633	ENSRNOP00000092332	Arsj-201	Arsj-201	arylsulfatase family, member J	RGD transcript name	MISC
        //ENSRNOG00000000001	ENSRNOT00000055633	ENSRNOP00000092332	NM_001047887	NM_001047887.1		RefSeq mRNA	DIRECT		50	100


        BufferedReader in = Utils.openReader(fixXrefFile());
        String header = in.readLine();
        String line;
        while( (line=in.readLine())!=null ) {

            String[] cols = line.split("[\\t]", -1);
            if( cols.length!=11 ) {
                continue;
            }
            String geneId = cols[0];
            EnsemblGene gene = genes.get(geneId);

            String xrefId = cols[3];
            String desc = cols[5];
            String dbName = cols[6];

            if( dbName.equals("RGD transcript name") ) {

                String geneSymbol = xrefId;
                // remove version from symbol
                int dashPos = geneSymbol.indexOf("-");
                if( dashPos>0 ) {
                    geneSymbol = geneSymbol.substring(0, dashPos);
                }

                gene.setGeneSymbol(geneSymbol);
                gene.setgene_description(desc);
            }
            else if( dbName.equals("RefSeq mRNA") ) {
                Set<String> refseqAccIds = gene.getRefseqAccIds();
                if( refseqAccIds==null ) {
                    refseqAccIds = new HashSet<>();
                    gene.setRefseqAccIds(refseqAccIds);
                }
                refseqAccIds.add(xrefId);
            }
        }
        in.close();

        resolveRefSeqAccIds(genes);
    }

    void resolveRefSeqAccIds(Map<String, EnsemblGene> geneMap) throws Exception {

        EnsemblDAO dao = new EnsemblDAO();
        List<EnsemblGene> genes = new ArrayList<>(geneMap.values());
        Collections.shuffle(genes);
        AtomicInteger conflictCount = new AtomicInteger();

        // map transcript acc ids to EntrezGene ids
        for( EnsemblGene g: genes ) {
            boolean conflict = false;
            int geneRgdId = 0;

            if( g.getRefseqAccIds()!=null ) {
                Set<Integer> geneRgdIds = new HashSet<>();
                for (String acc : g.getRefseqAccIds()) {
                    List<Transcript> trList = dao.getTranscriptsByAccId(acc);
                    for (Transcript tr : trList) {
                        geneRgdIds.add(tr.getGeneRgdId());
                    }
                }
                if( geneRgdIds.size()>1 ) {
                    // filter out inactive rgd ids
                    for( int incomingGeneRgdId: geneRgdIds ) {
                        RgdId id = dao.getRgdId(incomingGeneRgdId);
                        if( id.getObjectStatus().equals("ACTIVE") ) {
                            if( geneRgdId==0 ) {
                                geneRgdId = id.getRgdId();
                            } else {
                                if( geneRgdId!=incomingGeneRgdId ) {
                                    conflict = true;
                                    geneRgdId = 0;
                                    break;
                                }
                            }
                        }
                    }
                } else if( geneRgdIds.size()==1 ) {
                    geneRgdId = geneRgdIds.iterator().next();
                }

                if (geneRgdId != 0 && !conflict) {
                    g.setrgdid(Integer.toString(geneRgdId));
                }
                if (conflict) {
                    conflictCount.incrementAndGet();
                }
            }
        }

        Logger statuslog = Logger.getLogger("status");
        statuslog.info("resolveRefSeqAccIds(): conflicts: "+conflictCount.toString());
    }

    // Ensembl publishes xrefs file *without new lines*;
    // this fix adds new lines so we can easily parse the file
    String fixXrefFile() throws IOException {
        String input_fname = getXrefFile();
        String output_fname = "data/rapid/xref_fixed.tsv.gz";
        BufferedWriter out = Utils.openWriter(output_fname);
        BufferedReader in = Utils.openReader(input_fname);
        String line = in.readLine();
        System.out.println("line len="+line.length());

        int writePos = 0;
        int pos = 1;
        while( (pos=line.indexOf("ENSRNOG000", pos))>0 ) {
            out.write(line.substring(writePos, pos));
            out.write("\n");

            writePos=pos;
            pos++;
        }
        // write leftover
        out.write(line.substring(writePos));
        out.write("\n");
        out.close();
        in.close();

        return output_fname;
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

    public int getNcbiAssemblyMapKey() {
        return ncbiAssemblyMapKey;
    }

    public void setNcbiAssemblyMapKey(int ncbiAssemblyMapKey) {
        this.ncbiAssemblyMapKey = ncbiAssemblyMapKey;
    }

    public int getEnsemblAssemblyMapKey() {
        return ensemblAssemblyMapKey;
    }

    public void setEnsemblAssemblyMapKey(int ensemblAssemblyMapKey) {
        this.ensemblAssemblyMapKey = ensemblAssemblyMapKey;
    }
}
