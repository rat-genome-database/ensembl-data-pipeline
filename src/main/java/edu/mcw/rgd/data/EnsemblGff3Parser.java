package edu.mcw.rgd.data;

import edu.mcw.rgd.process.Utils;

import java.io.*;
import java.util.*;

/**
 * Parses Ensembl main-release GFF3 files (one species per file) into EnsemblGene objects and an
 * exon-level transcript file. The Entrez (NCBI) gene id is joined from the companion
 * '<species>.<assembly>.<release>.entrez.tsv.gz' file; the species xref (RGD/MGI/HGNC id) is taken
 * from the gene 'description' attribute. The produced EnsemblGene objects and transcript file match
 * what the BioMart path (Parser) produces -- verified field-by-field against prod BioMart output.
 *
 * Sample gene line (attributes in column 9, ';'-separated; the description's own ';' is %3B-encoded):
 * 1  ensembl  gene  78333971  78342685  .  +  .  ID=gene:ENSRNOG00000038600;Name=Dnaaf3;biotype=protein_coding;description=dynein%2C axonemal%2C assembly factor 3 [Source:RGD Symbol%3BAcc:2323487];gene_id=ENSRNOG00000038600;version=6
 */
public class EnsemblGff3Parser {

    private String genomeBuild;       // Ensembl assembly, f.e. 'GRCr8' -- must appear on the '#!genome-build' header line
    private String ensemblGenePrefix; // f.e. 'ENSRNOG'
    private String xrefAuthority;     // species xref authority in the description: 'RGD' | 'MGI' | 'HGNC' | null
    private int ncbiAssemblyMapKey;
    private int ensemblAssemblyMapKey;

    private String gff3File;
    private String entrezFile;

    public EnsemblGff3Parser() {
    }

    public List<EnsemblGene> parseGenes() throws Exception {

        Map<String, EnsemblGene> genes = parseGff3File();
        System.out.println("genes loaded from gff3 file: "+genes.size());

        // join NCBI (Entrez) gene ids from the companion entrez.tsv file
        parseEntrezFile(genes);

        return new ArrayList<>(genes.values());
    }

    Map<String, EnsemblGene> parseGff3File() throws IOException {

        Map<String, EnsemblGene> results = new HashMap<>();
        boolean genomeBuildVerified = false;

        BufferedReader in = Utils.openReader(getGff3File());
        String line;
        while( (line=in.readLine())!=null ) {

            // header parsing -- '#!genome-build  GRCr8'
            if( line.startsWith("#") ) {
                if( line.contains("genome-build") && line.contains(genomeBuild) ) {
                    genomeBuildVerified = true;
                }
                continue;
            }

            if( !genomeBuildVerified ) {
                System.out.println("assembly "+genomeBuild+" expected! not found in gff3 header");
                break;
            }

            // only gene lines for our species' gene-id prefix
            String[] cols = line.split("[\\t]", -1);
            String info = cols[8];
            if( !info.startsWith("ID=gene:"+getEnsemblGenePrefix()) ) {
                continue;
            }

            String geneId = attr(info, "gene_id=");
            String bioType = attr(info, "biotype=");
            if( geneId==null || bioType==null ) {
                System.out.println("unexpected gene line: "+line);
                continue;
            }

            // BioMart only loads genes that have a symbol -- mirror that
            String symbol = attr(info, "Name=");
            if( Utils.isStringEmpty(symbol) ) {
                continue;
            }

            EnsemblGene g = new EnsemblGene();
            g.setEnsemblGeneId(geneId);
            g.setChromosome(cols[0]);
            g.setStartPos(cols[3]);
            g.setStopPos(cols[4]);
            g.setStrand(cols[6]); // '+' / '-' -- same as BioMart path
            g.setGeneBioType(bioType.equals("protein_coding") ? "protein-coding" : bioType);
            g.setGeneSymbol(gff3Unescape(symbol));

            // gene name = description with the trailing '[Source:...]' removed (matches Parser.parseGene)
            String description = attr(info, "description=");
            g.setGeneName(geneNameFromDescription(description));

            // species xref id (rat rgd_id / mouse mgi_id / human hgnc_id) from the description Source
            g.setrgdid(extractSpeciesXref(description));

            results.put(geneId, g);
        }
        in.close();

        return results;
    }

    void parseEntrezFile( Map<String, EnsemblGene> genes ) throws Exception {

        if( Utils.isStringEmpty(getEntrezFile()) ) {
            return;
        }

        // gene_stable_id  transcript_stable_id  protein_stable_id  xref  db_name  info_type  ...
        // ENSRNOG00000009523  ENSRNOT00000012734  ENSRNOP00000012734  308003  EntrezGene  DEPENDENT  ...
        BufferedReader in = Utils.openReader(getEntrezFile());
        in.readLine(); // header
        String line;
        while( (line=in.readLine())!=null ) {
            String[] cols = line.split("[\\t]", -1);
            if( cols.length<5 || !cols[4].equals("EntrezGene") ) {
                continue;
            }
            EnsemblGene gene = genes.get(cols[0]);
            if( gene!=null ) {
                gene.setEntrezGeneId(cols[3]);
            }
        }
        in.close();
    }

    String generateTranscriptFile() throws IOException {

        // geneId to List of lines constituting a gene
        Map<String, List<String[]>> genes = new HashMap<>();
        boolean genomeBuildVerified = false;

        BufferedReader in = Utils.openReader(getGff3File());
        String geneId = null;
        List<String[]> geneLines = new ArrayList<>();

        String line;
        while( (line=in.readLine())!=null ) {
            // new gene boundary
            if( line.equals("###") ) {
                if( geneId!=null )
                    genes.put(geneId, geneLines);
                geneId = null;
                geneLines = new ArrayList<>();
                continue;
            }

            // header parsing
            if( line.startsWith("#") ) {
                if( line.contains("genome-build") && line.contains(genomeBuild) ) {
                    genomeBuildVerified = true;
                }
                continue;
            }

            if( !genomeBuildVerified ) {
                System.out.println("assembly "+genomeBuild+" expected! not found in gff3 header");
                break;
            }

            String[] cols = line.split("[\\t]", -1);
            String info = cols[8];
            if( info.startsWith("ID=gene:"+getEnsemblGenePrefix()) ) {
                geneId = attr(info, "gene_id=");
                continue;
            }

            geneLines.add(cols);
        }
        in.close();

        if( !genomeBuildVerified ) {
            return null;
        }

        return generateTranscriptFile(genes);
    }

    String generateTranscriptFile(Map<String, List<String[]>> genes) throws IOException {

        List<String> trLines = new ArrayList<>();

        genes.entrySet().parallelStream().forEach(e -> {
            List<String> lines = processGene(e.getKey(), e.getValue());
            synchronized (trLines) {
                trLines.addAll(lines);
            }
        });

        Collections.sort(trLines);

        String fname = "data/transcripts_"+genomeBuild+".txt.gz";
        BufferedWriter out = Utils.openWriter(fname);
        for (String line : trLines) {
            out.write(line);
        }
        out.close();

        return fname;
    }

    List<String> processGene(String geneId, List<String[]> lines) {

        // pass 1: transcripts
        Map<String, Tr> trs = new HashMap<>();
        for( String[] cols: lines ) {
            String info = cols[8];
            if( info.startsWith("ID=transcript:") ) {
                Tr tr = new Tr();
                tr.trId = attr(info, "transcript_id=");
                tr.trType = attr(info, "biotype=");
                String ver = attr(info, "version=");
                tr.trVer = ver!=null ? "."+ver : "";
                tr.chr = cols[0];
                tr.trStart = cols[3];
                tr.trStop = cols[4];
                tr.strand = cols[6].equals("+") ? "1" : "-1";
                trs.put(tr.trId, tr);
            }
        }

        // pass 2: exons
        for( String[] cols: lines ) {
            if( cols[2].equals("exon") ) {
                String info = cols[8];
                String trId = attr(info, "Parent=transcript:");
                Tr tr = trs.get(trId);
                if( tr==null ) continue;
                Exon exon = new Exon();
                exon.rank = attr(info, "rank=");
                exon.exonStart = Integer.parseInt(cols[3]);
                exon.exonStop = Integer.parseInt(cols[4]);
                tr.exons.add(exon);
            }
        }

        // pass 3: CDSs -> map coding region onto its exon
        for( String[] cols: lines ) {
            if( cols[2].equals("CDS") ) {
                String info = cols[8];
                String trId = attr(info, "Parent=transcript:");
                String proteinId = attr(info, "protein_id=");
                Tr tr = trs.get(trId);
                if( tr==null ) continue;
                int cdsStart = Integer.parseInt(cols[3]);
                int cdsStop = Integer.parseInt(cols[4]);
                if( proteinId!=null ) {
                    tr.proteinId = proteinId;
                }
                for( Exon e: tr.exons ) {
                    if( e.exonStart<=cdsStart && cdsStop<=e.exonStop ) {
                        e.cdsStart = cdsStart;
                        e.cdsStop = cdsStop;
                        break;
                    }
                }
            }
        }

        List<String> result = new ArrayList<>();
        for( Tr tr: trs.values() ) {
            String prefix = geneId+"\t"+tr.trId+"\t"+tr.trId+tr.trVer+"\t"+tr.chr+"\t"+tr.trStart+"\t"+tr.trStop+"\t"+tr.strand+"\t";
            String suffix = "\t"+Utils.defaultString(tr.proteinId)+"\t"+tr.trType+"\n";
            for( Exon e: tr.exons ) {
                String l = prefix+e.exonStart+"\t"+e.exonStop+"\t"+e.rank+"\t";
                if( e.cdsStart!=0 )
                    l += e.cdsStart;
                l += "\t";
                if( e.cdsStop!=0 )
                    l += e.cdsStop;
                l += suffix;
                result.add(l);
            }
        }
        return result;
    }

    // --- helpers -------------------------------------------------------------

    /// extract the value of a ';'-separated GFF3 attribute, f.e. attr(info, "gene_id=")
    static String attr(String info, String key) {
        for( String part: info.split("[\\;]") ) {
            if( part.startsWith(key) ) {
                return part.substring(key.length());
            }
        }
        return null;
    }

    /// gene name = description text before the trailing '[Source:...]', GFF3-unescaped
    static String geneNameFromDescription(String description) {
        if( Utils.isStringEmpty(description) ) {
            return "";
        }
        int bracketPos = description.indexOf('[');
        String name = bracketPos>0 ? description.substring(0, bracketPos) : description;
        return gff3Unescape(name).trim();
    }

    /// pull the species xref id from 'description', f.e. '...[Source:RGD Symbol%3BAcc:2323487]' -> '2323487'
    String extractSpeciesXref(String description) {
        if( Utils.isStringEmpty(xrefAuthority) || Utils.isStringEmpty(description) ) {
            return "0";
        }
        String marker = xrefAuthority+" Symbol%3BAcc:";
        int p = description.indexOf(marker);
        if( p<0 ) {
            return "0";
        }
        int start = p + marker.length();
        int end = start;
        while( end<description.length() && description.charAt(end)!=']' ) {
            end++;
        }
        String acc = description.substring(start, end).trim();
        return acc.isEmpty() ? "0" : acc;
    }

    /// decode GFF3 percent-escapes (%2C, %3B, %3D, %25, %26, ...) into characters
    static String gff3Unescape(String s) {
        if( s==null || s.indexOf('%')<0 ) {
            return s;
        }
        StringBuilder sb = new StringBuilder(s.length());
        for( int i=0; i<s.length(); i++ ) {
            char c = s.charAt(i);
            if( c=='%' && i+2<s.length() ) {
                try {
                    sb.append((char) Integer.parseInt(s.substring(i+1, i+3), 16));
                    i += 2;
                    continue;
                } catch( NumberFormatException ignore ) {
                }
            }
            sb.append(c);
        }
        return sb.toString();
    }

    // --- getters / setters ---------------------------------------------------

    public String getGenomeBuild() {
        return genomeBuild;
    }
    public void setGenomeBuild(String genomeBuild) {
        this.genomeBuild = genomeBuild;
    }

    public String getEnsemblGenePrefix() {
        return ensemblGenePrefix;
    }
    public void setEnsemblGenePrefix(String ensemblGenePrefix) {
        this.ensemblGenePrefix = ensemblGenePrefix;
    }

    public String getXrefAuthority() {
        return xrefAuthority;
    }
    public void setXrefAuthority(String xrefAuthority) {
        this.xrefAuthority = xrefAuthority;
    }

    public String getGff3File() {
        return gff3File;
    }
    public void setGff3File(String gff3File) {
        this.gff3File = gff3File;
    }

    public String getEntrezFile() {
        return entrezFile;
    }
    public void setEntrezFile(String entrezFile) {
        this.entrezFile = entrezFile;
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

    static class Exon {
        public int exonStart;
        public int exonStop;
        public String rank;
        public int cdsStart;
        public int cdsStop;
    }

    static class Tr {
        public String trId;
        public String trVer;
        public String chr;
        public String trStart;
        public String trStop;
        public String strand; // 1, -1
        public String proteinId = "";
        public String trType;
        public List<Exon> exons = new ArrayList<>();
    }
}
