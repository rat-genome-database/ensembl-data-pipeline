package edu.mcw.rgd.data;

import edu.mcw.rgd.process.MemoryMonitor;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;
import java.util.zip.GZIPInputStream;

/// loads PolyPhen-2 predictions for human variants from Ensembl variation databases into the POLYPHEN table.
/// <p>
/// Ensembl keeps one precomputed prediction matrix per protein sequence (every position x every amino acid),
/// keyed by the md5 of the protein sequence. RGD reference protein sequences of missense variant transcripts
/// are matched to Ensembl proteins by md5, so the match is exact and works for RefSeq and Ensembl transcripts alike.
/// A variant transcript is missense when its variant protein (VAR_AA, from the variant position to the first stop codon)
/// differs from the reference protein at the variant position only; rows of withdrawn transcripts are skipped.
/// <p>
/// Matrix format (Ensembl ProteinFunctionPredictionMatrix.pm): gzipped; 'VEP' header followed by a 16-bit
/// little-endian value for every protein position and every amino acid in order ACDEFGHIKLMNPQRSTVWY;
/// the top 2 bits are the prediction, the low 10 bits the score x 1000; 0xFFFF = no prediction.
/// <p>
/// All configured assemblies are processed in one run: GRCh37 and GRCh38 variants share variant RGD ids,
/// and POLYPHEN rows are matched to variant transcripts by (variant_rgd_id, transcript_rgd_id) only.
/// The POLYPHEN rows synced are all rows of the missense variants processed; rows of variants that are
/// no longer missense on any processed assembly are not seen, and are left as is.
public class PolyphenLoader {

    Logger log = LogManager.getLogger("status");

    private String version;
    private String model; // PolyPhen-2 classifier model: 'humdiv' or 'humvar'
    private List<PolyphenSource> sources;
    private EnsemblDataPuller dataPuller;

    EnsemblDAO dao = new EnsemblDAO();

    /// varTailMd5: first 64 bits of the md5 of the variant protein after the variant position (see EnsemblDAO)
    record VtRow(long variantRgdId, long transcriptRgdId, int pos, char refAa, char varAa, long seqKey, long varTailMd5) {}
    record SeqInfo(String md5, int length) {}

    static final int VAR_TAIL_LENGTH = 3899; // VAR_AA keeps up to 3900 amino acids, the changed one included

    static final String AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"; // matrix column order
    static final String[] PREDICTIONS = {"probably damaging", "possibly damaging", "benign", "unknown"};
    static final int NO_PREDICTION = 0xFFFF;
    static final String MATRIX_HEADER = "VEP";

    public void run() throws Exception {
        long time0 = System.currentTimeMillis();
        log.info(getVersion());
        log.info("  PolyPhen-2 model: "+getModel());
        if( !getModel().equals("humdiv") && !getModel().equals("humvar") ) {
            throw new Exception("Aborted: unknown PolyPhen-2 model '"+getModel()+"' -- must be 'humdiv' or 'humvar'");
        }

        MemoryMonitor memoryMonitor = new MemoryMonitor();
        memoryMonitor.start();
        try {
            java.util.Map<String, PolyphenRow> incoming = new HashMap<>();
            Set<Long> variantRgdIds = new HashSet<>(); // missense variants of all processed assemblies
            List<Integer> mapKeys = new ArrayList<>();
            for( PolyphenSource source: getSources() ) {
                mapKeys.add(source.getMapKey());
            }
            java.util.Map<Integer, List<VtRow>> vtRowsByMap = dao.getMissenseVariantTranscripts(mapKeys);
            for( PolyphenSource source: getSources() ) {
                List<VtRow> vtRows = vtRowsByMap.remove(source.getMapKey());
                processAssembly(source, vtRows==null ? Collections.emptyList() : vtRows, incoming, variantRgdIds);
            }
            syncWithDatabase(incoming, dao.getPolyphenRows(variantRgdIds));
        } finally {
            memoryMonitor.stop();
            log.info(memoryMonitor.getSummary());
        }

        log.info("POLYPHEN LOAD OK -- elapsed "+Utils.formatElapsedTime(time0, System.currentTimeMillis()));
        log.info("===");
    }

    void processAssembly(PolyphenSource source, List<VtRow> vtRows, java.util.Map<String, PolyphenRow> incoming, Set<Long> variantRgdIds) throws Exception {

        int mapKey = source.getMapKey();
        log.info("map_key "+mapKey+":");

        Set<Long> seqKeys = new HashSet<>();
        Set<Long> transcriptRgdIds = new HashSet<>();
        for( VtRow vt: vtRows ) {
            seqKeys.add(vt.seqKey());
            transcriptRgdIds.add(vt.transcriptRgdId());
            variantRgdIds.add(vt.variantRgdId());
        }
        java.util.Map<Long, String> sequences = dao.getProteinSequences(seqKeys);
        java.util.Map<Long, String[]> transcriptProteins = dao.getTranscriptProteins(transcriptRgdIds);
        log.info("  variant transcripts with a changed amino acid: "+Utils.formatThousands(vtRows.size())
            +", distinct reference proteins: "+Utils.formatThousands(sequences.size()));

        // reference proteins without the trailing stop codon '*': md5 and length, as in Ensembl's translation_md5
        java.util.Map<Long, SeqInfo> proteins = new HashMap<>();
        for( java.util.Map.Entry<Long, String> e: sequences.entrySet() ) {
            String protein = Utils.defaultString(e.getValue()).replaceAll("\\*+$", "");
            proteins.put(e.getKey(), new SeqInfo(md5(protein), protein.length()));
        }

        String prefix = "polyphen_"+mapKey+"_";
        int attribId = resolveAnalysisAttribId(dataPuller.downloadEnsemblFile(source.getAttribUrl(), prefix));
        Set<String> md5s = new HashSet<>();
        for( SeqInfo si: proteins.values() ) {
            md5s.add(si.md5());
        }
        java.util.Map<String, Integer> translationIds = loadTranslationIds(dataPuller.downloadEnsemblFile(source.getTranslationMd5Url(), prefix), md5s);
        log.info("  reference proteins matched to Ensembl proteins: "+Utils.formatThousands(translationIds.size()));

        // group the variant transcripts by Ensembl protein
        java.util.Map<Integer, List<VtRow>> rowsByTranslation = new HashMap<>();
        java.util.Map<Integer, Integer> proteinLengths = new HashMap<>();
        int stopCodon = 0, nonStandardAa = 0, noTranscript = 0, notSingleSubstitution = 0, noEnsemblProtein = 0;
        for( VtRow vt: vtRows ) {
            if( vt.refAa()=='*' || vt.varAa()=='*' ) {
                stopCodon++; // stop gained or stop lost: not scored by PolyPhen
                continue;
            }
            if( AMINO_ACIDS.indexOf(vt.refAa())<0 || AMINO_ACIDS.indexOf(vt.varAa())<0 ) {
                nonStandardAa++;
                continue;
            }
            if( !transcriptProteins.containsKey(vt.transcriptRgdId()) ) {
                noTranscript++; // stale VARIANT_TRANSCRIPT row of a withdrawn transcript
                continue;
            }
            if( !isSingleSubstitution(vt, sequences.get(vt.seqKey())) ) {
                notSingleSubstitution++;
                continue;
            }
            SeqInfo seq = proteins.get(vt.seqKey());
            Integer translationId = translationIds.get(seq.md5());
            if( translationId==null ) {
                noEnsemblProtein++;
                continue;
            }
            rowsByTranslation.computeIfAbsent(translationId, k -> new ArrayList<>()).add(vt);
            proteinLengths.put(translationId, seq.length());
        }
        vtRows = null;
        sequences = null;

        // stream the prediction matrices once, decoding only the cells of our variants
        int[] counts = new int[6]; // 0: with prediction, 1: no prediction in matrix, 2: position beyond protein end,
                                   // 3: proteins with malformed matrix, 4: rows on them, 5: rows already loaded for another assembly
        String predictionsFile = dataPuller.downloadEnsemblFile(source.getPredictionsUrl(), prefix);
        streamPredictionMatrices(predictionsFile, attribId, rowsByTranslation.keySet(), (translationId, matrix) -> {

            List<VtRow> rows = rowsByTranslation.remove(translationId);
            int length = proteinLengths.get(translationId);
            if( matrix.length != MATRIX_HEADER.length() + length*AMINO_ACIDS.length()*2
                || !new String(matrix, 0, MATRIX_HEADER.length()).equals(MATRIX_HEADER) ) {
                counts[3]++;
                counts[4] += rows.size();
                return;
            }

            for( VtRow vt: rows ) {
                if( vt.pos()<1 || vt.pos()>length ) {
                    counts[2]++;
                    continue;
                }
                int offset = MATRIX_HEADER.length() + ((vt.pos()-1)*AMINO_ACIDS.length() + AMINO_ACIDS.indexOf(vt.varAa()))*2;
                int value = (matrix[offset] & 0xFF) | ((matrix[offset+1] & 0xFF) << 8);
                if( value==NO_PREDICTION ) {
                    counts[1]++;
                    continue;
                }
                PolyphenRow r = createRow(vt, value, transcriptProteins.get(vt.transcriptRgdId()));
                if( incoming.putIfAbsent(r.getKey(), r)==null ) {
                    counts[0]++;
                } else {
                    counts[5]++;
                }
            }
        });

        int rowsNoMatrix = 0;
        for( List<VtRow> rows: rowsByTranslation.values() ) {
            rowsNoMatrix += rows.size();
        }

        log.info("  variant transcripts with PolyPhen prediction: "+Utils.formatThousands(counts[0]+counts[5])
            +(counts[5]>0 ? " (incl. "+Utils.formatThousands(counts[5])+" already loaded from another assembly)" : ""));
        log.info("  variant transcripts without prediction:");
        log.info("    transcript no longer in RGD: "+Utils.formatThousands(noTranscript));
        log.info("    not a single amino acid substitution (in-frame indel or several changes): "+Utils.formatThousands(notSingleSubstitution));
        log.info("    RGD reference protein not found in Ensembl (md5 mismatch): "+Utils.formatThousands(noEnsemblProtein));
        log.info("    Ensembl protein without PolyPhen matrix: "+Utils.formatThousands(rowsNoMatrix)
            +" (proteins: "+Utils.formatThousands(rowsByTranslation.size())+")");
        log.info("    PolyPhen matrix has no score for the substitution: "+Utils.formatThousands(counts[1]));
        log.info("    stop gained or stop lost (not scored by PolyPhen): "+Utils.formatThousands(stopCodon));
        if( nonStandardAa>0 ) {
            log.info("    other non-standard amino acid: "+Utils.formatThousands(nonStandardAa));
        }
        if( counts[2]>0 ) {
            log.info("    position beyond protein end: "+Utils.formatThousands(counts[2]));
        }
        if( counts[3]>0 ) {
            log.warn("    WARNING: malformed prediction matrix: "+Utils.formatThousands(counts[4])
                +" (proteins: "+Utils.formatThousands(counts[3])+")");
        }
    }

    /// VAR_AA is the variant protein from the variant position to its first stop codon: the variant is a single
    /// amino acid substitution when the reference has the reference amino acid at the variant position and
    /// the rest of the variant protein is the same as the rest of the reference protein (stop codon included)
    boolean isSingleSubstitution(VtRow vt, String refSequence) {
        if( refSequence==null || vt.pos()<1 || vt.pos()>refSequence.length() || refSequence.charAt(vt.pos()-1)!=vt.refAa() ) {
            return false;
        }
        String refTail = refSequence.substring(vt.pos(), Math.min(refSequence.length(), vt.pos()+VAR_TAIL_LENGTH));
        return md5Prefix(md5(refTail))==vt.varTailMd5();
    }

    static String md5(String s) {
        try {
            byte[] digest = java.security.MessageDigest.getInstance("MD5").digest(s.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
            return HexFormat.of().formatHex(digest);
        } catch( java.security.NoSuchAlgorithmException e ) {
            throw new RuntimeException(e);
        }
    }

    /// first 64 bits of an md5 hex string
    static long md5Prefix(String md5Hex) {
        return Long.parseUnsignedLong(md5Hex.substring(0, 16), 16);
    }

    PolyphenRow createRow(VtRow vt, int value, String[] transcriptProtein) {
        PolyphenRow r = new PolyphenRow();
        r.variantRgdId = vt.variantRgdId();
        r.transcriptRgdId = vt.transcriptRgdId();
        if( transcriptProtein!=null ) {
            r.proteinId = transcriptProtein[0];
            r.geneSymbol = transcriptProtein[1];
        }
        r.position = vt.pos();
        r.aa1 = String.valueOf(vt.refAa());
        r.aa2 = String.valueOf(vt.varAa());

        int prediction = value >> 14;
        r.prediction = PREDICTIONS[prediction];
        switch( prediction ) {
            case 0, 1 -> r.pph2Class = "deleterious";
            case 2 -> r.pph2Class = "neutral";
            default -> r.pph2Class = "none";
        }
        // score is stored x 1000; written without trailing zeros, like PolyPhen-2 reports: '0.96', '1', '0'
        r.pph2Prob = prediction==3 ? null : BigDecimal.valueOf(value & 0x3FF, 3).stripTrailingZeros().toPlainString();
        return r;
    }

    void syncWithDatabase(java.util.Map<String, PolyphenRow> incoming, Collection<PolyphenRow> inRgd) throws Exception {

        if( incoming.isEmpty() ) {
            throw new Exception("Aborted: no PolyPhen predictions loaded from Ensembl -- database left unchanged");
        }

        List<PolyphenRow> forInsert = new ArrayList<>();
        List<PolyphenRow> forUpdate = new ArrayList<>();
        List<PolyphenRow> forDelete = new ArrayList<>();
        int upToDate = 0;

        java.util.Map<String, PolyphenRow> rgdByKey = new HashMap<>();
        for( PolyphenRow r: inRgd ) {
            if( rgdByKey.putIfAbsent(r.getKey(), r)!=null ) {
                forDelete.add(r); // duplicate row
            }
        }

        for( PolyphenRow r: incoming.values() ) {
            PolyphenRow old = rgdByKey.remove(r.getKey());
            if( old==null ) {
                forInsert.add(r);
            } else if( old.hasSameData(r) ) {
                upToDate++;
            } else {
                r.polyphenId = old.polyphenId;
                forUpdate.add(r);
                log.debug("UPDATE OLD "+old.dump());
                log.debug("UPDATE NEW "+r.dump());
            }
        }
        forDelete.addAll(rgdByKey.values());
        for( PolyphenRow r: forDelete ) {
            log.debug("DELETE "+r.dump());
        }

        dao.insertPolyphenRows(forInsert);
        dao.updatePolyphenRows(forUpdate);
        dao.deletePolyphenRows(forDelete);

        log.info("POLYPHEN rows:");
        log.info("  inserted:   "+Utils.formatThousands(forInsert.size()));
        log.info("  updated:    "+Utils.formatThousands(forUpdate.size()));
        log.info("  deleted:    "+Utils.formatThousands(forDelete.size()));
        log.info("  up-to-date: "+Utils.formatThousands(upToDate));
    }

    /// attrib.txt: attrib_id, attrib_type_id, value -- the PolyPhen analyses are 'polyphen_humdiv' and 'polyphen_humvar'
    int resolveAnalysisAttribId(String attribFile) throws Exception {
        String analysis = "polyphen_"+getModel();
        Integer attribId = null;
        try( BufferedReader in = Utils.openReader(attribFile) ) {
            String line;
            while( (line=in.readLine())!=null ) {
                String[] cols = line.split("\t", -1);
                if( cols.length>=3 && cols[2].equals(analysis) ) {
                    if( attribId!=null ) {
                        throw new Exception("Aborted: multiple attrib ids for analysis "+analysis+" in "+attribFile);
                    }
                    attribId = Integer.parseInt(cols[0]);
                }
            }
        }
        if( attribId==null ) {
            throw new Exception("Aborted: no attrib id for analysis "+analysis+" in "+attribFile);
        }
        log.info("  analysis "+analysis+": attrib_id "+attribId);
        return attribId;
    }

    /// translation_md5.txt: translation_md5_id, translation_md5; returns md5 -> id for the given md5s only
    java.util.Map<String, Integer> loadTranslationIds(String translationMd5File, Set<String> md5s) throws Exception {
        java.util.Map<String, Integer> result = new HashMap<>();
        try( BufferedReader in = Utils.openReader(translationMd5File) ) {
            String line;
            while( (line=in.readLine())!=null ) {
                String[] cols = line.split("\t", -1);
                if( md5s.contains(cols[1]) ) {
                    result.put(cols[1], Integer.parseInt(cols[0]));
                }
            }
        }
        return result;
    }

    interface MatrixHandler {
        void handle(int translationId, byte[] matrix) throws Exception;
    }

    /// protein_function_predictions.txt is a MySQL tab dump: translation_md5_id, analysis_attrib_id, prediction_matrix;
    /// the binary matrix has backslash-escaped bytes (\0 for NUL; tab, newline and backslash are preceded by '\');
    /// only the wanted matrices are unescaped and gunzipped
    void streamPredictionMatrices(String predictionsFile, int attribId, Set<Integer> wantedTranslationIds, MatrixHandler handler) throws Exception {

        InputStream fileIn = new FileInputStream(predictionsFile);
        if( predictionsFile.endsWith(".gz") ) {
            fileIn = new GZIPInputStream(fileIn, 1<<16);
        }
        int matrices = 0;
        try( InputStream in = new BufferedInputStream(fileIn, 1<<20) ) {
            ByteArrayOutputStream blob = new ByteArrayOutputStream(1<<16);
            while( true ) {
                long translationId = readNumber(in);
                if( translationId<0 ) {
                    break; // end of file
                }
                long analysisId = readNumber(in);
                boolean wanted = analysisId==attribId && wantedTranslationIds.contains((int) translationId);
                blob.reset();
                readBlob(in, wanted ? blob : null);
                if( wanted && !(blob.size()==2 && blob.toByteArray()[0]=='\\') ) { // skip NULL matrix: \N
                    byte[] matrix = new GZIPInputStream(new ByteArrayInputStream(blob.toByteArray())).readAllBytes();
                    handler.handle((int) translationId, matrix);
                    matrices++;
                }
            }
        }
        log.info("  prediction matrices used: "+Utils.formatThousands(matrices));
    }

    /// reads a tab-terminated decimal number; returns -1 at the end of file
    long readNumber(InputStream in) throws Exception {
        long value = 0;
        int digits = 0;
        int b;
        while( (b=in.read())!='\t' ) {
            if( b<0 && digits==0 ) {
                return -1;
            }
            if( b<'0' || b>'9' ) {
                throw new Exception("unexpected byte "+b+" in a numeric field of prediction matrix file");
            }
            value = value*10 + (b-'0');
            digits++;
        }
        return value;
    }

    /// reads the last field of a row up to an unescaped newline (or end of file), unescaping it into 'out' if not null;
    /// a NULL field (\N) is kept as is
    void readBlob(InputStream in, ByteArrayOutputStream out) throws Exception {
        int b;
        while( (b=in.read())>=0 && b!='\n' ) {
            if( b=='\\' ) {
                int c = in.read();
                if( c<0 ) {
                    throw new Exception("unexpected end of prediction matrix file");
                }
                if( out!=null ) {
                    switch( c ) {
                        case '0' -> out.write(0);
                        case 'N' -> { out.write('\\'); out.write('N'); }
                        case 'n' -> out.write('\n');
                        case 't' -> out.write('\t');
                        case 'r' -> out.write('\r');
                        case 'b' -> out.write('\b');
                        case 'Z' -> out.write(26);
                        default -> out.write(c);
                    }
                }
            } else if( out!=null ) {
                out.write(b);
            }
        }
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getModel() {
        return model;
    }

    public void setModel(String model) {
        this.model = model;
    }

    public List<PolyphenSource> getSources() {
        return sources;
    }

    public void setSources(List<PolyphenSource> sources) {
        this.sources = sources;
    }

    public EnsemblDataPuller getDataPuller() {
        return dataPuller;
    }

    public void setDataPuller(EnsemblDataPuller dataPuller) {
        this.dataPuller = dataPuller;
    }
}
