package edu.mcw.rgd.data;

import java.io.*;

public class _XRefTool {

    public static void main(String[] args) throws IOException {
        String input_fname = "data/rapid/xref.tsv";
        String output_fname = "data/rapid/xref_fixed.tsv";
        BufferedWriter out = new BufferedWriter(new FileWriter(output_fname));
        BufferedReader in = new BufferedReader(new FileReader(input_fname));
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

        System.out.println("done");
    }
}
