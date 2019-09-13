package edu.mcw.rgd.data;

import java.sql.*;

/**
 * Created by sellanki on 8/13/2019.
 */
public class Databaseloading {
    public static Connection getConnection() {
        try {
            String JDBC_DRIVER = "oracle.jdbc.OracleDriver";
            String DB_URL = "jdbc:oracle:thin:@marcus.rgd.mcw.edu:1521:rgddev1";
            String USER = "dev_1";
            String PASS = "P3m2drefra8ucrat8a82abebn3w#19";
            Class.forName(JDBC_DRIVER);
            Connection conn = DriverManager.getConnection(DB_URL, USER, PASS);
            System.out.println("Connection Established");
            return conn;

        } catch (Exception e) {
            System.out.println("Connection not established");
            return null;
        }

    }

    public static void addgene(Connection conn, PreparedStatement preparedstatement, String ENSEMBL_GENE_ID, String CHROMOSE_NAME, String START_POSITION, String END_POSITION, String STRAND, String RGD_ID, String MAP_KEY, String ENTREZGENE_ID,String GENE_NAME,String GENE_BIOTYPE,String GENE_DESCRIPTION) throws SQLException {
        try {
            preparedstatement = conn.prepareStatement("insert into EXP(ENSEMBL_GENE_ID,CHROMOSE_NAME,START_POSITION,END_POSITION,STRAND,RGD_ID,MAP_KEY,ENTREZGENE_ID,GENE_NAME,GENE_BIOTYPE,GENE_DESCRIPTION) values(?,?,?,?,?,?,?,?,?,?,?)");
            preparedstatement.setString(1, ENSEMBL_GENE_ID);
            preparedstatement.setString(2, CHROMOSE_NAME);
            preparedstatement.setString(3, START_POSITION);
            preparedstatement.setString(4, END_POSITION);
            preparedstatement.setString(5, STRAND);
            preparedstatement.setString(6, RGD_ID);
            preparedstatement.setString(7, MAP_KEY);
            preparedstatement.setString(8, ENTREZGENE_ID);
            preparedstatement.setString(9, GENE_NAME);
            preparedstatement.setString(10, GENE_BIOTYPE);
            preparedstatement.setString(11, GENE_DESCRIPTION);
            preparedstatement.executeUpdate();
        } catch (SQLException e) {
            e.printStackTrace();
        } finally {
            try {
                if (preparedstatement != null)
                    preparedstatement.close();
            } catch (SQLException e) {
                e.printStackTrace();
            }
        }
    }


    public static int checkid(Connection conn, PreparedStatement preparedstatement, String RGD_ID) throws Exception

   /*{
           Connection conn = getConnection();
           String OBJECT_STATUS = null;
           Statement st = conn.createStatement();
           ResultSet rs = st.executeQuery("select RGD_ID,OBJECT_STATUS from RGD_IDS where RGD_ID=?");
           while (rs.next()) {
               RGD_ID = rs.getString("RGD_ID");
               OBJECT_STATUS = rs.getString("OBJECT_STATUS");
           }
           if (RGD_ID != "0" && OBJECT_STATUS != "ACTIVE")
               System.out.println("" + RGD_ID);
           System.out.println("" + OBJECT_STATUS);

   }*/ {
        String OBJECT_STATUS = "";
        try {
            preparedstatement = conn.prepareStatement(" select RGD_ID from RGD_IDS where RGD_ID=?");
            preparedstatement.setString(1, RGD_ID);
            ResultSet rs = preparedstatement.executeQuery();
            while (rs.next()) {
                if (RGD_ID.length() > 0)
                    return Integer.parseInt(RGD_ID);
                else
                    return 0;
            }

        } catch (SQLException e) {
            e.printStackTrace();
        } finally {
            try {
                if (preparedstatement != null)
                    preparedstatement.close();
            } catch (SQLException e) {
                e.printStackTrace();
            }
        }
return 0;
    }


    public void objectstatus(Connection conn,PreparedStatement preparedstatement,String RGD_ID)
    {
        try
        {
            preparedstatement = conn.prepareStatement(" select OBJECT_STATUS from RGD_IDS where RGD_ID=?");
            preparedstatement.executeUpdate();
        }
        catch (SQLException e) {
            e.printStackTrace();
        }

    }


}
