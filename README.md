# ensembl-data-pipeline
Load genes and transcripts from Ensembl into RGD.

The following species are processed: rat, mouse, human, dog, pig, bonobo, vervet, chinchilla and molerat.

For scaffold assemblies, GeneBank chromosome accessions are replaced with RefSeq accessions
in order to easily compare NCBI and Ensembl locus data.

Nomenclature events are generated only for genes:

    1. having GENE_SOURCE='Ensembl'
    2. having NOMEN_SOURCE either empty (null) or other than ('HGNC', 'MGI')

XDB_IDS (external database ids)

    1. genes: only Ensembl Gene Ids are loaded 
