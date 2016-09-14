A multicorn based foreign data wrapper for VCF 4.0 / 4.1 files using PyVCF.

Requirements
============

1. PostgreSQL 9.2 or later
2. Multicorn 1.3.2 or later 
3. PyVCF 0.6.8 or later
4. One VCF file per sample compressed with bgzip and tabix indexed

The filesystem layout must be like this:

* `<basedir>`
 * `<species>`  
  - `<sample><suffix>`
  - `<sample>.tbi`

This example is for data from the 1000 Genomes project, assuming a suffix of .vcf.gz:

* `<basedir>`
 * `<human>`  
  - `HG00099.vcf.gz`
  - `HG00099.tbi`
  - `HG00100.vcf.gz`
  - `HG00100.tbi`
  - `..more..`

Installation
============

1. Install the multicorn extension
2. sudo python setup.py install
3. Execute the following SQL commands in PostgreSQL:

```sql
CREATE EXTENSION multicorn;

CREATE SERVER multicorn_vcf FOREIGN DATA WRAPPER multicorn
OPTIONS (wrapper 'vcf_fdw.VCFForeignDataWrapper');

CREATE FOREIGN TABLE vcf.<species>
   (chrom text ,
    pos integer ,
    id text ,
    ref text ,
    alt text[] ,
    qual real ,
    heterozygosity real ,
    sample text ,
    species text ,
    info text ,
    depth integer ,
    partition integer ,
    genotype text ,
    filter text ,
    issnp boolean ,
    issv boolean ,
    isindel boolean ,
    ismonomorphic boolean ,
    isdeletion boolean ,
    issvprecise boolean ,
    istransition boolean ,
    source text )
   SERVER multicorn_vcf
   OPTIONS (basedir '<basedir>', species '<species>', suffix '<suffix>');
```

Basic Usage
============

CAUTION: Since tabix indexes only chromosome and position, query performance suffers for queries that do not specify chromosome and/or position in the WHERE clause!

The query below retrieves info and genotype for all samples in the vcf files within the specified region.

```sql
SELECT info, genotype FROM vcf.<species> WHERE chrom = '1' AND pos between 500000 AND 10000000;
```

The query below retrieves info and genotype for all named samples in the vcf files within the specified region.

```sql
SELECT info, genotype FROM vcf.<species> WHERE chrom = '1' AND pos between 500000 AND 10000000 AND sample IN ('HG00099','HG00100','HG00101','HG00230');
```

The query below retrieves info and genotype for all named chromosomes and all named samples in the vcf files within the specified region.

```sql
SELECT info, genotype FROM vcf.<species> WHERE chrom in ('1','X') AND pos between 500000 AND 10000000 AND sample IN ('HG00099','HG00100','HG00101','HG00230');
```

The query below retrieves info and genotype for all named samples in the vcf files within the specified region filtered by SNP.

```sql
SELECT info, genotype, issnp FROM vcf.<species> WHERE chrom = '1' AND pos between 500000 AND 10000000 AND sample IN ('HG00099','HG00100','HG00101','HG00230') AND issnp = TRUE;
```

Disclaimer
==========

This software comes without any warranty whatsoever. Use at your own risk.
