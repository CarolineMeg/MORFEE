# MORFEE

*MORFEE* (**Mutation on Open Reading FramE annotation**) is a tool that detects and annotates single nucleotide variants (SNVs) (i) creating canonical (ATG) and non-canonical Translation Initiation Sites (TIS); (ii) creating Stop codons (TAA, TAG and TGA), and ; (iii) deleting Stop codons in the 5’UTR (uTIS, uStop and New_uStop), the coding sequence (CDS, intTIS, intStop and New_intStop), and the 3’UTR (dTIS, dStop and New_dStop).

MORFEE is an R package that could be applied on any VCF file.

MORFEE algorithm is written in R language and can run on all operating
systems which have an R interpreter (including Linux, macOS and Windows).
MORFEE starts with a minimal VCF file (i.e. with at least the *chr*, *position*,
*reference allele* and *alternate allele* fields) as an input, but can also
work with already [ANNOVAR](http://annovar.openbioinformatics.org/)-annotated
([Wang et al., 2010](https://doi.org/10.1093/nar/gkq603)) VCF files.
MORFEE has some R packages dependencies which are available on
[CRAN](https://cran.r-project.org/) or
[Bioconductor](https://www.bioconductor.org/) repositories and uses the
[GENCODE](https://www.gencodegenes.org/) database.



# Installation 

```
library("devtools")
install_local("path_to_morfee")
```

# First step - Annotation with ANNOVAR

MORFEE takes as input a VCF file containing at least the chromosome, position, reference allele and alternate allele columns. An example file is available at inst/extdata/VCF_example.vcf.

In a first step, MORFEE reads the input VCF file and use ANNOVAR (that has to be installed beforehand) through the wrapper 
function vcf.annovar.annotation to annotate all variants. This step should be skipped if the input file has already been annotated. 
    
The necessary databases (refgene, ensgene, gwasCatalog, avsnp150, clinvar, gnomad, dbsnfp and contraint LoF metrics) for hg19 and hg38 can be downloaded with MORFEE using the following function :

```
download.annovar.db(hg=19)
```

The databases downloaded by this function are the following :
* refGene
* avsnp150
* gnomad211_genome for hg19 or gnomad312_genome for hg38
* clinvar_20221231
* dbnsfp42a
* gwasCatalog
* gnomad.v2.1.1.lof_metrics.by_gene 

If you want to download the databases individually, the lof metrics file can be downloaded there : https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz and the other databases can be downloaded with ANNOVAR using the following command :
```
annotate_variation.pl -downdb -webfrom annovar $NAME_OF_THE_DATABASE humandb -buildver hg38
```



Annotation can be done two ways :

1. Directly with ANNOVAR : 
For hg19 :

```
table_annovar.pl \
        inst/extdata/VCF_example.vcf \
        ${PATH_ANNOVAR_DB} \
        -buildver hg19 \
        -remove -protocol refGene,ensGene,gwasCatalog,avsnp150,clinvar_20221231,gnomad211_genome,dbnsfp42a \
        -operation gx,gx,r,f,f,f,f -nastring . -vcfinput \
        -xref ${PATH_ANNOVAR_DB}/morfee.lof_metrics.by_gene.txt
```

For hg38 :
```
table_annovar.pl \
        inst/extdata/VCF_example.vcf \
        ${PATH_ANNOVAR_DB} \
        -buildver hg38 \
        -remove -protocol refGene,ensGene,gwasCatalog,avsnp150,clinvar_20221231,gnomad312_genome,dbnsfp42a \
        -operation gx,gx,r,f,f,f,f -nastring . -vcfinput \
        -xref ${PATH_ANNOVAR_DB}/morfee.lof_metrics.by_gene.txt
```
        
2. Using MORFEE : 
Note : if the annotation is done through MORFEE, only the following version of the databases can be used : clinvar_20221231, dbnsfp42a, gnomad211_genome for hg19 and gnomad312_genome for hg38. To use different versions, ANNOVAR should be used directly.

```
library(VariantAnnotation)
library(MORFEE)

my_raw_vcf <- "inst/extdata/VCF_example.vcf"

# Annotation with ANNOVAR through the wrapper function
vcf.annovar.annotation(my_raw_vcf,hg=19,path_annovar_db)
```
The arguments path_annovar_db and path_to_annovar can be used to indicate respectively the path to ANNOVAR databases and to ANNOVAR software.





# Second step : MORFEE annotation
Once ANNOVAR annotation is done, MORFEE can be launched with the following commands :

```
library("VariantAnnotation")
library("MORFEE")

# Download database and create an R object used by MORFEE - by default GRCh37, can be changed with "GRCh=38" argument
MORFEE_DATA <- get.morfee.data()

# To avoid too much Internet traffic, save locally this object
save(MORFEE_DATA, file="MORFEE_DATA.RData") #file size should be about 850M
load("MORFEE_DATA.RData")

MY_VCF = "/path/to/annotated/file.vcf" 

# Read/Load the VCF into R
my_vcf <- readVcf(MY_VCF)

# Perform the MORFEE annotation
my_vcf_morfee <- morfee.annotation(my_vcf, MORFEE_DATA)

# Keep only variants which create a new TIS sequence, create a new STOP or delete a STOP
my_vcf_morfee_NEW_TIS <- my_vcf_morfee[[1]][!is.na(info(my_vcf_morfee[[1]])$MORFEE_uTIS) | !is.na(info(my_vcf_morfee[[1]])$MORFEE_uSTOP) | !is.na(info(my_vcf_morfee[[1]])$MORFEE_dTIS) | !is.na(info(my_vcf_morfee[[1]])$MORFEE_dSTOP) | !is.na(info(my_vcf_morfee[[1]])$MORFEE_intTIS)  | !is.na(info(my_vcf_morfee[[1]])$MORFEE_intSTOP) | !is.na(info(my_vcf_morfee[[1]])$MORFEE_intSTOP_creation) | !is.na(info(my_vcf_morfee[[1]])$MORFEE_uSTOP_creation) | !is.na(info(my_vcf_morfee[[1]])$MORFEE_dSTOP_creation)]

my_vcf_morfee_NEW_STOP <- my_vcf_morfee[[2]][!is.na(info(my_vcf_morfee[[2]])$MORFEE_uTIS) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_uSTOP) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_dTIS) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_dSTOP) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_intTIS) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_intSTOP) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_intSTOP_creation) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_uSTOP_creation) | !is.na(info(my_vcf_morfee[[2]])$MORFEE_dSTOP_creation)]

my_vcf_morfee_NEW_STOP_creation <- my_vcf_morfee[[3]][!is.na(info(my_vcf_morfee[[3]])$MORFEE_uTIS) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_uSTOP) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_dTIS) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_dSTOP) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_intTIS) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_intSTOP) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_intSTOP_creation) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_uSTOP_creation) | !is.na(info(my_vcf_morfee[[3]])$MORFEE_dSTOP_creation)]


# Write an XLSX file with the new MORFEE field/annotation, hg19 or hg38 must be specified with "19" or "38"
vcf_2_xlsx(my_vcf_morfee_NEW_TIS,my_vcf_morfee_NEW_STOP,my_vcf_morfee_NEW_STOP_creation, file="filename.xlsx",19)

```



# Output


Each output file includes the following general information:

orfSNVs_type: 
* u/d/intTIS*: creation of an upstream/downstream/internal to CDS TIS (ATG, CTG, TTG, GTG, AAG, ACG, AGG, ATC, ATA, ATT).
* u/d/intSTOP*: deletion of an upstream/downstream/ internal to CDS  STOP codon.
* new_u/d/intSTOP*: creation of an upstream/downstream/ internal to CDS STOP codon.

orfSNVs_frame: 
*in_frame* or *out_of_frame* according to whether the new TIS or the new or deleted STOP is in frame or not with the main ATG

type_of_generated_ORF: 
Possible values are *overlapping*, *not_overlapping* or *elongated_CDS* according to the type of the predicted upstream /downstream/internal ORF compared to the canonical CDS. 

Ratio_length_pred_obs: 
Ratio of the length of the predicted ORF over that of the CDS. In case of several transcripts, only one value is given, the one corresponding to the predicted ORF with length closest to that of the CDS.

NewAALength: 
Length (i.e. predicted number of amino acids) of the predicted peptide corresponding to the annotated ORF.

ORF_size:
Length (predicted number of nucleotides) of the predicted upORF. 

MORFEE_uTIS:
MORFEE VCF-style annotation for variants creating new upstream TIS. It summarizes all the information contained in other cases. 

MORFEE_uSTOP
MORFEE VCF-style annotation for variants deleting upstream STOP. It summarizes all the information contained in other cases. 


MORFEE_uSTOP_creation:
MORFEE VCF-style annotation for variants creating upstream STOP. It summarizes all the information contained in other cases. 

MORFEE_dTIS:
MORFEE VCF-style annotation for variants creating new downstream TIS. It summarizes all the information contained in other cases. 

MORFEE_dSTOP
MORFEE VCF-style annotation for variants deleting downstream STOP. It summarizes all the information contained in other cases. 

MORFEE_dSTOP_creation
MORFEE VCF-style annotation for variants creating downstream STOP. It summarizes all the information contained in other cases. 

MORFEE_intTIS
MORFEE VCF-style annotation for variants creating new exonic TIS. It summarizes all the information contained in other cases. 

MORFEE_intSTOP
MORFEE VCF-style annotation for variants deleting exonic STOP. It summarizes all the information contained in other cases. 

MORFEE_intSTOP_creation:
MORFEE VCF-style annotation for variants creating exonic STOP. It summarizes all the information contained in other cases. 

TIS_sequence:
Sequence of new TIS detected by MORFEE.

TIS_type:
Possible values for the identified predicted TIS are either canonical (ATG) or non-canonical.

Modification_type
Type of codon change induced by the variation.
Of note, deletion of the canonical TIS (ATG) could be found in this column under the label “canonical_TIS(ATG)_to_non_canonical_TIS(…)”

STOP_sequence: 
STOP codon sequence of the new predicted ORF or of the deleted STOP.

TIS_position:
Genomic position of the TIS associated with the newly predicted or altered ORF.

STOP_position
Genomic position of the STOP codon associated with the newly predicted ORF.

Kozak_sequence
Kozak sequence of the new predicted TIS.

Kozak_strength
Kozak strength of the new predicted TIS by taking into account the most conserved positions -3 and +4 regarding the first nucleotide of the TIS.

Kozak_score
Kozak score of the new predicted TIS as calculated by the TIS-predictor algorithm (KSS score).

TIS_c_pos
Coordinates of the first nucleotide of the TIS on the corresponding transcript(s).

STOP_c_pos
Genomic coordinates of the first nucleotide of the TIS.

Additional fields selected from ANNOVAR annotation (http://annovar.openbioinformatics.org/; Wang et al., 2010, https://doi.org/10.1093/nar/gkq603) include:

CLNDN: disease name from [ClinVar](https://www.ncbi.nlm.nih.gov/variation/).
CLNDISDB: database name and identifier for the disease name from ClinVar.
CLNSIG: clinical significance from ClinVar.
gnomad_AF: allele frequency from [gnomAD](https://gnomad.broadinstitute.org/).
gnomad_AF_popmax: maximum population allele frequency from gnomAD.
















