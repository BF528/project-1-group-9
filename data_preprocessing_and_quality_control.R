library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
"instructions at https://bf528.readthedocs.io/en/latest/content/projects/project_1_microarrays/project_1_microarrays.html"
## load the data file names
data_files<-list.files('/projectnb/bf528/users/group9/project1/samples/',full=TRUE)
affydata<-ReadAffy(filenames=data_files)

"what does rma actually do? normalize across probesets?"
normalized_expression<-rma(affydata)

"for computing Relative Log Expression (RLE) and Normalized Unscaled Standard Error (NUSE),
following this link:
https://www.bioconductor.org/packages/release/bioc/vignettes/affyPLM/inst/doc/QualityAssess.pdf"
###### fitPLM takes along time, and generates a large data file

"options used in fitPLM are:
normalize:  If TRUE background correct using RMA background correction
background:  If TRUE background correct using RMA background correction"
Probe_Set <- fitPLM(affydata) #this fits probe-level models to all of the probesets.
image(Probe_set,which=2, type="resids") #output can be visualized to check data quality


RLE(Probe_set,main="RLE for the cancer dataset")
NUSE(Probe_set,main="NUSE for the cancer dataset")

