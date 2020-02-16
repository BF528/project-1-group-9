library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(tidyverse)
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
normalize:  If TRUE performs quanitle normalization where each array is rank ordered 
(each array tests multiple genes and so gene expression is ordered), average values for 
most-expressed, 2nd most expressed etc are caclulated across arrays. Then the array readout is
replaced by what you would expect for a gene expression of that rank
background:  If TRUE background correct using RMA background correction. This corrects for the 
noise in the readout due to light from nearby cells in the microarray, and other such background noise
"
Probe_Set <- fitPLM(affydata,normalize=TRUE, background=TRUE) #this fits probe-level models to all of the probesets.
image(Probe_Set,which=2, type="resids") #output can be visualized to check data quality

##Plotting Relative Log expression and NUSE"
RLE(Probe_Set,main="RLE for the cancer dataset")
"indeed there are some arrays that have boxes not centered around zero, and with spread out boxes"
RLE_stats<-data.frame(RLE(Probe_Set, type="stats"))
"RLE_stats gives us file names as column and median as rows, so we transpose it
note that as.data.frame is required to transpose it. data.frame somehow reverses the transpose operation"
RLE_stats<-as.data.frame(t(RLE_stats) ) 
hist(RLE_stats$median, xlab="median RLE values")
"the histogram tells us which points have median RLE that are rather high."

NUSE(Probe_Set,main="NUSE for the cancer dataset")
NUSE_stats<-as.data.frame(t(NUSE(Probe_Set,type="stats")))
hist(NUSE_stats$median,xlab="median NUSE values")



metadata_file='/project/bf528/project_1/doc/proj_metadata.csv'
metadata<-read.csv(metadata_file)
expression_matrix<-exprs(normalized_expression) ### we need a matrix of expression times samples for ComBat
combat_output<-ComBat(dat=expression_matrix, batch=matrix(metadata$normalizationcombatbatch),  mod =matrix(metadata$normalizationcombatmod))
combat_output<-ComBat(dat=expression_matrix, batch=as.vector(metadata$normalizationcombatbatch),  mod =matrix(metadata$normalizationcombatmod))

