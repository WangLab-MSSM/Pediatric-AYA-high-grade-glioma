# Author: Francesca Petralia
# Affiliation: Icahn School of Medicine at Mount Sinai
#
# Purpose:
#   Run the BayesDeBulk immune deconvolution workflow used for the Figure 7
#   immune landscape analysis and write the posterior cell-fraction object to
#   the local output directory.

library(BayesDeBulk)
library(stringr)

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, mustWork = TRUE)) else getwd()
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -- load  Data
load(file.path(script_dir, "Data.rda"))

geneID<-pro[,1]
data<-pro[,-seq(1,4)]

index<-(rowSums(is.na(data))==0)
data<-data[index,]
geneID<-geneID[index]

replicates <- function(data,geneID){ 
  
  # -- When gene/protein replicates exist choose the one with highest interquantile range
  # -- data : (p x n) matrix with rowmanes indicating geneID/proteinID
  
  data<-data[!is.na(geneID),]; geneID<-geneID[!is.na(geneID)];
  dup<-duplicated(geneID)  
  geneDup<-unique(geneID[dup])
  dataNew<-data[is.na(match(geneID,geneDup)),]
  geneNew<-geneID[is.na(match(geneID,geneDup))]
  rownames(dataNew)<-geneNew
  
  
  for (j in 1:length(geneDup)){
    data.dup<-data[geneID==geneDup[j],]
    iqr<-apply(data.dup,1,function(x) {quantile(x,.75)-quantile(x,.25)})
    dataNew<-rbind(dataNew,data.dup[iqr==max(iqr),])
  }
  rownames(dataNew)<-c(geneNew,geneDup)
  
  return(dataNew)
}
pro<-replicates(data,geneID)


gene.rna<-rna[,1]
rna<-rna[,-seq(1,3)]
rownames(rna)<-gene.rna

# -- match RNA and proteome data
mg<-match(colnames(pro),colnames(rna))
rna.new<-rna[,mg[!is.na(mg)]]
pro.new<-pro[,!is.na(mg)]

mg<-match(colnames(pro),colnames(rna))
pro.new<-cbind(pro.new,pro[,is.na(mg)])
rna.new<-cbind(rna.new,matrix(NA,dim(rna.new)[1],sum(is.na(mg))))

colnames(rna.new)<-colnames(pro.new)

mg<-match(colnames(rna),colnames(pro))
rna.new<-cbind(rna.new,rna[,is.na(mg)])
pro.new<-cbind(pro.new,matrix(NA,dim(pro.new)[1],sum(is.na(mg))))
colnames(pro.new)<-colnames(rna.new)

rna<-rna.new
data<-pro.new
sampleid<-colnames(rna)
pro.id<-rownames(data)
gene.id<-rownames(rna)

# -- load cell type signatures for deconvolution

load(file = file.path(script_dir, "Gene_signature.rda"))

cell.type<-unique(c(index.matrix[,1],index.matrix[,2]))

# -- run bayesdebulk
n.iter<-5000
burn.in<-1000
k.fix<-length(cell.type)
data<-t(apply(data,1,function(x)(x-mean(x[!is.na(x)]))/sd(x[!is.na(x)])))
rna<-t(apply(rna,1,function(x)(x-mean(x[!is.na(x)]))/sd(x[!is.na(x)])))

gibbs<-BayesDeBulk(n.iter=n.iter,burn.in=burn.in,Y=list(data,rna),markers=index.matrix,prior=NULL) 

pi.post<-gibbs[[1]]

saveRDS(pi.post, file.path(output_dir, "Figure7_BayesDeBulk_pi_post.rds"))
