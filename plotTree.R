###################################################################
#############                  Script for tree            #########
###################################################################
###################################################################
## LOAD ANY EXTERIOR FUNCTIONS AND LIBRARIES
library("RColorBrewer")
library("grid")
source("functions/plot_genes.R")
source("functions/loadGeneticMap.R")
source("functions/line2user.R")
source("functions/compute.ld.R")
## 

## DEFINE PARAMETERS
chrom <- "16"
########################
## 00 READ IN DATA
## UNCOMMENT TO READ IN ONCE

# #Read in gen and samples files
# chr16.gen <- read.table (paste0("data/Kenya-2.5M_chr",chrom,":0-10000000_nature_2015_qc_set_merged.haps"), header=F)
# samples <- read.table (paste0("data/Kenya-2.5M_chr",chrom,":0-10000000_nature_2015_qc_set_merged.sample"), header=T)
# map <- read.table(paste0("data/genetic_map_chr",chrom,"_combined_b37.txt"),header = T, as.is = T)
# 
# 
# ## CONVERT TO SELSCAN FORMAT
# selscan.haps <- t(chr16.gen[,7:ncol(chr16.gen)])
# selscan.map <- chr16.gen[,c(1,3,4)]
# colnames(selscan.map) <- c("chrom","rsid","position")
# # need to interpolate genetic distance and this were copied from the cluster
# head(map)
# head(selscan.map)
# gen_dist <- approx(map$position, map$Genetic_Map.cM., selscan.map$position  )$y


hap.plot.beg <- 0
hap.plot.end <- 2e6
## gene/recrates plot region
gene.plot.beg <- 0
gene.plot.end <- 0.5e6

########################
## DEFINE REGION TO PLOT HAPLOTYPES
## here this is for all snps between 0 and 2Mb
hap.plot.region <- which(selscan.map$position>hap.plot.beg&selscan.map$position<=hap.plot.end)
## DEFINE REGION TO ZOOM IN FOR GENE PLOT AND RECRATES
## here this is for SNPs between 0 and 0.5Mb
gene.plot.region <- which(selscan.map$position>gene.plot.beg&selscan.map$position<=gene.plot.end) 

## make a phlogenetic tree of thal haplotypes
# plot individuals with alpha_thal
a.thal.haps <- selscan.haps[selscan.haps[,72]==1,hap.plot.region]
a.thal.samps <- as.character(samples$ID_1[-1])
a.thal.eth <- as.character(samples$curated_ethnicity[-1])
tmp <- c()
for(i in a.thal.samps) tmp <- c(tmp,i,i)
a.thal.samps <- tmp
tmp <- c()
for(i in a.thal.eth) tmp <- c(tmp,i,i)
a.thal.eth <- tmp


## LET'S TRY USING APE ...
## further options here:http://www.molecularecologist.com/2016/02/quick-and-dirty-tree-building-in-r/ 
# install.packages("ape")
# install.packages("phangorn")
# install.packages("seqinr")
library(ape)
library(phangorn)
library(seqinr)

## first we need to generate a distance matrix
thal.dat <- as.phyDat(a.thal.haps, type = "USER",
                      levels = c(0,1), names = a.thal.eth[selscan.haps[,72]==1])
thal_dist <- dist.ml(thal.dat)

## NOW TRY A FEW TREES
thal_NJ <- nj(thal_dist)
plot(thal_NJ)

thal_UPGMA <- upgma(thal_dist)
plot(thal_UPGMA)

## ethnic group of tip
new.eth <- a.thal.eth[as.numeric(gsub("V","",thal_UPGMA$tip.label))]
main.eths <- c("GIRIAMA","CHONYI","KAUMA")
new.eth[is.na(new.eth)] <- "OTHER"
new.eth[!new.eth%in%main.eths] <- "OTHER"

## MAKE A TABLE FOR TREE TIP LABELS
tip.col.palette <- c("hotpink","darkgreen","cornflowerblue","grey50")
tip.cols <- tip.col.palette[as.numeric(as.factor(new.eth))]
tree_labels <- cbind(new.eth,tip.cols)
colnames(tree_labels) <- c("eth","col2plot")
tree_labels <- as.data.frame(tree_labels)

tree2plot <- as.phylo(thal_UPGMA)
tree2plot$tip.label <- rep("-",length(new.eth)) 


## colour edges/
edge_cols <- c()
for(i in 1:nrow(tree2plot$edge))
{
  cl <- "#000000"
  edge_cols <- c(edge_cols,cl)
}

png("figures/AlphaThalHaplotypesUPGMATree.png",height = 2000,width =500, res = 150)
par(mar=c(1,3,1,1))
plot(tree2plot,lab4ut="axial",type="phylogram",xaxs="i",yaxs="i",
     edge.color=edge_cols,tip.color=as.character(tree_labels$col2plot),use.edge.length=F)
dev.off()

