###################################################################
#############                  Script for ld              #########
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
########################
## 01 COMPUTE LD STATISTICS
########################
## here this is for SNPs between 0 and 0.5Mb
gene.plot.region <- which(selscan.map$position>gene.plot.beg&selscan.map$position<=gene.plot.end) 
snp.pos.g <- selscan.map$position[gene.plot.region]
gene_region <- range(snp.pos.g)
ld.haps <- selscan.haps[,gene.plot.region]

## use Gavin's function
d <- compute.ld(t(ld.haps))

## we want to average LD across some increments
## lets try 100kb 
n <- 100
#n <- length(snp.pos.g)
pos.increments <- pretty(snp.pos.g,n)
ld.cols <- rev(heat.colors(n))

ld.subset <- cut(snp.pos.g, pos.increments)

r.mat <- d.mat <- matrix(0,nr = length(levels(ld.subset)),nc = length(levels(ld.subset)))
for(i in 1:nrow(r.mat))
{
  for(j in 1:ncol(r.mat))
  {
    if(j != i)
    {
      row.selection <- ld.subset==levels(ld.subset)[i]
      col.selection <- ld.subset==levels(ld.subset)[j]
      r <- max(d$R[row.selection,col.selection]^2)
      d2 <- max(d$D[row.selection,col.selection])
      if(is.infinite(r)) r <- NA
      if(is.infinite(d2)) d2 <- NA
      r.mat[i,j] <- r
      d.mat[i,j] <- d2
    }
  }
}

no.data <- which(is.na(r.mat[1,]))
diag(r.mat) <-  diag(d.mat) <-  NA
# # ld.mat <- r.mat
# # ld.mat[lower.tri(ld.mat)] <- d.mat[lower.tri(d.mat)]
r.mat[upper.tri(r.mat)] <- NA
d.mat[lower.tri(d.mat)] <- NA



# Make two square regions side by side
plot(0, type = "n", xlim=c(0,ncol(r.mat)),ylim=c(0,nrow(r.mat)))
image(1:ncol(r.mat),
      1:nrow(r.mat),
      r.mat,col = ld.cols,
      axes = F, xlab = "",ylab = "",
      add = F)
# ## ADD D MATRIX
image(1:ncol(r.mat),
      1:nrow(r.mat),
      d.mat,col = ld.cols,
      axes = F, xlab = "",ylab = "",
      add = T)


### BELOW IS A TRIAL VERSION, BUT DIFFICULT TO ADD THINGS TO
# ## plot R^2 matrix at bottom
# source("functions/plotRotatedMatrix.R)
# plot.mat <- r.mat#[-no.data,-no.data]
# plotRotatedMatrix(plot.mat, plot.new = T, col.palette = ld.cols)
# ## plot D on top
# # plotRotatedMatrix(d.mat[-no.data,-no.data], plot.new = F, col.palette = ld.cols)
# 
# ## ADD X AXIS
# d1 <- sqrt(nrow(plot.mat)^2 + ncol(plot.mat)^2)
# d2 <- 0.5*d1
# d2 - d1
# 
# x_at <- seq(0,d1,length = ncol(plot.mat))-0.5
# axis(1,at = x_at, pos = 0)
#   
# # ADD SCALES ## bit of a hack
# x1 <- 0.75
# x2 <- 0.95
# y1 <- 0.8
# y2 <- 0.95
# y3 <- 0.05
# y4 <- 0.2
# # R2 scale below
# par(new = T, fig = c(x1,x2,y3,y4))
# par(mar = c(5,1,4,1))
# scale.matrix <- matrix(seq(0,1,0.1))
# image(seq(0,1,0.1),1,scale.matrix,
#       col = ld.cols, axes = F,
#       xlab = "")
# x_at <- c(0,0.25,0.5,0.75,1)
# axis(1, x_at, labels = x_at )
# title(main = expression(italic(R)^2))
# box()
# # ## D scale on top
# # par(new = T, fig = c(x1,x2,y1,y2))
# # par(mar = c(5,1,4,1))
# # scale.matrix <- matrix(seq(0,1,0.1))
# # image(seq(0,1,0.1),1,scale.matrix,
# #       col = ld.cols, axes = F,
# #       xlab = "")
# # x_at <- c(0,0.25,0.5,0.75,1)
# # axis(1, x_at, labels = x_at )
# # title(main = expression(italic(D)))
# # box()
#   
  
    

plot.three.maps <- FALSE
leftmar <- 6
rightmar <- 1.5
axis.lab.cex <- 1
panel.label.cex <- 1.2
png("figures/AlphaThalLD.png",height = 800,width = 1000, res = 150)  

layout(matrix(c(3,2,1),3,1),heights = c(2,1,5))


###################################
## PLOT LD -- I'M USING R^2 HERE ##
par(mar=c(1,leftmar,1,rightmar))
plot.mat <- r.mat
plot.mat[lower.tri(plot.mat)] <- NA
## empty plot
plot(0,xlim=range(1,ncol(plot.mat)),ylim=c(nrow(plot.mat)/2,1),
     axes= F,type = "n", xlab = "", ylab = "", xaxs = "i",yaxs = "i")
title(ylab = paste0("Max LD in ",diff(pos.increments)[1]," bp window"),
      las = 0, line = 3, cex.lab = axis.lab.cex)
mtext(3,text = "C", adj = -0.1, cex = panel.label.cex, xpd = T, font = 2, line=-1)
## TURN R^2 MATRIX INTO A MATRIX OF COLOURS
x <- matrix(ld.cols[cut(r.mat,length(ld.cols))],nc=ncol(r.mat),byrow=T)

for(i in 1:nrow(x))
{
  for(j in 1:ncol(x))
  {
      r.col <- as.vector(x[i,j])
      if(!is.na(r.col))
      {
        
        # xl <- j - 0.5
        # yb <- i - 0.5
        # xr <- j + 0.5
        # yt <- i + 0.5
        # rect(xl,yb,xr,yt, col = r.col, border = NA)
        xl <- (i+j)/2
        yb <- (j-i)/2 
        xr <- xl + 1
        yt <- yb + 1
      
        polygon(list(x = c(xl,xl+0.5,xr,xr-0.5),y = c(yb+0.5,yt,yt-0.5,yb)),
                col = r.col, border = r.col, xpd =T)
      }
  }
}

scale.matrix <- matrix(seq(0,1,0.1))
y1 <- nrow(plot.mat)*0.4
x1 <- 0
for(i in scale.matrix)
{
  rec.col <- rev(heat.colors(10))[scale.matrix == i]
  rect(x1,y1,x1+1,y1+1, col = rec.col, border = rec.col)
  x1 <- x1 + 1
}
text(x1/2,y1-1,labels = expression(italic(R)^2),
     cex = axis.lab.cex)
text(x1-1,y1+1.5,labels = "1", cex = axis.lab.cex, xpd = T)
text(0,y1+1.5,labels = "0", cex = axis.lab.cex, xpd = T)
# x_at <- c(0,0.25,0.5,0.75,1)
# axis(1, x_at, labels = x_at )
# title(main = expression(italic(R)^2))
# box()

# ## ADD X AXIS
# x_lab <- pos.increments[-1]
# x_lab  <- pretty(x_lab)
# x_lab[1] <- pos.increments[-1][1]
# x_at <- (1:ncol(plot.mat))[pos.increments[-1]%in%x_lab]
# axis(3,at = x_at, labels = paste0(round(x_lab/1e6,2),"Mb"), las = 1)

###################################
## PLOT THE IMPORTANT GENES
genes <- load.genes( "data/UCSC_hg19_2015-08-18_refGene.tsv" )
genes <- genes[genes$cdsStartStat!="unk",]
genes <- genes[genes$txEnd<=range(snp.pos)[2],]
genes <- genes[grep("^HB",genes$name2),]
par(mar=c(0,leftmar,0,rightmar)) 
plot.genes( chromosome = chrom, region = range(pos.increments[-1]), genes,
            xaxt = "n" ,plot.ylab = "", label.cex = 1,
            xaxs = "i", yaxs = "i", axes = F)
title(ylab = "Genes", las = 0, line = 3, cex.lab = axis.lab.cex)
mtext(3,text = "B", adj = -0.1, cex = panel.label.cex, xpd = T, font = 2)
## plot axis?
# x_lab <- pos.increments[-1]
# x_lab  <- pretty(x_lab)
# x_lab[1] <- pos.increments[-1][1]
# axis(3, at = x_lab, labels = paste0(x_lab/1e6,"Mb"), line = 1)
# title(xlab = paste0("Position on chromosome ",chrom))
###################################
## PLOT THE RECRATES
filename <- "data/maps_b37.sqlite"
gm <- loadGeneticMap(chrom,filename=filename)  
par(mar=c(3,leftmar,2,rightmar))
plot(gm$position[gm$chromosome==chrom],gm$COMBINED_LD_rate[gm$chromosome==chrom]*1e6,
     xlim=range(snp.pos.g), type = "S", ylim = c(0,max(rec_y_at)),
     axes = F, xlab = "", xaxs = "i", ylab = "", yaxs = "i")# "Recom.  \nRate(cM/Mb)")
box()
title(ylab = "Recom.\nRate\n(cM/Mb)", las = 0, line = 3, adj = 0.4,cex.lab = axis.lab.cex, xpd = T)
axis(2, at = rec_y_at, labels = rec_y_at, las = 2, cex.axis = axis.lab.cex)

leg.2.plot <- 1
if(plot.three.maps ==  TRUE)
{
  ## PLOT AFRICA MAP [BLUE]
  points(gm$position[gm$chromosome==chrom],gm$YRI_LD_rate[gm$chromosome==chrom]*1e6,type = "S", col = "blue", lty = 2)
  ## PLOT AFRICA ENRICHED MAP [RED]
  points(gm$position[gm$chromosome==chrom],gm$African_Enriched_rate[gm$chromosome==chrom]*1e6,type = "S", col = "red", lty = 3)
  leg.2.plot <- 1:3
}

legend("topleft",xpd = T,inset=c(0,0.1),
       legend=c("HapMap combined","HapMap YRI", "Africa enriched")[leg.2.plot],x.intersp = 0.35,
       col = c("black","blue","red")[leg.2.plot], lty = c(1,2,3)[leg.2.plot], bty = "n", horiz = T, cex = axis.lab.cex)

## PLOT X AXIS
x_at <- pretty(snp.pos.g)
x_at[1] <- snp.pos.g[1]
axis(1,at=x_at,labels=paste0(round(x_at/1e6,2),"Mb"), cex.axis = axis.lab.cex, xpd = T)
title(xlab = paste0("Position on chromosome ",chrom), line = -6.5, cex.lab = axis.lab.cex)
mtext(3,text = "A", adj = -0.1, cex = panel.label.cex, xpd = T, font = 2)

dev.off()


## ADD GENES, SCALE, RECRATES, PLOT LABELS ....


