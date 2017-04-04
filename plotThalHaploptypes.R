###################################################################
#############         Script for plotting haplotypes      #########
###################################################################
###################################################################
## LOAD ANY EXTERIOR FUNCTIONS AND LIBRARIES
library("RSQLite") 
source("functions/plot_genes.R")
source("functions/loadGeneticMap.R")
source("functions/line2user.R")
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
## 01 PLOT
########################
## DEFINE SOME PLOTTING PARAMETERS
## haplotype colours
hap.cols <- c("white","black")
## colour for gene region highlight
gene.region.col <- "red"
## axis label size
axis.lab.cex <- 1
panel.label.cex <- 1.2
## haplotype plot region
hap.plot.beg <- 0
hap.plot.end <- 2e6
## gene/recrates plot region
gene.plot.beg <- 0
gene.plot.end <- 0.5e6
## y axix limits for recrates
rec_y_at <- c(0,5,10)

## should you plot all three recombination maps
plot.three.maps <- FALSE


########################
## DEFINE REGION TO PLOT HAPLOTYPES
## here this is for all snps between 0 and 2Mb
hap.plot.region <- which(selscan.map$position>hap.plot.beg&selscan.map$position<=hap.plot.end)
## DEFINE REGION TO ZOOM IN FOR GENE PLOT AND RECRATES
## here this is for SNPs between 0 and 0.5Mb
gene.plot.region <- which(selscan.map$position>gene.plot.beg&selscan.map$position<=gene.plot.end) 


png("figures/AlphaThalHaplotypes.png", width = 1000, height = 1200, res=150)

  layout(matrix(c(1,2,4,3),4,1, byrow=T),heights = c(4,4,2,2))
  snp.pos <- selscan.map$position[hap.plot.region]
  snp.pos.g <- selscan.map$position[gene.plot.region]
  plot_range <- range(snp.pos)
  gene_region <- range(snp.pos.g)
 
  par(mar=c(0,6,4,1))
  # plot individuals with alpha_thal
  a.thal.haps <- selscan.haps[selscan.haps[,72]==1,hap.plot.region]
  a.thal.haps <- a.thal.haps#[1:100,]
  image(snp.pos,
        1:nrow(a.thal.haps),
        t(a.thal.haps), col=hap.cols,
        axes = F, xlab = "", ylab = "",
        xaxs = "i", yaxs = "i")
  # add line
  abline(v=which(selscan.map$rsid[hap.plot.region]=="alphathal"),col="blue",lwd=3)
  abline(v = gene.plot.region, lwd = 2, xpd = F, col = "red")
  xat <- pretty(snp.pos)
  xat[1] <- selscan.map$position[1]
  axis(3, at = xat, labels = paste0(round(xat/1e6,1),"Mb"), xpd = T)
  segments(gene_region[1], nrow(a.thal.haps),gene_region[1],-500,col = gene.region.col, lwd = 2, xpd = T)
  segments(gene_region[2], nrow(a.thal.haps),gene_region[2],-500,col = gene.region.col, lwd = 2, xpd = T)
  title(ylab = expression(alpha~haplotypes))
  #mtext(3, text = paste0("Position on chromosome ", chrom), line = 2)
  mtext(3,text = "A", adj = -0.1, cex = panel.label.cex, xpd = T, font = 2)
  
  # plot individuals without alpha_thal
  n.thal.haps <- selscan.haps[selscan.haps[,72]==0,hap.plot.region]
  n.haps <- sample(1:nrow(n.thal.haps), size = nrow(a.thal.haps))
  n.thal.haps <- n.thal.haps[n.haps,]
  par(mar=c(3,6,1,1))
  image(snp.pos,
        1:nrow(n.thal.haps),
        t(n.thal.haps), col=hap.cols,
        axes = F, xlab = "", ylab = "",
        xaxs = "i", yaxs = "i")
  # add lines
  abline(v=which(selscan.map$rsid[hap.plot.region]=="alphathal"),col="blue",lwd=3)
  segments(gene_region[1], nrow(n.thal.haps)+500,gene_region[1],-500,col = gene.region.col, lwd = 2, xpd = T)
  segments(gene_region[2], nrow(n.thal.haps)+500,gene_region[2],-500,col = gene.region.col, lwd = 2, xpd = T)
  title(ylab = expression(-~haplotypes))
  mtext(3,text = "B", adj = -0.1, cex = panel.label.cex, xpd = T, font = 2)
  #######################################################
  ### Code to load genetic maps from Gav
  filename <- "data/maps_b37.sqlite"
  gm <- loadGeneticMap(chrom,filename=filename)  
  par(mar=c(3,6,0,1))
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
  axis(1,at=x_at,labels=paste0(round(x_at/1e6,1),"Mb"), cex.axis = axis.lab.cex, xpd = T)
  title(xlab = paste0("Position on chromosome ",chrom), line = 2, cex.lab = axis.lab.cex)
  mtext(3,text = "D", adj = -0.1, cex = panel.label.cex, xpd = T, font = 2)
  #######################################################
  ### PLOT GENES
  genes <- load.genes( "data/UCSC_hg19_2015-08-18_refGene.tsv" )
  genes <- genes[genes$cdsStartStat!="unk",]
  genes <- genes[genes$txEnd<=range(snp.pos)[2],]
  # genes <- genes[grep("^HB",genes$name2),]
  par(mar=c(1,6,1.5,1)) 
  plot.genes( chromosome = chrom, region = range(snp.pos.g), genes,
              xaxt = "n" ,plot.ylab = "", label.cex = 1,
              xaxs = "i", yaxs = "i", axes = F)
  title(ylab = "Genes", las = 0, line = 3, cex.lab = axis.lab.cex)
 
  ## PLOT CONNECTING LINES
  w <- which( genes$chromosome == chrom & genes$txEnd >= gene_region[1] & genes$txStart <= gene_region[2] )
  
  local.genes = genes[w,]
  local.genes$y = NA ;
  local.genes$y[1] = 1 ;
  if( nrow( local.genes ) > 1 ) {
    spacer = ( gene_region[2] - gene_region[1] ) / 10 ;
    maxes = ( local.genes[1,]$txEnd + spacer )
    for( i in 2:nrow( local.genes )) {
      for( l in 1:length( maxes )) {
        if( local.genes$txStart[i] >= maxes[l] ) {
          local.genes$y[i] = l ;
          maxes[l] = local.genes$txEnd[i] + spacer ;
          break ;
        }
      }
      if( is.na( local.genes$y[i] )) {
        maxes = c( maxes, local.genes$txEnd[i] + spacer )
        local.genes$y[i] = length( maxes ) ;
      }
    }
  }
  ymax <- max( 2, max(local.genes$y)+0.5 )
  
  ## add a plot on top 
  par(mar=c(1,6,1,1)) 
  par(new = T)
  plot(0,0,axes = F, ylim=c(0,ymax), type = "n",xlab = "", ylab = "",
       xlim=range(snp.pos), xaxs = "i", yaxs = "i")
  box(bty = "u", col = gene.region.col, lwd = 2)
  ## ADD LINES
  segments(plot_range[1], ymax, x1 = gene_region[1], y1 = line2user(1,3), col = gene.region.col ,xpd = T, lwd = 2)
  segments(plot_range[2], ymax, x1 = gene_region[2], y1 = line2user(1,3), col = gene.region.col ,xpd = T, lwd = 2)

  ## ## and lines around plot
  segments(plot_range[1], ymax,plot_range[1],0,col = gene.region.col, lwd = 2, xpd = T)
  segments(plot_range[2], ymax,plot_range[2],0,col = gene.region.col, lwd = 2, xpd = T)
  mtext(3,text = "C", adj = -0.1, font = 2, cex = panel.label.cex, xpd = T)
 dev.off()

