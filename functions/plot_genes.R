# load genes in tab-delimited refGene.txt-type format with headers
# file should have columns chrom, txStart, txEnd, name2, exonStarts, exonEnds
# This function add a 'chromosome' column that is chromosome formatted differently.
# By default load.genes excludes all but the longest transcript of each gene,
# use condense=FALSE to stop this behaviour.
load.genes <- function( filename, condense = TRUE ) {
    # Load genes
    gene <- read.delim(
	filename,
        header=TRUE,
        as.is=TRUE
    );
	gene <- gene[ -grep( "hap", gene$chrom ), ]
	
	if( condense ) {
	    gene <- gene[order(gene$txEnd - gene$txStart,decreasing=TRUE),];  #Get just longest transcript
	    gene <- gene[ !duplicated( gene$name2 ), ];
	}
    gene <- gene[ !is.na(gene$txStart), ];
	
    chromosome =  gsub( "^chr", "", gene$chrom )
    w1 = which( nchar( chromosome ) == 1 )
    chromosome[ w1 ] = sprintf( "0%s", chromosome[w1] )
    gene$chromosome = chromosome
    return( gene ) ;
}

# Return a data fram of exons for the given set of genes.
# genes should be a dataframe of the type returned by load.genes()
# Exons are found by parsing genes$exonStarts and genes$exonEnds.
get.exons <- function( genes ) {
	result = data.frame()
	for( i in 1:nrow( genes )) {
		gene = genes[i,]
		if( gene$exonCount > 0 ) {
			result = rbind(
				result,
				data.frame(
					name = gene$name,
					name2 = gene$name2,
					txStart = gene$txStart,
					txEnd = gene$txEnd,
					exonStart = as.integer( strsplit( gene$exonStarts, split = "," )[[1]] ),
					exonEnd = as.integer( strsplit( gene$exonEnds, split = "," )[[1]] )
				)
			)
		}
	}
	return( result )
}

# Plot a set of genes in a given region of the genome.
# All genes intersecting the region will be plotted, with exons denoted by rectangles
# and the gene name and transcription direction annotated.
# Arguments:
# - local.genes is a dataframe of genes e.g. as returned by load.genes()
# - chromosome refers to the 'chromosome' column of local.genes.
# - region is a vector of start, end positions.
# - height_in_inches is used to adjust gene and text height, making them (hopefully)
# suitably small when there are many genes in the region.
plot.genes <- function(
    chromosome,
    region,
    local.genes,
    height_in_inches = 1,
    exons = get.exons( local.genes ),
    vertical = FALSE,
    plot.ylab = "",
	colours = list(
		gene = "blue",
		arrows = rgb( 0.5, 0.5, 0.5, 0.7 ),
		exon = "blue",
		text = "black"
	),
	highlight.gene = NULL,
	highlight.gene.col = "red",
	switch.gene.side = FALSE,
	verbose = FALSE,
	label.cex = NULL,
	...
) {
	w = which( local.genes$chromosome == chromosome & local.genes$txEnd >= region[1] & local.genes$txStart <= region[2] )
	if( verbose ) {
		print(w)
	}
    if( length(w) > 0 ) {
		local.genes = local.genes[w,]
		local.genes$lty = 1
		local.genes$label.cex = label.cex
		wPseudoGene = which( local.genes$cdsStart == local.genes$cdsEnd )
		local.genes$lty[ wPseudoGene ] = 2
		local.genes$label.cex[ wPseudoGene ] = label.cex * 0.8
		stopifnot( nrow( local.genes ) > 0 )
        local.genes = local.genes[ order( local.genes$txStart ),, drop = FALSE ]
        local.genes$y = NA ;
        local.genes$y[1] = 1 ;
        if( nrow( local.genes ) > 1 ) {
            spacer = ( region[2] - region[1] ) / 10 ;
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
        
		if( is.null( label.cex )) {
			# We fit about 7 gene names per inch at cex = 1.  After that we need to start scaling.
			label.cex = min( 1, ( 7 * height_in_inches ) / max( local.genes$y ) )
		}
		
        exons$y = local.genes$y[ match( exons$name, local.genes$name )]
		if( verbose ) {
			print( exons )
		}
        height_of_genes = height_in_inches / max( local.genes$y )
        
        strands = c("+","-" ) ;
        level = c( 80, 70 ) ;
		if( vertical ) {
			plot( region[1], 0, pch = '', ylim = region, xlim = c( 0, max( 2, max(local.genes$y)+1 ) ), xlab = '', ylab = '', xaxt = 'n', xaxs = "i", ... ) ;
		} else {
#			plot( 0, region[1], pch = '', xlim = region, ylim = c( 0, max( 2, max(local.genes$y)+0.5 ) ), xlab = '', ylab = plot.ylab, yaxt = 'n', xaxs = "i", ... ) ;
		  plot( 0, region[1], pch = '', xlim = region, ylim = c( 0.5, max( 2, max(local.genes$y)+0.5 ) ), xlab = '', ylab = plot.ylab, yaxt = 'n', xaxs = "i", ... ) ;
		}

        arrow.sep = ( region[2] - region[1] ) / 200 ;
		arrow.voffset = 0
		
        relative.lengths = ( local.genes$txEnd - local.genes$txStart ) / ( region[2] - region[1] )

        local.genes$mark1 = ( 0.25 * local.genes$txStart + 0.75 * local.genes$txEnd ) ;
        local.genes$mark2 = ( 0.5 * local.genes$txStart + 0.5 * local.genes$txEnd ) ;
        local.genes$mark3 = ( 0.75 * local.genes$txStart + 0.25 * local.genes$txEnd ) ;
        local.genes$sign = 1 ;
        local.genes$sign[ which( local.genes$strand == '-' ) ] = -1 ;

		gene.height = 0.4
		exon.height = 0.25

		if( vertical ) {
	        segments(
	            y0 = local.genes$txStart, y1 = local.genes$txEnd,
	            x0 = local.genes$y, x1 = local.genes$y,
	            col = colours[['gene']],
				lty = local.genes$lty
	        )
	        segments(
	            y0 = local.genes$txStart, y1 = local.genes$txStart,
	            x0 = local.genes$y + gene.height, x1 = local.genes$y - gene.height,
	            col = colours[['gene']]
	        )
	        segments(
	            y0 = local.genes$txEnd, y1 = local.genes$txEnd,
	            x0 = local.genes$y + gene.height, x1 = local.genes$y - gene.height,
	            col = colours[['gene']]
	        )
		} else {
	        segments(
	            x0 = local.genes$txStart, x1 = local.genes$txEnd,
	            y0 = local.genes$y, y1 = local.genes$y,
	            col = colours[['gene']],
				lty = local.genes$lty
	        )
	        #segments(
	        #    x0 = local.genes$txStart, x1 = local.genes$txStart,
	        #    y0 = local.genes$y + gene.height, y1 = local.genes$y,
	        #    col = colours[['gene']]
	        #)
#	        segments(
#	            x0 = local.genes$txEnd, x1 = local.genes$txEnd,
#	            y0 = local.genes$y + gene.height, y1 = local.genes$y - gene.height,
#	            col = colours[['gene']]
#	        )
		
			w = which( local.genes$sign == 1 ) 
			arrow.top = exon.height + 0.2
			if( length(w) > 0 ) {
				segment.length = (region[2] - region[1])/120 ;
				segments( x0 = local.genes$txStart[w], x1 = local.genes$txStart[w], y0 = local.genes$y[w] - exon.height, y1 = local.genes$y[w] + arrow.top, col = colours[['gene']], xpd = NA )
		        segments(
		            x0 = c( local.genes$txStart[w], local.genes$txStart[w], local.genes$txStart[w] + segment.length/2, local.genes$txStart[w] + segment.length/2 ),
		            x1 = c( local.genes$txStart[w], local.genes$txStart[w] + segment.length, local.genes$txStart[w] + segment.length, local.genes$txStart[w] + segment.length ),
		            y0 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top + 0.05, local.genes$y[w] + arrow.top - 0.05 ),
		            y1 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top ),
		            col = colours[['gene']],
					xpd = NA
		        )
			}
			w = which( local.genes$sign == -1 ) 
			if( length(w) > 0 ) {
				segment.length = (region[2] - region[1])/120 ;
				segments( x0 = local.genes$txEnd[w], x1 = local.genes$txEnd[w], y0 = local.genes$y[w] - exon.height, y1 = local.genes$y[w] + arrow.top, col = colours[['gene']], xpd = NA )
		        segments(
		            x0 = c( local.genes$txEnd[w], local.genes$txEnd[w], local.genes$txEnd[w] - segment.length/2, local.genes$txEnd[w] - segment.length/2 ),
		            x1 = c( local.genes$txEnd[w], local.genes$txEnd[w] - segment.length, local.genes$txEnd[w] - segment.length, local.genes$txEnd[w] - segment.length ),
		            y0 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top + 0.05, local.genes$y[w] + arrow.top - 0.05 ),
		            y1 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top ),
		            col = colours[['gene']]
		        )
			}
		}
		
		
        wBigEnough = which( ( local.genes$txEnd - local.genes$txStart ) > arrow.sep * 2 ) ;
		if( length(wBigEnough) > 0 ) {
			arrows = data.frame()
			for( i in wBigEnough ) {
				arrows = rbind(
					arrows,
					data.frame(
						name2 = local.genes$name2[i],
						y = local.genes$y[i],
						x = seq( from = local.genes$txStart[i] + arrow.sep, to = local.genes$txEnd[i], by = arrow.sep ),
						sign = local.genes$sign[i]
					)
				)
			}
			
			if( vertical ) {
				segments(
					y0 = arrows$x, y1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					x0 = arrows$y + 0.2 + arrow.voffset, x1 = arrows$y + arrow.voffset,
					col = colours[['arrow']]
				)

				segments(
					y0 = arrows$x, y1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					x0 = arrows$y - 0.2 + arrow.voffset, x1 = arrows$y + arrow.voffset,
					col = colours[['arrow']]
				)

				wExon = which( exons$name2 %in% local.genes$name2[wBigEnough] )
				if( length(wExon) > 0 ) {
					print( exons )
					exons$exonMid = (exons$exonStart + exons$exonEnd)/2
					segmentLength = (region[2] - region[1])/100 ;
					exons$fixedExonStart = pmin( exons$exonStart, exons$exonMid - segmentLength, exon$exonMid + segmentLength ) ;
					exons$fixedExonEnd = pmax( exons$exonStart, exons$exonMid - segmentLength, exon$exonMid + segmentLength ) ;
					rect(
						ybottom = exons$fixedExonStart[wExon],
						ytop = exons$fixedExonEnd[wExon],
						xleft = exons$y[wExon] - 0.25,
						xright = exons$y[wExon] + 0.25,
						border = NA,
						col = colours[['exon']]
					)
				}
			} else {
				segments(
					x0 = arrows$x, x1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					y0 = arrows$y + 0.2 + arrow.voffset, y1 = arrows$y + arrow.voffset,
					col = colours[['arrow']]
				)

				segments(
					x0 = arrows$x, x1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
					y0 = arrows$y - 0.2 + arrow.voffset, y1 = arrows$y + arrow.voffset,
					col = colours[['arrow']]
				)

				wExon = which( exons$name2 %in% local.genes$name2[wBigEnough] )
				if( length(wExon) > 0 ) {
					exons$exonMid = (exons$exonStart + exons$exonEnd)/2
					segmentLength = (region[2] - region[1])/1000 ;
					exons$fixedExonStart = pmin( exons$exonStart, exons$exonMid - segmentLength, exons$exonMid + segmentLength ) ;
					exons$fixedExonEnd = pmax( exons$exonStart, exons$exonMid - segmentLength, exons$exonMid + segmentLength ) ;
					rect(
						xleft = exons$fixedExonStart[wExon],
						xright = exons$fixedExonEnd[wExon],
						ybottom = exons$y[wExon] - 0.25,
						ytop = exons$y[wExon] + 0.25,
						border = NA,
						col = colours[['exon']]
					)
				}
			}
		}
    gene.text.col = rep(colours[[ 'text' ]],length(local.genes$name2))
    if(!is.null(highlight.gene))
    {
      gene.text.col[which(local.genes$name2%in%highlight.gene)] <- highlight.gene.col
    }
        
		if( vertical ) {
			text( local.genes$y, local.genes$txEnd, label = local.genes$name2, adj = c( -0.1, 0.5 ), cex = label.cex, srt=90, col = gene.text.col, font = 3 )
		} else if(switch.gene.side){
		  text( local.genes$txStart, local.genes$y, label = local.genes$name2, adj = c( 1.1, 0.5 ), cex = label.cex, col = gene.text.col, font = 3 )
		} else {
			text( local.genes$txEnd, local.genes$y, label = local.genes$name2, adj = c( -0.1, 0.5 ), cex = label.cex, col = gene.text.col, font = 3 )
		}
    } else {
#        plot.new()
        plot( region[1], 0, pch = '', xlim = region, ylim = c( 0, 1 ), xlab = '', ylab = '', ...) ;
    }
}

# This is a wrapper around the inthinnerator executable.
# It takes a data frame and thins it so that no two picked SNPs are within the specified distance of
# each other.
# Arguments:
# - data: a data frame with (at least) chromosome, position, rsid, alleleA and alleleB columns
# - distance: the distance to thin by, can be specified in physical coords (bp kb Mb) or recombination
# distance (cM), or a combination as in '0.125cM+25kb'
# - rank.column (optional): the name of a column in the data to rank SNPs by.  Prepend the column name with a minus
# sign to rank in ascending order, e.g. rank.column = '-pvalue'.
# - genetic.map: path to genetic map files, in the format available with the 1000 Genomes Phase III release.
# to make sure and detect chromosome properly, use # as a wildcard in the filename for chromosome in these files.
# E.g. inthinnerate( mydata, rank.column = '-pvalue', )
inthinnerate <- function(
	data,
	distance = "0.125cM+25kb",
	rank.column = NULL,
	genetic.map = sprintf( "%s/reference/mathgen.stats.ox.ac.uk/impute/2014-10-07-1000G_Phase3/genetic_map_chr#_combined_b37.txt", Sys.getenv( 'MG_PROJECTDIR' ) ),
	genes = sprintf( "%s/reference/genome-mysql.cse.ucsc.edu/2015-08-18/UCSC_hg19_2015-08-18_refGene.tsv", Sys.getenv( 'MG_PROJECTDIR' ) ),
	excluded.regions = c(),
	extra.args = "",
	picked.only = TRUE
) {
    a = tempfile( pattern = c( "in", "out" ), fileext = ".csv" )
	excluded.region.opts = ""
	if( length( excluded.regions ) > 0 ) {
		excluded.region.opts = sprintf( "-excl-range %s", paste( excluded.regions, collapse = " " ) )
	}
	if( picked.only ) {
		extra.args = paste( extra.args, '-suppress-excluded' )
	}
	if( is.null( rank.column )) {
		write.csv( data[,c( "rsid", "rsid", "chromosome", "position", "alleleA", "alleleB" ), ], file = a[1], row.names = F, quote = F )
	    system( sprintf( 'inthinnerator_v2.0-dev -min-distance %s -map %s -genes %s -g %s -o %s %s %s', distance, genetic.map, genes, a[1], a[2], excluded.region.opts, extra.args ) )
	} else {
		write.csv( data[,c( "rsid", "rsid", "chromosome", "position", "alleleA", "alleleB", rank.column ), ], file = a[1], row.names = F, quote = F )
	    system( sprintf( 'inthinnerator_v2.0-dev -min-distance %s -map %s -genes %s -rank %s -rank-column -%s -g %s -o %s %s %s', distance, genetic.map, genes, a[1], rank.column, a[1], a[2], excluded.region.opts, extra.args ) )
	}
	thinned = read.table( sprintf( '%s.000', a[2] ), hea=T, comment = '#', as.is = T )
	thinned = thinned[ which( thinned$result == "picked" ), ]
	thinned$chromosome = fix_chromosome( thinned$chromosome )
	result = data[ match( paste( thinned$rsid, thinned$chromosome, thinned$position ), paste( data$rsid, data$chromosome, data$position ) ), ]
	result$inthinnerator_result = thinned$result
	result$inthinnerator_pick_index = thinned$pick_index
	result$cM_from_start_of_chromosome = thinned$cM_from_start_of_chromosome
	result$region_lower_bp = thinned$region_lower_bp
	result$region_upper_bp = thinned$region_upper_bp
	result$nearest_gene = sprintf( "%s (%dkb)", thinned$nearest_gene_in_region, thinned$distance_to_nearest_gene_in_region )
	result$nearest_gene[ is.na( thinned$nearest_gene )] = NA
	result$all_genes_in_region = thinned$all_genes_in_region
	return( result )
}

