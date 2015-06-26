#!/usr/bin/Rscript

## --------------------
## prints all arguments as msg
## --------------------
catMsg <- function( msg=c() ){	
	cat( MAIN_NAME, paste( msg, collapse="" ), "\n", sep='')
}


## ==================================================
## Start of analysis
## ==================================================
MAIN_NAME <- '[INFO] '
catMsg( "Starting QDNAseq-plot wrapper" )
catMsg( "Loading R libraries" )

## supress msg to allow R to finish with non-error msg
suppressWarnings( suppressMessages( library( QDNAseq, quietly = TRUE ) ) )

## only one param: the tmp config file
cmdLineArgs <- commandArgs(TRUE)
config      <- cmdLineArgs[1]

## sourcing the config file will load all input params
## many variables are imported via sourced "config"
source( config ) # outputPngPath, outputPdfPath, allOrOne, rdsFilePath
#cat( "ALL? ", allOrOne, sep='' )

## desparate tries to make png text scale well, damn you R...!
PLOT_RES  <- min( PLOT_WIDTH, PLOT_HEIGHT ) / 6.3 
PAR_SET   <- list( pch=22 )
systemUser <- system("whoami",T)
qdnaseqVersion <- packageDescription( "QDNAseq" )$Version
rVersion <- R.version.string
catMsg( c("QDNAseq version: ", qdnaseqVersion) )
catMsg( c( rVersion ) )

qdnaseqObject <- readRDS( rdsFilePath )
chromosomesToPlot <- unlist( strsplit( chromosomesToPlotString, ",") )

#cat( "CHROM: ", chromosomesToPlotString, "\n" )
#cat( "REGION: ", regionToPlotString, "\n" )
#cat( "What to plot: ", whatToPlot, "\n" )

## COPYNUMBER PLOT
sample <- SAMPLE_INDEX
png( outputPngPath, width=PLOT_WIDTH, height=PLOT_HEIGHT, res=PLOT_RES )
	par( PAR_SET )
	ylab_text <- "log2 read counts"
	if ( whatToPlot == 'everything' ){
		catMsg( c( "Plotting all data in object" ) )
		plot( qdnaseqObject[ ,sample ], ylab=ylab_text )
		abline( h=c(-4,-3,-2,-1,1,2,3,4), lty=1, lwd=0.5, col="grey" )
	} else if( whatToPlot == 'chromosomes' ){
		catMsg( c( "Plotting subset of chromosomes" ) )
		fdata <- qdnaseqObject@featureData@data
		idx_region <- which( fdata$chromosome %in% chromosomesToPlot )
		plot( qdnaseqObject[ idx_region, sample ], ylab=ylab_text )
		abline( h=c(-4,-3,-2,-1,1,2,3,4), lty=1, lwd=0.5, col="grey" )
	} else if( whatToPlot == 'region' ){
		regionC <- chrName
		regionS <- chrStart
		regionE <- chrEnd
		if ( regionS > regionE ) stop("Chosen start is > end")
		catMsg( c( "Plotting genomic region (chr=", regionC, " start=", regionS, " end=", regionE, ")" ) )
		
		fdata <- qdnaseqObject@featureData@data
		idx_region <- which( fdata$chromosome == regionC & fdata$start > regionS & fdata$end < regionE )
		cat( idx_region, "\n")
		
		plot( qdnaseqObject[ idx_region, sample ], doCalls=FALSE, ylab=ylab_text )
	}
	#mtext( "plotted in galaxy", 3 )
	
dev.off()


## all ok
q(status=0)
