#!/usr/bin/Rscript

## --------------------
## prints all arguments as msg
## --------------------
catMsg <- function( msg=c() ){	
	cat( MAIN_NAME, paste( msg, collapse="" ), "\n", sep='')
}
## --------------------
## return the location of this script
## --------------------
getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if( length(script.dir) == 0 ) stop("[ERR] Can't determine script dir: please call the script with Rscript\n")
    if( length(script.dir) > 1 ) stop("[ERR] Can't determine script dir: more than one '--file' argument detected\n")
    return(script.dir)
}
## --------------------
## Some html creation functions
## --------------------
htmlTableRow <- function( string_array=c() ){
	td_cells <- ''
	for ( i in string_array ){ 
		td_cells <- paste( td_cells, '<td>', i, '</td>', sep='' )	
	}
	return( paste( "<tr>", td_cells, "</tr>") )
}
htmlLink <- function( path, desc="LINK" ){
	return( paste( '<a href="', path, '">', desc, "</a>", sep='') )
}
## --------------------
## constructs a list with input bam file info
## --------------------
makeBamFileList <- function( paths, names ){	
	tmp <- list()
	l1 <- length(paths)
	l2 <- length(names)
	if ( l1 != l2 ) stop( "Unequal amount of bam-paths (", l1, ") and -names (", l2, ") in makeBamFileList!!!\n" )
	if ( l1 == 0 ){ return(tmp) } # empty list in debug mode

	for ( i in 1:length(paths) ){
		path <- paths[i]
		name <- names[i]
		file <- basename(path)

		tmp[[ file ]] <- name
		tmp[[ 'all_paths' ]] <- c( tmp[[ 'all_paths' ]], path )
		tmp[[ 'all_files' ]] <- c( tmp[[ 'all_files' ]], file )
		tmp[[ 'all_names' ]] <- c( tmp[[ 'all_names' ]], name )
	}
	return( tmp )
}

## --------------------
## copied code for extracting the regions by segment call status
## --------------------
fuseRegions <- function( obj, minRatio=0 ) {
	if ( ncol(obj) > 1 ) stop('Please specify which sample...')

	data <- data.frame( obj@featureData@data[,1:3], copynumber(obj), segmented(obj), check.names=FALSE, stringsAsFactors=FALSE)
	colnames( data ) <- c( "chr", "start", "end", "log2", "segmentval" )
	
	fused.data <- data.frame()
	curr.bin <- 1
	for ( chr in unique( data$chr ) ) {
		chr.data  <- data[ data$chr == chr, ]
		prev.bin  <- curr.bin
		prev.log2 <- chr.data[ 1, 'log2' ]
		prev.segm <- chr.data[ 1, 'segmentval' ]
		start     <- chr.data[ 1, 'start' ]

		if ( nrow(chr.data) > 1) {
			for ( i in 2:nrow(chr.data) ) {
				curr.bin  <- curr.bin + 1
				curr.segm <- chr.data[ i, 'segmentval']

				if ( curr.segm != prev.segm ) {
					fused.data <- rbind( fused.data, data.frame( chr=chr, start=start, end=chr.data[ i-1, 'end'], segmentval=round(prev.segm, digits=DECIMALS) ) )
					prev.segm <- curr.segm
					prev.bin  <- curr.bin
					start     <- chr.data[ i, 'start']
				}
			}
			fused.data <- rbind( fused.data, data.frame( chr=chr, start=start, end=chr.data[ i-1, 'end'], segmentval=round(prev.segm, digits=DECIMALS) ) )
		}else{
			fused.data <- rbind( fused.data, data.frame( chr=chr, start=start, end=chr.data[ i-1, 'end'], segmentval=round(prev.segm, digits=DECIMALS) ) )
		}
	}
	## remove regions with low amplitude
	fused.data <- fused.data[ abs(fused.data$segmentval) >= minRatio, ]
	fused.data
}

## DESC: takes the output of fuse.regions and outputs a txt file per sample
outputRegionsFromList <- function ( regionsList, outputBasename, outputDir="./", binSize, storeList ){
	if ( missing(regionsList) ) stop( 'Please provide regionsList...' )
	if ( missing(outputBasename) ) stop( 'Please provide outputBasename...' )
	if ( !is.list(regionsList) ) stop( 'Input not a list...?' )
	if ( length(regionsList) < 1 ) stop( 'List seems empty...?' )
	if ( file.exists( outputDir ) ) catMsg( c(" Using dir ", outputDir, " for output") )
	else dir.create( outputDir )

	## have to set R output options otherwise scientific method is used at some point
	options( "scipen"=100 )

	sampleCount <- length( regionsList )
	sampleNames <- names( regionsList )
	bedgraphColumns <- c( 'chr', 'start', 'end', 'segmentval' )
	
	catMsg( c( " There are ", sampleCount, " samples found in input list") )

	for ( sample in sampleNames ){		
		catMsg( c(" Working on sample ", sample ) )
		regionCount <- nrow( regionsList[[sample]] )
		
		outSampleBase   <- paste( outputBasename, '_', sample, '_', binSize, 'kbp', sep='')
		outBedgraphFile <- paste( outSampleBase, '.bedGraph', sep="" )
		outBedgraphPath <- paste( outputDir, '/', outBedgraphFile, sep="" )

		## ---------- BEDGRAPH ----------
		txt <- paste( "track type=bedGraph color=0,100,0 altColor=255,0,0 name=", sample," description=segmented_regions_from_QDNAseq_",binSize,"kbp\n", sep="")
		sink( outBedgraphPath )
			cat( txt )
		sink()
		write.table( regionsList[[sample]][,bedgraphColumns], outBedgraphPath, quote=F, sep="\t", row.names=F, append=T, col.names=F)
		#outFiles[[sample]] <- c( outBedgraphFile )
		storeList[[ paste( binSize, sample, 'bedgraph', sep="_")]] <- outBedgraphFile
	}
	return(storeList)
}

## ==================================================
## Unused but potential usefull code
## ==================================================
#@ a bit hacky galaxy way to allow an unknown number of output files based on param selection
#@ see: https://wiki.galaxyproject.org/Admin/Tools/Multiple%20Output%20Files
#historyName <- paste(binSize, 'kbp-IGV', sep="")
#igvFile <- paste( newFilePath, "/primary_", outputId, "_", historyName, "_visible_txt", sep="" )

## ==================================================
## Start of analysis
## ==================================================
MAIN_NAME <- '[INFO] '
catMsg( "Starting QDNAseq wrapper" )	
#catMsg( R.version.string )
catMsg( "Loading R libraries" )

## supress msg to allow R to finish with non-error msg
suppressWarnings( suppressMessages( library( QDNAseq, quietly = TRUE ) ) )
suppressWarnings( suppressMessages( library( CGHcall, quietly = TRUE ) ) )

## only one param: the tmp config file
cmdLineArgs <- commandArgs(TRUE)
config      <- cmdLineArgs[1]
TOOL_PATH   <- cmdLineArgs[2]
CSS_FILE  <- paste( TOOL_PATH, '/static/css/QDNAseq.css', sep="" )
DECIMALS  <- 3
WEB_LINK  <- 'http://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html'
PURE_CSS  <- 'http://yui.yahooapis.com/pure/0.5.0/pure-min.css'

## sourcing the config file will load all input params
## many variables are imported via sourced "config"
source( config )

## if calling requested we always need segmenting first as well
if ( doCall ){ doSegment <- TRUE }

## desparate tries to make png text scale well, damn you R...!
PLOT_RES  <- min( PLOT_WIDTH, PLOT_HEIGHT ) / 6.3 
PAR_SET   <- list( pch=22 )

systemUser <- system("whoami",T)
qdnaseqVersion <- packageDescription( "QDNAseq" )$Version
rVersion <- R.version.string
startTime <- Sys.time()
analysisStart <- as.character( startTime )
catMsg( c("QDNAseq version: ", qdnaseqVersion) )
catMsg( c( rVersion ) )

## get the comma separated list of chromosomes to exclude
excludeChrs <- unlist( strsplit( excludeChrsString, ",") )

## format binSizes back to integers because stupid galaxy doesn't do what I want
#print( binSizesString )
binSizes <- gsub( 'kb', '', binSizesString ) # remove the kb string to get integers
#print( binSizes )
binSizes <- gsub( 'bin', '', binSizes ) # remove the kb string to get integers
#print( binSizes )
binSizes <- as.numeric( unlist( strsplit( binSizes, ",") ) )
#print( binSizes )


## ------------------------
## DEBUG
if ( debug ){
	catMsg( c("Analysis run by user: ", systemUser ) )
	catMsg( c("DEBUG SessionInfo: " ) )
	print( sessionInfo() )
}
## /DEBUG
## ------------------------

## prepare output dir
if ( !file.exists( outputPath) ){
	dir.create( outputPath )
}

## copy source config file to output dir to include it in output zip
if ( inGalaxy ){
	file.copy( config, paste(outputPath, 'galaxyConfigFile.R', sep='/') )	
}

## setup bam filelist for easy retrieval later
fileList    <- makeBamFileList( bamsPaths, bamsNames )
bamCount    <- length( fileList[[ 'all_paths' ]] )

gzipOutputName <- paste( 'QDNAseqResults_', outputName, '.zip', sep='' )
gzipOutputPath <- paste( outputPath, '/', gzipOutputName, sep='')
htmlOutputName <- 'index.html'
htmlOutputPath <- paste( outputPath, '/', htmlOutputName, sep='')

plotted_images <- list() # to keep track of images for later linking
regions <- list() # will contain the segments

## ------------------------
## in case of debug just use inbuilt LGG data for speedup
if ( debug ){ 
	catMsg( c("Built in data only contains binsize 15kb so overriding chosen binSizes to single 15kb") )	
	binSizes <- c(15)
	bamsPaths  <- c( "BUILD_IN_DATA")
	bamsNames  <- c( "LGG150")
	fileList   <- makeBamFileList( bamsPaths, bamsNames )
	bamCount   <- length( fileList[[ 'all_paths' ]] )
}

for ( binSize in binSizes ){

	catMsg( c("Starting analysis for binSize: ", binSize) )	
	## ------------------------
	## construct output file-names and -paths
	## ------------------------
	rdsReadName <- paste( binSize, 'kbp_QDNAseqReadCounts.rds', sep='')
	rdsCopyName <- paste( binSize, 'kbp_QDNAseqCopyNumbers.rds', sep='')
	rdsSegmName <- paste( binSize, 'kbp_QDNAseqCopyNumbersSegmented.rds', sep='')
	rdsCallName <- paste( binSize, 'kbp_QDNAseqCopyNumbersCalled.rds', sep='')
	igvCopyName <- paste( binSize, 'kbp_QDNAseqCopyNumbers.igv', sep='')
	igvSegmName <- paste( binSize, 'kbp_QDNAseqCopyNumbersSegmented.igv', sep='')
	igvCallName <- paste( binSize, 'kbp_QDNAseqCopyNumbersCalled.igv', sep='')

	regiOutputName <- paste( binSize, 'kbp_QDNAseqRegions.rds', sep='')
	noiseImgName   <- paste( binSize, 'kbp_QDNAseqNoiseplot.png', sep='')
	
	rdsRegiPath <- paste( outputPath, '/', regiOutputName, sep='')
	noiseImgPath <- paste( outputPath, '/', noiseImgName, sep='')

	binAnnFile <- paste( TOOL_PATH, '/static/binannotation/', binSize, 'kbp_binAnnotations.rds', sep="" )
	if ( file.exists(binAnnFile) ){
		binAnnotations <- readRDS( binAnnFile )
		catMsg( c("Using local binAnnotations file" ) )
	}else{
		binAnnotations <- getBinAnnotations( binSize=binSize, type=experimentType )
	}

	## in case of debug just use inbuilt LGG data for speedup
	if ( debug ){
		data( LGG150 )
		readCounts <- LGG150
	}else{
		## provide bamnames because in galaxy everyting is called "dataset_###"
		readCounts <- binReadCounts( binAnnotations, bamfiles=fileList[[ 'all_paths' ]], bamnames=fileList[[ 'all_names' ]] )	
	}

	readCountsFiltered    <- applyFilters( readCounts, residual=TRUE, blacklist=filterBlacklistedBins, mappability=mappabilityCutoff, chromosomes=excludeChrs )
	readCountsFiltered    <- estimateCorrection( readCountsFiltered )
	copyNumbers           <- correctBins( readCountsFiltered )
	copyNumbersNormalized <- normalizeBins( copyNumbers )
	copyNumbersSmooth     <- smoothOutlierBins( copyNumbersNormalized )
	sampleNames           <- readCountsFiltered@phenoData@data$name

	## set file to output if output requested
	outputData  <- copyNumbersSmooth
	outputType  <- 'copynumber'
	outputLogT  <- TRUE
	rdsReadPath <- paste( outputPath, '/', rdsReadName, sep='')
	saveRDS( readCounts, rdsReadPath );
	rdsPath <- paste( outputPath, '/', rdsCopyName, sep='')
	igvPath <- paste( outputPath, '/', igvCopyName, sep='')

	## proceed with segmenting / calling if requested
	if ( doSegment ){
		copyNumbersSegmented  <- segmentBins( copyNumbersSmooth, undo.splits=undoSplits, undo.SD=undoSD, transformFun=c("sqrt") )
		copyNumbersSegmented  <- normalizeSegmentedBins( copyNumbersSegmented )
		outputData <- copyNumbersSegmented
		outputType <- 'segments'
		igvPath <- paste( outputPath, '/', igvSegmName, sep='')
		rdsPath <- paste( outputPath, '/', rdsSegmName, sep='')
	}
	if ( doCall ){
		copyNumbersCalled <- callBins(copyNumbersSegmented)
		outputData <- copyNumbersCalled
		outputType <- 'calls'
		outputLogT  <- FALSE # call values should not be transformed at output
		rdsPath <- paste( outputPath, '/', rdsCallName, sep='')
		igvPath <- paste( outputPath, '/', igvCallName, sep='')
	}

	## save the QDNAseq objects and tsv file of highest level (calls or segments)
	saveRDS( outputData, rdsPath );
	exportBins( outputData, file=igvPath, format="igv", type=outputType, logTransform=outputLogT )

	## also save objects for galaxy history output if requested
	if ( txt2history ){
		fileId <- paste('txt_', binSize, sep='')
		historyOutputPath <- historyOutputFiles[[ fileId ]]
		catMsg( c("About to export igv/txt file to history for ", binSize, "kbp bin") )
		exportBins( outputData, file=historyOutputPath, format="igv", type=outputType, logTransform=outputLogT )
	}
	if ( rds2history ){
		fileId <- paste('rds_', binSize, sep='')
		rdsHistoryOutputPath <- historyOutputFiles[[ fileId ]]
		catMsg( c("About to export rds file to history for ", binSize, "kbp bin") )
		saveRDS( outputData, file=rdsHistoryOutputPath )
	}

	## ------------------------
	## create output files
	## ------------------------
	png( noiseImgPath, width=PLOT_HEIGHT, height=PLOT_HEIGHT, res=PLOT_RES );
		par( PAR_SET )
		noisePlot( readCountsFiltered, main=paste( "Noise Plot ", binSize, "kbp", sep=''), col="darkgreen" )
	dev.off()

	binSize <- as.character( binSize ) # to avoid R using it as array index... (*#$^@ you R!)
	binSizeString <- paste( binSize, 'kbp', sep='')
	cgh <- makeCgh( outputData ) # needed for fuseRegions function

	for (i in 1:length(sampleNames) ){
		
		sample <- sampleNames[i]
		usedReads  <- readCountsFiltered@phenoData@data$used.reads[i]
		catMsg( c("Creating plots for sample: ", sample, " (", binSizeString, ")" ) )	

		type <- 'CopyNumbers'
		img_file <- paste( sample, '_', binSize, 'kbp_QDNAseq', type, '.png',  sep='')
		img_file_path <- paste( outputPath, '/', img_file, sep='' )

		## COPYNUMBER PLOT
		png( img_file_path, width=PLOT_WIDTH, height=PLOT_HEIGHT, res=PLOT_RES ); 
			par( PAR_SET )
			plot( copyNumbersSmooth[ ,sample ], main=paste(sample, ": CopyNumbers", sep="") ) 
			mtext( paste( "(", binSizeString, " bins)", sep=""), 3 )
			abline( h=c(-2,-1,1,2,3,4), lty=1, lwd=0.5, col="grey" )
		dev.off()
		
		plotted_images[[ paste(binSize, sample, type, sep="_" ) ]] <- img_file
		
		if ( doSegment ){
			type <- 'Segmented'
			img_file <- paste( sample, '_', binSize, 'kbp_QDNAseq', type, '.png',  sep='')
			img_file_path <- paste( outputPath, '/', img_file, sep='' )
			
			## COPYNUMBER + SEGMENTS PLOT
			png( img_file_path, width=PLOT_WIDTH, height=PLOT_HEIGHT, res=PLOT_RES ); 
				par( PAR_SET )
				plot( copyNumbersSegmented[ ,sample ], main=paste(sample, ": CopyNumbers (Segmented)", sep="") )
				mtext( paste( "(", binSizeString, " bins)", sep=""), 3 )
				abline( h=c(-2,-1,1,2,3,4), lty=1, lwd=0.5, col="grey" )
			dev.off()

			plotted_images[[ paste(binSize, sample, type, sep="_" ) ]] <- img_file

			## if segmented we can also retrieve the segment locations
			catMsg( c(" Fusing regions of sample: ", sample) )
			regions[[ sample ]] <- fuseRegions( cgh[, sample] )

			region_count <- nrow( data.frame( regions[[ sample ]] ) )
			catMsg( c( ' sample "', sample, '" has ', region_count, " regions" ) )
			plotted_images[[ paste(binSize, sample, 'region_count', sep="_" ) ]] <- region_count
		}

		if ( doCall ){
			type <- 'Called'
			img_file <- paste( sample, '_', binSize, 'kbp_QDNAseq', type, '.png',  sep='')
			img_file_path <- paste( outputPath, '/', img_file, sep='' )
			
			## COPYNUMBER + SEGMENTS + CALLS PLOT
			png( img_file_path, width=PLOT_WIDTH, height=PLOT_HEIGHT, res=PLOT_RES ); 
				par( PAR_SET )
				plot( copyNumbersCalled[ ,sample ], main=paste(sample, ": CopyNumbers (Segmented and Called)", sep="") )
				mtext( paste( "(", binSizeString, " bins)", sep=""), 3 )
				abline( h=c(-2,-1,1,2,3,4), lty=1, lwd=0.5, col="grey" )
			dev.off()

			plotted_images[[ paste(binSize, sample, type, sep="_" ) ]] <- img_file
		}
		
		## add USED read counts
		plotted_images[[ paste(binSize, sample, 'usedReads', sep="_" ) ]] <- usedReads
	}

	if ( doSegment ){
		saveRDS( regions, rdsRegiPath )
		plotted_images <- outputRegionsFromList( regions, outputBasename=outputName, outputDir=outputPath, binSize=binSize, storeList=plotted_images )	
	}
}# end bin


## ----- debug -----
#catMsg( "done" )
#q(status=0)
## ---- /debug -----


## ------------------------
## prepare output
## ------------------------
catMsg( "...zipping output")
zip_cmd <- paste( "zip -j", gzipOutputPath, paste(outputPath,'/*',sep='') ) ## -j is for removing dirs from the tree
system( zip_cmd )

## ------------------------
## get filesizes for report
## ------------------------
zippedSize <- paste( round( file.info( gzipOutputPath )[["size"]] / 1e+6, digits=2 ), 'MB' )
endTime <- Sys.time()
timeDiff <- format( round( endTime - startTime, 3 ) )
analysisEnd <- as.character( endTime )

## ------------------------
## creating html output to be linked to from the middle galaxy pane
## ------------------------
sink( file = htmlOutputPath, type = "output" )
		cat( "<html>\n")
		cat( "<head>\n")

			cat( "\t", '<title>QDNAseq Report | ', outputName,'</title>', "\n", sep='' )
			cat( "\t", '<link rel="stylesheet" href="', PURE_CSS, '">', "\n", sep='' )
			cat( "\t<style>\n", sep='')
				## include CSS into html file, makes it more portable
				cat( "\t\t", readLines( CSS_FILE ), sep="\n\t\t" )
				#cat( "\t\th1 {color:red;}", "\n")
			cat( "\n\t</style>\n" )
			
		cat( "\n</head>\n")
		cat( "\n<body>\n")

		cat( "<h1>QDNAseq Report</h1>", "\n")
		
		cat( '<h3 class="qdnaseq">About this analysis</h3>', "\n")
		cat( '<p>This page provides access to all results. To have a local copy of this report just download the <a href="', gzipOutputName, '" class="button">zipfile</a> with all output (', zippedSize, ')</p>', "\n", sep='')		
		
		## ------------------------
		## table with general info
		## ------------------------
		cat( '<h3 class="qdnaseq">Settings</h3><p>', "\n")
		cat( '<table class="pure-table pure-table-striped"><tbody>' )
			cat( htmlTableRow( c( "AnalysisName", outputName ) ) )
			cat( htmlTableRow( c( "AnalysisStart", analysisStart ) ) )
			cat( htmlTableRow( c( "AnalysisEnd", analysisEnd ) ) )
			cat( htmlTableRow( c( "AnalysisTime", timeDiff ) ) )
			cat( htmlTableRow( c( "BinSizes (kbp)", paste(binSizes,collapse=", ") ) ) )
			cat( htmlTableRow( c( "R info", rVersion ) ) )
			cat( htmlTableRow( c( "QDNAseq info", qdnaseqVersion ) ) )
			
			sampleStrings <- c()
			for ( galaxyName in fileList[[ 'all_files' ]] ){
				sampleName <- fileList[[ galaxyName ]]
				sampleStrings <- c( sampleStrings, paste( galaxyName, ' (', sampleName, ')', sep='' ) )
			}
			cat( htmlTableRow( c( "InputBams", paste( sampleStrings, collapse=", ") ) ) )

		cat( "</tbody></table></p>", "\n")
		
		## ------------------------
		## list with links to all output files
		## ------------------------
		cat( '<h3 class="qdnaseq">Output files</h3><p>', "\n")
		cat( '<p>This table contains output files that can be used for local downstream analysis with the bioconductor QDNAseq package. For each bin-size / data-level there is a R data structure file with data of all samples. See ', htmlLink( WEB_LINK, 'the bioconductor QDNAseq documentation' ), ' for more information on how to work with these files.</p>', "\n", sep='')
		cat( '<table class="pure-table pure-table-striped">', "\n" )
		cat( '<thead><th>Type</th>', as.vector( mapply( paste, "<th>", binSizes, "kbp</th>", sep="" ) ),'</thead>', "\n" )
		cat( "<tbody>", "\n")
			files <- list()
			#fileTypes <- c( 'ReadCounts.rds', 'CopyNumbers.rds' )
			fileTypes <- c( 'ReadCounts.rds' )
			if ( doCall ){ 
				fileTypes <- c( fileTypes, 'CopyNumbersCalled.rds') 
			}else if ( doSegment ){ 
				fileTypes <- c( fileTypes, 'CopyNumbersSegmented.rds') 
			}else { 
				fileTypes <- c( fileTypes, 'CopyNumbers.rds') 
			}

			for ( fileType in fileTypes ){
				fileNames <- mapply( paste, binSizes, paste( 'kbp_QDNAseq', fileType, sep=''), sep='')
				fileLinks <- mapply( htmlLink, fileNames, paste( binSizes, "kbp", sep="" ) )
				cat( htmlTableRow( c( fileType, fileLinks ) ) )	
			}
		cat( "\n</tbody></table></p>", "\n")

		## ------------------------
		## table with links to files	
		## ------------------------
		ratio <- PLOT_WIDTH / PLOT_HEIGHT
		width <- 960; height <- width / ratio ## bigger img
		width_t <- 100; height_t <- 40 ## thumb img

		cat( '<h3 class="qdnaseq">Results: overview</h3><p>', "\n")
		cat( '<p>This table contains the visual results of the copy number aberration analysis. You can click on an image to jump to the larger version. If segmentation was performed as well the number of segments is shown and a file with genomic regions can be downloaded (just remember to inspect the results carefully as this is, together with optional calling afterwards, a more experimental type of analysis).</p>', "\n", sep='')
		plots_html <- ''

		colspan <- 1
		binHeader <- "<th>Image</th>"
		if ( doSegment ){ # extra column with segment info
			colspan <- 2 
			binHeader <- "<th>Image</th><th>Segments</th>"
		} 
		cat( '<table class="pure-table pure-table-striped">', "\n" )
		cat( '<thead><tr><th></th><th></th>', as.vector( mapply( paste, "<th colspan=\"", colspan,"\">", binSizes, "kbp</th>", sep="" ) ), '</tr></thead>' )
		cat( '<thead><tr><th>Sample / File</th><th>Reads</th>', rep( binHeader, length(binSizes) ), '</tr></thead>' )
		cat( '<tbody>' )

			for ( bam_file in bamsNames ){
				
				usedReads <- plotted_images[[ paste(binSize, bam_file, 'usedReads', sep="_" ) ]]
				usedReads <- format( as.integer(usedReads), digits=4, decimal.mark=".", big.mark="," )
				htmlRow <- paste( '<tr><td>', bam_file, '</td><td>', usedReads, '</td>', sep='' )

				for ( binSize in binSizes ){
					
					## add thumbnails to table with links to anchors on html page
					copy_img <- plotted_images[[ paste(binSize, bam_file, 'CopyNumbers', sep="_" ) ]]
					html_copy_thumb <- htmlLink( path=paste('#', copy_img, sep=''), paste('<img src="',copy_img,'" alt="', bam_file, '" width="', width_t, '" height="', height_t, '">', sep='') )
					html_copy_img <- htmlLink( path=copy_img, paste('<img id="', copy_img,'" src="',copy_img,'" alt="',bam_file, '" width="', width, '" height="', height, '">', sep='') )
					html_segm_img <- ''
					html_call_img <- ''
					html_bedGraph <- ''
					region_count <- ''
					htmlRow <- paste( htmlRow, '<td>', html_copy_thumb, '</td>' )

					if ( doSegment ){
						segm_img <- plotted_images[[ paste(binSize, bam_file, 'Segmented', sep="_" ) ]]
						region_count <- plotted_images[[ paste(binSize, bam_file, 'region_count', sep="_" ) ]]

						html_bedGraph <- htmlLink( path=plotted_images[[ paste(binSize, bam_file, 'bedgraph', sep="_" ) ]], 'bedGraph' )
						html_segm_img <- htmlLink( path=segm_img, paste('<img id="', segm_img,'" src="', segm_img,'" alt="', bam_file, '" width="', width, '" height="', height,'">', sep='') )
						htmlRow <- paste( htmlRow, '<td>', region_count, ' (', html_bedGraph, ')</td>', sep="" )
					}
					if ( doCall ){
						call_img <- plotted_images[[ paste(binSize, bam_file, 'Called', sep="_" ) ]]
						html_call_img <- htmlLink( path=call_img, paste('<img id="', call_img,'" src="', call_img,'" alt="', bam_file, '" width="', width, '" height="', height,'">', sep='') )
					}
					plots_html <- paste( plots_html, html_copy_img, "\n", html_segm_img, "\n", html_call_img, "\n<br \\>\n", sep='' )	
				}
				plots_html <- paste( plots_html, "\n<hr \\>\n", sep='' )
				## add info to overview table, including small thumbnails
				htmlRow <- paste( htmlRow, '</tr>', sep='' )
				cat( htmlRow, "\n" )				
			}
		cat( "</tbody></table></p>", "\n")
		
		## ------------------------
		## section with various output shown
		## ------------------------
		cat( '<h3 class="qdnaseq">Results: Sample plots</h3><p>', "\n")
		## now include (large) images in html page
		cat( plots_html, "\n")
		cat( "\n</p></body>\n")
		cat( "\n</html>\n")
sink()

## ------------------------
## creating main html output for galaxy history
## ------------------------
if ( inGalaxy ){ # dont create when running outside Galaxy
	sink( file = outputHtml, type = "output" )
			
		cat( "<head>", "\n")
			cat( "\t", '<link rel="stylesheet" href="', PURE_CSS, '">', "\n", sep='' )

			cat( "<style>", "\n")
				## include CSS directly into html file
				cat( paste( "\t", '/* the css here originates from ', CSS_FILE,' */', "\n") )
				cat( paste( "\t", readLines( CSS_FILE, n = -1)), sep="\n" )
			cat( "</style>", "\n")
		cat( "</head>", "\n")

		cat( '<h1>QDNAseq Results (', outputName,')</h1>', "\n", sep="")
		cat( '<p>Explore <a href="', htmlOutputName, '" class="button">the results</a> directly within galaxy</p>', "\n", sep="")
		cat( '<p>Or download a <a href="', gzipOutputName, '" class="button">zipfile</a> with all output (', zippedSize, ')</p>', "\n", sep="" )

	sink()
}

## ------------------------
## create final zip and quit with status 0 to tell galaxy all was fine
## ------------------------
catMsg( "zipping all output")
system( paste( "zip -j ", gzipOutputPath, paste(outputPath,'/', htmlOutputName, sep='') ) )
catMsg( "done" )
q(status=0)
