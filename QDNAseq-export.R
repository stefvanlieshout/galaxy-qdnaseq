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
catMsg( "Starting QDNAseq-export wrapper" )

## supress msg to allow R to finish with non-error msg
catMsg( "Loading R libraries" )
suppressWarnings( suppressMessages( library( QDNAseq, quietly = TRUE ) ) )
suppressWarnings( suppressMessages( library( CGHcall, quietly = TRUE ) ) )

## only one param: the tmp config file
cmdLineArgs <- commandArgs(TRUE)
config      <- cmdLineArgs[1]

## sourcing the config file will load all input params
## many variables are imported via sourced "config"
source( config ) # outputFile, outputName, outputFormat, sampleIndex, filterBins

systemUser <- system("whoami",T)
qdnaseqVersion <- packageDescription( "QDNAseq" )$Version
rVersion <- R.version.string
rPath <- R.home()
catMsg( c("QDNAseq version ", qdnaseqVersion) )
catMsg( c( rVersion ) )
qdnaseqObject <- readRDS( rdsFilePath )
logTransform <- TRUE

sampleNames <- sampleNames( qdnaseqObject )
elements <- assayDataElementNames(qdnaseqObject)
element <- dataLevel
if (dataLevel == "segments") element <- "segmented"
if (element == "calls") logTransform <- FALSE

## sanity checks
if ( ! element %in% elements ) stop( paste( "Data-level \"", element, "\" not present in object", sep='') ) 
if ( outputFormat == "bed" ){
	if ( sampleIndex == "None") stop("Bed option requires sample index")
	if ( sampleIndex > length(sampleNames) ) stop("Chosen sample index not present in object")
	qdnaseqObject <- qdnaseqObject[ ,sampleIndex]
}

## output
exportBins( qdnaseqObject, 
	file=outputFile, 
	format=outputFormat, 
	filter=filterBins, 
	type=dataLevel, 
	logTransform=logTransform 
)

## tell galaxy all seems ok
q(status=0)
