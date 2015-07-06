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
suppressWarnings( suppressMessages( library( CGHregions, quietly = TRUE ) ) )

## only one param: the tmp config file
cmdLineArgs <- commandArgs(TRUE)
config      <- cmdLineArgs[1]

## sourcing the config file will load all input params
## many variables are imported via sourced "config"
source( config ) # outputFile, outputName

systemUser <- system("whoami",T)
qdnaseqVersion <- packageDescription( "QDNAseq" )$Version
rVersion <- R.version.string
rPath <- R.home()
catMsg( c("QDNAseq version ", qdnaseqVersion) )
catMsg( c( rVersion ) )

qdnaseqObject <- readRDS( rdsFilePath )
sampleNames <- sampleNames( qdnaseqObject )

## sanity checks
elements <- assayDataElementNames(qdnaseqObject)
if ( ! "calls" %in% elements ) stop( "No calls present in object, regioning with CGHregions only work on called data" ) 
if ( length(sampleNames) < 2 ) stop( "Object contains too few samples, regioning with CGHregions only works with at least two samples" ) 

## analysis
cgh <- makeCgh( qdnaseqObject, filter=TRUE )
regions <- CGHregions( cgh, averror=0.00001 )
outputData <- cbind( regions@featureData@data, regions@assayData$regions )

## output
write.table( outputData, outputFile, sep="\t", quote=FALSE, row.names=FALSE )

## tell galaxy all seems ok
q(status=0)
