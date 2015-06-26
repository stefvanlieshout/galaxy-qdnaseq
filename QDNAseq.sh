#!/bin/sh

## This script creates a config file for use with the QDNAseq R wrapper
## For use outside Galaxy

QDNASEQ_R_SCRIPT='/ccagc/home/stef/code/QDNASEQ/QDNAseq.R'
R_LOCATION='/ccagc/lib/R/R-3.1.0/bin/'
OUTPUT_DIR='QDNASEQ'
SCRIPT_DIR=$(dirname $0)
SCRIPT_NAME=$(basename $0)

binSizesString='1000,100,15'
outputName=`basename $PWD`
doSegment='TRUE'
debug='FALSE'

usage()
{
cat<<EOF
  This script creates a QDNAseq wrapper config.
  
  Experiment name and/or output dir are optional
    $SCRIPT_NAME [-n myExperiment] [-o outputDir] *.bam

  Debug mode for working with build in LGG150 data
    $SCRIPT_NAME -d

  OPTIONS:
    -n [s]  provide analysis name [$outputName]
    -o [s]  output dir [$OUTPUT_DIR]
    -d      debug mode (only create config)
    
EOF
}

PARAM_COUNT=0
if [ $# -eq 0 ] || [ "$1" == "-h" ]
then
	usage
	exit 1
else
	while getopts "n:o:d" opt; do
	  case $opt in
	    n)
	      #echo "$opt was triggered, Parameter: $OPTARG" >&2
	      outputName=$OPTARG
	      PARAM_COUNT=`bc <<< $PARAM_COUNT+2`
	      ;;
	    o)
	      OUTPUT_DIR=$OPTARG
	      PARAM_COUNT=`bc <<< $PARAM_COUNT+2`
	      ;;
	    d)
	      debug='TRUE'
	      PARAM_COUNT=`bc <<< $PARAM_COUNT+1`
	      ;;
	    \?)
	      echo "Invalid option: -$OPTARG" >&2
	      exit 1
	      ;;
	    :)
	      echo "Option -$OPTARG requires an argument." >&2
	      exit 1
	      ;;
	  esac
	done
fi

## set param dependent variables
printf '%s\n' "[INFO] --- QDNAseq wrapper config creator ---"
CONFIG_FILE=$OUTPUT_DIR'/qdnaseqConfig.R'

## remove the params from input string leaving only the BAMs
shift $PARAM_COUNT

## sanity checks
if [ -e $CONFIG_FILE ]; then
	echo "[ERR] Config file ($CONFIG_FILE) already exists, exit.."; 
	exit 0
fi
if [ -d $OUTPUT_DIR ]; then
	echo "[ERR] Output dir ($OUTPUT_DIR) already exists, exit.."; 
	exit 0
else
	mkdir $OUTPUT_DIR
fi

## variables setup
CONFIG_TXT="## ==========
## QDNAseq pipeline config file
## ==========

## ----------
inGalaxy <- FALSE
binSizesString <- '$binSizesString'
experimentType <- 'SR50'
outputName <- '$outputName'

## ----------
outputHtml  <- 'galaxyIndex.html'
outputId    <- NA
newFilePath <- '$OUTPUT_DIR'
outputPath  <- '$OUTPUT_DIR'
doSegment   <- as.logical( $doSegment )
debug       <- as.logical( $debug )
undoSD      <- as.double( 1.0 )
binAnnotations <- ''

## ----------
filterBlacklistedBins <- as.logical( 'TRUE' )
mappabilityCutoff <- as.integer( 0 )
undoSplits <- 'sdundo'
doOutputCopynumbersIgv <- FALSE

## ----------
PLOT_WIDTH <- as.integer( 1440 )
PLOT_HEIGHT <- as.integer( 720 )
excludeChrsString <- 'X,Y'

## ----------
bamsPaths <- c()
bamsNames <- c()
"

## create config file
printf '%s\n' "$CONFIG_TXT" > $CONFIG_FILE

## add bam files to config file
BAMS=$@
for bam_path in $BAMS; do
	
	if [ ! -e $bam_path  ]; then
		echo "[ERR] Bam file does not exist ($bam_path), exit.."; 
		exit 0
	fi

	bam_name=`basename "$bam_path"`
	printf '%s\n' "bamsPaths <- c( bamsPaths, \"$bam_path\" )" >> $CONFIG_FILE
	printf '%s\n' "bamsNames <- c( bamsNames, \"$bam_name\" )" >> $CONFIG_FILE
	printf '%s\n' "[INFO] BAM: $bam_name"
done

## done
printf '%s\n' "[INFO] OUTPUT DIR:  $OUTPUT_DIR"
printf '%s\n' "[INFO] OUTPUT NAME: $outputName"
printf '%s\n' "[INFO] CONFIG FILE can be found in $CONFIG_FILE"
printf '%s\n' "[INFO] Run with: $R_LOCATION/Rscript $QDNASEQ_R_SCRIPT $CONFIG_FILE"
printf '%s\n' "[INFO] Done"