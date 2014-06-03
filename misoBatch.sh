#!/bin/bash

usage()
{
cat << EOF
usage: $0 [options] annotationFile sampleFile

OPTIONS:
  -t tmpDir	Temporal directory, default to tmp/
  -o outDir	Output directory, default is ./
  -h      	Show this message

This script analyze a group of samples specfied in sampleFile using miso. Sample
File format is a csv with the following header

groupName,treatment,sampleName,filePath

EOF
}

# Default values
TMP_DIR=tmp
OUT_DIR=.
LOG_DIR=log 
DBFILE=misodb.sqlite

# Argument processing
while getopts "ho:t:" OPTION
do
  case $OPTION in
    t)  TMP_DIR=$OPTARG;;
    o)  OUT_DIR=$OPTARG;;
    h)  usage; exit 0 ;;
    # Unknown option. No need for an error, getopts informs
    # the user itself.
    \?) exit 1;;
  esac
done

# Positional arguments
shift $(( OPTIND - 1 ));
ANNOTATION_FILE=$1
SAMPLE_GRUOPS=$2


# validation
if [[ -z $ANNOTATION_FILE ]] || [[ -z $SAMPLE_GRUOPS ]] 
then
	 echo "Error: missing annotation or samples configuration file"	
     usage
     exit 1
fi

# Annotation processing
ANNOTATION_INDEX=$TMP_DIR/$ANNOTATION_FILE.index

if [ -e $ANNOTATION_INDEX ];  then 
	echo "Using current index"
else
	echo "Creating index for annotations"
	index_gff.py --index $ANNOTATION_FILE $ANNOTATION_INDEX
fi

# sample processing
# File format:
# groupName,treatment,replica,sampleName,filePath
CURRENT_GROUP=""
SAMPLES_COUNT=0

mkdir $OUT_DIR $OUT_DIR/genes-psi $OUT_DIR/comparisons $OUT_DIR/summary  $OUT_DIR/plots $OUT_DIR/logs $OUT_DIR/logs/plots &> /dev/null

while IFS=',' read groupName treatment sampleName filePath; do
	if [ $groupName = "groupName" ]; then 
	    # create sashami conf file
    	sed -e "s/MISO_PREFIX/$OUT_DIR\/genes-psi\//" conf/sashimi_plot_template.txt > $OUT_DIR/sashami_settings.txt
    	rm $OUT_DIR/annotationIndex 
    	ln -s ../$ANNOTATION_INDEX $OUT_DIR/annotationIndex
		continue; # skip header
	fi
	
	# generate psi values and summary for genes 
	if [ ! -e $OUT_DIR/genes-psi/$sampleName/ ]; then
    	echo "Generating PSI values for $sampleName"
		READ_LENGTH=`samtools view $filePath | awk '{print $10}' | head -1  | wc -c`
		READ_LENGTH=$(( $READ_LENGTH - 1 ))
		run_events_analysis.py --compute-genes-psi $ANNOTATION_INDEX $filePath --output-dir $OUT_DIR/genes-psi/$sampleName --read-len $READ_LENGTH --settings-filename=conf/miso.conf &> $OUT_DIR/logs/$sampleName-psi.log
		run_miso.py --summarize-samples $OUT_DIR/genes-psi/$sampleName $OUT_DIR/summary/$sampleName  --settings-filename=conf/miso.conf &> $OUT_DIR/logs/$sampleName-summary.log
	fi
	
	# format and filter summay for upload to database
	FILTERED_SAMPLE_FILE=$OUT_DIR/filtered/${sampleName}.tsv
	#if [ ! -e $FILTERED_SAMPLE_FILE ]; then	
	    echo "Filtering $sampleName"
	#	./filter_isoforms.py --psi=0.0 --filter $OUT_DIR/summary/$sampleName/summary/${sampleName}.miso_summary --output $FILTERED_SAMPLE_FILE
		./misodb import $DBFILE ${OUT_DIR,,} psi $FILTERED_SAMPLE_FILE ${sampleName}
	#fi

    # Process comparations
	if [ -n "$COMPARATION_SAMPLE" -a "$CURRENT_GROUP" == "$groupName" ] ; then

		echo "Comparation between $COMPARATION_SAMPLE and $sampleName"
		COMPARATION_NAME=${COMPARATION_SAMPLE}_vs_${sampleName}
		if [ ! -e $OUT_DIR/comparisons/$COMPARATION_NAME ]; then 
			run_miso.py  --settings-filename=conf/miso.conf --compare-samples $OUT_DIR/genes-psi/$COMPARATION_SAMPLE/ $OUT_DIR/genes-psi/$sampleName/ $OUT_DIR/comparisons/ &> $OUT_DIR/logs/$COMPARATION_NAME.log
		fi	
		
		# Filter comparations and convert format
		FILTERED_FILE=$OUT_DIR/filtered/${COMPARATION_NAME}.tsv
#		if [ ! -e $FILTERED_FILE ]; then 
			echo "Filtering differencial events or isoforms"
			#./filter_isoforms.py --filter $OUT_DIR/comparisons/${COMPARATION_NAME}/bayes-factors/${COMPARATION_NAME}.miso_bf --delta-psi 0.0 --bayes-factor 0 --output ${FILTERED_FILE}			
			./misodb import $DBFILE ${OUT_DIR,,} diff $FILTERED_FILE ${COMPARATION_SAMPLE} ${sampleName}
#		fi	
        
        # Plot events
		ARGS_FILE=$TMP_DIR/argsFile.tmp    
		#echo "Plotting events"
		#echo "" > $ARGS_FILE 
		#while IFS=$'\t' read event_name sample1_posterior_mean sample1_ci_low sample1_ci_high sample2_posterior_mean sample2_ci_low sample2_ci_high diff bayes_factor isoforms sample1_counts sample1_assigned_counts sample2_counts sample2_assigned_counts chrom strand mRNA_starts mRNA_ends; do
		#	echo " --plot-event $event_name  $ANNOTATION_INDEX conf/sashimi_plot_settings.txt --output-dir $OUT_DIR/plots/ &> $OUT_DIR/logs/plots/$COMPARATION_NAME-$event_name.log " >> $ARGS_FILE
		#done < $FILTERED_FILE
		#xargs -n 1 --max-procs=10 --arg-file=$ARGS_FILE ./plot.py
		#rm $ARGS_FILE
	
		COMPARATION_SAMPLE=""
	else
		COMPARATION_SAMPLE="$sampleName"
	fi
	CURRENT_GROUP="$groupName"

done < $SAMPLE_GRUOPS
