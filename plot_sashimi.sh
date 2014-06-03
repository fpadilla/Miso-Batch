#!/bin/bash

for event in `cat $1`; do

	./plot.py --plot-event ${event} ./tmp/Homo_sapiens.GRCh37.65_fixed.gff.index/ conf/sashimi_plot_settings.txt --output-dir plotBF20/

done
