#!/bin/bash
Usage(){
	echo "Usage : $0 -i sample.pep.fa -o . -n sampe.pep -t <thread_count>
from pep.fa to treefile "
        echo "-i   Input file"
        echo "-o   directory"
        echo "-n   Output filename pre"
        echo "-t   <thread_count>"
        exit 1
}

while getopts "hi:o:n:t:" opt
do
	case $opt in
		"h")
			Usage
			exit
			;;
		"i")
			input="$OPTARG"
			;;
		"o")
			outdir="$OPTARG"
			;;
		"n")
			outname="$OPTARG"
			;;
		"t")
			nt="$OPTARG"
			;;
	esac 
done


muscle -in ${input} -out ${outdir}/${outname}.mfa
iqtree2 -s ${outdir}/${outname}.mfa -m MFP -bb 1000 -nt ${nt} -pre ${outdir}/${outname}.iqtree
ls ${outdir} |grep ${outname} |grep -v "treefile" | xargs rm  
