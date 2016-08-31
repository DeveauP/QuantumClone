#!/bin/bash
PROGNAME=$(basename $0)
########################
#
# Basic stats
#
########################

#number of clones in the samples
# def refers to default value, i.e. value used when the parameter is not tested
# range refers to the different values that will be used 

#### Defaults
number_of_tests=50

def_clones=4
def_mut=100
def_ploidy="AB"
def_samp=2
def_dpth=100
def_ctm=0
### Ranges

range_clones="2 4 6 7 8 9 10"  #"2 4 6 8 10"
range_mutations="50 75 100 150 200"
range_ploidy="AB AAB AABB 2"
range_samples="1 2 3 4 5"
range_depth="50 100 200 500 1000"
range_contamination="0 0.1 0.2 0.3 0.5 0.7"

purity=1

### For debugging:
### see range clones after ~ line 150
#number_of_tests=1
#def_mut=30
##### end debugging
########################

function usage
{
	###
	#	Function to display usage message
	#	No arguments
	###
	echo "Usage: ${PROGNAME} [-h | --help] | [-p | --parameter] "
	echo "Commandline should look like this:"
	echo "${PROGNAME} -p=contamination"
}

function helptext
{
	###
	#	Function to display help message
	#	No argument
	###

	cat <<- _EOF_
	
	${PROGNAME}
	This program benchmark QuantumClone, sciClone, and pyClone, testing a single parameter
	
	$(usage)
	
	Options:
	
	-h, --help	Display this help message
	-p, --parameter	The parameter to test, one of: clones, mutations, samples, depth, contamination, CNA
			
	_EOF_

}


################################################################################################################
#
#	Program starts here
#
################################################################################################################

#trap term_exit TERM HUP
#trap int_exit INT

################################################################################################################
# Command line arguments: checking input is correct							       #
################################################################################################################

p_flag=0

while [ "$1" != "" ]; do

	PARAM=`echo $1 | awk -F= '{print $1}'`
	VALUE=`echo $1 | awk -F= '{print $2}'`
	
	case $PARAM in
		-p | --parameter)
			param=$VALUE
			p_flag=1
			;;
		-h | --help)	
			helptext
			exit
			;;
		*)	usage
			exit 1
	esac
	shift
done

ex_flag=0
if [ $p_flag = 0 ]; then
	err_msg="\t\t${PROGNAME}: parameter not specified"
	echo -e ${err_msg} >&2
	ex_flag=1
fi

if [[ 'clones mutations samples depth contamination CNA' != *$param* ]]; then
 	err_msg="\t\tParameter must be one of: clones, mutations, samples, depth, contamination, CNA.\n\t\tCheck provided parameter."
	echo -e ${err_msg} >&2
	ex_flag=1

fi

echo $param

if [ $ex_flag = 1 ]; then
	exit
fi

input=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source ${input}/config


echo $input
################################################################################################################
# CREATING DIRECTORIES AND R ANALYSES									       #
################################################################################################################

cd $output
if [ $param = clones ];then
	idx=$range_clones
	for p in $range_clones;do
		mkdir -p $output/$param/$p
		cd $output/$param/$p
		$R < $input/Rscript.R $p $def_mut $def_ploidy $def_samp $def_dpth $def_ctm $number_of_tests --no-save
	done
elif [ $param = mutations ];then

	idx=$range_mutations
	for p in $range_mutations;do

		mkdir -p $output/$param/$p
		cd $output/$param/$p
		$R < $input/Rscript.R $def_clones $p $def_ploidy $def_samp $def_dpth $def_ctm $number_of_tests --no-save
	done
elif [ $param = samples ];then
	idx=$range_samples
	for p in $range_samples;do
		mkdir -p $output/$param/$p
		cd $output/$param/$p
		$R < $input/Rscript.R $def_clones $def_mut $def_ploidy $p $def_dpth $def_ctm $number_of_tests --no-save
	done
elif [ $param = depth ];then
	idx=$range_depth
	for p in $range_depth;do
		mkdir -p $output/$param/$p
		cd $output/$param/$p
		$R < $input/Rscript.R $def_clones $def_mut $def_ploidy $def_samp $p $def_ctm $number_of_tests --no-save
	done
elif [ $param = contamination ];then
	idx=$range_contamination
	for p in $range_contamination ; do
		mkdir -p $output/$param/$p
		cd $output/$param/$p
		$R < $input/Rscript.R $def_clones $def_mut $def_ploidy $def_samp $def_dpth $p $number_of_tests --no-save
	done
elif [ $param = CNA ];then
	idx=$range_ploidy
	for p in $idx ; do
		mkdir -p $output/$param/$p
		cd $output/$param/$p
		$R < $input/Rscript.R $def_clones $def_mut $p $def_samp $def_dpth $def_ctm $number_of_tests --no-save
	done
fi


#SUBDIRS=$(ls -d $output/$param/)

echo $SUBDIRS

### sourcing python!
#export PYTHONPATH=$pythonpath:$pyclone

################################################################################################################
# pyClone ANALYSES											       #
################################################################################################################
################################################################################################################
# CHECKED TILL HERE											       #
################################################################################################################

if [ $param != "contamination" ]; then
	SUBDIRS=$output/$param/
	for p in $idx;do
		> ${SUBDIRS}/$p/time_PyClone
		if [ "$param" == "samples" ]; then 
			seqp=$(seq 1 $p)
		else
			seqp=$(seq 1 $def_samp)
		fi
		for test in $(seq 1 $number_of_tests);do
			SUB=${SUBDIRS}/$p/test$test/
			mkdir -p  ${SUB}/YAML
			STARTTIME=$(date +%s)
			sed "s|WORKING_DIR|$SUB|g" $pyclone_config > ${SUB}/config.yaml
			for sp in $seqp;do
				echo "sp = $sp"
				echo "  SAMPLE_$sp:" >> ${SUB}/config.yaml
				echo "    mutations_file: $SUB/YAML/sample${sp}.yaml" >> ${SUB}/config.yaml
				echo "" >> ${SUB}/config.yaml
				echo "    tumour_content:" >> ${SUB}/config.yaml
				echo "      value: $purity" >> ${SUB}/config.yaml
				echo "" >> ${SUB}/config.yaml
				echo "    error_rate: 0.001" >> ${SUB}/config.yaml
				$pyclone build_mutations_file ${SUB}/sample${sp}.tsv ${SUB}/YAML/sample${sp}.yaml 
			done
			$pyclone run_analysis --config_file ${SUB}/config.yaml #--in_file ${SUB}cluster.txt
			$pyclone build_table --config_file ${SUB}/config.yaml --table_type old_style --out_file ${SUB}cluster.txt
			ENDTIME=$(date +%s)
			ELAPSED=$(($ENDTIME-$STARTTIME))
			echo "$ELAPSED" >> ${SUBDIRS}/$p/time_PyClone
			#R < $input/post_process.R ${output}/$param/$p/ $param $p --no-save
			#rm ${output}/$param/$p/Rplots.pdf	
		done
	done
	#$R < $input/post_process.R $SUBDIRS --no-save
else ## conta
	SUBDIRS=$output/$param/

	for p in $idx;do
		echo $p
		> ${SUBDIRS}/$p/time_PyClone
		seqp=$(seq 1 $def_samp)
		purity=$(echo 1.0 - $p | bc)
		echo $purity

		for test in $(seq 1 $number_of_tests);do

			SUB=${SUBDIRS}/$p/test$test/

			mkdir -p  ${SUB}/YAML
			STARTTIME=$(date +%s)
			sed "s|WORKING_DIR|$SUB|g" $pyclone_config > ${SUB}/config.yaml
			for sp in $seqp;do
				echo "sp = $sp"
				echo "  SAMPLE_$sp:" >> ${SUB}/config.yaml
				echo "    mutations_file: $SUB/YAML/sample${sp}.yaml" >> ${SUB}/config.yaml
				echo "" >> ${SUB}/config.yaml
				echo "    tumour_content:" >> ${SUB}/config.yaml
				echo "      value: $purity" >> ${SUB}/config.yaml
				echo "" >> ${SUB}/config.yaml
				echo "    error_rate: 0.001" >> ${SUB}/config.yaml
				$pyclone build_mutations_file --in_file ${SUB}/sample${sp}.tsv --out_file ${SUB}/YAML/sample${sp}.yaml 
			done

			$pyclone run_analysis --config_file ${SUB}/config.yaml > ${SUB}/analysis.log #--in_file ${SUB}cluster.txt
			$pyclone build_table --config_file ${SUB}/config.yaml --table_type old_style --out_file ${SUB}cluster.txt > ${SUB}/table.log
			ENDTIME=$(date +%s)
			ELAPSED=$(($ENDTIME-$STARTTIME))
			echo "$ELAPSED" >> ${SUBDIRS}/$p/time_PyClone	
			#R < $input/post_process.R ${output}/$param/test$test/ $test $param --no-save
			#rm ${output}/$param/test$test/Rplots.pdf
		done
	done
	#$R < $input/post_process.R $SUBDIRS --no-save
fi
