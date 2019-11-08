#!/bin/bash

############################################################################
# This script assumes that HISAT2, SAMtools, Stringtie are installed       #
# and that this script is downloaded and chmod executable #		   #
# Furthermore, this script uses ballgown as ultimate analysis so initial   #
# "stringtie" call will be without the "-e", so opt to not use "-e" in     #
# initial "stringtie" assembly but it will be used in abundance estimate   #
# It also assumes that a reference gene files will be used for 		   #
# mapping and assemblying. And that sample names follow the convention of  #
# sample names in JHU HISAT2 example names. If not then this               # 
# script should be modified to reflect those choices.                      # 
############################################################################
# This script starts with mapping via HISAT2,
# Sorting based on genomic position and conversion by  SAMtools, 
# Transcriptome assembly by Stringtie, and merging for downstream analysis
# Generation of gene and transcript count matrices by PrepDEanalysis.py

#### Additionally, and this is CRITICAL as some scripts depend on this structure:
# It assumes that the following datasets and directories are already created
# This is the reference annotation genome for assembly:       "chrX_data/genes/chrX.gtf"
# The genome index files are located at and have prefix :     "chrX_data/indexes/chrX_tran"
# The  samples are located at :        		              "chrX_data/samples"

# Testing git update remotely

################################# Optional ################################
##### Can download the chrX_data folder from JHU and untar it
##### At least the indexes and reference gene annotation can be downloaded
##### The samples subdirectory is the reason for such 2.2GB file

#wget -c -r ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
#2019-09-23 16:13:44 (40.5 MB/s) - ‘ftp.ccb.jhu.edu/pub/RNAseq_protocol/.listing’ saved [733]
#FINISHED --2019-09-23 17:50:23--
#Total wall clock time: 1h 36m 41s
#Downloaded: 1 files, 2.0G in 1h 36m 38s (356 KB/s)



# original expression was: 
#read -p "Are all folders index, genes, samples in proper place. Continue (y/n)? " -n 1 -r
# The (-n 1 -r) argument to (read -p) command is to not depend on user hitting enter

read -p "Question 1: Are all the indexes, genes, and samples folders with their respective contents in their proper place? (y/n)?"

echo    # (optional) move to a new line


if [[ $REPLY =~ ^[Yy]$ ]]
then
echo "Great. Next question..."
#echo "Will proceed with script"

elif [[ $REPLY =~ ^[Nn]$ ]]
then
	echo "got to this line"
	read -p "Should I download the JHU sample index, ref. genes, and samples(it may take more than an hour)?(y/n)"
	
	if [[ $REPLY =~ ^[Yy]$ ]]
	then
		curl -L "ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz" > chrX_data.tar.gz
		# The download took almost 2 hours
		# % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
		#                                 Dload  Upload   Total   Spent    Left  Speed
		#100 2015M  100 2015M    0     0   339k      0  1:41:16  1:41:16 --:--:--  404k


		tar -xzvf chrX_data.tar.gz

	elif [[ $REPLY =~ ^[Nn]$ ]]
	then  echo "Then it will not be downloaded"

	fi


#exit 1
else
echo "Will proceed with script assuming things are in order"
fi








#################################### Some checks before starting script #######


###### To generate a non-redundant sample names txt file for convenience if user didn't supply it!!!
###### Thus user can supply the text file name as argument during script call
# Use an if statement to determine if users has supplied a sample_list file
# test if user supplied input argument
# https://stackoverflow.com/questions/6482377/check-existence-of-input-argument-in-a-bash-shell-script

if [ $# -eq 0 ]
  then
        ls -1 chrX_data/samples/ | cat > samplenames.txt

        # this assumes that suffix of filenames contain "_chrX..."
        sed -n -i -e 's/_chrX.*$//p' samplenames.txt #> samplenames.txt

        # inelegant but works. can't do it in place or use same name.
        sort -u samplenames.txt > samplenamestemp.txt
        mv samplenamestemp.txt samplenames.txt

        files=samplenames.txt

	echo "You didn't supply a text file of sample names so script will scan samples folder and construct this text file for you. It will be called $files."
	#echo "$# args supplied. No text file of sample names so will create based on \n filenames supplied in <chrX_data/samples/*> folder"
	echo "This script also assumes that sample names are in format <basename_chrX_n.fastq.gz> where n=null or 1 and 2 for paired end reads"

else

#	files=("$@") # This also works
	files=$1 # $0 is command/script name
        echo "You provided $# text file(s) of sample names called $files !!!"

fi

echo # blank lin

# The following didn't work
#sed -n -e 's/_chrX.*$//p' samplenames.txt
#sort -u samplenames.txt > samplenames.txt

#read -p "Were able to check file input and user input without conflict. Should we stop the script here?"
#if [[ $REPLY =~ ^[Yy]$ ]]
#then exit 1
#elif [[ $REPLY =~ ^[Nn]$ ]]
#then #exit 1
#echo "will proceed"
#fi



############ Create directories for output ############
# make map directory for HISAT2 if it doesn't exist
[ -d map ] || mkdir map
#  also make assembly directory for output of Stringtie if it doesn't exist
[ -d assembly ] || mkdir assembly





## A conditional to see if user has made appropriate folders with samples etc.
# https://stackoverflow.com/questions/1885525/how-do-i-prompt-a-user-for-confirmation-in-bash-script

#read -p "Are all folders index, genes, samples in proper place. Continue (y/n)? " -n 1 -r

echo "There are 3 QUESTIONS to answer before script executes. TYPE response and hit ENTER"

echo # blank line

sleep 1



#### conditional to see if indexes, genes, and samples directories exist and have content

#if [ -z "$(ls -A chrX_data/indexes/)" ]; then echo "directory exists but is empty"; else echo "it exists and is not empty"; fi

if [ -z "$(ls -A chrX_data/indexes/)" ]; 
then echo "chrX_data/indexes/ directory doesn't exist or is empty. Please create and populate it with chromosome index files" 
elif [ -z "$(ls -A chrX_data/genes/)" ]; 
then echo "chrX_data/genes/ directory doesn't exist or is empty. Please create and populate it with reference genes files"
elif [ -z "$(ls -A chrX_data/samples/)" ]; 
then echo "chrX_data/samples/ directory doesn't exist or is empty. Please create and populate it with your sample files"
else
   echo #"chromosome indexes, reference genes, and samples directories exist and are not Empty" # not useful message
fi

# original expression was: 
#read -p "Are all folders index, genes, samples in proper place. Continue (y/n)? " -n 1 -r
# The (-n 1 -r) argument to (read -p) command is to not depend on user hitting enter

read -p "Question 1: Last check; are all the indexes, genes, and samples folders with their respective contents in their proper place? Continue (y/n)? " -t 10

echo    # (optional) move to a new line


if [[ $REPLY =~ ^[Yy]$ ]]
then
echo "Great. Next question..."
#echo "Will proceed with script"

elif [[ $REPLY =~ ^[Nn]$ ]]
then exit 1
else
echo "Will proceed with script assuming things are in order"
fi


	#http://www.daniloaz.com/en/differences-between-physical-cpu-vs-logical-cpu-vs-core-vs-thread-vs-socket/
        # To get logicals, which includes systems with hyperthreading (not exactly equivalent to a core processing power)
if [[ $(uname) == "Linux" ]]; 
then echo "You have a linux system" && cpu_cores=$(lscpu -p | egrep -v '^#' | wc -l) || echo "NO Problem; Will use an older command to get cpu cores" && cpu_cores=$(getconf _NPROCESSORS_ONLN) ;
elif [[ $(uname) == "Darwin" ]]; 
then echo "You have a macOS" && cpu_cores=$(sysctl -n hw.logicalcpu_max) || echo "NO Problem; Will use an older command to get cpu cores" && cpu_cores=$(getconf _NPROCESSORS_ONLN) ;
else  echo "Not macOS or Linux system so this script will not work";
fi

read -p "Question 2: How many CPU Cores do you want to use? Your system has $cpu_cores cores. Give integer value(from 1 to ...):"  user_cores
#echo # (optional) move to a new line
if [[ -z "$user_cores" || ! "$user_cores" =~ ^[0-9]+$ ]]
then
echo "You didn't supply the number of cores so script will use default value of 1"
user_cores=1
fi

echo "The script will use $user_cores cores"
echo # blank line


read -p "Question 3: Stringtie with (-e) option removes potentially novel transcripts and false positive transcripts. Useful for DEG analysis ex. DESEQ2. However, for ballgown it is typically omitted. Should Stringtie run with (-e) option? (y/n)"
echo    # (optional) move to a new line



# The essence of this argument is to generate a value or leave it empty and 
# the logic in "stringtie" argument will determine if e is empty or contains string etc.
if [[ $REPLY =~ ^[Yy]$ ]]
then
# The essence of this argument is to generate a value or leave it empty and the logic in "stringtie" determine if e is empty or contains string etc.
e="yes"

echo "Great! Will proceed with script using stringtie WITH the (-e) ...."
 # 
else e=""
echo "Great! Will proceed with script using stringtie WITHOUT the (-e) ...."
fi



sleep 1

#echo "Great! Will proceed with script ...."
echo # blank line
echo # blank 









#############################################################################################################

				# HISAT2 #

# $ hisat2 --version
#Projects/RNAseqPreprocess/HISAT2/hisat2-2.1.0/hisat2-align-s version 2.1.0
# 64-bit
# Built on login-node03
# Wed Jun  7 15:53:42 EDT 2017
# Compiler: gcc version 4.8.2 (GCC) 
# Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY
# Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}




# Here is the alignment/mapping command
#$ hisat2 -p 4 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188428_chrX_1.fastq.gz -2 chrX_data/samples/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam
# https://ccb.jhu.edu/software/hisat2/manual.shtml


# This script, for simplicity's sake, assumes that the chromosome indexes and samples are in 
# their respective folders("indexes" and "samples") within the "chrX_data" parent directory 
# and that a "map" subdirectory exists in this current directory from which the script is run
# also in this case this version of script will request the a text file containing the sample 
# names(one sample name per line) in current directory where this script lives.

# make map directory for HISAT2 if it doesn't exist
#[ -d map ] || mkdir map
#  also make assembly directory for output of Stringtie if it doesn't exist
#[ -d assembly ] || mkdir assembly

# For case when User provides sample names

# files=("$@")
# for samplenames in $(cat $files);


# For case when script above generates the sample names based on 
# file name structures in "chrX_data/samples/" folder

for samplenames in $(cat $files)
do 
	# This creates an array
	var1=($(grep -l "$samplenames" chrX_data/samples/* ))
	#var1=($(grep -HRl "$samplenames" chrX_data/samples/* | cat ))

	##var1=($(grep -lf "$samplenames" chrX_data/samples/* )) # not work
	##var1=($(grep -ef "$samplenames" chrX_data/samples/* | cat )) # not work
	#grep -HRl
	##var1=($(echo chrX_data/samples/* | grep -l "$samplenames" )) # not work

        # This stores the matched files into variable not just their name
	#var1=$(grep -l "$samplenames" chrX_data/samples/* )  # this creates an array of length "1" although it prints two elements
	var2=$(grep -l "$samplenames" chrX_data/samples/* | cat | wc -l)
	#var2=$(grep -c "$samplenames" chrX_data/samples/*)

	if [[ $var2 -eq 2 ]]; 
	then 
#		echo "$var1"
#		declare -a FL
#		for f in $var1; #$(cat $var1);
#		do 	FL[$#FL[@]}+1]=$(echo "$f");
#			echo ${FL[1]}
#		done
#	hisat2 -p 4 --dta -x chrX_data/indexes/chrX_tran -1 ${var1[0]} -2 ${var1[1]} -S map/$samplenames\_chrX.sam
	hisat2 -p $user_cores --dta -x chrX_data/indexes/chrX_tran -1 ${var1[0]} -2 ${var1[1]} -S map/$samplenames\_chrX.sam
	#This prints out length of array
#        echo "${var1[@]}" "${#var1[@]}" "This is paired end"
	echo "This is first sample of pair ${var1[0]}" "This is second sample of pair ${var1[1]}" "the length of array is ${#var1[@]}"
	echo "This is paired end";
#	echo $(wc -l "$var1[2]") #${#var1};
	# else # initially I used an else statement but elif is better
	# technically it should be either single end($var2=1) or paired end ($var1=2), but to catch all possibilities
	elif [[ $var -eq 1 ]]; 
	then 
	hisat2 -p $user_cores --dta -x chrX_data/indexes/chrX_tran -U ${var1[0]} -S map/$samplenames\_chrX.sam
	echo "This is not paired end and sample is ${var1[0]}";
	else echo "Neither single end or paired samples found"
	fi
	

#	declare -a FILELIST
#for f in *; do 
    #FILELIST[length_of_FILELIST + 1]=filename
#    FILELIST[${#FILELIST[@]}+1]=$(echo "$f");
#done
	#[ -z "$var1" ] && echo "Empty" || echo "Not empty"
	#var2=$(cat $var1 | wc -l ) 

#	echo ${var1[1]}

#	[[ $var1 =~ OK$ ]]

	
#	printf "%s\n"  "This is $var1"
#	var2=$(wc -l "$var1")
#	echo "$var1" | wc -l
#	printf "%s\n" "$var2"
	echo "$var2"
#	echo $(wc -l $var1)
#	var2=$(wc -l "$var1")
#	echo $var2
#	if [[ $var2 -eq 2 ]]; 
#        then echo "This is paired end";

	#then echo $var1 | cut -d" " -f2
	#then set -- $var1 
        #echo "$1"
        #echo "$2"; 

#        else echo "NOT paired end"; 
#        fi




	#echo $(wc -l $var1)

done


################################################################################################# 

# $ samtools --version
# samtools 1.9
# Using htslib 1.9
# Copyright (C) 2018 Genome Research Ltd.

# This part of script will automate the run of the samtools sorting and 
# conversion of ".sam" files into ".bam" in the "map" directory

# general syntax of samtools the "-@ 4" is the thread number
# $ samtools sort -@ 4 -o map/ERR188428_chrX.bam map/ERR188428_chrX.sam


# This script is for samtools analysis automation. Basic syntax is
# $ samtools sort -@ 4 -o map/ERR188428_chrX.bam map/ERR188428_chrX.sam

# This works in listing files in directory
#$ ls map/*.sam | cat

for samfiles in  $( ls map/*.sam | cat );
do
	# This outputs ex. "map/ERR188044_chrX.sam"
        echo "$samfiles"
	# This is a temp file name for samtools command. sed with "-i"  will update the file
	new="$( echo "$samfiles" | sed 's/sam/bam/g' )"
        echo "$new"

#	new=$( sed 's/sam/bam $samfiles') # needed a g
	samtools sort -@ $user_cores -o $new $samfiles 
done


############ remove the sam files ##############

rm map/*.sam



#################################################################################################

# $ stringtie --version
# 2.0.1



# This script to automate the Stringtie analysis, especially for many sample case
# The general syntax is
# stringtie map/ERR188044_chrX.bam -l ERR188044 -p 4 -e -G chrX_data/genes/chrX.gtf -o assembly/ERR188044_chrX.gtf


#for files in $( ls -1 testdir2/*.bam | tr '\n' '\0' | xargs -0 -n 1 basename | cat );
# This if statement structure is to determine if user specified to run Stringtie with the use of "-e" option
if [ -z "$e" ];
        then    echo  " Then therefore no (-e) option is used"
	for bamfiles in $( ls -1 map/*.bam | cat );

	do
        prefxfile=$( echo $bamfiles | tr '\n' '\0' | xargs -0 -n 1 basename )
#       the " tr '\n' '\0' " is not required here. It converts 
#       newlines to nulls so "xargs" can work with whitespaces in names
#       prefix="$( echo "$bamfiles" | sed 's/txt/bam/g' )"
        echo "$bamfiles"
        # The prefix argument assumes that the filename ends in "_chrX.bam" and replaces it with empty string
        prefix=$(echo $prefxfile | sed 's/_chrX.bam//g')
        outputname=$( echo $prefix\_chrX.gtf)

	stringtie $bamfiles -l $prefix -p $user_cores -G chrX_data/genes/chrX.gtf  -o assembly/$outputname
	done

else    echo "$samplenames the stringtie option is $e"

	for bamfiles in $( ls -1 map/*.bam | cat );

        do
        	prefxfile=$( echo $bamfiles | tr '\n' '\0' | xargs -0 -n 1 basename )
#       	the " tr '\n' '\0' " is not required here. It converts 
#       	newlines to nulls so "xargs" can work with whitespaces in names
#       	prefix="$( echo "$bamfiles" | sed 's/txt/bam/g' )"
        	echo "$bamfiles"
        	# The prefix argument assumes that the filename ends in "_chrX.bam" and replaces it with empty string
        	prefix=$(echo $prefxfile | sed 's/_chrX.bam//g')
        	outputname=$( echo $prefix\_chrX.gtf)

                stringtie $bamfiles -l $prefix -p $user_cores -e -G chrX_data/genes/chrX.gtf  -o assembly/$outputname
	done

fi





#for files in $( ls -1 map/*.bam | cat );

#do
#        file=$( echo $files | tr '\n' '\0' | xargs -0 -n 1 basename )
#       the " tr '\n' '\0' " is not required here. It converts 
#	newlines to nulls so "xargs" can work with whitespaces in names
#       prefix="$( echo "$files" | sed 's/txt/bam/g' )"
#        echo "$files"
#	# The prefix argument assumes that the filename ends in "_chrX.bam" and replaces it with empty string
#        prefix=$(echo $file | sed 's/_chrX.bam//g')
#        outputname=$( echo $prefix\_chrX.gtf)

#	if [ -z "$e" ];
#        then    echo  " Then therefore no "-e" option is used"

#		stringtie $files -l $prefix -p $user_cores -G chrX_data/genes/chrX.gtf  -o assembly/$outputname

#        else    echo "$samplenames the stringtie option is $e"
#		stringtie $files -l $prefix -p $user_cores -e -G chrX_data/genes/chrX.gtf  -o assembly/$outputname
#        fi

#done


############################### In this version of Script #################################

# I will add the generation of a merglist which will be given to 
# Stringtie command, which will output a file(...generate a non-redundant set of 
# transcripts observed in all the RNA-Seq samples assembled previously) for 
# PrepDEanalysis.py script to generate counts for DESeq2 analysis.
# Also consider saving the stdout output of "hisat2..." command, which 
# displays the quality report of alignment process. Typically it is not saved.


########## First step is to generate a text file containing the "assembly/samples_n" name
# Initially, I automatically generated the mergelist.txt
#ls -1 assembly/* | cat > mergelist.txt
# but 
# I decided to check if user already has mergelist.txt in directory. 
# This assumes that the file if it exists is in the current directory


mergefilel=$(ls | grep -e "mergelist.txt")
if [ -z "$mergefilel"  ];
        then ls -1 assembly/* | cat > mergelist.txt
else echo "A mergelist.txt file has been found in current directory. So it will not be automatically generated"
fi


########### To generate a specific subset of samples for mergelist ########
#  Then manually generate the text file  in the format

#	assembly/samplex_chrX.gtf

# merge all transcripts from the different samples
# stringtie --merge -p 4 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt

# Using stringtie to merge transcripts
stringtie --merge -p $user_cores -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf ./mergelist.txt


################################ PrepDEanalysis.py ########################################

	###############  Generating the text file for PrepDEanalysis.py ##################



for samplepath in $( less mergelist.txt | cat );
do
        file=$(echo "$samplepath" | tr '\n' '\0' | xargs -0 -n 1 basename) 
        # the " tr '\n' '\0' " is not required. It converts newlines to nulls so "xargs" can work with whitespaces in names
        prefix=$(echo "$file" | sed 's/_chrX.gtf//g')

        echo "$prefix" "$samplepath" >> sample_lst.txt
done



#       # This will append the prefix to original line text and concatenate to a file. 
#       # Using "sed" will add the same prefix to each line. For simplicity no need for conditional for match
######### echo "$prefix" "$files" >> sample_lst.txt # a single ">" will overwrite ">>" will append
######### Apparently the appending at "done" is more efficient than at "echo..." command
#       # https://stackoverflow.com/questions/9980425/loop-and-append-write-to-the-same-file-without-overwriting

#done > sample_lst.txt



# To get to this point it took
# RNAseqAlignment_tut4$ time ./hisat_sam_stringtie_merge_prepDEanalysispy_v3.sh
# real	7m8.520s
# user	15m59.331s
# sys	1m20.399s
# The second run of script gave
# real	7m42.852s
# user	16m4.658s
# sys	1m20.640s



##################################### Generate files using Stringtie for ballgown analysis #####################

############# This step will estimate transcript abundance for ballgown analysis ###################
# Basic syntax
#stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf map/ERR188044_chrX.bam

# This file created in the beginning has the sample names without directories or extensions
#samplenames.txt

# Although the "-e" option is omitted in transcript assembly step above, here it will be used


[ -d ballgown ] || mkdir ballgown


#if [ -z "$e" ];

#then
#	echo "No (-e) option used"
#	for samplenames in $( less samplenames.txt | cat );

#	do
	#       cpu_cores=$(lscpu -p | egrep -v '^#' | wc -l)
	#       mappath=$(less map/* | egrep -v '$samplenames' | cat)

#        	mappath=$(ls -1 map/* | grep -e '$samplenames')
#        	stringtie -B -p $user_cores -G stringtie_merged.gtf -o ballgown/$samplenames/$samplenames\_chrX.gtf $mappath
#else 


echo "In this step the (-e) option is used for downstream ballgown analysis in R environment"
for samplenames in $( less samplenames.txt | cat );

do
        #       cpu_cores=$(lscpu -p | egrep -v '^#' | wc -l)
        #       mappath=$(less map/* | egrep -v '$samplenames' | cat)
#	mappath=$(ls -1 map/* | grep -e '$samplenames') # quoting the $samplenames was incorrect
	mappath=$(ls -1 map/* | grep -e $samplenames) 
	stringtie -e -B -p $user_cores -G stringtie_merged.gtf -o ballgown/$samplenames/$samplenames\_chrX.gtf $mappath

done

#fi

#stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf map/ERR188044_chrX.bam





	#################### Actual run of script ##################

#$ python2 ../PrepForDESEQ2/prepDEanalysis.py -i sample_lst.txt

# python2 ./prepDEanalysis.py -i sample_lst.txt

#$ time python2 ./prepDEanalysis.py -i sample_lst.tx
# This part took
# real	0m0.515s
# user	0m0.459s
# sys	0m0.044s

