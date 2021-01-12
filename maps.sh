#!/bin/bash
python_path=python #should have pysam, pybedtools installed. bedtools, samtools should be in the path
Rscript_path=Rscript
cwd="/hpc/home/yx157/apps/maps_bin/bin"  ## maps script folder !
###################################################################
# setnames=("H1")  
# dataset_name='H1' ## samplename !!! 
# macs2_filepath="./H1.bed"   ### atac peaks!!!
# outdir="h1_hg38_10k"   ### output folder !!!
# feather_output_symlink="./h1_reads"   ## short long range folder !!!
# organism="hg38" ## set the reference !!! # mm9 mm10 hg19 hg38
# bin_size=10000  ## bin size!!!
# genomic_feat_filepath="/datacommons/ydiaolab/genome_ref/maps_feature/hg38_cviqi/hg38/CviQI/F_GC_M_CviQI_10Kb_el.txt"  ## the backgound!!!

dataset_name=$1 ## samplename !!! 
macs2_filepath=$2   ### atac peaks!!!
outdir=$3   ### output folder !!!
feather_output_symlink=$4  ## short long range folder !!!
organism=$5 ## set the reference !!! # mm9 mm10 hg19 hg38
bin_size=$6  ## bin size!!!
genomic_feat_filepath=$7
##################################################################
# feather=0 #start from feather or run only MAPS
# maps=1
# number_of_datasets=1

#dataset_name="naive_R1"
# fastq_format=".fastq"
# fastq_dir="/home/jurici/MAPS/examples/test_set1"

binning_range=100000000
fdr=2 # this is used just for labeling. do not change
filter_file="None"
length_cutoff=1000
model="pospoisson" #"negbinom"
sex_chroms_to_process=""
####################################################################
###SET THE VARIABLES AT THIS PORTION ONLY IF 
### number_of_datasets > 1 (merging exisitng datasets)
### specify as many datasets as required
####################################################################
#...
##################################################################
###SET THESE VARIABLES ONLY IF FEATHER = 0 AND YOU WANT TO RUN
###USING A SPECIFIC FEATHER OUTPUT RATHER THAN $datasetname_Current
###################################################################



##################################################################

# DATE=`date '+%Y%m%d_%H%M%S'`
DATE=''
#####Armen:
# fastq1=$fastq_dir/$dataset_name"_R1"$fastq_format
# fastq2=$fastq_dir/$dataset_name"_R2"$fastq_format

feather_output=$outdir"/feather_output/"$dataset_name"_"$DATE
if [ "$feather_output_symlink" == "" ]; then
feather_output_symlink=$outdir"/feather_output/"$dataset_name"_current"
fi

resolution=$(bc <<< "$bin_size/1000")
per_chr='True' # set this to zero if you don't want per chromosome output bed and bedpe files

feather_logfile=$feather_output"/"$dataset_name".feather.log"

resolution=$(bc <<< "$bin_size/1000")

if [ $organism == "mm10" ]; then
chr_count=19
elif [ $organism == "mm9" ]; then
chr_count=19
elif [ $organism == "hg19" ]; then
chr_count=22
elif [ $organism == "hg38" ]; then
chr_count=22
fi


####Ivan:"
sex_chroms_to_process="NA"
sex_chroms=""

long_bedpe_dir=$feather_output_symlink"/"
short_bed_dir=$feather_output_symlink"/"

maps_output=$outdir"/MAPS_output/"$dataset_name"/"
# maps_output_symlink=$outdir"/MAPS_output/"$dataset_name"_current"
# genomic_feat_filepath="/home/jurici/MAPS/MAPS_data_files/"$organism"/genomic_features/"$genomic_features_filename


# if [ $maps -eq 1 ]; then
mkdir -p $maps_output
echo "$dataset_name $maps_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output"
## create the parameter table
# $python_path $cwd/MAPS/make_maps_runfile.py $dataset_name $maps_output $macs2_filepath $genomic_feat_filepath $long_bedpe_dir $short_bed_dir $bin_size $chr_count $maps_output $sex_chroms_to_process --BINNING_RANGE $binning_range
# echo "first"
# $python_path $cwd/MAPS/MAPS_from_bedpe.py $maps_output"maps_"$dataset_name".maps"
# echo "second"
$Rscript_path $cwd/MAPS/MAPS_regression_and_peak_caller.r $maps_output $dataset_name"."$resolution"k" $bin_size $chr_count$sex_chroms $filter_file $model 
$Rscript_path $cwd/MAPS/MAPS_peak_formatting.r $maps_output $dataset_name"."$resolution"k" $fdr $bin_size
echo "third"
cp "$(readlink -f $0)" $maps_output"/execution_script_copy"
# chmod 777 $maps_output
# ln -sfn $maps_output $maps_output_symlink
# fi



