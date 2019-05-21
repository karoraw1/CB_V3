#!/bin/bash

#SBATCH
#SBATCH --job-name=@SID@_trim
#SBATCH --time=6:00:00
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --exclusive
#SBATCH --partition=parallel
#SBATCH --mail-type=END
#SBATCH --mail-user=karoraw1@jhu.edu
#SBATCH --error=trim_@SID@.err
#SBATCH --output=trim_@SID@.out

module load R/3.5.1

# these will be filled in automatically
SCRIPTS_=~/scratch/CB_V3/otu_scripts/utility_scripts
TSTAT=~/scratch/CB_V3/otu_data/TrimParams.csv
TRIM_DATA=~/scratch/CB_V3/otu_data/trim_stats
SEQ_ID=@SID@
BASE_OUT=~/work/sprehei1/Keith_Files/Processed_data_group
DEMUX_DIR=$BASE_OUT/$SEQ_ID/FASTQ/Demux
TRIM_DIR=$BASE_OUT/$SEQ_ID/FASTQ/Trim

mkdir -p $TRIM_DIR $TRIM_DATA/FASTQC_Summaries $TRIM_DATA/Trim_Plots

Rscript $SCRIPTS_/FilterNTrim.R $SEQ_ID $DEMUX_DIR $TSTAT $TRIM_DIR

source activate otu_caller
source deactivate
source activate otu_caller
echo $(python --version)

python $SCRIPTS_/randomSampleofLibs.py $TRIM_DIR _F_filt.fastq
python $SCRIPTS_/randomSampleofLibs.py $TRIM_DIR _R_filt.fastq
cat $TRIM_DIR/*_F_filt.fastq.sample > $BASE_OUT/$SEQ_ID/FASTQ/RandomSampleT.R1.fastq
cat $TRIM_DIR/*_R_filt.fastq.sample > $BASE_OUT/$SEQ_ID/FASTQ/RandomSampleT.R2.fastq
rm $TRIM_DIR/*_filt.fastq.sample

mkdir $BASE_OUT/$SEQ_ID/FASTQ/pre-QC_report;
export PERL5LIB=~/scratch/miniconda3/envs/otu_caller/bin
fastqc -q -t 24 -o $BASE_OUT/$SEQ_ID/FASTQ/pre-QC_report -f fastq $BASE_OUT/$SEQ_ID/FASTQ/RandomSampleT.R1.fastq
mv $BASE_OUT/$SEQ_ID/FASTQ/pre-QC_report/RandomSampleT.R1_fastqc.html $TRIM_DATA/FASTQC_Summaries/${SEQ_ID}_fastqc_report.html
rm -rf FASTQ/pre-QC_report

module load R/3.5.1
Rscript $SCRIPTS_/quality_plot.R $BASE_OUT/$SEQ_ID FASTQ RandomSampleT ${SEQ_ID}_Trim
mv $BASE_OUT/$SEQ_ID/FASTQ/${SEQ_ID}_Trim_R1_and_R2_quals.png $TRIM_DATA/Trim_Plots
