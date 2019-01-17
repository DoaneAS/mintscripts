#!/bin/bash -l
#$ -N Mintathena
#$ -j y
#$ -m e
#$ -cwd
#$ -l athena=true
#$ -M ashley.doane@gmail.com
#$ -l h_rt=22:00:00
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -R y
#$ -o /home/asd2007/joblogs


## This script takes a folder with
## mintchip fastq files and demuxes
## per sample barcodes


while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -f|--folder)
            path="$2"
            shift # past argument
            ;;
        -g|--gtf)
            gtf_path="$2"
            ;;
        *)
            # unknown option
            ;;
    esac
    shift
done


# Uses job array for each Sample in the folder
file=$(ls ${path} | tail -n +${SGE_TASK_ID}| head -1)

#change directory to node directory
cd $TMPDIR

echo "File : $file Read from $path\n"

#Obtain the name of the Sample
Sample=$(basename "$file")
Sample=${Sample%%.*}

rsync -r -v -a $path/$file/* ./


R1=Sample.R1.fastq.gz
R2=Sample.R2.fastq.gz

# prefix is chip target

mkdir -p logs 

barcode_splitter.py --bcfile $PWD/bcftab.txt \
                    --idxread 1 --mismatches 2 $R1 $R2 2> logs/$R1_$R2.detailed.log \
                    > logs/$PREFIX_$R1_$R2.log

mkdir -p
rsync -r -a -v $TMPDIR/ $path/${Sample}


