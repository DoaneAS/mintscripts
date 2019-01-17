#!/bin/bash -l
#$ -N Mintathena
#$ -j y
#$ -m e
#$ -cwd
#$ -l athena=true
#$ -l h_rt=24:00:00
#$ -pe smp 4-6
#$ -l h_vmem=4G
#$ -o /home/asd2007/joblogs


## This script takes mint chip fastq files after sample demux
## and aligns, calls peaks, and generates read density tracks.
## Several methods are used for peak calling and tracks

while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -f|--folder)
            FOLDERPATH="$2"
            shift # past argument
            ;;
        -g|--genome)
            GENOME="$2"
            shift # past argument
            ;;
        -t|--trim)
            TRIM=YES
            shift # past argument
            ;;
        --align)
            BWA=1
            ;;
        *)
            # unknown option
            ;;
    esac
    shift # past argument or value
done
echo FOLDER PATH   = "${FOLDERPATH}"
echo GENOME     = "${GENOME}"
echo TRIM READS    = "${TRIM}"
echo BWA ALN    = "${BWA}"
if [[ -n $1 ]]; then
    echo "Last line of FILEPATH specified as non-opt/last argument:"
    tail -1 $1
fi



# Uses job array for each Sample in the folder
FILEPATH=$(ls ${FOLDERPATH} | tail -n +${SGE_TASK_ID}| head -1)

#change directory to node directory
cd $TMPDIR

echo "File : $FILEPATH Read from $FOLDERPATH\n"

#Obtain the name of the Sample
Sample=$(basename "$FILEPATH")
Sample=${Sample%%.*}


## copy sample files to cluster local scratch

rsync -r -v -a $FOLDERPATH/$FILEPATH/* ./


## rsync symlinks for input file
rsync -L $FOLDERPATH/$FILEPATH/* ./ 

mkdir -p ${Sample}

echo "Processing  $Sample ..."

#Figuring out the reference genome


if [[ $ATHENA == 1 ]] ; then
    REFDIR="/athena/elementolab/scratch/asd2007/reference"
    ANNOTDIR="/athena/elementolab/scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.athena.conf"
    export PATH="/home/asd2007/anaconda2/bin:$PATH"
else
    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.conf"
fi

# get genome args
if [[ $GENOME == "hg19" ]] ; then
    REF="/athena/elementolab/scratch/asd2007/reference/hg19/bwa_index/male.hg19.fa"
    #REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	  REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed"
    #BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg19.chrom.sizes"
    cp $chrsz $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/reference/hg19/seq/male.hg19.fa"
elif [[ $GENOME == "hg38" ]] ; then
    echo "genome is ${GENOME}"
    DNASE_BED="${ANNOTDIR}/${GENOME}/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
	  REF="/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  REFbt2="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  bwt2_idx="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    RG="hg38"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    chrsz=/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes
    seq=/athena/elementolab/scratch/asd2007/reference/hg38/seq
    gensz=hs
    REF_FASTA=/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
elif [[ $GENOME == "mm10" ]]; then
    genome=mm10
	  gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Mus_UCSC_ref.gtf"
	  REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex"
    REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index"
    BLACK="/athena/elementolab/scratch/asd2007/reference/mm10-blacklist.bed"
    RG="mm10"
    SPEC="mm"
    REFGen="/athena/elementolab/scratch/asd2007/bin/bcbio/genomes/Mmusculus/mm10/seq/"
    fetchChromSizes mm10 > mm10.chrom.sizes
    chrsz=`pwd`/mm10.chrom.sizes
    rsync -av /home/asd2007/Scripts/picardmetrics.Mouse.conf ./
fi



echo "----------bwa-mem aligning-------------"
cat *read-1.fastq.gz > ${Sample}.R1.fastq.gz
cat *read-2.fastq.gz > ${Sample}.R2.fastq.gz

bwa mem -t ${NSLOTS} -M ${REF} $TMPDIR/${Sample}.R1.fastq.gz $TMPDIR/${Sample}.R2.fastq.gz | samtools view -bS - >  $TMPDIR/${Sample}/${Sample}.bam


samtools sort ${Sample}/${Sample}.bam -o ${Sample}/${Sample}.sorted.bam

samtools index ${Sample}/${Sample}.sorted.bam

mv ${Sample}/${Sample}.sorted.bam ${Sample}/${Sample}.bam

mv ${Sample}/${Sample}.sorted.bam.bai ${Sample}/${Sample}.bam.bai


##bwa mem -t 2 /athena/elementolab/scratch/asd2007/local/share/bcbio/genomes/Mmusculus/mm10/bwa/mm10.fa
#Smc3_het_NBC_7_AdBC07-read-1.fastq.gz Smc3_het_NBC_7_AdBC07-read-2.fastq.gz | sam

echo "----------Processing alignment and filtering for duplicates and mitochondrial mapping reads-----------------"

mkdir ${Sample}/outputMetrics


## read filtering for mint chip alignments.  outputs ${Sample}.sorted.nodup.bam
processMintBam.sh ${Sample}/${Sample}.bam


#mv $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam 

samtools index ${Sample}/${Sample}.sorted.nodup.bam

sambamba sort --memory-limit 20GB -n \
         -t ${NSLOTS} --out $TMPDIR/${Sample}/${Sample}.nsorted.nodup.bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam



samtools flagstat ${Sample}/${Sample}.bam > ${Sample}/outputMetrics/${Sample}.flagstat.log 


p1="${TMPDIR}/${Sample}/${Sample}.sorted.nodup.bam"

o1=$(echo ${p1} | sed -r 's/\.bam$/.tagAlign.gz/g')

o2=$(echo ${p1} | sed -r 's/\.bam$/.chipEncode.tagAlign.gz/g')

o3=$(echo ${p1} | sed -r 's/\.bam$/.bedpe.gz/g')

#o4=$(echo ${p1} | sed -r 's/\.bam$/.tn5.tagAlign.gz/g')


bamToBed -i ${p1} | awk 'BEGIN{OFS="\t"} $6=="+" { $2=$2+4; $3=$3 ; $4="N" ; print $0} $6=="-"{ $2=$2; $3=$3-5; $4="N" ; print $0}' | gzip -c > ${o4}

bedtools bamtobed -i ${p1} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${o2}

bedtools bamtobed -bedpe -mate1 -i ${p1} | gzip -c > ${o3}

#shifted_tag    = "$prefix.tn5.tagAlign.gz"

#zcat ${o2} | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0 }' | gzip -c > ${o4}



macs2 callpeak -t ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.bam -f BAMPE -n ${Sample}/${Sample}.broad --broad -g hs --keep-dup all --broad-cutoff 0.1 --bdg --SPMR

macs2 callpeak -t  $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam -f BAMPE -n $TMPDIR/${Sample}/${Sample}.nps.broad -g hs --nomodel --shift 37 --extsize 73 --keep-dup all --broad --broad-cutoff 0.1 --bdg --SPMR



#callPeaks.sh $TMPDIR/${Sample}/${Sample}.nodup.bedpe.gz ${Sample} $SPEC


rsync -r -a -v $TMPDIR/${Sample} $FOLDERPATH/${Sample}


mkdir -p $TMPDIR/${Sample}/bamPE


macs2 callpeak -t ${TMPDIR}/${Sample}/${Sample}.sorted.nodup.bam -f BAMPE -n $TMPDIR/${Sample}/bamPE/${Sample}.bamPE -g ${SPEC} --nomodel --keep-dup all --broad --broad-cutoff 0.1  --bdg --SPMR 


## get FE over background

getFC.sh $TMPDIR/${Sample}/bamPE/${Sample}.bamPE ${chrsz}


rsync -r -a -v $TMPDIR/${Sample}/bamPE $FOLDERPATH/${Sample}/



#### deeptools bamCoverage additional tracks
##
##bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam --binSize 1 \
##            --outFileFormat bigwig --MNase \
##            --minFragmentLength 70 \
##            --maxFragmentLength 300 \
##            --normalizeUsingRPKM \
##            -o $TMPDIR/${Sample}/${Sample}.mnase.70.300.bw --numberOfProcessors ${NSLOTS}
##
##
##bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam --binSize 5 \
##            --outFileFormat bigwig --smoothLength 150 \
##            --normalizeUsingRPKM \
##            -o $TMPDIR/${Sample}/${Sample}.smooth150.center.extend.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}
##
##bamCoverage --bam $TMPDIR/${Sample}/${Sample}.sorted.nodup.bam --binSize 20 \
##            --outFileFormat bigwig --smoothLength 300 \
##            --normalizeUsingRPKM \
##            -o $TMPDIR/${Sample}/${Sample}.smooth300.center.extend.bw --centerReads --extendReads --numberOfProcessors ${NSLOTS}
##
##
##

rsync -av $TMPDIR/${Sample} $FOLDERPATH/${Sample}
#adjustedBed="/home/ole2001/PROGRAMS/SOFT/bedtools2/bin/slopBed -i $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.black.bed -g sizes -l 75 -r -75 -s"
