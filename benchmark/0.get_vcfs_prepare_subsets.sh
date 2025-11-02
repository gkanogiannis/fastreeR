#!/bin/bash

# Download raw 100k Genomes vcf files for 1-22 chromosomes
mkdir -p ${PWD}/vcf && pushd ${PWD}/vcf
for i in {1..22}; do
	wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
	wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
done
popd

# create incremental merged vcfs, chromosomes 1 with 2, 1 with 2 and 3, 1 with 2,3 amd 4, etc
pushd ${PWD}/vcf
for ((i=1; i<=22; i++)); do
    output=""
    line=""
    for ((j=1; j<=i; j++)); do
        line+="ALL.chr${j}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz "
    done
    output+="${line%" "}"$'\n'
    bcftools concat -Oz -o chr1-${i}.vcf.gz $output 
done
popd


# fixed 91 samples, varying number of SNPs (1 to 10 million)
# use file chr1-2.vcf.gz as it has enough snps (at least 10 million)
mkdir -p ${PWD}/vcf/subsets && pushd ${PWD}/vcf/subsets
SAMPLES=91
cp ../chr1-2.vcf.gz .
VCF=chr1-2.vcf.gz
bcftools index -f ${VCF}
SAMPLE_LIST=$(bcftools query -l $VCF | shuf | head -n $SAMPLES | paste -sd, -)
echo ${SAMPLE_LIST} > samples_S${SAMPLES}.list
bcftools view -s $SAMPLE_LIST -Oz -o chr1-2_S${SAMPLES}.vcf.gz $VCF
bcftools index -f chr1-2_S${SAMPLES}.vcf.gz
VCF=chr1-2_S${SAMPLES}.vcf.gz
for (( SNP_COUNT=1000000; SNP_COUNT<=10000000; SNP_COUNT+=1000000 )); do
    bcftools view -s $SAMPLE_LIST -H $VCF | head -n $SNP_COUNT | \
    awk '{print $1"\t"$2}' > regions_M`echo "$((SNP_COUNT / 1000000))"`.txt
    bcftools view -s $SAMPLE_LIST -R regions_M`echo "$((SNP_COUNT / 1000000))"`.txt -Oz -o M`echo "$((SNP_COUNT / 1000000))"`_S${SAMPLES}.vcf.gz $VCF
done
popd

# fixed 2 million SNPs, varying number of samples (100 to 2500, step 100)
# use file chr1-2.vcf.gz as it has enough samples (at least 91)
pushd ${PWD}/vcf/subsets
SNP_COUNT=2000000
VCF=../chr1-2.vcf.gz
ALL_SAMPLES=($(bcftools query -l $VCF))
for SAMPLE_COUNT in $(seq 100 100 2500); do
    SAMPLE_LIST=$(printf "%s\n" "${ALL_SAMPLES[@]}" | shuf | head -n $SAMPLE_COUNT | paste -sd, -)
    echo ${SAMPLE_LIST} > samples_${SAMPLE_COUNT}.txt
    bcftools view -S samples_${SAMPLE_COUNT}.txt -R regions_M2.txt -Oz -o M2_S${SAMPLE_COUNT}.vcf.gz $VCF
done
popd
