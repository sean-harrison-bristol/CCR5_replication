#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=00:12:00:00
#PBS -N GRS_creation
cd $PBS_O_WORKDIR
 
bed="/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/raw_downloaded/genotype_calls/ukb_cal_chr3_v2.bed"
bim="/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/raw_downloaded/genotype_calls/ukb_snp_chr3_v2.bim"
fam="/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/directly_genotyped/released/2017-07-04/data/raw_downloaded/data.fam"
plink --bed ${bed} --bim ${bim} --fam ${fam} --from-bp 45373456 --to-bp 47373456 --chr 3 --export A