

# Retriveing the Assignment_1 .fq.gz files from Dr. Lukens' scratch to my own scratch
cp -r Assignment_1 /home/mirza/scratch

# Unzipping the .fq.gz files using a for loop
for i in ./*.gz; do gunzip "$i"; done

# Retrieving the most recent version of the fish genome
curl -o GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gasterosteus_aculeatus/reference/GCF_016920845.1_GAculeatus_UGA_version5/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna.gz

# Unzipping the reference genome
gunzip GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna.gz

# Viewing the reference genome
head GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna

# Loading the bwa module and indexing the fish genome for downstream alignment
module load bwa
bwa index GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna 

# Loop through each .fq file and convert it to a .sam file for alignment
for i in ./*.fq; do
    name="$(basename "$i")"
    bwa mem GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna "$i" > "$name.sam"
done

# Further convert the .sam files to .bam files for further downstream analysis
for i in ./*.sam; do
        name="$(basename "$i" .sam)"   #This will remove the .sam extension
        samtools view -b "$i" | samtools sort -o "$name.sorted.bam"
done

# This code changes the '.fq.sorted.bam' to just .bam for naming convention purposes
for i in ./*.fq.sorted.bam; do mv "$i" "${i/.fq.sorted/}"; done

# Move all .bam files to the directory 'sorted_bam_files' for better organization
mkdir sorted_bam_files
for i in ./*.bam; do mv "$i" sorted_bam_files/; done

# Create a population map for use in Gstacks
cd sorted_bam_files
ls *.bam | sed 's/\.bam//' | awk -F '_' '{print $0 "\t" $1}' > popmap.tsv

# Load the appropriate modules for running Gstacks
module load StdEnv/2020
module load stacks/2.64

# Using the population map, Run Gstacks on all bam files and output to directory g_stacks_output
mkdir g_stacks_output
gstacks -I . -M popmap.tsv -O ./g_stacks_output # the '.' indicates the current directory we are in (so, sorted_bam_files where all the BAM files are lcoated)

# Run the populations command in order to convert Gstacks outputs into .vcf files for PCA
cd g_stacks_output
populations -P . -O ./populations_output -M ../popmap.tsv --vcf

# Convert the .vcf file into PLINK format
cd populations_output
mkdir pca_data
vcftools --vcf ./populations.snps.vcf --plink --out pca_data

# Convert PLINK files to binary format for PCA analysis
plink --file pca_data --make-bed --out pca_plink

# turn the pca_plink to the pca_results
plink --bfile pca_plink --pca 10 --out pca_results

#These pca_result files, along with the popmap, will be used to make the PCA plot