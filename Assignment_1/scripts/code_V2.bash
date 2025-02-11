

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
gstacks -I ./sorted_bam_files -M popmap.tsv -O ./g_stacks_output

# Run the populations command in order to convert Gstacks outputs into .vcf files for PCA
cd g_stacks_output
mkdir populations_output

#Why these parameters?
    # -r 0.8 : This ensures that a locus must have at least 80% of the samples covered to be included in the analysis - ensuring I have a high porportion of samples with coverage for each locus. This is chosen because it is a balanced choice - excludes loci with extremely poor coverage while still including a significant number of samples.
    # --min-maf 0.05 : Filters out very rare variants that might not be informative for population-level studies. In this case, I am ensuring that only loci with at least 5% frequency of the minor allele are included - This number was chosen as low-frequency variants (<5%) are often seen as noise or sequencing errors, and they can introduce false positives especially if working with large datasets.
    # --fstats - This is an important metric because it measures the level of genetic differentiation between populations. It compares the allele frequencies between populations and shows if there is a significnat divergence at each locus. Overall, this will help in visualizing and highlighting genetically divergent loci, and in this case, show regions of hte genome that show significant differentiation between populations.

populations -P ./g_stacks_output -M ../popmap.tsv --vcf -r 0.8 --min-maf 0.05 --fstats -O ./populations_output/


# Convert the .vcf file into PLINK format
cd populations_output
mkdir pca_data
vcftools --vcf ./populations.snps.vcf --plink --out ./pca_data/

# Convert PLINK files to binary format for PCA analysis
plink --file pca_data --make-bed --out pca_plink

# turn the pca_plink to the pca_results
plink --bfile pca_plink --pca 10 --out pca_results

#These pca_result files, along with the popmap, will be used to make the PCA plot