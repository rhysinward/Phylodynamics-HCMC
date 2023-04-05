# Phylodynamics-HCMC

This respoitory contains all the code and data used for the analysis presented in 'Cryptic transmission and re-emergence of Cosmopolitan genotype of Dengue Virus Serotype 2 within Ho Chi Minh City and Southern Vietnam' - INSERT LINK TO PRE-PRINT WHEN UPLOADED


# Step by step process of pipeline 

## Step 1: Download genomic and metadata from Genbank

- Download data for DENV-2 from [GenBank](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Dengue%20virus%20type%202,%20taxid:11060&utm_source=data-hub)
- Here we downloaded both the fasta file and metadata for the sequences - N.B. we are able to do some additional filtering on the interface if needed
  - Make sure when you download the metadata you add the Country column

## Step 2: Process and clean metadata and output fasta

- Load the sequences and metadata and using the R script Processing_sequences_and_metadata.R to clean the metadata and output fasta sequences with correct naming 
- This code can be found here Code/Processing_sequences_and_metadata.R

## Step 3: Genome Detective 

- Due to errors in the labelling of DENV serotypes within the GenBank metadata we need to confirm serotype using [Genome Detective](https://www.genomedetective.com/app/typingtool/dengue/)
  - Note that genome detective has a maximum of 2000 sequences it can run for free
  - Code to process the outputs of Genomic Detective and produce the final datasets can be found at Code/genome_detective_processing.R
- This is the final processing step

## Step 3: Installing nextstrain and nextalign locally

#### Step 3a (if you're using the ubuntu virtual box)
We need to install some basic software such as [conda](https://docs.conda.io/en/latest/miniconda.html) and curl first:
```
sudo apt --assume-yes install curl
sudo apt --assume-yes install gcc
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
close the terminal & reopen it -- you'll now be in the default conda environment ("base").

#### Step 3b (if you're using MacOS)
You should have `curl` installed, but may not have `conda` -- run `which conda` to check!
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

#### Step 3c Install nextstrain, nextalign and activate environment
```
conda install -c bioconda nextalign
curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml
conda env create -f nextstrain.yml
conda activate nextstrain
npm install --global auspice
```

## Step 4: Alignment 

We aligned our sequences using MAFFT](https://mafft.cbrc.jp/alignment/software/) which is the defult of [nextstrain](https://nextstrain.org/). This can be done by :
```
conda activate nextstrain
augur align --sequences data/dengue_denv2_unaligned.fasta --reference-sequence config/KF955363.gb --output results/aligned.fasta --fill-gaps --remove-reference
```


Alternatively, [Nextalign](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextalign-cli.html) can also be used:
```
mkdir results
nextalign run \
--input-ref=Data/reference_denv1.fasta \
--output-fasta=results/dengue_denv1_aligned.fasta \
Data/dengue_denv1_unaligned.fasta
```

```
nextalign run \
--input-ref=Data/reference_denv2.fasta \
--output-fasta=results/dengue_denv2_aligned.fasta \
Data/dengue_denv2_unaligned.fasta
```
  - Ensure you are in the correct directory through the cd command 
  - Ensure you have changed the names of the fasta files or in the command if needed

After alignment we can run this code Code/Split_to_E_gene_and_WG.R to split dataset into E gene and WG, trim, and remove alignments with a large number of N's. We choose a threshold of 95% completeness.

## Step 5: Tree Building 

- We can generate a maximum-liklihood tree (ML-tree) using [IQ-tree](http://www.iqtree.org/) with the following command:

ML tree with boostrapping
```
nohup iqtree2 -m TIM2+F+R4 -s denv2_E_gene_cosmo_cleaned.fasta -bb 1000  -alrt 1000
```
ML tree without boostrapping
```
nohup iqtree2 -m TIM2+F+R4 -s denv2_E_gene_cosmo_cleaned.fasta
```
## Step 6: Quality Control

After we have produced an ML-tree we want to visuliase the output in [TempEst](http://tree.bio.ed.ac.uk/software/tempest/) and remove any problimatic sequences. 

- Make sure to download the output of the TempEst and use this code found Code/remove_bad_sequences_tempest.R to remove problamatic sequences
- Repeat step 5
- After step 5 is repeated, again visulise within TempEst to ensure all probalamtic sequences are removed

## Step 7: Time-scaled trees

- A time-scaled ML-tree can be created using [TreeTime](https://treetime.biozentrum.unibas.ch/)
- Additional metadata needs to be created using the R code found here Code/metadata_treetime.R
- NB to get confidence intervals for dates need to use non-bootstrapped tree

```
treetime --tree denv2_E_gene_cosmo_cleaned.fasta.treefile --aln denv2_E_gene_cosmo_cleaned.fasta --dates metadata_denv2_cosmo.csv  --clock-filter 4 --confidence
```
## Step 8: Mugration anlysis

N.b. non-bootstrapped tree needs to be used

```
treetime mugration --tree denv2_E_gene_cosmo_cleaned.fasta.treefile --states metadata_denv2_cosmo.csv --attribute country
```
