# genetic-evolutionary-mycobiome
Code for analyses and figure generation associated with the manuscript Van Syoc et al, under review.

Code is currently updated for submitted version of manuscript. Restricted dbGaP data for the Human Microbiome Project cannot be publicly shared and scripts using this information are redacted as necessary. The ITS2 sequencing for the Human Microbiome Project and hominid datasets are publicly available on SRA and the enclosed scripts show how to access and download this data. 

**Folders:**   

"amplicon-processing" - basic processing steps for ITS amplicon data for the Human Microbiome Project and hominid datasets  
1. `get-runinfo.sh`: inputs an SRA/EBI accession number and pulls basic project information from SRA and tests the connection to download read files  
2. `fastq-dump.sh`: inputs an SRA/EBI accession number and pulls reads from SRA, does fastqc and multiqc  
3. `ITS_qualitycontrol_hominid_HMP.R`: run with slightly different parameters for the HMP and hominid datasets; inputs raw sequencing files, then removes N's, primer sequences, and does basic quality filtering using the dada2 platform  
4. `ITS_vsearch_hominid_HMP.sh`: inputs a directory of filtered reads, then runs the VSEARCH pipeline to merge and dereplicate reads, then clusters OTUs and detects de novo chimeras.  
5. `vsearch-to-phyloseq_rarefy_filter_hominid.R`: for the hominid dataset, puts together the VSEARCH output into phyloseq and runs a rarefaction and filtering pipeline; output is the finalized phyloseq object for downstream analyses. Note: the HMP version is in the 'gwas' folder to follow a slightly different method.  


"gwas" - PLINK scripts to run GWAS and R scripts to run MWAS on the HMP dataset  
1. `plink-code.sh`: contains the PLINK code to do quality filtering on the HMP variants, run PCA for genetic stratification, and has the basic structure of the PLINK GLM  
2. `ITS-vsearch-to-phyloseq_prepGWAS.R`: for the HMP dataset, puts together the VSEARCH output into phyloseq and runs a rarefaction and filtering pipeline, then collapses taxonomic bins, CLR-transforms, and filters for 30% prevalence and gets into the PLINK GWAS structure. Outputs are phyloseq objects for downstream analysis and PLINK 'phenotype' file 
3. `get-sig-results.R`: combines PLINK output files, calculates genome-wide significance for each fungal taxa, writes out significant FAVs and preps for input to SNPNexus  
4. `GTEx-analysis.R`: Queries the GTEx lookup table to match b37 to b38 ID's, then queries the list of significant variant-eQTL associations for our FAVs  
5. `stats_Kazachstania.R`: runs basic stats for Kazachstania; mostly retracted since this data is restricted on dbGaP  
6. `phewas_IEUGWAS.R`: tests the two heart-associated FAVs against other coronary disease GWAS cohorts using the IEUGWAS platform  
7. `mendelian-randomization.R`: uses the TwoSampleMR package to conduct Two-Sample Mendelian Randomization analysis on the heart-associated FAVs  
8. Subfolder `wgs-scripts`: bash scripts to process HMP shotgun metagenomics data  
8.1 `match-ITS-to-WGS.R`: get sample IDs for ITS and matches them to sample IDs for WGS (should be so simple, but it is not!)  
8.2 `get-runinfo.sh`: similarly to above, takes PRJNA acession number and gets basic info  
8.3 `fastq-dump.sh`: gets reads from SRA  
8.4 `hmp-qualitycontrol.sh`: basic structure is taken from the original HMP bioinformatics  https://github.com/genome/genome/blob/master/lib/perl/Genome/Site/TGI/Hmp/HmpSraProcess/trimBWAstyle.usingBam.pl with a substitution made to use CD-HIT instead of Picard to remove duplicate reads (either the function is deprecated or I just can't find the right documentation)  
8.5 `humann-array.sh`: runs humann3.0 on the quality filtered reads  
8.6 `humann-regrouptables.sh`: input is humann3 default output, then regroups in various ways based on the Kolde 2018 analyses  
9. `mwas.R`: tests each fungal taxa and prokaryotic gene pathway in a linear model  

"phylosymbiosis" - scripts used to compare topological congruency of hominid and mycobiomes. Note; script to create random trees is located at, https://github.com/awbrooks19/phylosymbiosis (written in Python 2)
1. `topological_congruency.R`: runs TC tests at the OTU and Genus level; basic structure comes from Brooks et al 2016  

"codivergence" - scripts used to compare fungal exact sequence variants from each OTU to hominid phylogeny  
1. `get_genus_ESVs.R`: from phyloseq object, get exact sequence variants from the sequencing data for each genus and write out to fasta file  
2. `run-vsearch-derep.sh`: dereplicate exact sequence variants by length and identifier (decreases alignment size)  
3. `mafft-alignments.sh`: align sequences with MAFFT and trim with clipkit  
4. `run-seqkit-rename.sh`: uses seqkit to assign unique identifiers to all sequences; required by downstream tree building softwares  
5. `make_NJ_trees.R`: make neighbor-joining trees using the ape R package  
6. `replace_neg_branchlengths.R`: quick utility script to replace negative branch lengths; these are a side effect of neighbor-joining trees but can cause downstream issues  
7. `get_OTU_subtrees.R`: from each genus-level tree, subset to the OTUs that have representation from all four hominids.  
8. `run-PACO.R`: R script that is submitted as an HPC job to run PACo in a parallel loop (PACo is very slow on bigger trees)  
9. `run-parafit.R`: run Parafit on each OTU tree  
10. `results_paco-parafit.R`: combines results from PACo and Parafit, adjusts for multiple comparisons  
11. `determine_randomchance.R`: shuffles tip labels from the fungal OTU trees to create 100 random trees, then run both PACo and Parafit to determine the odds of observing cophylogeny given random chance  

"figure-generation" - R scripts that perform analyses and generate main and supplementary figures (most figures paneled in Illustrator but otherwise fully made in R)  
`fig1.R`  
`fig2.R`  
`fig3_and_figS8_codiversification.R`  
`figS1_QQplots.R`
`figS2_linkage_disequilibrium.R`  
`figS3_FAVannotation.R`
`figS4_fungi-eqtl-abundance.R`  
`figS5_kazachstania_abund-prev.R`  
`figS6_mwas.R`  
`figS7_hominid-alphadiv.R`  

"R" - miscellaneous helper functions and scripts  
1. `fx_myCollapse.R`: custom function to collapse count abundance at each taxonomic level 
2. `fx_myPvalue.R`: custom function to calculate probability of observing a given topological congruency score given random chance  
3. `colors.R`: sets a color palette for the hominid dataset for all figures  


