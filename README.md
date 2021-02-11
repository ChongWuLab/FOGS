# FOGS

Transcriptome-wide association studies (TWAS) have been recently applied to successfully identify many novel genes associated with complex traits. While appealing, TWAS tend to identify multiple significant genes per locus, and many of them may not be causal due to confounding through linkage disequilibrium (LD) among SNPs. Here we introduce a powerful fine-mapping method called **FOGS** that prioritizes putative causal genes by accounting for local LD. We apply a weighted adaptive test with eQTL-derived weights to maintain high power across various scenarios.  

In this software, we implement **FOGS** as well as one alternate method [**FOCUS**](https://github.com/bogdanlab/focus).



Please cite the following manuscript for using this software:

>  Wu, C., & Pan, W. (2020). A powerful fine-mapping method for transcriptome-wide association studies. *Human genetics*, *139*(2), 199-213.




### Updates

1. Version 1.0: the preliminary release
2. Version 2.0: the standard alone release. It should take less than one hour to learn and configure the software.



## Outline

1. [Installation](#Installation)
2. [Typical analysis and output](#Analysis)
3. [Command-line parameters](#Command)
4. [FAQ](#FAQ)



## <a name="Installation"></a>Installation

- Download and unpackage the FOGS package from GitHub. Download through [this link](https://github.com/ChongWuLab/FOGS/releases/tag/2.0) or by the following commands:

  ~~~
  wget https://github.com/ChongWuLab/FOGS/releases/download/2.0/FOGS.zip
  ~~~

- Install the required packages (in R, run the following)

```R
list.of.packages <- c("data.table","optparse","Rcpp","RcppArmadillo","mvtnorm","BEDMatrix","bigmemory","dplyr","mvnfast")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages([,"Package"])]
if(length(new.packages)) install.packages(new.packages)
   
# install GenomicRanges
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
                                   
```





## <a name="Analysis"></a>Typical analysis and output

The FOGS analysis takes pre-computed Gene expression prediction models, LD reference panel, and GWAS summary data to prioritize putative causal genes. This example assumes you have set up the required environment and data, as illustrated in the previous section. 

To help users better use our software, we provided [a detailed pipeline]() for running FOGS with [COVID19-hg GWAS meta-analyses round 5 data](https://www.covid19hg.org/results/). *We will provide this pipeline within two weeks.*

### Input 1: GWAS summary statistics

We write a wrapping code (munge.R) and try to support all publically available GWAS summary data. Please run the following code to reformat GWAS summary data.

```R
source("munge.R")
sumstats = "/gpfs/research/chongwu/shared/summary_statistics/COVID19/release5/processed/ANA_B2_eur_V5.txt"

data = munge_sumstat(sumstats)
write.table(data,"processed_data.txt",col.names=TRUE,row.names=FALSE,quote=FALSE)
```

The processed data should look like:

| CHR  | POS    | SNP         | A1   | A2   | beta   | se     | N       | Z      |
| :--: | ------ | ----------- | ---- | ---- | ------ | ------ | ------- | ------ |
|  1   | 777232 | rs112618790 | T    | C    | 0.0135 | 0.0402 | 1175143 | 0.336  |
|  1   | 791853 | rs6684487   | A    | G    | 0.0298 | 0.0409 | 1172527 | 0.7274 |

**Note:** The performance of FOGS depends on the density of summary-level data. We highly recommend running FOGS with raw summary-level data. Pre-process steps such as pruning and restricting to top SNPs may harm the performance.

### Input 2: Gene expression prediction models

The format of prediction models are quite different for the different TWAS methods. To expand the usefulness of FOGS, we only requires the prediction models with the following format:

| rsid       | gene            | weight   | ref_allele | eff_allele |
| ---------- | --------------- | -------- | ---------- | ---------- |
| rs16861623 | ENSG00000000457 | -0.00015 | C          | A          |
| rs857633   | ENSG00000000457 | 0.0035   | C          | T          |



This format can be easily obtained. For example, for the PrediXcan type prediction models, we can use the following codes to obtain this:

```R
library(RSQLite)

weights = "/gpfs/research/chongwu/shared/TWAS_JTI/UTMOST_Lung.db"
weightsave = "Lung_weights.txt"

sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,dbname = weights)

dbListTables(db)
weights =  dbReadTable(db, "weights")

write.table(weights,weightssave,col.names=TRUE,row.names=FALSE,quote=FALSE)
```

*Note: See the pipeline for details.*



### Running the FOGS

After we prepared the data, we can run CMO via the following single line.

```
Rscript FOGS.R --refld /gpfs/research/chongwu/shared/1000Genomes/1000G.EUR.ALLSNP.QC.CHR --outd /gpfs/research/chongwu/Chong/MWAS/Finemapping/Blood --loci /gpfs/research/chongwu/shared/LDetect_LD_regions/EUR/ --weights /gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/Lung_weights.txt --genelist /gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/Lung_gene_list.txt --sumstat /gpfs/research/chongwu/Chong/Application/COVID19/Finemapping/processed_data.txt --saveprefix Blood --chr_id 21 --locus_id 3
```



### Output: Gene-disease association

The results are stored in a user-defined output file. For illustration, we explain the meaning of each entry in the first two lines of the output.

| Col. num. | Column name | Value           | Explanations                                   |
| --------- | ----------- | --------------- | ---------------------------------------------- |
| 1         | CHR         | 21              | Chromosome ID                                  |
| 2         | ID          | ENSG00000159110 | Ensemble ID                                    |
| 2         | P0          | 34602206        | Gene start                                     |
| 3         | P1          | 34637980        | Gene end                                       |
| 4         | n.SNP       | 98              | Number of SNPs in the prediction models        |
| 5         | n.condSNP   | 81              | Number of SNPs that FOGS conditioned on        |
| 6         | FOGS-aSPU   | 0.0003          | P value for FOGS (the underlying test is aSPU) |
| 7         | TWAS        | 0.00017         | P value for TWAS                               |
| 8         | Focus       | 0.99            | Posterior probability of FOCUS                 |
| 9         | Runtime(s)  | 583.16          | Running time in seconds                        |



## <a name="Command"></a>Command-line parameters

### FOGS.R

| Flag         | Usage                                                        | Default  |
| :----------- | ------------------------------------------------------------ | -------- |
| --refld      | LD reference panel directory (in plink bim/fam/bed format; chromosome specific) | Required |
| --outd       | The output directory                                         | Required |
| --loci       | The directory for the LDetect, which can be downloaded by [this link](https://bitbucket.org/nygcresearch/ldetect-data/src/master/). | Required |
| --weights    | The prediction model weights file (Input 2 described above)  | Required |
| --genelist   | The corresponding gene list for the weights files. It should contain the following four columns:left, right,chr, and gene. | Required |
| --sumstat    | Summary statistics (needs to be processed by munge.R; Input 1 described above) | Required |
| --saveprefix | The prefix for the output files                              | Required |
| --chr_id     | Chromosome ID for the locus of interest                      | Required |
| --locus_id   | Locus ID for the locus of interest                           | Required |

## <a name="Analysis"></a>FAQ

If you have questions, please submit an issue. We will summarize commonly asked questions here. 





## License

Maintainer: [Chong Wu](http://wuchong.org/index.html) (cwu3@fsu.edu)

[MIT](http://opensource.org/licenses/MIT)

Copyright (c) 2013-present, Chong Wu (cwu3@fsu.edu)