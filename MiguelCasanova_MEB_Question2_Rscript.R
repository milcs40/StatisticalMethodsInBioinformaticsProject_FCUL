# Metodos Estatisticos em Bioinformatica
# Masters' Degree in Bioinformatics and Computational Biology
# FCUL - 2020/2021
# Jun 2021
# author: Miguel Casanova, fc24475

############################################################################################################################
########################## Let's start by loading all our required libraries ###############################################
############################################################################################################################
library("pheatmap")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library ("magrittr")
library("vsn")
library("NMF")
library("org.Hs.eg.db")
library("grDevices")
library("DESeq2")
library("ggrepel")
library("data.table")
library("EDASeq")
library("biomaRt")
library("plyr")
library("limma")
library("edgeR")
library("DESeq2")
library("ggvenn")
library("UpSetR")
library("baySeq")

############################################################################################################################
########################## Preparing the counts table and gene annotations #################################################
############################################################################################################################
# Set working directory
getwd() # shows the directory where R is currently looking for files and saving files to
setwd("G:/My Drive/Mestrado_BCG/4th_Semester/MetodosEstatisticosBioinformatica/Projects/Project_II/featureCounts/") # You can change the working directory

# Load relevant datafiles and metadata (in this case, I create a matrix)
countData <- read.table(file = "YuviaUnique/YuviaUnique_genes_gencode.counts", header = T, sep = "\t")
row.names(countData) = countData$Geneid # Give the Geneid as the row.names for the table.

# Remove all the columns that don't have counts
countData_simple = countData[ , -c(1:6)]

# Change the names for something a little bit more intelligible.
names(countData_simple) = c("ASD_1", "ASD_2", "ASD_3", "TCLAB_1", "TCLAB_2", "TCLAB_3")
# Check the first 8 lines.
head(countData_simple, 20)

# There are lots of genes with almost no reads. There are several ways to remove genes with low reads,
# that are statistically not robust to give confident DE results. 
# One of the ways of filtering, is simply getting rid of the genes that, across all samples, do not have at least 50 reads mapped. 
# We could consider increasing this number, as this is rather low.
#countData_simple <- countData_simple[rowSums(countData_simple) > 50, ]

##############################################
# Creating a metadata table for iPSCs
##############################################
# Next, we create a data frame with information about each sample, genotype and other relevant information.
sample_info <- data.frame(phenotype = c("Angelman", "Angelman", "Angelman", "Control", "Control", "Control"),
                          mutation_type = c("Megadeletion", "Megadeletion", "Megadeletion", "None", "None", "None"),
                          xist_status = c("+", "+", "+", "-", "-", "-"),
                          cell_line = c("AS-D", "AS-D", "AS-D", "F0002.1A.13", "F0002.1A.13", "F0002.1A.13"),
                          age_donor = c("Adult", "Adult", "Adult", "Adult", "Adult", "Adult"),
                          sex = c("Female", "Female", "Female", "Female", "Female", "Female"),
                          passage = c("P20", "P20", "P20", "P30-40", "P30-40", "P30-40"),
                          source = c("Fibroblasts", "Fibroblasts", "Fibroblasts", "Fibroblasts", "Fibroblasts", "Fibroblasts"),
                          method = c("Lentiviral", "Lentiviral", "Lentiviral", "Retroviral", "Retroviral", "Retroviral"),
                          type = c("iPSC", "iPSC", "iPSC", "iPSC", "iPSC", "iPSC"),
                          row.names = names(countData_simple))

# Even though we have lot's of information that we could use in our DEG models, we will only use the "phenotype" information.
# As such, we will define the "phenotype" as a factor.
sample_info$phenotype <- factor(sample_info$phenotype)
# The levels of the condition (phenotype, in this case), are important to determine the order of the comparison.
str(sample_info$phenotype)
# Levels in a factor are ordered alphanumerically by default, 
# but re-specification of the reference can be carried out using the relevel function in the stats package.
# Set "Control" as the first-level-factor. This is not absolutelly required. We can always get the information about it later, through contrast matrices.
sample_info$phenotype = relevel(sample_info$phenotype, "Control")
sample_info

#########################################################################################################################
######################################### Obtaining Gene annotations ####################################################
#########################################################################################################################
# Get names of all genes:
geneSymbols <- unique(rownames(countData_simple))
geneSymbols <- gsub("\\.[0-9]*$", "", geneSymbols)
head(geneSymbols)

# Use biomaRt to retrieve gene annotations from the genes in the count table.
# Let's first check the list of attributes in the datamart from Ensembl, containing all gene information.
# We will use this, to check which attributes we are going to collect for the genes we are working with.
#View(listAttributes(useDataset(dataset = "hsapiens_gene_ensembl",
#                               mart = useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")))) 

# For this project, I am sticking with genes that are annotated to the main chromosomes (and not mitochondrial or 
# other incomplete/unmapped scaffolds)
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", mirror = "www")
chrom=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, "X", "Y")

annotLookupChr <- getBM(
  mart = mart,
  attributes = c("external_gene_name", "ensembl_gene_id", "hgnc_symbol",
                 "chromosome_name", "gene_biotype", "ensembl_transcript_id",
                 "transcript_is_canonical", "transcript_length", 
                 "percentage_gene_gc_content"),
  filters = c("external_gene_name", "chromosome_name"),
  values = list(geneSymbols,chrom),
  uniqueRows = TRUE)

head(annotLookupChr)

# Several genes, have multiple transcript IDs. Biomart, has a canonical transcript annotation. 
# Let's Remove all entries that are non-canonical (have a NA in the "transcript_is_canonical" column)
cleanAnnotLookupChr <- na.omit(annotLookupChr)
head(cleanAnnotLookupChr)

# We can also decide to keep the transcripts with the largest "transcript_length".
# For this, we can use ddply:
annotLookupChrMax <- ddply(annotLookupChr, ~external_gene_name, function(d)d[which.max(d$transcript_length),])
head(annotLookupChrMax)

##############################################
# Using EDASeq
##############################################
# We can use the EDASeq library, to retrieve the GC and lenght of genes. For this, we first need to get 
# the EnsemblID for the genes (that we did above).

# geneEnsemblIDs <- unique(annotLookupChr$ensembl_gene_id[!is.na(annotLookupChr$ensembl_gene_id)])
# #geneEnsemblIDs <- geneSymbols
# 
# out <- getGeneLengthAndGCContent(
#   id = geneEnsemblIDs,
#   org = 'hsa',
#   mode = c('biomart', 'org.db'))
# 
# View(out[order(row.names(out)), ])

##############################################
# Creating gene annotation files
##############################################
# All of the above objects, take a considerable space in memory, and can take a while to run.
# For comodity, the retrieved annotations can be stored as .tsv files.

# write.table(as.data.frame(annotLookupChr), file = "annotLookupChr.tsv",
#              sep = "\t", quote = FALSE, row.names = TRUE)
# 
# write.table(as.data.frame(cleanAnnotLookupChr), file = "cleanAnnotLookupChr.tsv",
#             sep = "\t", quote = FALSE, row.names = TRUE)
# 
# write.table(as.data.frame(annotLookupChrMax), file = "annotLookupChrMax.tsv",
#             sep = "\t", quote = FALSE, row.names = TRUE)
# 
# write.table(as.data.frame(out), file = "out.tsv",
#             sep = "\t", quote = FALSE, row.names = TRUE)

#########################################################################################################################
############################################ Filtering counts ###########################################################
#########################################################################################################################
# For this exercise, we will create a single count table that will be used for all the different normalizations 
# and methodologies.

# As we have previously discussed, we can filter low read counts, in several ways. 
# Here. we start by filtering genes that have low counts. For this, we will only keep genes with an average read
# greater than 10 across all lanes (replicates - columns).
head(countData_simple)

filter <- apply(countData_simple,1,function(x) mean(x)>10)
head(filter)

table(filter)

# Next, we will keep the gene for which we have length and GC-content information.
geneGC <- annotLookupChrMax$percentage_gene_gc_content/100
names(geneGC) <- annotLookupChrMax$external_gene_name 
head(geneGC)

geneLength <- annotLookupChrMax$transcript_length
names(geneLength) <- annotLookupChrMax$external_gene_name 
head(geneLength)

# For doing so, we will intersect the original countData, filtered for count number,
# with genes that have information about the GC content (which are the same that have length information).
common <- intersect(names(geneGC),rownames(countData_simple[filter,]))
head(common)
length(common)

# Finally, we have our final count data, after all filtering has been applied.
countsFiltered <- as.matrix(countData_simple[common,])

############################################################################################################################
########################## Quality assessment of the libraries and pre-processing ##########################################
############################################################################################################################
# With this object created, we can start by looking at the individual libraries to see if we can identify a problem with them.
# As plotting raw number of counts will be uninformative, we will log2-transform the countsFiltered.
countsFilteredLog2 = log2(countsFiltered + 1)

# Next, we will plot the log2 read counts, to have a feeling about the dispersion of the data.
# For plotting, ggplot2 will be used.
# ggplot2 uses dataframes, for ploting. Let's prepare the dataframe for efficient plotting:
countsFilteredLog2 = as.data.frame(countsFilteredLog2)
df_raw = reshape2::melt(countsFilteredLog2, variable.name = "Samples", value.name = "count") #melt reshapes the matrix
# Create a column with the "ASD" and "CT" information for the different samples.
df_raw$Condition <- ifelse(df_raw$Samples == "ASD_1" | df_raw$Samples == "ASD_2" | df_raw$Samples == "ASD_3",
                       "ASD", "WT")

ggplot(df_raw, aes(x = Samples, y = count, fill = Condition)) + 
  theme_bw() +
  geom_point(position = "jitter", alpha = 0.1, size = 0.3, aes(color = Condition)) + 
  geom_boxplot(fill = NA) + xlab("") +
  ylab(expression(log[2](count + 1))) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Let's also visualize the count distribution in density plots
ggplot(df_raw, aes(x = count, colour = Samples, fill = Samples)) +
  theme_bw() +
  ylim(c(0,0.17)) +
  geom_density(alpha = 0.2, size = 0.5) +
  facet_wrap(~ Samples) +
  theme(legend.position = "top") +
  xlab(expression(log[2](count + 1)))

# The same can be represented with histograms.
ggplot(df_raw, aes(x = count, colour = Samples, fill = Samples)) +
  theme_bw() +
  geom_histogram(alpha = 0.2, size = 0.5, bins = 50) +
  #geom_density(alpha = 0.2, size = 0.5) +
  facet_wrap(~ Samples) +
  theme(legend.position = "top") +
  xlab(expression(log[2](count + 1)))

#######################################################################
# PCA Plot???
#######################################################################


############################################################################################################################
####################################### EDASeq normalization ###############################################################
############################################################################################################################
# We next create a feature data, to store GC-content and transcript length information
feature <- data.frame(gc=geneGC,length=geneLength)
head(feature)

# Finally, we will create a EDASeq ExpressionSet object, that will store:
# - Raw count table
# - Feature table, with Gc% and length information
# - Metadata table, with the phenotype of each sample.
data <- newSeqExpressionSet(
  counts = as.matrix(countData_simple[common,]), # Filtered count table
  featureData = feature[common,], # Information about GC contente and Length of filtered genes.
  phenoData = sample_info) # Metadata information

# Let's visualize the EDASeq object we have just created
data
# The associated counts table
head(counts(data))
# The metadata information
pData(data)
# And finally, the gene feature data, with GC content and gene length
head(fData(data))

# Let's now proceed with normalization. EDASeq allows for normalization, considering two types of effects on gene-level counts:
# 1- within-lane gene specific effects (as GC-contente or gene length);
# 2- Effects related to between-lane distributional differences (such as sequencing depth).
# Let's consider GC-content. Let's first visualize the dependence of gene-level counts on GC-content.
par(mfrow=c(1,1))
biasPlot(data, "gc", log=TRUE, ylim = c(3,8))
# We can see that there is a bias towards lower GC-content genes. 
# Let's tackle this, by doing a gene-level normalization, taking GC-content into consideration.
# It has been argued that its better to leave the count data unchanged to preserve the sampling properties and, instead,
# use an offset for normalization purposes in the statistical model. For this, we simply use the offset argument in both normalization functions. 
dataWithin <- withinLaneNormalization(data,"gc", which = "full", offset = TRUE)
# Next, we normalize for sequencing depth.
EDASeqNorm <- betweenLaneNormalization(dataWithin, which = "full", offset = TRUE)
# Let's visualize if this had any effect on the gene-level count dependence on GC-content.
biasPlot(EDASeqNorm, "gc", log=TRUE, ylim = c(3,8))
# As we can see, the bias is now almost eliminated.

# Let's have a look at the normalized read counts and the count distribution, after normalization.
countsEDASeqNorm <- normCounts(EDASeqNorm)
countsEDASeqNormLog2 = log2(countsEDASeqNorm + 1)

# ggplot2 uses dataframes, for ploting. Let's prepare the dataframe for efficient plotting:
countsEDASeqNormLog2 = as.data.frame(countsEDASeqNormLog2)
df_EDASeq = reshape2::melt(countsEDASeqNormLog2, variable.name = "Samples", value.name = "count") #melt reshapes the matrix
# Create a column with the "ASD" and "CT" information for the different samples.
df_EDASeq$Condition <- ifelse(df_EDASeq$Samples == "ASD_1" | df_EDASeq$Samples == "ASD_2" | df_EDASeq$Samples == "ASD_3",
                       "ASD", "WT")

ggplot(df_EDASeq, aes(x = Samples, y = count, fill = Condition)) + 
  theme_bw() +
  geom_point(position = "jitter", alpha = 0.1, size = 0.3, aes(color = Condition)) + 
  geom_boxplot(fill = NA) + xlab("") +
  ylab(expression(log[2](count + 1))) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################################################################################
######################### Limma Normalization using the EdgeR package ######################################################
############################################################################################################################
# edgeR can use the same count table we have previously defined and our metadata table.
# Use the DGEList() function to convert the count matrix into the edgeR object
edgeR_DGElist <- DGEList(counts = countsFiltered,
                         group = sample_info$phenotype)

# Check results
head(edgeR_DGElist$counts)
edgeR_DGElist$samples

# Calculate normalization factors using trimmed mean M-values (TMM).
edgeR_DGElist <- calcNormFactors(edgeR_DGElist, method = "TMM")
edgeR_DGElist$samples

# Let's have a look at the normalized read counts and the count distribution, after normalization.
# EdgeR provides the pseudo-counts as a CPM value (counts per milion). For this, we run the cpm funcion, 
# and use log argument, to have log2CPM counts.
# As edgeR works with CPMs, we will also create CPM pseudocounts of raw counts.
countsRawLog2CPM <- cpm(edgeR_DGElist, log=TRUE, normalized.lib.sizes = F)
countsEdgeRNormLog2CPM <- cpm(edgeR_DGElist, log=TRUE) 

# Pseudocounts - Raw CPMs
# ggplot2 uses dataframes, for ploting. Let's prepare the dataframe for efficient plotting:
countsRawLog2CPM = as.data.frame(countsRawLog2CPM)
df_edgeR_raw = reshape2::melt(countsRawLog2CPM, variable.name = "Samples", value.name = "count") #melt reshapes the matrix
# Create a column with the "ASD" and "CT" information for the different samples.
df_edgeR_raw$Condition <- ifelse(df_edgeR_raw$Samples == "ASD_1" | df_edgeR_raw$Samples == "ASD_2" | df_edgeR_raw$Samples == "ASD_3",
                                 "ASD", "WT")

ggplot(df_edgeR_raw, aes(x = Samples, y = count, fill = Condition)) + 
  theme_bw() +
  geom_point(position = "jitter", alpha = 0.1, size = 0.3, aes(color = Condition)) + 
  geom_boxplot(fill = NA) + xlab("") +
  ylab(expression(log[2](CPM))) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Pseudocounts - Normalized CPMs
# ggplot2 uses dataframes, for ploting. Let's prepare the dataframe for efficient plotting:
countsEdgeRNormLog2CPM = as.data.frame(countsEdgeRNormLog2CPM)
df_edgeR_TMM = reshape2::melt(countsEdgeRNormLog2CPM, variable.name = "Samples", value.name = "count") #melt reshapes the matrix
# Create a column with the "ASD" and "CT" information for the different samples.
df_edgeR_TMM$Condition <- ifelse(df_edgeR_TMM$Samples == "ASD_1" | df_edgeR_TMM$Samples == "ASD_2" | df_edgeR_TMM$Samples == "ASD_3",
                       "ASD", "WT")

ggplot(df_edgeR_TMM, aes(x = Samples, y = count, fill = Condition)) + 
  theme_bw() +
  geom_point(position = "jitter", alpha = 0.1, size = 0.3, aes(color = Condition)) + 
  geom_boxplot(fill = NA) + xlab("") +
  ylab(expression(log[2](CPM))) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################################################################################
############################################# DESeq2 normalization #########################################################
############################################################################################################################
# For the DESeq2 normalization, let's start by creating the DESeq2 object for all further analysis.
# What you put on the design matrix is very important. You need to know exactly what effect each feature has on the data.
# Also, you need to be careful with not overfitting your design matrix, into the specificity of your data. 
# As you provide more variables, you increase the risk of overfitting.
dds = DESeqDataSetFromMatrix(countData = countsFiltered,
                             colData = sample_info,
                             design = ~ phenotype)
dds

# DESeq2 does normalization, by estimating a size factor for each sample. Each gene count for that sample, is then multiplied by the factor. 
# To add the size factor to the dds, we do the following:
dds = estimateSizeFactors(dds)
sizeFactors(dds)

# Let's have a look at the normalized read counts and the count distribution, after normalization.
countsDESeq2Norm = counts(dds, normalized = TRUE) # Extract normalized counts
countsDESeq2NormLog2 = log2(countsDESeq2Norm + 1) # Convert to log2-scale for visualization, as we did before.

# ggplot2 uses dataframes, for ploting. Let's prepare the dataframe for efficient plotting:
countsDESeq2NormLog2 = as.data.frame(countsDESeq2NormLog2)
df_DESeq2 = reshape2::melt(countsDESeq2NormLog2, variable.name = "Samples", value.name = "count") #melt reshapes the matrix
# Create a column with the "ASD" and "CT" information for the different samples.
df_DESeq2$Condition <- ifelse(df_DESeq2$Samples == "ASD_1" | df_DESeq2$Samples == "ASD_2" | df_DESeq2$Samples == "ASD_3",
                                 "ASD", "WT")

ggplot(df_DESeq2, aes(x = Samples, y = count, fill = Condition)) + 
  theme_bw() +
  geom_point(position = "jitter", alpha = 0.1, size = 0.3, aes(color = Condition)) + 
  geom_boxplot(fill = NA) + xlab("") +
  ylab(expression(log[2](count + 1))) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################################################################################
########################################### Testing Normalization ##########################################################
############################################################################################################################
# We can visualize all the normalization methodologies employed, and how do they influence the gene count distribution/dispersion.
df_raw$method <- rep("Raw counts", nrow(df_raw))  
head(df_raw)

df_EDASeq$method <- rep("EDASeq", nrow(df_EDASeq))  
head(df_EDASeq)

df_DESeq2$method <- rep("DESeq2", nrow(df_EDASeq))  
head(df_DESeq2)

df_edgeR_raw$method <- rep("edgeR raw", nrow(df_EDASeq))  
head(df_edgeR_raw)

df_edgeR_TMM$method <- rep("edgeR TMM", nrow(df_EDASeq))  
head(df_edgeR_TMM)

df_allnorm <- rbind(df_raw, df_EDASeq, df_DESeq2, df_edgeR_raw, df_edgeR_TMM)
df_allnorm$method <- factor(df_allnorm$method, levels = c("Raw counts", "EDASeq", "DESeq2", "edgeR raw", "edgeR TMM"))

ggplot(df_allnorm, aes(x=Samples, y=count)) +
  theme_bw() +
  geom_point(position = "jitter", alpha = 0.02, size = 0.01, aes(color = Condition)) + 
  geom_boxplot(fill = NA) + xlab("") +
  ggtitle("Boxplots of normalized pseudo counts\nfor all samples by normalization methods") +
  facet_grid(. ~ method) +
  guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.8))) +
  ylab(expression(log[2] ~ (normalized ~ count + 1))) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

############################################################################################################################
############################# edgeR Differential Expression after EDASeq normalization #####################################
############################################################################################################################
# edgeR can use the EDASeq count table (with offsets) and the metadata table.
# Use the DGEList() function to convert the count matrix into the edgeR object
# For our example, we will use as a group for the DE analysis, the phenotype of the samples ("Angelman" or "Control").
EDASeqNorm_edgeR_DGElist <- DGEList(counts = counts(EDASeqNorm),
                                    group = sample_info$phenotype)
# We will add the offset object obtained from the EDASeq normalization to the edgeR object,
# to be used for the statistical analysis for differential expression.
EDASeqNorm_edgeR_DGElist$offset <- -offst(EDASeqNorm)

# Check results
head(EDASeqNorm_edgeR_DGElist$counts)
head(EDASeqNorm_edgeR_DGElist$offset)
EDASeqNorm_edgeR_DGElist$samples
# Notice that here, the normalization factors haven't been corrected. This is due to the use of the offset argument, 
# when performing the EDASeq normalization. edgeR will use the raw counts, together with the offset, to maintain the original
# sampling distribution of read counts.

# A very important step, on differential expression analysis, is the definition of the design matrix.
# We will use the same design matrix for all of our methods, namely on using the phenotype information as the 
design <- model.matrix(~ sample_info$phenotype)
rownames(design) <- colnames(EDASeqNorm_edgeR_DGElist)

# To calculate the differential expression in a gene-wise manner, edgeR first estimates the dispersion.
# It then tests whether the observed gene counts fit the respective negative binomial model.
# Estimate the dispersion for all read counts across all samples (common and tagwise dispersion)
EDASeqNorm_edgeR_DGElist <- estimateDisp(EDASeqNorm_edgeR_DGElist, design)
#EDASeqNorm_edgeR_DGElist <- estimateGLMRobustDisp(EDASeqNorm_edgeR_DGElist, design)

plotBCV(EDASeqNorm_edgeR_DGElist)

# Fit a negative binomial model
edgeR_fit <- glmFit(EDASeqNorm_edgeR_DGElist, design)

# Perform the testing for every gene using the neg. binomial model
edgeR_LRT  <- glmLRT(edgeR_fit)

# Obtain the log2 fold changes
DGE_results_edgeR <- topTags(edgeR_LRT, n = Inf, # to retrieve all gene
                             sort.by = "P", adjust.method = "BH")

head(DGE_results_edgeR, 10)

############################################################################################################################
########################## Limma Differential Expression after edgeR/limma normalization ###################################
############################################################################################################################
# Limma uses the edgeR_DGElist_filtered and same design matrix
# Transform the count data to log2cpm, and estimate mean-variance relatioship.
# This will be used to compute the weights for each count, that can be used
# to be able to apply linear models to the read counts
design <- model.matrix (~ sample_info$phenotype)
rownames(design) <- colnames(edgeR_DGElist)
# We will use the voom function to transform the data, for linear modeling.
# Normally, for RNA-seq, the normalization methods used for microarray data aren't used.
# We can add a normalization method to the voom function, keeping in mind that this should only be used 
# if the data is very noise. Which isn't the case.
voomTransformed <- voom(edgeR_DGElist, design, plot = TRUE)

# What is voom doing?
  
# - Counts are transformed to log2 counts per million reads (CPM), where "per million reads" is defined based on 
#   the normalization factors we calculated earlier
# - A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
# - A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot)
# - The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.

# Fit a linear model using weighted least squares for each gene.
voom_fit <- lmFit(voomTransformed, design = design)

# Apply Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or 
# smaller than those from other genes towards the average standard error)
# This also allows computing moderated t-statistics, moderated F-statistics, and log-odds of differential expression
voom_fit <- eBayes(voom_fit)

plotSA(voom_fit, main="Final model: Mean-variance trend")


# Extract DGE list
# Check how the coefficient is names
colnames(design)
DGE_results_limma <- topTable(voom_fit, coef = "sample_info$phenotypeAngelman",
                              number = Inf, adjust.method = "BH", sort.by = "P")

head(DGE_results_limma, 10)

############################################################################################################################
############################## DESeq2 Differential Expression after DESeq2 normalization ###################################
############################################################################################################################
# As the DESeq2 has already been specified with an experimental design, when we first created the DESeqDataSet, we can now just run the differential expression pipeline
# to the raw counts with a single call to the function DESeq
dds = DESeq(dds) # Run the differential expression pipeline on the raw counts. ALWAYS raw counts!!!
dds

plotDispEsts(dds)

# Extract the results, either using a defined contrast or letting the program choose.
# Using the contrast matrix, tries to analyze the differential gene expression between the conditions provided, while taking all information into consideration.
# This is the advantage of using these models, as it confers an added efficiency of predicting DEGs that are really associated with the different conditions provided in the contrast matrix.
DGE_results_DESeq2 = results(dds,  independentFiltering = TRUE) 
# Using independentFiltering = FALSE, the results obtained are not independently filtered, using the mean of normalized counts as a filter statistic.
# This shouldn't be used as a standard approach, as it might give less robust results.
mcols(DGE_results_DESeq2, use.names = TRUE)
summary(DGE_results_DESeq2)

head(DGE_results_DESeq2[order(DGE_results_DESeq2$padj),], 10)

############################################################################################################################
############################## Bayseq differential expression after edgeR normalization ####################################
############################################################################################################################
# baySeq uses empirical Bayesian methods to estimate the posterior likelihoods of each of a set
# of models that define patterns of differential expression for each row. This approach begins
# by considering a distribution for the row defined by a set of underlying parameters for which
# some prior distribution exists. By estimating this prior distribution from the data, we are able
# to assess, for a given model about the relatedness of our underlying parameters for multiple
# libraries, the posterior likelihood of the model.

# To define our model, we start by defining it's parameters. Namely, the replicates:
replicates <- as.factor(c("Angelman","Angelman","Angelman",
                          "Control","Control","Control")) 
# And the groups that will be tested. For this, we assume that certain genes won't change ("NDE"),
# while some will change ("DE"). For the genes that don't change, we assume that the expression has no
# dependance of the groups and, as such, all samples are considered within the same group.
# For genes which expression will depend on the group, we define two levels, one for control (1), and
# one for Angelman (2).
groups <- list(NDE = factor(rep(1,6)),
               DE = factor(c("2","2","2","1","1","1")))

# We can now create our countData object.
CD <- new("countData",
          data = countsFiltered,
          replicates = replicates, 
          groups = groups)

# Similarly to edgeR and DESeq2, baySeq assumes a negative binomial distribution for gene counts.
# Let's set a negative binomial density function in the countData object.
densityFunction(CD) <- nbinomDensity

# For normalization, baySeq uses the library size and normalizes following different estimates.
# Let's estimate library size for countData object, using TMM normalization (similar to edgeR).
libsizes(CD) <- getLibsizes(CD, estimationType= "edgeR")
head(libsizes(CD))

# We need to finally provide the gene names as an annotation to the countData object.
geneNames <- rownames(countsFiltered)
CD@annotation <- as.data.frame(geneNames)

# baySeq can use multiple CPU threads, increasing the performance of the algorithm.
# For this, we use the following function:
cl <- makeCluster(6, "PSOCK")

# We first estimate an empirical distribution on the parameters of the Negative Binomial distribution by bootstrapping from the data, 
# taking individual counts and finding the quasilikelihood parameters for a Negative Binomial distribution. 
# For efficiently estimating the empirical posterior distribution, a large enough sample size has to be used (10e5, in this case).
#CD <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl, verbose = T)
CD <- getPriors(CD, samplesize = 10^5, cl = cl, verbose = T)

plotPriors(CD, group= "NDE")

# We then acquire posterior likelihoods, estimating the proportions of differentially expressed
# counts.
#CD <- getLikelihoods.NB(CD, pET = 'BIC', nullData = T, cl = cl, verbose = T)
CD <- getLikelihoods(CD, pET = 'BIC', cl = cl, verbose = T)

# To finalize, we can ask for the top candidates for differential expression using the topCounts function.
# Here, we ask for all of the results for all of the genes (with Inf)
DGE_results_Bayseq <- topCounts(CD, group = "DE", number = Inf, normaliseData = TRUE)
head(DGE_results_Bayseq, 10)

# We next make the results table more easy to understand and to add the log2FC.
DGE_results_Bayseq[c("log2FC")] <-NA
DGE_results_Bayseq[c("Control","Angelman")] <-NA
rownames(DGE_results_Bayseq) <- DGE_results_Bayseq$geneNames

DGE_results_Bayseq$Control = rowMeans(DGE_results_Bayseq[,c(5:7)])+1
DGE_results_Bayseq$Angelman = rowMeans(DGE_results_Bayseq[,c(2:4)])+1
DGE_results_Bayseq$log2FC = log2((DGE_results_Bayseq$Angelman/DGE_results_Bayseq$Control))
head(DGE_results_Bayseq, 10)

###########################################################################################################################
############################## RankProd differential expression after limma normalization #################################
###########################################################################################################################
# library(RankProd)
# 
# cl <- c(1,1,1,0,0,0)
# fc <- RankProducts(as.matrix(voomTransformed), cl, logged = TRUE, na.rm = TRUE, gene.names = rownames(voomTransformed) ,
#                  plot = FALSE, rand = NULL, calculateProduct = TRUE, MinNumOfValidPairs = NA,
#                  RandomPairs = NA, huge = FALSE, fast = TRUE, tail.time = 0.05)
# 
# fc
# met3 <- topGene(fc, cutoff = NULL, method="pfp",logged=TRUE,
#               logbase=2, gene.names = rownames(voomTransformed), num.gene = length(rownames(voomTransformed)))
# upregulatedGenes <- as.data.frame(met3$Table1)
# downregulatedGenes <- as.data.frame(met3$Table2)
# downregulatedGenes$`FC:(class1/class2)` <- downregulatedGenes$`FC:(class1/class2)`*-1
# DGE_results_RankProd <- rbind(upregulatedGenes, downregulatedGenes)
# 
# plotRP(fc, cutoff = 0.05)

############################################################################################################################
########################################## Comparing the different methods #################################################
############################################################################################################################

#########################################
# Venn diagram and Upset Plots for DEGs #
#########################################
# Let's get a list of significantly misregulated genes  
DE_list <- list(edgeR = c(rownames(subset(DGE_results_edgeR$table, FDR <= 0.05))),
                deseq2 = c(rownames(subset(DGE_results_DESeq2, padj <= 0.05))), 
                limma = c(rownames(subset(DGE_results_limma, adj.P.Val <= 0.05))),
                bayseq = c(rownames(subset(DGE_results_Bayseq, FDR.DE <= 0.05))))

# The code below, filters gene more stringently, by applying a log2FC threshold of 2.
# DE_list <- list(edgeR = c(rownames(subset(DGE_results_edgeR$table, FDR <= 0.05 & logFC >= 2 | FDR <= 0.05 & logFC <= -2))),
#                 deseq2 = c(rownames(subset(DGE_results_DESeq2, padj <= 0.05 & log2FoldChange >= 2 | padj <= 0.05 & log2FoldChange <= -2))), 
#                 limma = c(rownames(subset(DGE_results_limma, adj.P.Val <= 0.05 & logFC >= 2 | adj.P.Val <= 0.05 & logFC <= -2))),
#                 bayseq = c(rownames(subset(DGE_results_Bayseq, FDR.DE <= 0.05 & log2FC >= 2 | FDR.DE <= 0.05 & log2FC <= -2))))

ggvenn(DE_list, 
       stroke_linetype = 3, 
       stroke_size = 0.5,
       set_name_color = "black",
       set_name_size = 5,
       fill_color = c("coral", "darkgoldenrod1", "lightblue", "darkseagreen"),
       fill_alpha = 0.4,
       show_percentage = T,
       digits = 1)

# Comparing DEGs between methods
DE_gns <- fromList(DE_list)
upset(DE_gns, 
      order.by = "freq", 
      text.scale = 1.5,
      sets.bar.color=c("coral", "darkgoldenrod1", "lightblue", "darkseagreen"))

#########################################
# Correlation plots for logFC for DEGS ##
#########################################
# Correlation of logFC for genes found DE in all four methods
table_DE_methods <- table(stack(DE_list))
DE_gns_all <- row.names(table_DE_methods[rowSums(table_DE_methods) == 4, ]) # Extract names of rows

# Make data frame of fold change values
DE_fc <- data.frame(edgeR = DGE_results_edgeR[DE_gns_all,]$table$logFC,
                    limma = DGE_results_limma[DE_gns_all,]$logFC,
                    deseq2 = DGE_results_DESeq2[DE_gns_all,]$log2FoldChange,
                    bayseq = DGE_results_Bayseq[DE_gns_all,]$log2FC,
                    row.names = DE_gns_all)
pairs(DE_fc)

#########################################
########## Heatplots with DEGs ##########
#########################################
pheatmap(as.matrix(DE_fc),
         cluster_rows = T,
         show_rownames = F,
         border_color = NA)

# For the cover image.
# pheatmap(as.matrix(DE_fc),
#          cluster_rows = T,
#          show_rownames = F,
#          show_colnames = F,
#          border_color = NA,
#          treeheight_row = 0,
#          treeheight_col = 0,
#          legend = F)

#########################################
############# Volcano Plots #############
#########################################

#############
### edgeR ###
#############

# Volcano plot with gene name information and information about gene LFC and statistics
# We first create a columns with gene names that will be used to label the 20 most significant DEGs.
results_ordered_edgeR <- data.frame(DGE_results_edgeR)
results_ordered_edgeR <- results_ordered_edgeR[order(results_ordered_edgeR$FDR), ]
results_ordered_edgeR$genelabels = ""
results_ordered_edgeR$genelabels[1:20] <- rownames(results_ordered_edgeR)[1:20]
results_ordered_edgeR["UBE3A", "genelabels"] <- "UBE3A"

# Next, we crate another column to include the LFC and statistical information for the genes.
results_ordered_edgeR$Class <- "NS"
results_ordered_edgeR$Class[which(results_ordered_edgeR$logFC > 2 | results_ordered_edgeR$logFC < -2)] <- "logFC > |2|"
results_ordered_edgeR$Class[which(results_ordered_edgeR$FDR < 0.05)] <- "FDR < 0.05"
results_ordered_edgeR$Class[which(results_ordered_edgeR$FDR < 0.05 & (results_ordered_edgeR$logFC > 2 | results_ordered_edgeR$logFC < -2))] <- "FDR < 0.05 & logFC > |2|"

# Now, we only have to plot the volcano plot.
ggplot(results_ordered_edgeR) +
  theme_bw() +
  geom_point(aes(x = logFC, y = -log10(FDR), color = Class), size = 2, alpha = .4) +
  scale_colour_manual(breaks = c("NS", "logFC > |2|", "FDR < 0.05", "FDR < 0.05 & logFC > |2|"),
                      values = c("gray50", "seagreen", "royalblue", "red3")) + # Add a manual scale color (and order), for the elements in column "Class".
  geom_text_repel(max.overlaps = 100000,  
                  force = 10,
                  aes(x = logFC, y = -log10(FDR), label = ifelse(genelabels != "", rownames(results_ordered_edgeR),""))) + # Add gene names to most significant DEGs
  ggtitle("Volcano plot for edgeR\nFDR < 0.05") +
  xlab("log2 fold change") + 
  ylab("-log10 FDR") +
  scale_y_continuous(trans = "log1p") +
  #scale_y_continuous(limits = c(0,700)) +
  #scale_x_continuous(limits = c(-10.5, 10.5)) +
  geom_hline(yintercept = -log10(0.05), colour = "firebrick", linetype = 4) + # Add line to intercept -log10(0.05).
  geom_vline(xintercept = c(-2, 2), colour = "firebrick", linetype = 4) + # Add lines to intercept the + and - 2 log2FC.
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + # Over-ride aesthetics for figure legend.
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white",
                                         size = 0.5, 
                                         linetype = 3, 
                                         colour ="black"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12)) + # Add several theme features.
  annotate("rect", xmin = c(-2, 2), xmax = c(-Inf, Inf), ymin = -log10(0.05), ymax = Inf, alpha = .1, fill = "mistyrose") # Add a rectangle to shade the area of DEGs.

#############
### limma ###
#############

# Volcano plot with gene name information and information about gene LFC and statistics
# We first create a columns with gene names that will be used to label the 20 most significant DEGs.
results_ordered_limma <- data.frame(DGE_results_limma)
results_ordered_limma <- results_ordered_limma[order(results_ordered_limma$adj.P.Val), ]
results_ordered_limma$genelabels = ""
results_ordered_limma$genelabels[1:20] <- rownames(results_ordered_limma)[1:20]
results_ordered_limma["UBE3A", "genelabels"] <- "UBE3A"

# Next, we crate another column to include the LFC and statistical information for the genes.
results_ordered_limma$Class <- "NS"
results_ordered_limma$Class[which(results_ordered_limma$logFC > 2 | results_ordered_limma$logFC < -2)] <- "logFC > |2|"
results_ordered_limma$Class[which(results_ordered_limma$adj.P.Val < 0.05)] <- "Padj < 0.05"
results_ordered_limma$Class[which(results_ordered_limma$adj.P.Val < 0.05 & (results_ordered_limma$logFC > 2 | results_ordered_limma$logFC < -2))] <- "Padj < 0.05 & logFC > |2|"

# Now, we only have to plot the volcano plot.
ggplot(results_ordered_limma) +
  theme_bw() +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), color = Class), size = 2, alpha = .4) +
  scale_colour_manual(breaks = c("NS", "logFC > |2|", "Padj < 0.05", "Padj < 0.05 & logFC > |2|"),
                      values = c("gray50", "seagreen", "royalblue", "red3")) + # Add a manual scale color (and order), for the elements in column "Class".
  geom_text_repel(max.overlaps = 100000,  
                  force = 20,
                  aes(x = logFC, y = -log10(adj.P.Val), label = ifelse(genelabels != "", rownames(results_ordered_limma),""))) + # Add gene names to most significant DEGs
  ggtitle("Volcano plot for limma\nPadj < 0.05") +
  xlab("log2 fold change") + 
  ylab("-log10 Padj") +
  scale_y_continuous(trans = "log1p") +
  #scale_y_continuous(limits = c(0,700)) +
  #scale_x_continuous(limits = c(-10.5, 10.5)) +
  geom_hline(yintercept = -log10(0.05), colour = "firebrick", linetype = 4) + # Add line to intercept -log10(0.05).
  geom_vline(xintercept = c(-2, 2), colour = "firebrick", linetype = 4) + # Add lines to intercept the + and - 2 log2FC.
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + # Over-ride aesthetics for figure legend.
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white",
                                         size = 0.5, 
                                         linetype = 3, 
                                         colour ="black"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12)) + # Add several theme features.
  annotate("rect", xmin = c(-2, 2), xmax = c(-Inf, Inf), ymin = -log10(0.05), ymax = Inf, alpha = .1, fill = "mistyrose") # Add a rectangle to shade the area of DEGs.


#############
## DESeq2 ###
#############

# Volcano plot with gene name information and information about gene LFC and statistics
# We first create a columns with gene names that will be used to label the 20 most significant DEGs.
results_ordered_DESeq2 <- data.frame(DGE_results_DESeq2)
results_ordered_DESeq2 <- results_ordered_DESeq2[order(results_ordered_DESeq2$padj), ]
results_ordered_DESeq2$genelabels = ""
results_ordered_DESeq2$genelabels[1:20] <- rownames(results_ordered_DESeq2)[1:20]
results_ordered_DESeq2["UBE3A", "genelabels"] <- "UBE3A"

# Next, we crate another column to include the LFC and statistical information for the genes.
results_ordered_DESeq2$Class <- "NS"
results_ordered_DESeq2$Class[which(results_ordered_DESeq2$log2FoldChange > 2 | results_ordered_DESeq2$log2FoldChange < -2)] <- "logFC > |2|"
results_ordered_DESeq2$Class[which(results_ordered_DESeq2$padj < 0.05)] <- "Padj < 0.05"
results_ordered_DESeq2$Class[which(results_ordered_DESeq2$padj < 0.05 & (results_ordered_DESeq2$log2FoldChange > 2 | results_ordered_DESeq2$log2FoldChange < -2))] <- "Padj < 0.05 & logFC > |2|"

# Now, we only have to plot the volcano plot.
ggplot(results_ordered_DESeq2) +
  theme_bw() +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = Class), size = 2, alpha = .4) +
  scale_colour_manual(breaks = c("NS", "logFC > |2|", "Padj < 0.05", "Padj < 0.05 & logFC > |2|"),
                      values = c("gray50", "seagreen", "royalblue", "firebrick")) + # Add a manual scale color (and order), for the elements in column "Class".
  geom_text_repel(max.overlaps = 100000,  
                  force = 10,
                  aes(x = log2FoldChange, y = -log10(padj), label = ifelse(genelabels != "", rownames(results_ordered_DESeq2),""))) + # Add gene names to most significant DEGs
  ggtitle("Volcano plot for DESeq2\nPadj < 0.05") +
  xlab("log2 fold change") + 
  ylab("-log10 Padj") +
  scale_y_continuous(trans = "log1p") +
  #scale_y_continuous(limits = c(0,700)) +
  #scale_x_continuous(limits = c(-10.5, 10.5)) +
  geom_hline(yintercept = -log10(0.05), colour = "firebrick", linetype = 4) + # Add line to intercept -log10(0.05).
  geom_vline(xintercept = c(-2, 2), colour = "firebrick", linetype = 4) + # Add lines to intercept the + and - 2 log2FC.
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + # Over-ride aesthetics for figure legend.
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white",
                                         size = 0.5, 
                                         linetype = 3, 
                                         colour ="black"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12)) + # Add several theme features.
  annotate("rect", xmin = c(-2, 2), xmax = c(-Inf, Inf), ymin = -log10(0.05), ymax = Inf, alpha = .1, fill = "mistyrose") # Add a rectangle to shade the area of DEGs.


#############
### baySeq ##
#############

# Volcano plot with gene name information and information about gene LFC and statistics
# We first create a columns with gene names that will be used to label the 20 most significant DEGs.
results_ordered_baySeq <- data.frame(DGE_results_Bayseq)
results_ordered_baySeq <- results_ordered_baySeq[order(results_ordered_baySeq$FDR.DE), ]
results_ordered_baySeq$genelabels = ""
results_ordered_baySeq$genelabels[1:20] <- rownames(results_ordered_baySeq)[1:20]
results_ordered_baySeq["UBE3A", "genelabels"] <- "UBE3A"

# Next, we crate another column to include the LFC and statistical information for the genes.
results_ordered_baySeq$Class <- "NS"
which(abs(results_ordered_baySeq$log2FC) > 2)
results_ordered_baySeq$Class[which(abs(results_ordered_baySeq$log2FC) > 2)] <- "logFC > |2|"
results_ordered_baySeq$Class[which(results_ordered_baySeq$FDR.DE < 0.05)] <- "FDR < 0.05"
results_ordered_baySeq$Class[which(results_ordered_baySeq$FDR.DE < 0.05 & (abs(results_ordered_baySeq$log2FC) > 2))] <- "FDR < 0.05 & logFC > |2|"

# Now, we only have to plot the volcano plot.
ggplot(results_ordered_baySeq) +
  theme_bw() +
  geom_point(aes(x = log2FC, y = -log10(FDR.DE), color = Class), size = 2, alpha = .4) +
  scale_colour_manual(breaks = c("NS", "logFC > |2|", "FDR < 0.05", "FDR < 0.05 & logFC > |2|"),
                      values = c("gray50", "seagreen", "royalblue", "red3")) + # Add a manual scale color (and order), for the elements in column "Class".
  geom_text_repel(max.overlaps = 100000,  
                  force = 10,
                  aes(x = log2FC, y = -log10(FDR.DE), label = ifelse(genelabels != "", rownames(results_ordered_baySeq),""))) + # Add gene names to most significant DEGs
  ggtitle("Volcano plot for baySeq\nFDR < 0.05") +
  xlab("log2 fold change") + 
  ylab("-log10 FDR") +
  scale_y_continuous(trans = "log1p") +
  #scale_y_continuous(limits = c(0,700)) +
  #scale_x_continuous(limits = c(-10.5, 10.5)) +
  geom_hline(yintercept = -log10(0.05), colour = "firebrick", linetype = 4) + # Add line to intercept -log10(0.05).
  geom_vline(xintercept = c(-2, 2), colour = "firebrick", linetype = 4) + # Add lines to intercept the + and - 2 log2FC.
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 0.8))) + # Over-ride aesthetics for figure legend.
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white",
                                         size = 0.5, 
                                         linetype = 3, 
                                         colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12)) + # Add several theme features.
  annotate("rect", xmin = c(-2, 2), xmax = c(-Inf, Inf), ymin = -log10(0.05), ymax = Inf, alpha = .1, fill = "mistyrose") # Add a rectangle to shade the area of DEGs.


#######################################################################
# Heatplots
#######################################################################
# Let's start by create a dataframe with annotations.
annotation_dataframe = as.data.frame(sample_info[ , c("phenotype", "mutation_type", "xist_status", "sex")])
colnames(annotation_dataframe) = c("Phenotype", "Mutation Type", "XIST status", "Sex")
rownames(annotation_dataframe) = rownames(sample_info)

#### edgeR
# Heatplot for the 30 most significantly differentially expressed genes.
# Heatplot of normalized readcounts
results_ordered <- DGE_results_edgeR$table[head(order(DGE_results_edgeR$table$FDR), 30), ]
significant_genes <- subset(results_ordered, FDR < 0.05)
significant_genes <- rbind(significant_genes, DGE_results_edgeR$table["UBE3A", ])
normalized_significant_genes <- countsEDASeqNormLog2[rownames(significant_genes), ]
normalized_significant_genes <- normalized_significant_genes[,c(4, 5, 6, 1, 2, 3)]

pheatmap(normalized_significant_genes,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         border_color = NA,
         scale = "row", # using scale = "row" allows a better visualization, by plotting z-scores
         # (gene by gene basis, by subtracting the mean and then dividing by STDev)
         annotation_col = annotation_dataframe,
         fontsize_col =10,
         fontsize_row = 8,
         cutree_rows = 2, 
         cutree_cols = 2,
         main = "Top 30 differentially expressed genes\nedgeR FDR")

#### limma
# Heatplot for the 30 most significantly differentially expressed genes.
# Heatplot of normalized readcounts
results_ordered <- DGE_results_limma[head(order(DGE_results_limma$adj.P.Val), 30), ]
significant_genes <- subset(results_ordered, adj.P.Val < 0.05)
significant_genes <- rbind(significant_genes, DGE_results_limma["UBE3A", ])
normalized_significant_genes <- countsEdgeRNormLog2CPM[rownames(significant_genes), ]
normalized_significant_genes <- normalized_significant_genes[,c(4, 5, 6, 1, 2, 3)]

pheatmap(normalized_significant_genes,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         border_color = NA,
         scale = "row", # using scale = "row" allows a better visualization, by plotting z-scores
         # (gene by gene basis, by subtracting the mean and then dividing by STDev)
         annotation_col = annotation_dataframe,
         fontsize_col =10,
         fontsize_row = 8,
         cutree_rows = 2, 
         cutree_cols = 2,
         main = "Top 30 differentially expressed genes\nlimma padj")


#### DESeq2
# Heatplot for the 30 most significantly differentially expressed genes.
# Heatplot of normalized readcounts
results_ordered <- DGE_results_DESeq2[head(order(DGE_results_DESeq2$padj), 30), ]
significant_genes <- subset(results_ordered, padj < 0.05)
significant_genes <- rbind(significant_genes, DGE_results_DESeq2["UBE3A", ])
normalized_significant_genes <- countsDESeq2NormLog2[rownames(significant_genes), ]
normalized_significant_genes <- normalized_significant_genes[,c(4, 5, 6, 1, 2, 3)]

pheatmap(normalized_significant_genes,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         border_color = NA,
         scale = "row", # using scale = "row" allows a better visualization, by plotting z-scores
         # (gene by gene basis, by subtracting the mean and then dividing by STDev)
         annotation_col = annotation_dataframe,
         fontsize_col =10,
         fontsize_row = 8,
         cutree_rows = 2, 
         cutree_cols = 2,
         main = "Top 30 differentially expressed genes\nDESeq2 padj")

#### baySeq
# Heatplot for the 30 most significantly differentially expressed genes.
# Heatplot of normalized readcounts
results_ordered <- DGE_results_Bayseq[head(order(DGE_results_Bayseq$FDR.DE), 30), ]
significant_genes <- subset(results_ordered, FDR.DE < 0.05)
significant_genes <- rbind(significant_genes, DGE_results_Bayseq["UBE3A", ])
normalized_significant_genes <- CD@data[significant_genes$geneNames, ]
#normalized_significant_genes <- normalized_significant_genes[,c(4, 5, 6, 1, 2, 3)]

pheatmap(normalized_significant_genes,
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         border_color = NA,
         scale = "row", # using scale = "row" allows a better visualization, by plotting z-scores
         # (gene by gene basis, by subtracting the mean and then dividing by STDev)
         annotation_col = annotation_dataframe,
         fontsize_col =10,
         fontsize_row = 8,
         cutree_rows = 2, 
         cutree_cols = 2,
         main = "Top 30 differentially expressed genes\nbaySeq FDR")
