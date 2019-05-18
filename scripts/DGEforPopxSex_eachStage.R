setwd("~/Dropbox/1_Research/DungBeetles/OtaurusPopGenomics/DGE")
library("DESeq2")

library("ggplot2")

library(pheatmap)
library(RColorBrewer)

############################ SUBSETTING THE DATA ON STAGE AND TESTING FOR SEX*POP INTERACTION


# Test that the sexes are doing different things in each population, subset data for each developmental stage, import, re-estimate dispersion from each, test for interaction between sex and population separately for each stage.

countsTable <- read.delim('allcountsdataRN_noIT_AD4.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_noIT_AD4.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ sex + population + sex:population)
dds <- DESeq(dds, parallel=T)

dds <- dds[ rowSums(counts(dds)) > 1, ]   ###  Note - I may have forgotten to do this step in the above ~group interaction tests....
dim(dds)
[1] 15953    12

resultsNames(dds)
# [1] "Intercept"           "sex_M_vs_F"          "population_WA_vs_NC" "sexM.populationWA" 

res <- results(dds, name="sexM.populationWA", alpha=0.05)
res <- res[order(res$padj),]
head(res,21)
# log2 fold change (MAP): sexM.populationWA 
# Wald test p-value: sexM.populationWA 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
# OTAU006803-RA   23.50450      2.3633938 0.3684764  6.413962 1.417855e-10 1.794153e-06
# OTAU008166-RA   31.47753      2.1682044 0.3730642  5.811881 6.177488e-09 3.908497e-05
# OTAU013892-RA  485.50809      1.8338054 0.3201012  5.728830 1.011259e-08 4.265489e-05
# OTAU017172-RA   53.32943      2.1273265 0.3796359  5.603597 2.099482e-08 6.641712e-05
# OTAU003421-RA 2164.48886      0.9849516 0.1785216  5.517268 3.443094e-08 8.488777e-05
# OTAU003558-RA  247.79398      2.0073206 0.3656488  5.489750 4.025025e-08 8.488777e-05

# log2 fold change (MLE): sexM.populationWA 
# Wald test p-value: sexM.populationWA 
# DataFrame with 21 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   OTAU006803-RA    23.50450       3.915918 0.6399810  6.118803 9.428082e-10 1.205475e-05
# OTAU008166-RA    31.47753       3.712071 0.6485891  5.723302 1.044735e-08 4.452660e-05
# OTAU013892-RA   485.50809       2.385467 0.4159557  5.734905 9.756705e-09 4.452660e-05
# OTAU017172-RA    53.32943       4.191088 0.7464842  5.614436 1.972043e-08 6.303634e-05
# OTAU003421-RA  2164.48886       1.046271 0.1888935  5.538943 3.043022e-08 7.781615e-05
# ...                   ...            ...       ...       ...          ...          ...
# OTAU003656-RA  219.164112      2.3971621 0.5845103  4.101146 4.111085e-05   0.03092019
# OTAU010507-RA  799.453863      1.5665320 0.3842488  4.076869 4.564615e-05   0.03242398
# OTAU007185-RA  296.240117     -2.0854949 0.5207343 -4.004911 6.204071e-05   0.04010661
# OTAU011105-RA    9.934481      3.7340531 0.9329816  4.002279 6.273520e-05   0.04010661
# OTAU004599-RA 7352.467987      0.8550012 0.2169823  3.940418 8.133966e-05   0.04952423

summary(res)
# out of 15953 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 32, 0.2% 
# LFC < 0 (down)   : 6, 0.038% 
# outliers [1]     : 123, 0.77% 
# low counts [2]   : 3176, 20% 
# (mean count < 4.8)

# out of 15953 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 20, 0.13% 
# LFC < 0 (down)   : 1, 0.0063% 
# outliers [1]     : 74, 0.46% 
# low counts [2]   : 3093, 19% 
# (mean count < 4)

d <- plotCounts(dds, gene="OTAU007185-RA", intgroup=(c("population","sex")), returnData=TRUE)
d
p <- ggplot(d, aes(x= sex, y=count, colour = population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

d <- plotCounts(dds, gene="OTAU006803-RA", intgroup=(c("population","sex")), returnData=TRUE)
d
p <- ggplot(d, aes(x= sex, y=count, colour = population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
p

# These look like great/real interactions!!!

d <- plotCounts(dds, gene="OTAU006803-RA", intgroup=(c("population","sex")), returnData=TRUE)
d
p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100)) + scale_shape_manual(values=c(21,24))
p

p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + ylim(0,50) + scale_shape_manual(values=c(21,24))
p


ggsave("dot_plot-OTAU006803_AD4_SxPop_top.png", p, width=8, height=4, dpi=300)

d <- plotCounts(dds, gene="OTAU007185-RA", intgroup=(c("population","sex")), returnData=TRUE)
d
p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + ylim(0,750) + scale_shape_manual(values=c(21,24))
p


############### Make heatmap for AD4 genes with significant interaction

sig <- res[which(res$padj < 0.05), ]
dim(sig) # 21 significant here for interaction between sex and pop at AD4 stage


## pull out only results with padj <0.05 
sig_df<-as.data.frame(sig)
sig_df$Row.names<-rownames(sig_df)
dim(sig_df)
head(sig_df)

genesOfInterest<-c(sig_df$Row.names)
head(genesOfInterest)
length(genesOfInterest) #check; yes, 21 total

vsd <- vst(dds, blind=FALSE)
str(vsd)

m <- assay(vsd)[genesOfInterest, ]  # pulls out a matrix m of the normalized counts for each of the sig genes for all samples; can specifiy which columns/samples if wanting to plot just some
head(m)

colnames(vsd)
# [1] "NC_AD4_F1_" "NC_AD4_F2_" "NC_AD4_F3_" "NC_AD4_M1_" "NC_AD4_M2_" "NC_AD4_M3_" "WA_AD4_F1_" "WA_AD4_F2_" "WA_AD4_F3_" "WA_AD4_M1_"
# [11] "WA_AD4_M2_" "WA_AD4_M3_"


pheatmap(m, scale="row", cluster_rows=T, cluster_cols=F,show_colnames=T, col= rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), labels_col=c("NC_F1", "NC_F2", "NC_F3", "NC_M1", "NC_M2", "NC_M3", "WA_F1", "WA_F2", "WA_F3", "WA_M1","WA_M2", "WA_M3")) # leave samples in order but cluster by gene patterns (rows)


############## output to save and for enrichment analyses.

# write.csv(res, file = "DGE_sexBypopNCvsWA_AD4.csv", row.names = T, quote = F)


###################### NOW FOR L3L

countsTable <- read.delim('allcountsdataRN_noIT_L3L.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_noIT_L3L.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ sex + population + sex:population)
dds <- DESeq(dds, parallel=T)

dds <- dds[ rowSums(counts(dds)) > 1, ]   
dim(dds)
[1] 15790    12

resultsNames(dds)
# [1] "Intercept"           "sex_M_vs_F"          "population_WA_vs_NC" "sexM.populationWA" 


res <- results(dds, name="sexM.populationWA", alpha=0.05)
res <- res[order(res$padj),]
head(res)
# log2 fold change (MLE): sexM.populationWA 
# Wald test p-value: sexM.populationWA 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>   <numeric>
#   OTAU015610-RA  7.466924      -9.285562 1.7814710 -5.212301 1.865130e-07 0.002929932
# OTAU016477-RA 85.662158      -4.310282 0.9101607 -4.735738 2.182595e-06 0.017143195
# OTAU017106-RA 38.375238      -6.999490 1.7603464 -3.976200 7.002524e-05 0.366675498
# OTAU010970-RA 20.992169      -3.076200 0.7877515 -3.905039 9.421037e-05 0.369987656
# OTAU004804-RA  7.919264      -7.594537 1.9928582 -3.810877 1.384749e-04 0.419414578
# OTAU011031-RA 58.666858      -1.723069 0.4564772 -3.774710 1.601940e-04 0.419414578


summary(res)
# out of 15790 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 0, 0% 
# LFC < 0 (down)   : 2, 0.013% 
# outliers [1]     : 81, 0.51% 
# low counts [2]   : 0, 0% 
# (mean count < 0)

d <- plotCounts(dds, gene="OTAU015610-RA", intgroup=(c("population","sex")), returnData=TRUE)
d
p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100)) + scale_shape_manual(values=c(21,24))
p

p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + ylim(0,50) + scale_shape_manual(values=c(21,24))
p


ggsave("dot_plot-OTAU015610_L3L_SxPop_top.png", p, width=8, height=4, dpi=300)

############### Make heatmap for L3L genes with significant interaction

sig <- res[which(res$padj < 0.05), ]
dim(sig) # 2 significant here for interaction between sex and pop at AD4 stage


## pull out only results with padj <0.05 
sig_df<-as.data.frame(sig)
sig_df$Row.names<-rownames(sig_df)
dim(sig_df)
head(sig_df)

genesOfInterest<-c(sig_df$Row.names)
head(genesOfInterest)
length(genesOfInterest) #check; yes, 2 total

vsd <- vst(dds, blind=FALSE)
str(vsd)

m <- assay(vsd)[genesOfInterest, ]  # pulls out a matrix m of the normalized counts for each of the sig genes for all samples; can specifiy which columns/samples if wanting to plot just some
head(m)

colnames(vsd)
# [1] "NC_AD4_F1_" "NC_AD4_F2_" "NC_AD4_F3_" "NC_AD4_M1_" "NC_AD4_M2_" "NC_AD4_M3_" "WA_AD4_F1_" "WA_AD4_F2_" "WA_AD4_F3_" "WA_AD4_M1_"
# [11] "WA_AD4_M2_" "WA_AD4_M3_"

mat_scaled = t(apply(m, 1, scale))  
head(mat_scaled) 
dim(mat_scaled)

pheatmap(mat_scaled, cluster_rows=T, cluster_cols=F,show_colnames=T, labels_col=c("NC_F1", "NC_F2", "NC_F3", "NC_M1", "NC_M2", "NC_M3", "WA_F1", "WA_F2", "WA_F3", "WA_M1","WA_M2", "WA_M3")) # leave samples in order but cluster by gene patterns (rows)


pheatmap(m, scale="row", cluster_rows=T, cluster_cols=F,show_colnames=T, col= rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), labels_col=c("NC_F1", "NC_F2", "NC_F3", "NC_M1", "NC_M2", "NC_M3", "WA_F1", "WA_F2", "WA_F3", "WA_M1","WA_M2", "WA_M3")) # leave samples in order but cluster by gene patterns (rows)






###################### NOW FOR PP1

countsTable <- read.delim('allcountsdataRN_noIT_PP1.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_noIT_PP1.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ sex + population + sex:population)
dds <- dds[ rowSums(counts(dds)) > 1, ]  
dim(dds)
# [1] 16138    12

dds <- DESeq(dds, parallel=T)


resultsNames(dds)
# [1] "Intercept"           "sex_M_vs_F"          "population_WA_vs_NC" "sexM.populationWA" 

res <- results(dds, name="sexM.populationWA", alpha=0.05)
res <- res[order(res$padj),]
head(res)
# log2 fold change (MLE): sexM.populationWA 
# Wald test p-value: sexM.populationWA 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   OTAU014765-RA 362.040864      -1.413773 0.2589802 -5.459002 4.788187e-08 0.0007658226  # looks real
# OTAU000111-RA 791.785859      -1.479710 0.2998556 -4.934743 8.025633e-07 0.0032090493  # looks real
# OTAU009178-RA 163.879301      -4.204571 0.8389846 -5.011500 5.400742e-07 0.0032090493  # looks real
# OTAU009605-RA   8.393102      30.000000 6.0128769  4.989292 6.060093e-07 0.0032090493  # looks like crap
# OTAU014528-RA  71.009842      -8.485299 1.8232606 -4.653914 3.256922e-06 0.0104182411  # looks a little strange (2 outliers in NC males)
# OTAU005791-RA  41.763200       2.820404 0.6518636  4.326678 1.513752e-05 0.0403515861  # looks real

summary(res)

# out of 16138 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 2, 0.012% 
# LFC < 0 (down)   : 4, 0.025% 
# outliers [1]     : 144, 0.89% 
# low counts [2]   : 0, 0% 
# (mean count < 0)

d <- plotCounts(dds, gene="OTAU014528-RA", intgroup=(c("population","sex")), returnData=TRUE)
d
p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100)) + scale_shape_manual(values=c(21,24))
p

p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + ylim(0,500) + scale_shape_manual(values=c(21,24))
p


ggsave("dot_plot-OTAU014765_PP1_SxPop_top.png", p, width=8, height=4, dpi=300)

############### Make heatmap for genes with significant interaction

sig <- res[which(res$padj < 0.05), ]
dim(sig) # 6 significant here for interaction between sex and pop at AD4 stage


## pull out only results with padj <0.05 
sig_df<-as.data.frame(sig)
sig_df$Row.names<-rownames(sig_df)
dim(sig_df)
head(sig_df)

c("OTAU014765-RA","OTAU000111-RA", "OTAU009178-RA", "OTAU005791-RA")

genesOfInterest<-c(sig_df$Row.names)
head(genesOfInterest)
length(genesOfInterest) #check; yes, 6 total

vsd <- vst(dds, blind=FALSE)
str(vsd)

m <- assay(vsd)[c("OTAU014765-RA","OTAU000111-RA", "OTAU009178-RA", "OTAU005791-RA"), ]  # pulls out a matrix m of the normalized counts for each of the sig genes for all samples; can specifiy which columns/samples if wanting to plot just some
head(m)


pheatmap(m, scale="row", cluster_rows=T, cluster_cols=F,show_colnames=T, col= rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), labels_col=c("NC_F1", "NC_F2", "NC_F3", "NC_M1", "NC_M2", "NC_M3", "WA_F1", "WA_F2", "WA_F3", "WA_M1","WA_M2", "WA_M3")) # leave samples in order but cluster by gene patterns (rows)





###################### NOW FOR PD1

countsTable <- read.delim('allcountsdataRN_noIT_PD1.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_noIT_PD1.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ sex + population + sex:population)
dds <- dds[ rowSums(counts(dds)) > 1, ]  
dim(dds)
#[1] 16154    12


dds <- DESeq(dds, parallel=T)


resultsNames(dds)
# [1] "Intercept"           "sex_M_vs_F"          "population_WA_vs_NC" "sexM.populationWA" 

res <- results(dds, name="sexM.populationWA", alpha=0.05)
res <- res[order(res$padj),]
head(res)
# log2 fold change (MLE): sexM.populationWA 
# Wald test p-value: sexM.populationWA 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue       padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>  <numeric>
#   OTAU012633-RA  54.96144      -6.418280 1.2932156 -4.963039 6.939860e-07 0.01108365
# OTAU001840-RA 254.16288       2.174578 0.4875971  4.459785 8.204192e-06 0.06551458
# OTAU007383-RA 320.91363       1.785536 0.4285082  4.166866 3.088159e-05 0.16440328
# OTAU000450-RA 111.18757       4.045255 1.0277097  3.936185 8.278724e-05 0.16527438
# OTAU002935-RA 146.63567      -1.322743 0.3358953 -3.937962 8.217651e-05 0.16527438
# OTAU003937-RA  22.22253      -5.592099 1.4112999 -3.962374 7.420803e-05 0.16527438


summary(res)

# out of 16154 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 0, 0% 
# LFC < 0 (down)   : 1, 0.0062% 
# outliers [1]     : 183, 1.1% 
# low counts [2]   : 0, 0% 
# (mean count < 0)


d <- plotCounts(dds, gene="OTAU012633-RA", intgroup=(c("population","sex")), returnData=TRUE)
d
p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100)) + scale_shape_manual(values=c(21,24))
p

p <- ggplot(d, aes(x= population, y=count, shape = sex, colour = population, fill=population)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + ylim(0,350) + scale_shape_manual(values=c(21,24))
p


ggsave("dot_plot-OTAU012633_PD1_SxPop_top.png", p, width=8, height=4, dpi=300)

############### Make heatmap for genes with significant interaction

sig <- res[which(res$padj < 0.05), ]
dim(sig) # 1 significant here for interaction between sex and pop at AD4 stage


## pull out only results with padj <0.05 
sig_df<-as.data.frame(sig)
sig_df$Row.names<-rownames(sig_df)
dim(sig_df)
head(sig_df)


genesOfInterest<-c(sig_df$Row.names)
head(genesOfInterest)
length(genesOfInterest) #check; yes, 2 total

vsd <- vst(dds, blind=FALSE)
str(vsd)

m <- assay(vsd)[genesOfInterest, ]  # pulls out a matrix m of the normalized counts for each of the sig genes for all samples; can specifiy which columns/samples if wanting to plot just some
dim(m)



genesOfInterest<-c(sig_df$Row.names,"OTAU012633-RA") # add a gene so there are two to be able to make the heatmap
head(genesOfInterest)
length(genesOfInterest) #check; yes, 1 total

vsd <- vst(dds, blind=FALSE)

m <- assay(vsd)[genesOfInterest, ]  # pulls out a matrix m of the normalized counts for each of the sig genes for all samples; can specifiy which columns/samples if wanting to plot just some
head(m)

pheatmap(m, scale="row", cluster_rows=T, cluster_cols=F,show_colnames=T, col= rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), labels_col=c("NC_F1", "NC_F2", "NC_F3", "NC_M1", "NC_M2", "NC_M3", "WA_F1", "WA_F2", "WA_F3", "WA_M1","WA_M2", "WA_M3")) # leave samples in order but cluster by gene patterns (rows)


logd<-log(d$count,10)
logd
logd<-as.data.frame(logd)
logd<-t(logd)

pheatmap(logd, cluster_cols=F,show_colnames=T, col= rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), labels_col=c("NC_F1", "NC_F2", "NC_F3", "NC_M1", "NC_M2", "NC_M3", "WA_F1", "WA_F2", "WA_F3", "WA_M1","WA_M2", "WA_M3")) # leave samples in order but cluster by gene patterns (rows)


### Try merging significant genes from each stage
