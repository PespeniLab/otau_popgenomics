###### Code to be able to calculate population branch length statistic (PBS) from 
###### pairwise FST values of the 3 beetle populations.

library(tidyverse)
library(ggplot2)
library(ggpubr)
library("ggExtra")

setwd("~/Dropbox/1_Research/DungBeetles/OtaurusPopGenomics/fst")

#Estimating Fst between all populations

fstWAvsNC <- read.table("WAvsNC_Fst.weir.fst.txt", header=TRUE)

fstWAvsIT <- read.table("WAvsIT_Fst.weir.fst", header=TRUE)

fstNCvsIT <- read.table("NCvsIT_Fst.weir.fst", header=TRUE)

str(fstWAvsNC)
str(fstWAvsIT)
str(fstNCvsIT)

# make a fourth column of SNP by combining CHROM and POS

fstWAvsNC$SNP = paste(fstWAvsNC$CHROM,fstWAvsNC$POS)
fstWAvsIT$SNP = paste(fstWAvsIT$CHROM,fstWAvsIT$POS)
fstNCvsIT$SNP = paste(fstNCvsIT$CHROM,fstNCvsIT$POS)

# merge the dataframes based on SNP column

WAmerge <- merge(fstWAvsNC,fstWAvsIT,by.x='SNP',by.y='SNP')
str(WAmerge)
AllMerge <- merge(WAmerge,fstNCvsIT,by.x='SNP',by.y='SNP')
str(AllMerge)

# improve col names

AllMerge$FstWAvsNC <- AllMerge$WEIR_AND_COCKERHAM_FST.x
AllMerge$fstWAvsIT <- AllMerge$WEIR_AND_COCKERHAM_FST.y
AllMerge$fstNCvsIT <- AllMerge$WEIR_AND_COCKERHAM_FST

dim.data.frame(AllMerge)
head(AllMerge)

#newTable <- select(AllMerge$SNP, AllMerge$FstWAvsNC, AllMerge$fstWAvsIT, AllMerge$fstNCvsIT)

# Calculate divergence time between pops, T = -log(1 - FST) (and check with plots)

AllMerge$T_WAvsNC = -log(1 - AllMerge$FstWAvsNC)
plot(AllMerge$FstWAvsNC,AllMerge$T_WAvsNC)


AllMerge$T_WAvsIT = -log(1 - AllMerge$fstWAvsIT)
plot(AllMerge$fstWAvsIT,AllMerge$T_WAvsIT)

AllMerge$T_NCvsIT = -log(1 - AllMerge$fstNCvsIT)
plot(AllMerge$fstNCvsIT,AllMerge$T_NCvsIT)


# Then calc PBS for WA, PBS_WA = (T_WAvsNC + T_WAvsIT - T_NCvsIT)/2 ; for each pop

AllMerge$PBS_WA = (AllMerge$T_WAvsNC + AllMerge$T_WAvsIT - AllMerge$T_NCvsIT)/2
plot(AllMerge$FstWAvsNC,AllMerge$PBS_WA)

AllMerge$PBS_NC = (AllMerge$T_WAvsNC + AllMerge$T_NCvsIT - AllMerge$T_WAvsIT)/2
plot(AllMerge$FstWAvsNC,AllMerge$PBS_NC)

AllMerge$PBS_IT = (AllMerge$T_WAvsIT + AllMerge$T_NCvsIT - AllMerge$T_WAvsNC)/2
plot(AllMerge$FstWAvsNC,AllMerge$PBS_IT)

wilcox.test(PBS_WA, PBS_NC)
# Wilcoxon rank sum test with continuity correction
# 
# data:  PBS_WA and PBS_NC
# W = 1.3521e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(PBS_WA, PBS_IT)
# Wilcoxon rank sum test with continuity correction
# 
# data:  PBS_WA and PBS_IT
# W = 1.5876e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(PBS_NC, PBS_IT)
# Wilcoxon rank sum test with continuity correction
# 
# data:  PBS_NC and PBS_IT
# W = 1.7932e+10, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0


######## Pairwise scatter plots of PBS statistic, color points by Fst value for pop pair!

AllMerge$PBS_IT <- PBS_IT
AllMerge$PBS_WA <- PBS_WA
AllMerge$PBS_NC <- PBS_NC
head(AllMerge)
head(PBS_WA)



# density plot

library(reshape2)
x <- subset.data.frame(AllMerge[,c("PBS_IT","PBS_WA","PBS_NC")])
data<- melt(x)

pbs_density <- ggdensity(data, "value", fill = "variable", alpha = 0.6, xlim=c(-0.1,0.1), xlab="Population branch statistic",
                   palette = c("green","red","blue"), legend.title="", legend=c(0.8,0.75))

##################
png("multipanel_PBS_scatter4.png", height=175, width=175, units="mm", res=300)

WA_NC <- ggscatter(AllMerge, x = "PBS_WA", y = "PBS_NC", color="FstWAvsNC", xlab="AU population branch statistic", 
                   ylab="US population branch statistic",
                   legend.title="", alpha = 0.6, size=3, legend=c(0.8,0.75)) + scale_color_gradient(low="blue", high="red") +
  xlim(-0.5,2.3) + ylim(-0.5,2.3) 

IT_NC <- ggscatter(AllMerge, x = "PBS_IT", y = "PBS_NC", color="fstNCvsIT", xlab="IT population branch statistic", 
                   ylab="US population branch statistic",
                   legend.title="", alpha = 0.6, size=3, legend=c(0.8,0.75)) + scale_color_gradient(low="blue", high="red") +
  xlim(-0.5,2.3) + ylim(-0.5,2.3) 

IT_WA <- ggscatter(AllMerge, x = "PBS_IT", y = "PBS_WA", color="fstWAvsIT", xlab="IT population branch statistic", 
                   ylab="AU population branch statistic",
                   legend.title="", alpha = 0.6, size=3, legend=c(0.8,0.75)) + scale_color_gradient(low="blue", high="red") +
  xlim(-0.5,2.3) + ylim(-0.5,2.3) 

ggarrange(WA_NC, IT_NC, IT_WA, pbs_density, labels = c("A", "B", "C", "D"),nrow = 2, ncol=2)

dev.off()
