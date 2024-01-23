##################################################################
########### Data analysis for Renwick et al. 2024 Nutrients
##################################################################


##############################
############################## Section 1: Load libraries and import data
##############################

## Clear work space
rm(list=ls())
graphics.off()

## Set working directory
setwd("C:/")

## Load libraries
library(tidyverse)
library(ggplot2)
library(eoffice)
library(phyloseq)
library(ggpubr)
library(mblm)

## Read in data and sample metadata
hmo <- read.csv("Renwick et al. 2024 Nutrients HMO profiling.csv")
met <- read.csv("Renwick et al. 2024 Nutrients metadata.csv")

## Format metadata
met$Secretor_status <- gsub("1", "Secretor", met$Secretor_status, fixed=TRUE)
met$Secretor_status <- gsub("0", "Non-secretor", met$Secretor_status, fixed=TRUE)
met$Pregnancy <- gsub("1", "first", met$Pregnancy, fixed=TRUE)
met$Pregnancy <- gsub("2", "subsequent", met$Pregnancy, fixed=TRUE)
met$Pregnancy <- gsub("3", "subsequent", met$Pregnancy, fixed=TRUE)
met <- transform(met , Participant = as.character(Participant))

##############################
############################## Section 2: Calculate difference in days of lactation (first and subsequent delivery (t2 - t1)) and time in between donations for each participant 
##############################

## Difference in days of lactation
dif <- subset(met, select = c(Participant, Pregnancy, Duration_of_Lactation.days.))
dif2 <- spread(dif, Pregnancy, Duration_of_Lactation.days.)
dif2$Difference_lac <- dif2$subsequent - dif2$first

met$Difference_lac <- dif2$Difference_lac[match(met$Participant, dif2$Participant)]

## Difference in time of donation
dif <- subset(met, select = c(Participant, Pregnancy, Collection_date))
dif2 <- spread(dif, Pregnancy, Collection_date)
dif2$Difference_don <- ceiling(difftime(dif2$subsequent, dif2$first, units = "days"))
dif2$Difference_don <- gsub("days", "", dif2$Difference_don, fixed=TRUE)
dif2$Difference_don = as.numeric(as.character(dif2$Difference_don))

met$Difference_don <- dif2$Difference_don[match(met$Participant, dif2$Participant)]

##############################
############################## Section 3: Create absolute and relative abundance plots
##############################

## Format data further
hmo2 <- subset(hmo, select = -c(SUM, Fuc, Sia))
hmo2 <- gather(hmo2, HMO, Concentration, -Sample)

hmo2$HMO <- gsub("LNFP.III", "LNFP3", hmo2$HMO, fixed=TRUE)
hmo2$HMO <- gsub("LNFP.II", "LNFP2", hmo2$HMO, fixed=TRUE)
hmo2$HMO <- gsub("LNFP.I", "LNFP1", hmo2$HMO, fixed=TRUE)
hmo2$HMO <- gsub("X2FL", "2'FL", hmo2$HMO, fixed=TRUE)
hmo2$HMO <- gsub("X3FL", "3FL", hmo2$HMO, fixed=TRUE)
hmo2$HMO <- gsub("X3SL", "3'SL", hmo2$HMO, fixed=TRUE)
hmo2$HMO <- gsub("X6SL", "6'SL", hmo2$HMO, fixed=TRUE)

hmo2$Participant <- met$Participant[match(hmo2$Sample, met$Sample)]
hmo2$Participant_n <- met$Participant_n[match(hmo2$Sample, met$Sample)]
hmo2$Pregnancy <- met$Pregnancy[match(hmo2$Sample, met$Sample)]
hmo2$Secretor_status <- met$Secretor_status[match(hmo2$Sample, met$Sample)]
hmo2$Difference_lac <- met$Difference_lac[match(hmo2$Sample, met$Sample)]


colours <- c("#FFC0CB","#FF1493","#DA70D6","#D8BFD8","#800080","#7B68EE","#336699","#00BFFF","#7FFFD4","#008B8B","#006400","#6B8E23","#FFFF00", "#FF8C00","#FF4500","#CD5C5C","#DC143C","#800000","#808080")

p <- ggplot(hmo2, aes(fill = HMO, y = Concentration, x = Pregnancy)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Secretor_status + Participant_n) +
  theme_bw() +
  scale_fill_manual(values = colours) +
  ylab("Absolute Abundance (nmol/mL)") +
  scale_x_discrete(labels = hmo2$Difference_lac) +
  theme(axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) +
  labs(fill = "HMO") 

p

topptx(file = "Parity.pptx", append = TRUE, width = 12, height = 5, units = "cm")


## Plot relative abundance
q <- ggplot(hmo2, aes(fill=HMO, y=Concentration, x=Pregnancy)) +
  geom_bar(stat="identity", position = "fill") +
  facet_grid(~Secretor_status + Participant_n) +
  theme_bw() +
  scale_fill_manual(values = colours) +
  ylab("Relative Abundance") +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    label = c("0", "25%", "50%", "75%", "100%")) +
  scale_x_discrete(labels = hmo2$Difference_lac) +
  theme(axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) +
  labs(fill = "HMO") 

q

topptx(file = "Parity.pptx", append = TRUE, width = 12, height = 5, units = "cm")



##############################
############################## Section 4: Calculate Bray-Curtis dissimilarity between paired samples of each participant
##############################

hmo4 <- hmo
hmo4$Participant <- met$Participant[match(hmo4$Sample, met$Sample)]
hmo4 <- subset(hmo4, select = -c(Sample, SUM, Fuc, Sia))

## Calculate Bray-Curtis dissimilarity using a loop
Participants <- unique(hmo4$Participant)

bc <- data.frame(matrix(ncol = 2, nrow = 34))
colnames(bc) <- c("Participant", "BC")
bc$Participant <- Participants

for (a in 1:nrow(bc)) {
  hmo5 <- subset(hmo4, hmo4$Participant == bc[a,1])
  hmo5$Participant <- NULL
  bc3 <- sum(apply(hmo5, 2, function(x) abs(max(x)-min(x)))) / sum(rowSums(hmo5))
  bc[a,2] <- bc3
}


bc$Difference_lac <- met$Difference_lac[match(bc$Participant, met$Participant)]

r <- ggplot(data = bc, aes(x = Difference_lac, y = BC)) +
  geom_point(size = 7) +
  geom_smooth(col = 'black') +
  theme_bw() +
  theme(axis.title = element_text(size = (25)),
        axis.text = element_text(size = (23))) +
  labs(y = "Bray-Curtis dissimilarity", x = "Difference in time of lactation (t2 - t1) (days)") +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  geom_hline(yintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  scale_y_continuous(limits = c(0,0.6))

r

topptx(file = "Parity.pptx", append = TRUE, width = 6, height = 5, units = "cm")


## Investigating a linear relationship and correlation coefficient

bc$Difference_abs <- abs(bc$Difference_lac)

r2 <- ggplot(data = bc, aes(x = Difference_abs, y = BC)) +
  geom_point(size = 7) +
  geom_smooth(method = glm, col = 'black') +
  theme_bw() +
  theme(axis.title = element_text(size = (25)),
        axis.text = element_text(size = (23))) +
  labs(y = "Bray-Curtis dissimilarity", x = "Absolute difference in time of lactation (t2 - t1) (days)") +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  geom_hline(yintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  scale_y_continuous(limits = c(0,0.6)) +
  stat_cor(method="spearman", size = 8)

r2

topptx(file = "Parity.pptx", append = TRUE, width = 6, height = 5, units = "cm")



##############################
############################## Section 5: Calculate Bray-Curtis dissimilarity between paired samples of randomized participants within the same secretor status
##############################

hmo7 <- subset(hmo, select = -c(SUM, Fuc, Sia))
hmo7$Participant <- met$Participant[match(hmo7$Sample, met$Sample)]
hmo7$Pregnancy <- met$Pregnancy[match(hmo7$Sample, met$Sample)]
hmo7$Secretor_status <- met$Secretor_status[match(hmo7$Sample, met$Sample)]
hmo7$Difference_lac <- met$Difference_lac[match(hmo7$Sample, met$Sample)]

bc_final <- data.frame(matrix(ncol = 4))
colnames(bc_final) <- c("Participant", "BC", "Difference", "Secretor_status")

Status <- unique(hmo7$Secretor_status)

for (s in 1:length(Status)) {
  
  hmo8 <- subset(hmo7, hmo7$Secretor_status == Status[s])
  
  Participants2 <- unique(hmo8$Participant)
  
  bc2 <- data.frame(matrix(ncol = 4, nrow = (length(Participants2))))
  colnames(bc2) <- c("Participant", "BC", "Difference", "Secretor_status")
  bc2$Participant <- Participants2
  
  for (a in 1:nrow(bc2)) {
    hmo9 <- subset(hmo8, hmo8$Participant == bc2[a,1])
    
    hmo10 <- subset(hmo9, hmo9$Pregnancy == "first")
    hmo11 <- subset(hmo8, hmo8$Pregnancy == "subsequent")
    hmo11 <- hmo11[!grepl(hmo8[a,1], hmo11$Participant),]
    hmo11 <- sample_n(hmo11, 1)
    
    hmo12 <- rbind(hmo10, hmo11)
    bc2[a,3] <- hmo11$Difference
    hmo12 <- subset(hmo12, select = -c(Participant, Sample, Pregnancy, Secretor_status))
    
    bc3 <- sum(apply(hmo12, 2, function(x) abs(max(x)-min(x)))) / sum(rowSums(hmo12))
    bc2[a,2] <- bc3
    bc2$Secretor_status <- Status[s]
  }
  
  bc_final <- rbind(bc_final, bc2)
  
}

s <- ggplot(data = bc_final, aes(x = Difference, y = BC)) +
  geom_point(size = 7) +
  geom_smooth(col = 'black') +
  theme_bw() +
  theme(axis.title = element_text(size = (25)),
        axis.text = element_text(size = (23))) +
  labs(y = "Bray-Curtis dissimilarity", x = "Difference in time of lactation (t2 - t1) (days)") +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  geom_hline(yintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  scale_y_continuous(limits = c(0,0.6))

s

topptx(file = "Parity.pptx", append = TRUE, width = 6, height = 5, units = "cm")


## Investigating a linear relationship and correlation coefficient

bc_final$Difference_abs <- abs(bc_final$Difference)

s2 <- ggplot(data = bc_final, aes(x = Difference_abs, y = BC)) +
  geom_point(size = 7) +
  geom_smooth(method = glm, col = 'black') +
  theme_bw() +
  theme(axis.title = element_text(size = (25)),
        axis.text = element_text(size = (23))) +
  labs(y = "Bray-Curtis dissimilarity", x = "Absolute difference in time of lactation (t2 - t1) (days)") +
  geom_vline(xintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  geom_hline(yintercept = 0, linetype="dashed", alpha = 0.6, colour = "black") +
  scale_y_continuous(limits = c(0,0.6)) +
  stat_cor(method="spearman", size = 8)

s2

topptx(file = "Parity.pptx", append = TRUE, width = 6, height = 5, units = "cm")



##############################
############################## Section 6: Evaluating individual HMOs
##############################

hmo5 <- hmo 
hmo5 <- gather(hmo5, HMO, Concentration, -Sample)

hmo5$HMO <- gsub("LNFP.III", "LNFP3", hmo5$HMO, fixed=TRUE)
hmo5$HMO <- gsub("LNFP.II", "LNFP2", hmo5$HMO, fixed=TRUE)
hmo5$HMO <- gsub("LNFP.I", "LNFP1", hmo5$HMO, fixed=TRUE)
hmo5$HMO <- gsub("X2FL", "2'FL", hmo5$HMO, fixed=TRUE)
hmo5$HMO <- gsub("X3FL", "3FL", hmo5$HMO, fixed=TRUE)
hmo5$HMO <- gsub("X3SL", "3'SL", hmo5$HMO, fixed=TRUE)
hmo5$HMO <- gsub("X6SL", "6'SL", hmo5$HMO, fixed=TRUE)

hmo5$Participant <- met$Participant[match(hmo5$Sample, met$Sample)]
hmo5$Pregnancy <- met$Pregnancy[match(hmo5$Sample, met$Sample)]
hmo5$Secretor_status <- met$Secretor_status[match(hmo5$Sample, met$Sample)]
hmo5$Difference_lac <- met$Difference_lac[match(hmo5$Sample, met$Sample)]
hmo5$Duration_of_Lac <- met$Duration_of_Lactation.days.[match(hmo5$Sample, met$Sample)]

hmo5$Concentration2 <- hmo5$Concentration/1000

hmo5$HMO <- factor(hmo5$HMO, levels = c("2'FL", "3FL", "DFLac", "LNFP1", "LNFP2", "LNFP3", "Fuc", "3'SL", "6'SL", "LSTb", "LSTc", "DSLNT", "DSLNH", "Sia", "LNT", "LNnT", "LNH", "DFLNT", "FLNH", "DFLNH", "FDSLNH", "SUM"), ordered = TRUE)

t <- ggplot(data = hmo5, aes(x = Duration_of_Lac, y = Concentration2, group = Pregnancy)) +
  geom_smooth(aes(color = Pregnancy, fill = Pregnancy)) +
  facet_wrap(~HMO, scales = "free", ncol = 7) +
  theme_bw() +
  theme(axis.title = element_text(size = (18)),
        axis.text = element_text(size = (10)),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  labs(y = "Concentration (umol/mL)", x = "Time of lactation (days)") +
  scale_color_manual(values=c("black","orange")) +
  scale_fill_manual(values=c("black","orange"))

t

topptx(file = "Parity.pptx", append = TRUE, width = 15, height = 6, units = "cm")


## Consider percent change
ind <- subset(hmo5, select = c(Participant, Secretor_status, HMO, Difference_lac, Pregnancy, Concentration))
ind <- spread(ind, Pregnancy, Concentration)
ind$Change <- ind$subsequent - ind$first
ind$Change2 <- (ind$Change/ind$first)*100

ind$HMO <- factor(ind$HMO, levels = c("2'FL", "3FL", "DFLac", "LNFP1", "LNFP2", "LNFP3", "Fuc", "3'SL", "6'SL", "LSTb", "LSTc", "DSLNT", "DSLNH", "Sia", "LNT", "LNnT", "LNH", "DFLNT", "FLNH", "DFLNH", "FDSLNH", "SUM"), ordered = TRUE)

u <- ggplot(data = ind, aes(x = Difference_lac, y = Change2, colour = Secretor_status)) +
  geom_point(size = 3) +
  geom_smooth(col = 'black') +
  theme_bw() +
  theme(axis.title = element_text(size = (18)),
        axis.text = element_text(size = (10)),       
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12)) +
  labs(y = "Percent change in concentration [(t2 -  t1)/t1]", x = "Difference in time of lactation (t2 - t1) (days)") +
  facet_wrap(~HMO, scales = "free", ncol = 7) +
  geom_vline(xintercept = 0, linetype ="dashed", alpha = 0.6, colour = "black") +
  geom_hline(yintercept = 0, linetype ="dashed", alpha = 0.6, colour = "black")

u

topptx(file = "Parity.pptx", append = TRUE, width = 15, height = 6, units = "cm")




##############################
############################## Section 7: Create principal coordinate analysis (PCoA) based on Bray-Curtis dissimilarities of the HMO composition between pregnancies
##############################

## Add category column to metadata to discriminate samples donated within 30 days post delivery
met$Category[met$Difference_lac %in% c("-8","3","9","-25","-7","-27","-30","-19","-11")] <- "within a month"
met$Category[!(met$Difference_lac %in% c("-8","3","9","-25","-7","-27","-30","-19","-11"))] <- "longer than a month"
met$Category <- factor(met$Category, levels = c("within a month", "longer than a month"), ordered = TRUE)

## Create phyloseq object

## Prepare "otu" data
otu_mat <- subset(hmo, select = -c(SUM, Fuc, Sia))
otu_mat <- as.data.frame(t(otu_mat))
colnames(otu_mat) <- otu_mat[1,]
otu_mat <- otu_mat[-c(1),]
otu_mat[,1:ncol(otu_mat)] <- sapply(otu_mat[,1:ncol(otu_mat)],as.numeric)

## Prepare "taxonomy" data
tax_mat <- as.data.frame(rownames(otu_mat))
names(tax_mat)[names(tax_mat) == "rownames(otu_mat)"] <- "HMO"
row.names(tax_mat) <- tax_mat$HMO

## Prepare metadata 
samples_df <- met
samples_df <- transform(samples_df , Participant = as.character(Participant))
samples_df$Participant <- factor(samples_df$Participant, levels = c("Participant_1","Participant_2","Participant_3","Participant_4","Participant_5","Participant_6","Participant_7","Participant_8","Participant_9","Participant_10","Participant_11","Participant_12","Participant_13","Participant_14","Participant_15","Participant_16","Participant_17","Participant_18","Participant_19","Participant_20","Participant_21","Participant_22","Participant_23","Participant_24","Participant_25","Participant_26","Participant_27","Participant_28","Participant_29","Participant_30","Participant_31","Participant_32","Participant_33","Participant_34"))
samples_df <- tibble::column_to_rownames(samples_df, "Sample")

## Create phyloseq object
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
samples_df <- as.data.frame(samples_df)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

ps <- phyloseq(OTU, TAX, samples)

## Perform ordination
ps.ord <- ordinate(ps, "PCoA", "bray")

## Plot PCoA
colours <- c("#FFC0CB","#DA70D6","#FF00FF","#D8BFD8","#800080","#9400D3","#7B68EE","#483D8B","#4169E1","#87CEFA","#1E90FF","#00BFFF","#6495ED","#7FFFD4","#008B8B","#3CB371","#00FA9A","#90EE90","#32CD32","#006400","#6B8E23","#9ACD32","#FFFF00","#DAA520","#FFD700","#FF8C00","#FF4500","#FFA07A","#CD5C5C","#FF6347","#FF0000","#DC143C","#800000","#000000","#808080")

v <- plot_ordination(ps, ps.ord, type = "samples", color = "Participant", shape = "Category") +
  geom_point(size = 7) +
  scale_colour_manual(values = colours) +
  scale_shape_manual(values=c(16,10)) +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 13))

v

topptx(file = "Parity.pptx", append = TRUE, width = 11, height = 5, units = "cm")


##############################
############################## Section 8: Assessing the impact of factors using multivariate regression
##############################

## Transfer metadata details to bc dataframe and check normality
bc2 <- bc 

bc2$Infant_sex_combinations <- met$Infant_sex_combinations[match(bc2$Participant, met$Participant)]
bc2$Difference_don <- met$Difference_don[match(bc2$Participant, met$Participant)]
bc2$Secretor_status <- met$Secretor_status[match(bc2$Participant, met$Participant)]
bc2$Difference_abs <- abs(bc2$Difference_lac)

## Normalize variables 
shapiro.test(bc2$BC)
bc2$BC <- log(bc2$BC)
shapiro.test(bc2$BC)

shapiro.test(bc2$Difference_don)
shapiro.test(bc2$Difference_lac)

shapiro.test(bc2$Difference_abs)
bc2$Difference_abs <- log(bc2$Difference_abs)
shapiro.test(bc2$Difference_abs)

## Check distribution of HMOs
ind2 <- subset(ind, select = c(Participant, HMO, Change2))
ind2 <- spread(ind2, HMO, Change2)
#ind2 <- ind2[!(ind2$Participant %in% c("Participant_32")),]
ind3 <- ind2
ind3$Participant <- NULL
apply(ind3,2,shapiro.test)

## Positive shift data (to overcome negative values) and normalize by log transformation
ind4 <- ind3 + (min(ind3)*(-1)) + 1

ind4$`2'FL` <- log(ind4$`2'FL`)
ind4$DFLac <- log(ind4$DFLac)
ind4$LNFP1 <- log(ind4$LNFP1)
ind4$`3'SL`  <- log(ind4$`3'SL`)
ind4$`6'SL`  <- log(ind4$`6'SL`)
ind4$LSTb  <- log(ind4$LSTb)
ind4$LSTc  <- log(ind4$LSTc)
ind4$DSLNT  <- log(ind4$DSLNT)
ind4$DSLNH  <- log(ind4$DSLNH)
ind4$Sia  <- log(ind4$Sia) 
ind4$LNT  <- log(ind4$LNT)
ind4$LNnT  <- log(ind4$LNnT)
ind4$LNH  <- log(ind4$LNH)
ind4$FLNH  <- log(ind4$FLNH)

apply(ind4,2,shapiro.test)

## Combine variables into a single dataframe
ind4$Participant <- ind2$Participant
reg <- merge(bc2, ind4, by = "Participant")

## Perform multivariate regression
regression <- lm(cbind(BC,`2'FL`,`3FL`,DFLac,LNFP1,LNFP2,LNFP3,Fuc,`3'SL`,`6'SL`,LSTb,LSTc,DSLNT,DSLNH,Sia,LNT,LNnT,LNH,DFLNT,FLNH,DFLNH,SUM) ~ Secretor_status + Infant_sex_combinations + Difference_abs + Difference_don, data = reg)
summary(regression)

## Perform non-parametric regression on FDSLNH
regression2 <- mblm(FDSLNH ~ Difference_abs, data = reg)
summary(regression2)
regression3 <- mblm(FDSLNH ~ Difference_don, data = reg)
summary(regression3)

## Kruskal Wallis tests on FDSLNH
kruskal.test(FDSLNH ~ Secretor_status, data = reg)
kruskal.test(FDSLNH ~ Infant_sex_combinations, data = reg)


## Generate plots for significant comparisons
bc3 <- bc

bc3$Infant_sex_combinations <- met$Infant_sex_combinations[match(bc3$Participant, met$Participant)]
bc3$Difference_don <- met$Difference_don[match(bc3$Participant, met$Participant)]
bc3$Secretor_status <- met$Secretor_status[match(bc3$Participant, met$Participant)]
bc3$Difference_abs <- abs(bc3$Difference_lac)

plots <- merge(bc3, ind2, by = "Participant")

w <- ggplot(data = plots, aes(x = Difference_abs, y = BC)) +
  geom_smooth(method = lm, col = 'black') +
  geom_point(size = 7) +
  theme_bw() +
  theme(axis.title = element_text(size = (25)),
        axis.text = element_text(size = (23))) +
  labs(y = "Bray-Curtis dissimilarity", x = "Difference in time of lactation (t2 - t1) (days)")

w

topptx(file = "Parity.pptx", append = TRUE, width = 6, height = 5, units = "cm")


x <- ggplot(data = plots, aes(x = Difference_don, y = (FDSLNH)/1000)) +
  geom_smooth(method = lm, col = 'black') +
  geom_point(size = 7) +
  theme_bw() +
  theme(axis.title = element_text(size = (25)),
        axis.text = element_text(size = (23))) +
  labs(y = "FDSLNH (umol/mL)", x = "Time in between collections (days)")

x

topptx(file = "Parity.pptx", append = TRUE, width = 6, height = 5, units = "cm")

plots3 <- hmo2
plots3 <- subset(plots3, plots3$HMO == "LNFP3")
plots3$Infant_sex_combinations <- met$Infant_sex_combinations[match(plots3$Participant, met$Participant)]
plots3 <- plots3[!(is.na(plots3$Infant_sex_combinations)), ]

y <- ggplot(data = plots3, aes(x = Pregnancy, y = Concentration/1000, group = Participant)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_wrap(~Infant_sex_combinations, ncol = 4) +
  theme_bw() +
  theme(axis.title = element_text(size = (18)),
        axis.text = element_text(size = (18)),
        strip.text = element_text(size = 18)) +
  labs(y = "LNFP3 (umol/mL)", x = "Delivery") 

y

topptx(file = "Parity.pptx", append = TRUE, width = 9, height = 5, units = "cm")


colours <- c("#FFC0CB","#DA70D6","#FF00FF","#D8BFD8","#800080","#9400D3","#7B68EE","#483D8B","#4169E1","#87CEFA","#1E90FF","#00BFFF","#6495ED","#7FFFD4","#008B8B","#3CB371","#00FA9A","#90EE90","#32CD32","#006400","#6B8E23","#9ACD32","#FFFF00","#DAA520","#FFD700","#FF8C00","#FF4500","#FFA07A","#CD5C5C","#FF6347","#FF0000","#DC143C","#800000","#000000","#808080")

met$Participant <- factor(met$Participant, levels = c("Participant_1","Participant_2","Participant_3","Participant_4","Participant_5","Participant_6","Participant_7","Participant_8","Participant_9","Participant_10","Participant_11","Participant_12","Participant_13","Participant_14","Participant_15","Participant_16","Participant_17","Participant_18","Participant_19","Participant_20","Participant_21","Participant_22","Participant_23","Participant_24","Participant_25","Participant_26","Participant_27","Participant_28","Participant_29","Participant_30","Participant_31","Participant_32","Participant_33","Participant_34"))

z <- ggplot(met, aes(x = Duration_of_Lactation.days., y = Participant, col = Pregnancy, group = Participant)) +
  geom_point(aes(Duration_of_Lactation.days., Participant, col = Participant, shape = Pregnancy), size = 3) +
  theme_bw() + 
  scale_colour_manual(values = colours) +
  scale_shape_manual(values = c(16,1)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) +
  labs(x = "Time of Lactation (days)") +
  scale_x_continuous(breaks = c(0,100,200,300,400,500,600,700))

z

topptx(file = "Parity.pptx", append = TRUE, width = 8, height = 5, units = "cm")

