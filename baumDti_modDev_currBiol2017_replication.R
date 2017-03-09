### This script will run group-level analyses using structural network measures derived in Baum et al., Current Biology 2017 
### Participation coefficient, Global efficiency, and Edge betweenness centrality were calculated using functions from the Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/
### Modularity index was calculated using the genlouvain.m function described here: http://netwiki.amath.unc.edu/GenLouvain/GenLouvain

###############################
### Load Relevant Libraries ###
###############################
require(ggplot2)
require(Hmisc)
require(MASS)
require(mgcv)
require(ppcor)
require(stringr)
require(visreg)
require(parallel)
require(multilevel)
require(stats)
require(Formula)
require(lavaan)


## Set working directory 
setwd("/data/joy/BBL/projects/pncBaumDti/Modular_Development_paper/CurrBiol_replication/Network_measures")

## Load in csv with demographics and network measures
LTN_n882_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Modular_Development_paper/CurrBiol_replication/demog/LTN_n882_demog_netMeasures_20170306.csv")

## Define appropriate variables as factors/numerical values
LTN_n882_df$Sex <- as.factor(LTN_n882_df$Sex)
LTN_n882_df$Race <- as.factor(LTN_n882_df$Race)
LTN_n882_df$F1_Complex_Reasoning_Efficiency <- as.numeric(as.character(LTN_n882_df$F1_Complex_Reasoning_Efficiency))
LTN_n882_df$F2_Memory_Efficiency <- as.numeric(as.character(LTN_n882_df$F2_Memory_Efficiency))
LTN_n882_df$F3_Executive_Efficiency <- as.numeric(as.character(LTN_n882_df$F3_Executive_Efficiency))
LTN_n882_df$F4_Social_Cognition_Efficiency <- as.numeric(as.character(LTN_n882_df$F4_Social_Cognition_Efficiency))

## Define quadratic and cubic age terms (demeaned)
LTN_n882_df$ageSq_demeaned<-(LTN_n882_df$age_in_yrs-mean(LTN_n882_df$age_in_yrs))^2
LTN_n882_df$ageCub_demeaned <-(LTN_n882_df$age_in_yrs-mean(LTN_n882_df$age_in_yrs))^3

###############################################################################################
### Define cognitive measures dataframe after removing subjects with missing cognitive data ###
###############################################################################################

Cog_n880_df <- LTN_n882_df[!is.na(LTN_n882_df$F3_Executive_Efficiency),]

################
### FIGURE 1 ###
################
hist(LTN_n882_df$age_in_yrs,col="gray")
ExecEff_Age_gam <- gam(F3_Executive_Efficiency ~ s(age_in_yrs,k=4) + Sex, fx=TRUE, data = Cog_n880_df)
visreg(ExecEff_Age_gam,"age_in_yrs",xlab="Age (years)", ylab="Executive Efficiency (Z)")

#######################################################################
### FIGURE 3: Age Effects on Yeo 7-system Participation Coefficient ###
#######################################################################

## Figure 3a
mean_Yeo7system_PC_gam <- gam(FA_Yeo7system_wholebrain_partCoeff_scale125 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
visreg(mean_Yeo7system_PC_gam,"age_in_yrs", xlab="Age (years)", ylab="Mean Participation Coefficient", ylim=c(0.76,0.79))

## Calculate Partial r coefficient for Age effect
covs <- cbind(LTN_n882_df$ageSq_demeaned,LTN_n882_df$Sex, LTN_n882_df$meanRELrms, LTN_n882_df$FA_Total_Network_Strength_scale125)
pcor.test(LTN_n882_df$FA_Yeo7system_wholebrain_partCoeff_scale125, LTN_n882_df$age_in_yrs, covs)

## Run GAMs for Module-specific participation coefficient measures (Figure 3B)
Visual_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_Visual ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
Visual_partCoeff_Age_pval <- summary(Visual_partCoeff_gam)$s.table[1,4]

Somatomotor_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_Somatomotor ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
Somatomotor_partCoeff_Age_pval <- summary(Somatomotor_partCoeff_gam)$s.table[1,4]

DorsalAtt_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_DorsalAttention ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
DorsalAtt_partCoeff_Age_pval <- summary(DorsalAtt_partCoeff_gam)$s.table[1,4]

VentralAtt_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_VentralAttention ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
VentralAtt_partCoeff_Age_pval <- summary(VentralAtt_partCoeff_gam)$s.table[1,4]

Limbic_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_Limbic ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
Limbic_partCoeff_Age_pval <- summary(Limbic_partCoeff_gam)$s.table[1,4]

Frontoparietal_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_Frontoparietal ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
Frontoparietal_partCoeff_Age_pval <- summary(Frontoparietal_partCoeff_gam)$s.table[1,4]

Default_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_Default ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
Default_partCoeff_Age_pval <- summary(Default_partCoeff_gam)$s.table[1,4]

Subcortical_partCoeff_gam <- gam(Avg_FA_partCoeff_Yeo7system_Subcortical ~ FA_Total_Network_Strength_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex, fx=TRUE, data = LTN_n882_df)
Subcortical_partCoeff_Age_pval <- summary(Subcortical_partCoeff_gam)$s.table[1,4]

Module_specific_Age_pvals <- c(Visual_partCoeff_Age_pval, Somatomotor_partCoeff_Age_pval, DorsalAtt_partCoeff_Age_pval, VentralAtt_partCoeff_Age_pval, Limbic_partCoeff_Age_pval, Frontoparietal_partCoeff_Age_pval, Default_partCoeff_Age_pval, Subcortical_partCoeff_Age_pval)

## Calculate Z-scores for spline age pvalues (for each module) -- used in 
Module_specific_Age_Zscores <- qnorm(Module_specific_Age_pvals, lower.tail=FALSE)
Module_specific_Age_Zscores <- as.data.frame(Module_specific_Age_Zscores)
Module_specific_Age_Zscores$System_idx <- c(1:8)
# Make dummy variable according to strength of Age effect (1=strongest,8=weakest)
Module_specific_Age_Zscores$Effect <- base::rank(-Module_specific_Age_Zscores$Module_specific_Age_Zscores)
## Fig3B Barplot
Fig3B_ModSeg_Yeo7system_barplot <- ggplot(data=Module_specific_Age_Zscores, aes(Effect, Module_specific_Age_Zscores)) + geom_bar(fill=c("dark blue","firebrick4","skyblue","orange","forestgreen","khaki1", "darkgray","purple"), col="black",stat="identity",position=position_dodge())
Fig3B_ModSeg_Yeo7system_barplot + scale_x_continuous(breaks=1:8, labels=c("Somatomotor","Default","Ventral Attention","Frontoparietal","DorsalAttention","Limbic","Subcortical","Visual"))+ theme(axis.text = element_text(size= 14)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border= element_rect(colour="black",fill=FALSE, size=1))
    
## Figure 3c -- Read in regional PC values for each subject
Regional_FA_Yeo7system_participationCoefficient_scale125 <- read.csv("full_Participation_coeff_Yeo7system_FA_detWB_1mill_10_400mm_end2end_n882.txt",header=FALSE)

# Merge with dataframe containing demographics and covariates
full_FA_Yeo7system_participationCoefficient_scale125 <- cbind(Regional_FA_Yeo7system_participationCoefficient_scale125,LTN_n882_df)

# Set categorical variables as factors
full_FA_Yeo7system_participationCoefficient_scale125$Sex <- as.factor(full_FA_Yeo7system_participationCoefficient_scale125$Sex)

###############################################
### Run node-wise GAM estimating Age Effect ###
###############################################

covariates=" ~  s(age_in_yrs,k=4) + Sex + meanRELrms + FA_Total_Network_Strength_scale125"    
m <- mclapply(names(full_FA_Yeo7system_participationCoefficient_scale125[,1:234]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=2)
NodeWise_YeoPC_GAM_Age_pvals <- mclapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=full_FA_Yeo7system_participationCoefficient_scale125, REML=T))$s.table[1,4]},mc.cores=2)
NodeWise_YeoPC_GAM_Age_pvals <- as.data.frame(NodeWise_YeoPC_GAM_Age_pvals)
NodeWise_YeoPC_GAM_Age_pvals <- t(NodeWise_YeoPC_GAM_Age_pvals)
NodeWise_YeoPC_GAM_Age_pvals <- as.data.frame(NodeWise_YeoPC_GAM_Age_pvals)
colnames(NodeWise_YeoPC_GAM_Age_pvals) <- "NodeWise_YeoPC_GAM_Age_pvals"
NodeWise_YeoPC_GAM_Age_pvals$Node_index <- 1:234

## Apply linear model to get t statistic --> use to determine direction of GAM effects 
NodeWise_YeoPC_GAM_Age_pvals$lm_Tvals <- lapply(full_FA_Yeo7system_participationCoefficient_scale125[,1:234],function(x) summary(lm(x ~ full_FA_Yeo7system_participationCoefficient_scale125$age_in_yrs + full_FA_Yeo7system_participationCoefficient_scale125$ageSq_demeaned + full_FA_Yeo7system_participationCoefficient_scale125$meanRELrms  + full_FA_Yeo7system_participationCoefficient_scale125$Sex + full_FA_Yeo7system_participationCoefficient_scale125$FA_Total_Network_Strength_scale125))$coefficients[2,3])

## Calculate Z-scores for spline age p-values
NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore <- 0
NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore <- qnorm(NodeWise_YeoPC_GAM_Age_pvals$NodeWise_YeoPC_GAM_Age_pvals,lower.tail=FALSE)

## Set Z-score sign to positive/negative based on T-value ###
NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore[which(NodeWise_YeoPC_GAM_Age_pvals$lm_Tvals < 0)] <- -(NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore[which(NodeWise_YeoPC_GAM_Age_pvals$lm_Tvals < 0)])

# NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore <- abs(NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore)
# NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore[which(NodeWise_YeoPC_GAM_Age_pvals$lm_Tvals < 0)] <- -(NodeWise_YeoPC_GAM_Age_pvals$AgeEffect_Zscore)
## FDR correction
FDRcorr_NodeWise_YeoPC_GAM_Age_pvals <- p.adjust(NodeWise_YeoPC_GAM_Age_pvals$NodeWise_YeoPC_GAM_Age_pvals,method="fdr")
FDRcorr_NodeWise_YeoPC_GAM_Age_pvals  <- cbind(NodeWise_YeoPC_GAM_Age_pvals,FDRcorr_NodeWise_YeoPC_GAM_Age_pvals)

## Define subset of regions showing significant Age effects on Yeo PC
sig_FDRcorr_NodeWise_YeoPC_GAM_Age_pvals <- FDRcorr_NodeWise_YeoPC_GAM_Age_pvals[which(FDRcorr_NodeWise_YeoPC_GAM_Age_pvals$FDRcorr_NodeWise_YeoPC_GAM_Age_pvals <.05),]

#################################################
### Edgewise GAM: Age effect on Edge Strength ###
#################################################

## Read in edge vectorization file (contains all unique edges for all subjects)
FA_EdgeStrength <- read.csv("new_EdgeVec_FA_detWB_1mill_10_400mm_end2end_scale125_n882.txt",header = FALSE)
FA_EdgeStrength <- as.data.frame(FA_EdgeStrength)

# Merge demographics and edge strength data
full_FA_EdgeStrength  <- cbind(FA_EdgeStrength, LTN_n882_df)
full_FA_EdgeStrength$Sex <- as.factor(full_FA_EdgeStrength$Sex)

## Read in Within/Between module index 
Yeo_7system_withinBetween_index <- read.table("Yeo_7system_Lausanne234_withinBetween_index.txt",header=FALSE)
Yeo_7system_withinBetween_index <- as.data.frame(Yeo_7system_withinBetween_index)
colnames(Yeo_7system_withinBetween_index) <- "Yeo_7system_withinBetween_index"

## Mean Normalized Edge Betweenness Centrality (for defining Hub Edges)
mean_subj_Normalized_edgeBetweenness <- read.csv("mean_Normalized_0to1_edge_betweenness_n882.txt",header=FALSE)
mean_subj_Normalized_edgeBetweenness  <- as.data.frame(mean_subj_Normalized_edgeBetweenness)
colnames(mean_subj_Normalized_edgeBetweenness) <- "mean_subj_Normalized_edgeBetweenness"

## Apply GAM across each unique edge -- extract spline Age effect
covariates= "~ s(age_in_yrs,k=4) + Sex + meanRELrms + FA_Total_Network_Strength_scale125"    
m <- mclapply(names(full_FA_EdgeStrength[,1:27261]), function(x) {as.formula(paste(x, covariates, sep=""))},mc.cores=4)
EdgeStrength_GAM_Age_pvals <- lapply(m, function(x) { summary(gam(formula = x, fx=TRUE,data=full_FA_EdgeStrength, REML=T))$s.table[1,4]})
EdgeStrength_GAM_Age_pvals <- as.data.frame(EdgeStrength_GAM_Age_pvals)
EdgeStrength_GAM_Age_pvals <- t(EdgeStrength_GAM_Age_pvals)
EdgeStrength_GAM_Age_pvals<- as.data.frame(EdgeStrength_GAM_Age_pvals)
colnames(EdgeStrength_GAM_Age_pvals) <- "EdgeStrength_GAM_Age_pvals"
EdgeStrength_GAM_Age_pvals$Edge_index <- 1:27261

## Apply linear model to get t statistic --> use to determine direction of GAM effects 
EdgeStrength_GAM_Age_pvals$lm_Tvals <- lapply(full_FA_EdgeStrength[,1:27261],function(x) summary(lm(x ~ full_FA_EdgeStrength$age_in_yrs + full_FA_EdgeStrength$ageSq_demeaned + full_FA_EdgeStrength$meanRELrms  + full_FA_EdgeStrength$Sex + full_FA_EdgeStrength$FA_Total_Network_Strength_scale125))$coefficients[2,3])

## Calculate Z-scores for spline age p-values
EdgeStrength_GAM_Age_pvals$Age_Zscore <- 0
EdgeStrength_GAM_Age_pvals$Age_Zscore <- qnorm(EdgeStrength_GAM_Age_pvals$EdgeStrength_GAM_Age_pvals, lower.tail=FALSE)

## Merge results with Edge betweenness centrality data and within-module index
EdgeStrength_GAM_Age_pvals <- cbind(EdgeStrength_GAM_Age_pvals,Yeo_7system_withinBetween_index, mean_subj_Normalized_edgeBetweenness)

## FDR Correction
FDRcorr_EdgeStrength_GAM_Age_pvals <- p.adjust(EdgeStrength_GAM_Age_pvals$EdgeStrength_GAM_Age_pvals,method="fdr")
FDRcorr_EdgeStrength_GAM_Age_pvals <- cbind(EdgeStrength_GAM_Age_pvals,FDRcorr_EdgeStrength_GAM_Age_pvals)

## Bonferonni correction
BonfCorr_EdgeStrength_GAM_Age_pvals <- p.adjust(EdgeStrength_GAM_Age_pvals$EdgeStrength_GAM_Age_pvals,method="bonferroni")

## Merge vectors of FDR and Bonferonni-corrected pvals
FDRcorr_EdgeStrength_GAM_Age_pvals <- cbind(FDRcorr_EdgeStrength_GAM_Age_pvals, BonfCorr_EdgeStrength_GAM_Age_pvals)

## Identify edges with significant age effects
sig_FDR_EdgeStrength_Age_pvals <- FDRcorr_EdgeStrength_GAM_Age_pvals[which(FDRcorr_EdgeStrength_GAM_Age_pvals$FDRcorr_EdgeStrength_GAM_Age_pvals <.05),]
dim(sig_FDR_EdgeStrength_Age_pvals)

No_Effect_EdgeStrength_Age_pvals <- FDRcorr_EdgeStrength_GAM_Age_pvals[which(FDRcorr_EdgeStrength_GAM_Age_pvals$FDRcorr_EdgeStrength_GAM_Age_pvals >= 0.05),]
dim(No_Effect_EdgeStrength_Age_pvals)

Positive_sig_FDR_EdgeStrength_Age_pvals <- sig_FDR_EdgeStrength_Age_pvals[which(sig_FDR_EdgeStrength_Age_pvals$lm_Tvals > 0),]
dim(Positive_sig_FDR_EdgeStrength_Age_pvals)

Negative_sig_FDR_EdgeStrength_Age_pvals <- sig_FDR_EdgeStrength_Age_pvals[which(sig_FDR_EdgeStrength_Age_pvals$lm_Tvals < 0),]
dim(Negative_sig_FDR_EdgeStrength_Age_pvals)

## Define index for edges showing significant positive Age effects
FDRcorr_EdgeStrength_GAM_Age_pvals$FDR_sigPos_Index <- 0
FDRcorr_EdgeStrength_GAM_Age_pvals$FDR_sigPos_Index[which(FDRcorr_EdgeStrength_GAM_Age_pvals$FDRcorr_EdgeStrength_GAM_Age_pvals < 0.05 & FDRcorr_EdgeStrength_GAM_Age_pvals$lm_Tvals > 0)] <- 1

## Define index for edges showing significant positive Age effects
FDRcorr_EdgeStrength_GAM_Age_pvals$FDR_sigNeg_Index <- 0
FDRcorr_EdgeStrength_GAM_Age_pvals$FDR_sigNeg_Index[which(FDRcorr_EdgeStrength_GAM_Age_pvals$FDRcorr_EdgeStrength_GAM_Age_pvals < 0.05 & FDRcorr_EdgeStrength_GAM_Age_pvals$lm_Tvals < 0)] <- 1

## Bonferroni sigPos Index (for BrainNet Renderings)
FDRcorr_EdgeStrength_GAM_Age_pvals$Bonf_sigPos_Index <- 0
FDRcorr_EdgeStrength_GAM_Age_pvals$Bonf_sigPos_Index[which(FDRcorr_EdgeStrength_GAM_Age_pvals$BonfCorr_EdgeStrength_GAM_Age_pvals < 0.05 & FDRcorr_EdgeStrength_GAM_Age_pvals$lm_Tvals > 0)] <- 1

########################################################################################
### Define Hub Edges based on top quartile of normalized Edge Betweenness Centrality ###
########################################################################################
hubThreshold <- quantile(FDRcorr_EdgeStrength_GAM_Age_pvals$mean_subj_Normalized_edgeBetweenness,.75)
FDRcorr_EdgeStrength_GAM_Age_pvals$Normalized_Hub_index <- 0
FDRcorr_EdgeStrength_GAM_Age_pvals$Normalized_Hub_index[which(FDRcorr_EdgeStrength_GAM_Age_pvals$mean_subj_Normalized_edgeBetweenness >= hubThreshold)] <- 1

## All Bonferonni-corrected hub edges that strengthen with age
Bonf_sigPos_Norm_Hub_Edges_df <- as.data.frame(FDRcorr_EdgeStrength_GAM_Age_pvals[which(FDRcorr_EdgeStrength_GAM_Age_pvals$Normalized_Hub_index==1 & FDRcorr_EdgeStrength_GAM_Age_pvals$Bonf_sigPos_Index==1),]) 
dim(Bonf_sigPos_Norm_Hub_Edges_df)

## All FDR-corrected hub edges that strengthen with age
sigPos_Norm_Hub_Edges_df <- as.data.frame(FDRcorr_EdgeStrength_GAM_Age_pvals[which(FDRcorr_EdgeStrength_GAM_Age_pvals$Normalized_Hub_index==1 & FDRcorr_EdgeStrength_GAM_Age_pvals$FDR_sigPos_Index==1),]) 
dim(sigPos_Norm_Hub_Edges_df)

## All FDR-corrected Within-module hub edges that strengthen with age
sigPos_Norm_within_Hub_Edges_df <- as.data.frame(FDRcorr_EdgeStrength_GAM_Age_pvals[which(FDRcorr_EdgeStrength_GAM_Age_pvals$Yeo_7system_withinBetween_index==0 & FDRcorr_EdgeStrength_GAM_Age_pvals$Normalized_Hub_index==1 & FDRcorr_EdgeStrength_GAM_Age_pvals$FDR_sigPos_Index==1),])
dim(sigPos_Norm_within_Hub_Edges_df)

## All FDR-corrected Between-module hub edges that strengthen with age
sigPos_Norm_between_Hub_Edges_df <- as.data.frame(FDRcorr_EdgeStrength_GAM_Age_pvals[which(FDRcorr_EdgeStrength_GAM_Age_pvals$Yeo_7system_withinBetween_index==1 & FDRcorr_EdgeStrength_GAM_Age_pvals$Normalized_Hub_index==1 & FDRcorr_EdgeStrength_GAM_Age_pvals$FDR_sigPos_Index==1),])
dim(sigPos_Norm_between_Hub_Edges_df)

####################################################################
### FIGURE 4: Avg. Within-Module and Between-Module Connectivity ###
####################################################################

## Figure 4A
Avg_Yeo_WithinMod_gam <- gam(FA_Yeo7system_Avg_Within_Module_connectivity_strength ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
summary(Avg_Yeo_WithinMod_gam)$s.table[1,4]
visreg(Avg_Yeo_WithinMod_gam,"age_in_yrs")

## Figure 4B
Avg_Yeo_BetweenMod_gam <- gam(FA_Yeo7system_Avg_Between_Module_connectivity_strength ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
summary(Avg_Yeo_BetweenMod_gam)$s.table[1,4]
visreg(Avg_Yeo_BetweenMod_gam,"age_in_yrs")

## Figure 4D

# Subset non-zero edges
nzEdges_df <- FDRcorr_EdgeStrength_GAM_Age_pvals[!is.na(FDRcorr_EdgeStrength_GAM_Age_pvals$EdgeStrength_GAM_Age_pvals),]

## Within- and Between-Module edges that Strengthen with Age
sigPos_within <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==0 & nzEdges_df$FDR_sigPos_Index==1),]))
tot_within <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==0),]))
perc_strengthen_within <- as.numeric(sigPos_within / tot_within) * 100)

sigPos_between <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==1 & nzEdges_df$FDR_sigPos_Index==1),]))
tot_between <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==1),]))
perc_strengthen_between <- as.numeric((sigPos_between / tot_between) * 100)

## Within- and Between-Module edges that Weaken with Age
sigNeg_within <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==0 & nzEdges_df$FDR_sigNeg_Index==1),]))
tot_within <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==0),]))
perc_weaken_within <- as.numeric((sigNeg_within / tot_within) * 100)

sigNeg_between <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==1 & nzEdges_df$FDR_sigNeg_Index==1),]))
tot_between <- as.numeric(nrow(nzEdges_df[which(nzEdges_df$Yeo_7system_withinBetween_index==1),]))
perc_weaken_between <- as.numeric((sigNeg_between / tot_between) * 100)

## Create dataframe with percentages of edges that significantly strengthen or weaken with age
ageEffect <- c("Strengthen", "Strengthen", "Weaken", "Weaken")
withinBetween <- c("Within", "Between", "Within", "Between")
percent_edges <- c(perc_strengthen_within, perc_strengthen_between, perc_weaken_within, perc_weaken_between)
Fig4D_df <- as.data.frame(cbind(ageEffect,withinBetween,percent_edges))

## Reorder factors for plotting
Fig4D_df$withinBetween  <- factor(Fig4D_df$withinBetween, levels = rev(levels(Fig4D_df$withinBetween)))
Fig4D_df$percent_edges <- as.numeric(as.character(Fig4D_df$percent_edges))

## Create barplot for Figure 4D
percSigEdges_Yeo7system_plot <- ggplot(data=Fig4D_df, aes(ageEffect,as.numeric(percent_edges),fill=withinBetween)) + geom_bar(stat="identity",position=position_dodge())
percSigEdges_Yeo7system_plot + theme(axis.text = element_text(size= 14)) + scale_fill_brewer(palette = "Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border= element_rect(colour="black",fill=FALSE, size=1))

#############################################
### FIGURE 5: Methodological Replications ###
#############################################

## Figure 5a
YeoQ_scale125_gam <- gam(FA_Yeo7system_subject_Q_scale125 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
summary(YeoQ_scale125_gam)$s.table[1,4]
visreg(YeoQ_scale125_gam,"age_in_yrs",xlab="Age (years)",ylab="Subject-specific Q Index")

## Figure 5b	
StructPart_gamma2.5_wholebrain_PC_gam <- gam(Wholebrain_Avg_Participation_coeff_FA_n882_StructPart_gamma_2.5_scale125 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
summary(StructPart_gamma2.5_wholebrain_PC_gam)$s.table[1,4]
visreg(StructPart_gamma2.5_wholebrain_PC_gam,"age_in_yrs", xlab="Age (years)",ylab="Mean PC - Group StructPart")

## Figure 5c	
StructPart_Rastko_subjectQ_gamma2.5_gam <- gam(FA_SubjectLevel_Q_RastkoConsensus_gamma_2.5_scale125 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
summary(StructPart_Rastko_subjectQ_gamma2.5_gam)$s.table[1,4]
visreg(StructPart_Rastko_subjectQ_gamma2.5_gam,"age_in_yrs",xlab="Age (years)",ylab="Group StructPart Q Index")

## Figure 5d
FA_Yeo7system_meanPC_scale250_gam <- gam(FA_Yeo7system_wholebrain_partCoeff_scale250 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale250,fx=TRUE, data = LTN_n882_df)
summary(FA_Yeo7system_meanPC_scale250_gam)$s.table[1,4]
visreg(FA_Yeo7system_meanPC_scale250_gam,"age_in_yrs", xlab="Age (years)",ylab="Mean PC - Yeo 463-region")

## Figure 5e
rawSC_Yeo7system_meanPC_scale125_gam <- gam(rawSC_Yeo7system_wholebrain_partCoeff_scale125 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + rawSC_Total_NetworkStrength_scale125,fx=TRUE, data = LTN_n882_df)
summary(rawSC_Yeo7system_meanPC_scale125_gam)$s.table[1,4]
visreg(rawSC_Yeo7system_meanPC_scale125_gam,"age_in_yrs",xlab="Age (years)",ylab="Mean PC - Streamline Count")

## Figure 5f
volNormSC_Yeo7system_meanPC_scale125_gam <- gam(volNormSC_Yeo7system_wholebrain_partCoeff_scale125 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + volNormSC_Avg_NodeStrength_scale125,fx=TRUE, data = LTN_n882_df)
summary(volNormSC_Yeo7system_meanPC_scale125_gam)$s.table[1,4]
visreg(volNormSC_Yeo7system_meanPC_scale125_gam,"age_in_yrs",xlab="Age (years)",ylab="Mean PC - Streamline Density")

##################################################
### Probabilistic Replication - Figure 5 (G-I) ###
##################################################
PTx_n878_df <- read.csv("/data/joy/BBL/projects/pncBaumDti/Modular_Development_paper/CurrBiol_replication/demog/PTx_n878_ModDev_df_20170221.csv")

############################
### Probabilistic Yeo PC ###
############################
Integrated_rawSC <- gam(PTx_integrated_Araw_Yeo_mean_PartCoeff_groupThresh_totalStrengthNorm_scale125 ~ s(age_in_yrs,k=4) + Sex + meanRELrms, data = PTx_n878_df)
v1 <- visreg(Integrated_rawSC,"age_in_yrs",plot=FALSE)

Integrated_volNorm <- gam(PTx_integrated_AvolNorm_Yeo_mean_PartCoeff_groupThresh_totalStrengthNorm_scale125 ~ s(age_in_yrs,k=4) + Sex + meanRELrms, data = PTx_n878_df)
v2 <- visreg(Integrated_volNorm ,"age_in_yrs",plot=FALSE)

Integrated_Aprop <- gam(PTx_integrated_Aprop_Yeo_mean_PartCoeff_groupThresh_totalStrengthNorm_scale125 ~ s(age_in_yrs,k=4) + Sex + meanRELrms, data = PTx_n878_df)
v3 <- visreg(Integrated_Aprop,"age_in_yrs",plot=FALSE)

## Figures 5G-I
par(mfrow=c(1,3))
plot(v1,xlab="Age (years)",ylab="Prob. Integrated Yeo PC - Streamline Count",cex.axis=1.3)
plot(v2,xlab="Age (years)",ylab="Prob. Integrated Yeo PC - Streamline Density",cex.axis=1.3)
plot(v3,xlab="Age (years)",ylab="Prob. Integrated Yeo PC - Connectivity Probability",cex.axis=1.3)

## Age pvals for Yeo PartCoeff
Integrated_rawSC_YeoPC_agePval <- summary(gam(PTx_integrated_Araw_Yeo_mean_PartCoeff_groupThresh_totalStrengthNorm_scale125 ~ s(age_in_yrs,k=4) + Sex + meanRELrms,fx=TRUE, data = PTx_n878_df))$s.table[1,4]
Integrated_volNorm_YeoPC_agePval <- summary(gam(PTx_integrated_AvolNorm_Yeo_mean_PartCoeff_groupThresh_totalStrengthNorm_scale125 ~ s(age_in_yrs,k=4) + Sex + meanRELrms,fx=TRUE, data = PTx_n878_df))$s.table[1,4]
Integrated_Aprop_YeoPC_agePval <- summary(gam(PTx_integrated_Aprop_Yeo_mean_PartCoeff_groupThresh_totalStrengthNorm_scale125 ~ s(age_in_yrs,k=4) + Sex + meanRELrms,fx=TRUE, data = PTx_n878_df))$s.table[1,4]

## Reset graphical parameters for visualization
par(mfrow=c(1,1))

###########################################
### FIGURE 6: Modularity and Efficiency ###
###########################################

## Figure 6A: Global Efficiency over Age
FA_GlobalEfficiency_gam <- gam(FA_GlobalEfficiency_scale125 ~ s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
summary(FA_GlobalEfficiency_gam)
visreg(FA_GlobalEfficiency_gam,"age_in_yrs",xlab="Age (years)",ylab="Global Efficiency")

## Figer 6B: Global Efficiency over Modular Segregation
FA_GlobalEff_ModSeg_gam <- gam(FA_GlobalEfficiency_scale125 ~ FA_Yeo7system_wholebrain_partCoeff_scale125 + s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125,fx=TRUE, data = LTN_n882_df)
summary(FA_GlobalEff_ModSeg_gam)
visreg(FA_GlobalEff_ModSeg_gam,"FA_Yeo7system_wholebrain_partCoeff_scale125", xlim=c(0.75,0.79),xlab="Mean Participation Coefficient",ylab="Global Efficiency")

# Fig 6E
avg_Yeo7_within_EdgeStrength_gam <- gam(FA_GlobalEfficiency_scale125  ~ new_k4_Yeo7system_FA_Avg_sigAge_withinModule_edgeStrength_scale125 +  s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125, fx=TRUE, data = LTN_n882_df)
summary(avg_Yeo7_within_EdgeStrength_gam)$p.table[2,4]
visreg(avg_Yeo7_within_EdgeStrength_gam, "new_k4_Yeo7system_FA_Avg_sigAge_withinModule_edgeStrength_scale125", xlab="Avg Within-Module Connectivity",ylab="Global Efficiency", ylim=c(0.23,0.265))

## Fig 6F 
avg_Yeo7_between_EdgeStrength_gam <- gam(FA_GlobalEfficiency_scale125 ~ new_k4_Yeo7system_FA_Avg_sigAge_BetweenModule_edgeStrength_scale125 +  s(age_in_yrs,k=4) + meanRELrms + Sex + FA_Total_Network_Strength_scale125, fx=TRUE, data = LTN_n882_df)
summary(avg_Yeo7_between_EdgeStrength_gam)$p.table[2,4]
visreg(avg_Yeo7_between_EdgeStrength_gam, "new_k4_Yeo7system_FA_Avg_sigAge_BetweenModule_edgeStrength_scale125", xlab="Avg Between-Module Connectivity",ylab="Global Efficiency", ylim=c(0.23,0.265))

####################################
### FIGURE 7 - Cognitive Effects ###
####################################

## System-specific Executive Effeciency effects on Avg PC
Visual_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_Visual ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
Visual_Exec_pval <- summary(Visual_Exec_gam)$p.table[2,4]

Somatomotor_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_Somatomotor ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
Somatomotor_Exec_pval <- summary(Somatomotor_Exec_gam)$p.table[2,4]

DorsAtt_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_DorsalAttention ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
DorsAtt_Exec_pval <- summary(DorsAtt_Exec_gam)$p.table[2,4]

VentAtt_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_VentralAttention ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
VentAtt_Exec_pval <- summary(VentAtt_Exec_gam)$p.table[2,4]

Limbic_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_Limbic ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
Limbic_Exec_pval <- summary(Limbic_Exec_gam)$p.table[2,4]

FP_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_Frontoparietal ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
FP_Exec_pval <- summary(FP_Exec_gam)$p.table[2,4]

DMN_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_Default ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
DMN_Exec_pval <- summary(DMN_Exec_gam)$p.table[2,4]

Subcortical_Exec_gam <- gam(Avg_FA_partCoeff_Yeo7system_Subcortical ~ s(age_in_yrs,k=4) + F3_Executive_Efficiency + FA_Total_Network_Strength_scale125 + meanRELrms + Sex, fx=TRUE, data=Cog_n880_df)
Subcortical_Exec_pval <- summary(Subcortical_Exec_gam)$p.table[2,4]

## Calculate Z-scores for Cognitive Effects
ModularSeg_Exec_Zscores <- as.data.frame(c(Visual_Exec_pval,Somatomotor_Exec_pval, DorsAtt_Exec_pval, VentAtt_Exec_pval, Limbic_Exec_pval, FP_Exec_pval, DMN_Exec_pval, Subcortical_Exec_pval))
colnames(ModularSeg_Exec_Zscores) <- "pvals"
ModularSeg_Exec_Zscores$Zscores <- abs(qnorm(ModularSeg_Exec_Zscores$pvals, lower.tail=FALSE))
ModularSeg_Exec_Zscores$System_names <- c("Visual","Somatomotor","DorsalAtt","VentralAtt","Limbic", "Frontoparietal","Default","Subcortical")

## FIGURE 7A
Fig7A_ModSeg_Exec_barplot <- ggplot(data=ModularSeg_Exec_Zscores, aes(System_names, Zscores)) + geom_bar(fill=c("orange", "firebrick4", "purple", "forestgreen", "skyblue", "darkgray","dark blue", "khaki1"), col="black",stat="identity",position=position_dodge()) + scale_x_discrete(limits=c("Frontoparietal","Default","Visual","DorsalAtt","VentralAtt","Subcortical","Somatomotor","Limbic")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border= element_rect(colour="black",fill=FALSE, size=1))
Fig7A_ModSeg_Exec_barplot

#############################################
### Cognitive mediation analysis (Fig 7B) ###
#############################################

## Regress covariates out
Age_lm <- lm(age_in_yrs ~ Sex + meanRELrms + FA_Total_Network_Strength_scale125, data=Cog_n880_df)
Age_resid <- resid(Age_lm)

ExecEff_lm <- lm(F3_Executive_Efficiency ~ Sex + meanRELrms + FA_Total_Network_Strength_scale125, data=Cog_n880_df)
ExecEff_resid <- resid(ExecEff_lm)

ModSeg_lm <- lm(FA_Yeo7system_wholebrain_partCoeff_scale125  ~ Sex + meanRELrms + FA_Total_Network_Strength_scale125, data=Cog_n880_df)
ModSeg_resid <- resid(ModSeg_lm)

## Standardize independent (X), dependent (Y), and mediating (M) variables
X <- as.data.frame(scale(Age_resid))
Y <- as.data.frame(scale(ExecEff_resid)) 
M <- as.data.frame(scale(ModSeg_resid))

Data <- data.frame(X=X, Y=Y, M=M)
Data <- data.frame(cbind(X,Y,M))
colnames(Data) <- c("X", "Y", "M")

model <- ' # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
fit_sem <- sem(model, data = Data, se="bootstrap", bootstrap=10000)
summary(fit_sem, fit.measures=TRUE, standardize=TRUE, rsquare=TRUE)

## Calculate bootstrapped confidence intervals for the indirect (c') effect
boot.fit <- parameterEstimates(fit_sem, boot.ci.type="perc",level=0.95, ci=TRUE,standardized = TRUE)
boot.fit
