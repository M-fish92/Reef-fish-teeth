#Ecomorphology of coral reef fishes
#written by Michalis Mihalitsis
#-------------------------------------------------------------------------------------
rm(list=ls()) #clears R workspace (helps prevent errors from old data)
#-------------------------------------------------------------------------------------
setwd("add personal path")
library(readxl)
#load packages
library(stats)
library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(broom.mixed)
library(ggfortify) #for model diagnostics
library(sjPlot)    #for outputs
library(knitr)     #for kable
library(effects)   #for partial effects plots
library(ggeffects) #for effects plots in ggplot
library(emmeans)   #for estimating marginal means
library(MASS)      #for glm.nb
library(MuMIn)     #for AICc
library(tidyverse) #for data wrangling
library(DHARMa)   #for residuals and diagnostics
library(lme4)      #for glmer
library(glmmTMB)    #for glmmTMB
library(performance) #for diagnostic plots
library(see)         #for diagnostic plots
library(ggrepel)     #to avoid overlap of text in plots
library(ggConvexHull) #to add convex hulls in ggplot
library(vegan)
library(ape)
##############################################################################################
teeth = read.csv('Teeth_measurements2_final.csv',header = T, stringsAsFactors=T)

### Next, we import the dataset from Burns et al. (2024) that contains SL, gape size, Trophic ID etc of our specimens,
#and we merge that with our teeth dataset based on the ID number of each specimen.
jaws = read.csv('Burns_dataset.csv', header = T, stringsAsFactors = T)
jaws = jaws[,c(1,5,6,15,17,19,27,29,44)]
teeth2 = merge(teeth,jaws,by='ID',all.x=F)
###

### Next, we standardise the tooth variables, i.e., we divide measurements
#by the SL of each individual fish
teeth3=sweep(teeth2[,c(5:72)],1,teeth2[,79],'/')
teeth3= teeth3 %>% mutate(ID=teeth2$ID,Species=teeth2$Species,Trophic=teeth2$Trophic,
                          Function=teeth2$Function)

### Next, we summarise measurements for each species
teeth3 = teeth3 %>% group_by(Species,Trophic,Function) %>% 
  summarise(across(everything(), mean), .groups = 'drop')

### Next, we summarise measurements across the same variable and jaw systems
teeth3$a <- rowMeans(teeth3[,4:13], na.rm = T)
teeth3$b <- rowMeans(teeth3[,14:15], na.rm = T)
teeth3$c <- rowMeans(teeth3[,16:17], na.rm = T)
teeth3$d <- rowMeans(teeth3[,18:23], na.rm = T)
teeth3$e <- rowMeans(teeth3[,24:29], na.rm = T)
teeth3$f <- rowMeans(teeth3[,30:35], na.rm = T)
teeth3$g <- rowMeans(teeth3[,36:41], na.rm = T)
teeth3$h <- rowMeans(teeth3[,42:47], na.rm = T)
teeth3$i <- rowMeans(teeth3[,48:53], na.rm = T)
teeth3$j <- rowMeans(teeth3[,54:59], na.rm = T)
teeth3$k <- rowMeans(teeth3[,60:65], na.rm = T)
teeth3$l <- rowMeans(teeth3[,66:71], na.rm = T)
teeth3$m <- apply(teeth3[,18:23],1,var,na.rm=T)
teeth3$n <- apply(teeth3[,4:13],1,var,na.rm=T)
teeth3$o <- teeth3$a/teeth3$f
teeth3$p <- teeth3$d/teeth3$f
teeth3$q <- teeth3$f/teeth3$g
teeth3$r <- teeth3$d/teeth3$e
teeth3$s <- teeth3$h/teeth3$i
teeth3$t <- teeth3$l/teeth3$k
teeth3$u <- teeth3$j/teeth3$k

################################################################################
#Next we make a dataset that has all multicuspid variables, as we will need that later
compdat = teeth3 %>% select(Species,Trophic,Function,h,i,j,k,l,s,t,u)
### next, we remove all variables related to multicuspid teeth
teeth3 = teeth3[,-c(42:50)]
### Next we remove all columns with NAs
teeth3 <- teeth3 %>% select(where(~ !anyNA(.x)))
### Next, we do a Principal Component Analysis
test = prcomp(teeth3[,c(16:28)],scale. = T,center = T)
biplot = biplot(test,var.axes=T,cex=0.8,xlabs=1:nrow(teeth3),
                pch = 16,xlab = "PC 1", ylab = "PC 2")
summary(test)
###
library(ggplot2)
library(ggfortify)
library(ggforce)
library(reshape2)
library(corrplot)
###
cor_matrix <- cor(teeth3[,c(16:28)])
corrplot(cor_matrix)
# Create a biplot using ggplot2 and ggfortify
autoplot(test, data = teeth3,label = TRUE, label.label = "Species", label.size=3,label.colour='Function',
         loadings = F, loadings.label = F,size=4, loadings.label.size = 3)+
  geom_point(aes(colour=Function),size=4)+
  scale_color_brewer(palette = "Set2")+
  theme_classic()
autoplot(test, data = teeth3,label = F,label.colour='Function',colour='Function',
         loadings = F, loadings.label = F,size=4, loadings.label.size = 3,frame=T,frame.type='convex')+
  geom_point(aes(colour=Function),size=4)+
  scale_color_brewer(palette = "Set2")+
  theme_classic()
# Next, we add the variables of diet and bite force as well.
preytype = read.csv('finaldat.csv',header = T, stringsAsFactors=T)
teeth3 = merge(teeth3,preytype,by='Species',all.x=F)
###
autoplot(test, data = teeth3,label = F,label.colour='Function2',colour='Function2',
         loadings = F, loadings.label = F,size=4, loadings.label.size = 3,frame=T,frame.type='convex')+
  geom_point(aes(colour=Function2),size=4)+
  scale_color_manual(values = c("high" = "black",'medium'= '#d1995d','low'='#9b455b'))+
  theme_classic()
autoplot(test, data = teeth3,label = F,label.colour='Diet',colour='Diet',
         loadings = F, loadings.label = F,size=4, loadings.label.size = 3,frame=T,frame.type='convex')+
  geom_point(aes(colour=Diet),size=4)+
  scale_color_manual(values = c('GC'= '#f9ab64','MIS'='#f9cb64','PK'='#f97564','OM'='#bcbcbc','MIB'='#67c1bc',
                                'HD'='#67a6c1','SI'='#6781c1'))+
  scale_fill_manual(values = c('GC'= '#f9ab64','MIS'='#f9cb64','PK'='#f97564','OM'='#bcbcbc','MIB'='#67c1bc',
                                'HD'='#67a6c1','SI'='#6781c1'))+
  theme_classic()
###
library("FactoMineR")
library("factoextra")
res.pca <- PCA(teeth3[,c(16:28)], graph = T)
########################################################################################################
#### COMPLEX TEETH ONLY  ###############################################################################
########################################################################################################
compdat = merge(compdat,preytype,by='Species',all.x=F)
compdat = na.omit(compdat)
test = prcomp(compdat[,c(4:11)],scale. = T,center = T)
biplot = biplot(test,var.axes=T,cex=0.8,xlabs=1:nrow(compdat),
                pch = 16,xlab = "PC 1", ylab = "PC 2")
summary(test)
autoplot(test, data = compdat,label = TRUE, label.label = "Species", label.size=3,label.colour='Function',
         loadings = F, loadings.label = F,size=4, loadings.label.size = 3)+
  geom_point(aes(colour=Function),size=4)+
  theme_classic()
autoplot(test, data = compdat,
         loadings = T, loadings.label = T,size=4, loadings.label.size = 3)+
  theme_classic()
autoplot(test, data = compdat,label = TRUE, label.label = "Species", label.size=3,label.colour='Function2',
         colour='Function2',loadings = F, loadings.label = F,size=4, loadings.label.size = 3,frame=T,
         frame.type='convex')+
  geom_point(aes(colour=Function2),size=4)+
  scale_color_manual(values = c("high" = "black",'medium'= '#d1995d','low'='#9b455b'))+
  theme_classic()
res.pca <- PCA(compdat[,c(4:11)], graph = T)
###############################################################################
########     NEXT WE INCORPORATE PHYLOGENY INTO OUR ANALYSIS      #############
###############################################################################
library(phytools)
library(ape)
library(ggtree)
library(tidyverse)
library(treeio)
library(geiger)
tr = read.tree('actinopt_12k_treePL.tre')
jaws = read.csv('Burns_dataset.csv', header = T, stringsAsFactors = T)
jaws = jaws %>% select(Species,Tree_Name,out_lever_lower_jaw,closing_in_lever,
                       opening_in_lever,Body_Mass,Adductor_Mass,Protrusion_diff_abs,SL)
jaws$Protrusion_std <- jaws$Protrusion_diff_abs/jaws$SL
jaws$adductor_std <- jaws$Adductor_Mass/jaws$Body_Mass
jaws$jaw_closing_lever_ratio <- jaws$closing_in_lever/jaws$out_lever_lower_jaw
jaws$velocity_advantage <- jaws$out_lever_lower_jaw/jaws$closing_in_lever
jaws = jaws %>% group_by(Species,Tree_Name) %>% 
  summarise(across(everything(), mean, na.rm=T), .groups = 'drop')
teeth3 = left_join(teeth3,jaws,by='Species')
teeth3$Tree_Name <- jaws$Tree_Name[match(teeth3$Species, jaws$Species)]
species <- teeth3$Tree_Name
rownames(teeth3) <- teeth3$Tree_Name
species_to_drop <- setdiff(tr$tip.label, species)
tr = drop.tip(tr,species_to_drop)
name.check(tr,teeth3)
############################
library(phylolm)
test = prcomp(teeth3[,c(16:28)],scale. = T,center = T)
teeth3$pc1_values = test$x[,1]
#Define the models of evolution you want to test
model_types <- c("BM", "lambda", "OUrandomRoot", "kappa", "delta", "EB")
#Run phylolm across all types using a loop
fitted_models <- lapply(model_types, function(m) {
  phylolm(adductor_std ~ pc1_values, data = teeth3, phy = tr, model = m)})
names(fitted_models) <- model_types
#Extract AIC and Log-Likelihood into a comparison table
aic_table <- data.frame(Evolutionary_Model = names(fitted_models),AIC    = sapply(fitted_models, AIC),
  LogLik = sapply(fitted_models, function(x) x$logLik))
# Sort from best (lowest AIC) to worst
aic_table <- aic_table[order(aic_table$AIC), ]
print(aic_table)
best_model <- phylolm(adductor_std ~ pc1_values, data = teeth3, phy = tr, model = "lambda")
summary(best_model)
############################ THE SAME CODE IS USED FOR WHEN TESTING THE OTHER VARIABLES SUCH AS 
# JAW LEVER RATIO AND PROTRUSION. JUST CHANGE OUT THE NAME OF THE VARIABLE AND RERUN THAT CHUNK OF CODE.
#####   NEXT WE MAKE PLOTS FOR THE REGRESSIONS.
#Create a grid of X values covering your data range to get a smooth line
x_range <- range(teeth3$pc1_values, na.rm = TRUE)
prediction_grid <- data.frame(pc1_values = seq(from = x_range[1], to = x_range[2], length.out = 100))
#Predict the Y values using your phylolm model coefficients
prediction_grid$predicted_Y <- predict(best_model, newdata = prediction_grid)
#
ggplot() +
  geom_point(data = teeth3, aes(x = pc1_values, y = adductor_std)) +
  geom_line(data = prediction_grid, aes(x = pc1_values, y = predicted_Y)) +
  theme_classic() 
#############################################################################################
##### NEXT WE INCORPORATE PHYLOGENY IN THE MULTIVARIATE ANALYSES
#############################################################################################
library(geomorph)
feeding_mode <- teeth3$Feeding_mode
names(feeding_mode) <- teeth3$Tree_Name
bite_force <- teeth3$Function2
names(bite_force) <- teeth3$Tree_Name
Diet <- teeth3$Diet
names(Diet) <- teeth3$Tree_Name
#
Y_matrix <- as.matrix(teeth3[, 16:28])
#
gdf <- geomorph.data.frame(Y = Y_matrix, Feeding_Mode=feeding_mode,
                           Bite_Force=bite_force,Diet=Diet, phy = tr)
# Test the effect of the different factors on your multivariate traits (Y)
anova_feeding_mode <- procD.pgls(f1 = Y ~ Feeding_Mode+Diet+Bite_Force, phy = tr,
                                 data = gdf, iter = 10000,SS.type   = "III")
summary(anova_feeding_mode)
paircomps=pairwise(anova_feeding_mode,groups=gdf$Bite_Force)
summary(paircomps)
disp <- morphol.disparity(anova_feeding_mode, groups = ~ feeding_mode, iter = 10000, seed = 42)
disp <- morphol.disparity(anova_feeding_mode, groups = ~ Diet, iter = 10000, seed = 42)
disp <- morphol.disparity(anova_feeding_mode, groups = ~ bite_force, iter = 10000, seed = 42)
summary(disp)
#############################################################################################
##### NEXT WE MAKE A TREE FOR COMPLEX TEETH AND DO THE SAME ANALYSIS
#############################################################################################
compdat$Tree_Name <- jaws$Tree_Name[match(compdat$Species, jaws$Species)]
Tree_Name <- compdat$Tree_Name
rownames(compdat) <- compdat$Tree_Name
Tree_Name_to_drop <- setdiff(tr$tip.label, Tree_Name)
tr_compdat = drop.tip(tr,Tree_Name_to_drop)
#######
feeding_mode <- compdat$Feeding_mode
names(feeding_mode) <- compdat$Tree_Name
bite_force <- compdat$Function2
names(bite_force) <- compdat$Tree_Name
Diet <- compdat$Diet
names(Diet) <- compdat$Tree_Name
#######
Y_matrix <- as.matrix(compdat[, 4:11])
gdf <- geomorph.data.frame(Y = Y_matrix, Feeding_Mode=feeding_mode,
                           Bite_Force=bite_force,Diet=Diet, phy = tr_compdat)
anova_feeding_mode <- procD.pgls(f1 = Y ~ Feeding_Mode+Diet+Bite_Force, phy = tr_compdat,
                                 data = gdf, iter = 10000,SS.type   = "III")
anova_feeding_mode <- procD.pgls(f1 = Y ~ Feeding_Mode, phy = tr_compdat,
                                 data = gdf, iter = 10000,SS.type   = "III")
anova_feeding_mode <- procD.pgls(f1 = Y ~ Diet, phy = tr_compdat,
                                 data = gdf, iter = 10000,SS.type   = "III")
anova_feeding_mode <- procD.pgls(f1 = Y ~ Bite_Force, phy = tr_compdat,
                                 data = gdf, iter = 10000,SS.type   = "III")
summary(anova_feeding_mode)
disp <- morphol.disparity(anova_feeding_mode, groups = ~ bite_force, iter = 10000, seed = 42)
disp <- morphol.disparity(anova_feeding_mode, groups = ~ Diet, iter = 10000, seed = 42)
disp <- morphol.disparity(anova_feeding_mode, groups = ~ feeding_mode, iter = 10000, seed = 42)
summary(disp)
#################################################################################
#########      NEXT, WE RUN THE PHYLOLM MODELS FOR COMPLEX TEETH
#################################################################################
compdat = left_join(compdat,jaws,by='Species')
test = prcomp(compdat[,c(4:11)],scale. = T,center = T)
compdat$pc1_values = test$x[,1]
compdat$Tree_Name <- jaws$Tree_Name[match(compdat$Species, jaws$Species)]
Tree_Name <- compdat$Tree_Name
rownames(compdat) <- compdat$Tree_Name
Tree_Name_to_drop <- setdiff(tr_compdat$tip.label, Tree_Name)
tr_compdat = drop.tip(tr_compdat,Tree_Name_to_drop)
name.check(tr_compdat,compdat)
#Define the models of evolution you want to test
model_types <- c("BM", "lambda", "OUrandomRoot", "kappa", "delta", "EB")
#Run phylolm across all types using a loop
fitted_models <- lapply(model_types, function(m) {
  phylolm(jaw_closing_lever_ratio ~ pc1_values, data = compdat, phy = tr_compdat, model = m)})
names(fitted_models) <- model_types
#Extract AIC and Log-Likelihood into a comparison table
aic_table <- data.frame(
  Evolutionary_Model = names(fitted_models),AIC    = sapply(fitted_models, AIC),
  LogLik = sapply(fitted_models, function(x) x$logLik))
# Sort from best (lowest AIC) to worst
aic_table <- aic_table[order(aic_table$AIC), ]
print(aic_table)
best_model <- phylolm(jaw_closing_lever_ratio ~ pc1_values, data = compdat, phy = tr_compdat, model = "OUrandomRoot")
summary(best_model)
#Create a grid of X values covering your data range to get a smooth line
x_range <- range(compdat$pc1_values, na.rm = TRUE)
prediction_grid <- data.frame(pc1_values = seq(from = x_range[1], to = x_range[2], length.out = 100))
#Predict the Y values using your phylolm model coefficients
prediction_grid$predicted_Y <- predict(best_model, newdata = prediction_grid)
#
ggplot() +
  geom_point(data = compdat, aes(x = pc1_values, y = jaw_closing_lever_ratio)) +
  # Draw the phylogenetic regression line using our prediction grid
  geom_line(data = prediction_grid, aes(x = pc1_values, y = predicted_Y)) +
  theme_classic() 



