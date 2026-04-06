#Read data Sets
setwd("/Users/selenekalra/Desktop/UofE/Dissertation/selenedissertation")

#Install Packages
install.packages("MASS")
install.packages("factoextra")
install.packages("ggplot2")
install.packages("vegan")
install.packages("broom")
install.packages("ggpubr")


#Read Packages
library(MASS)
library(factoextra)
library(ggplot2)
library(vegan)
library(broom)
library(stats)
library(ggpubr)

#import data
df <- read.csv('Master_Kelp_data.csv', header=TRUE, sep=",")

#compact summary
str(masterdata)

#quick names
age <- masterdata$age_d
Season <- masterdata$season
mg <- masterdata$Alginate.M.G.Ratio..area.based.
mannitol <- masterdata$percent_Mannitol
thick <-masterdata$relative_thickness
length <- masterdata$length_cm
wave <- masterdata$Hs_Weekly.average

#---------------------------Compound profile groups-----------------------------
df3 <- read.csv("4_group_check.csv")
table(df$Group)

mg <-df3$M.G
mani <-df3$Mannitol
manu <- df3$Man.URA

#Relationship with M/G and compound groups
ggplot(df, aes(Group, mg)) +
  geom_point(shape = 1, size=3)+
  xlab("Development Stage and Season") +
  ylab("M/G Alginate Ratio (Area-Based)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))

mean <- aggregate(mg ~ Group, data = df3, FUN = mean)
##ANOVA - significant
group.aov <-aov(mg~Group, data=df3)
summary(group.aov)
# normality - normal
shapiro.test(group.aov$residuals)
TukeyHSD(mg.aov)

# variance - bartlett for normal distribution
bartlett.test(mg~Group, data=df3)
TukeyHSD(mani.aov)


#Relationship with mannitol and compound groups
ggplot(df, aes(Group, mani)) +
  geom_point(shape = 1, size=3)+
  xlab("Development Stage and Season") +
  ylab("M/G Alginate Ratio (Area-Based)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))
mean <- aggregate(mani ~ Group, data = df, FUN = mean)
##ANOVA - significant
mani.aov<-aov(mani~Group, data =df3)
mani.aov
# normality - normal
shapiro.test(group.aov$residuals)
TukeyHSD(mani.aov)

#Relationship with mannuronic acid and compound groups
m <- ggplot(df3, aes(Group, MU))

m+geom_boxplot(aes(color=Group))+
  scale_colour_manual(values = c(
    "Spring" = "royalblue4",
    "Summer" = "tomato"
  )) +
  labs(x="Season of Growth", y="Mannuronic Acid Concentration (%)")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))

aggregate(MU ~ Group, data = df3, FUN = mean)


#--------------------------- Hypotheses --------------------------------------
#PCA Set-Up:
#For standardisation
install.packages("deco")
library(deco) #not for this version of R
install.packages("FactoMineR")
library(FactoMineR)

#import data
setwd("/Users/selenekalra/Desktop/UofE/Dissertation/Data Collection/R")
df2check <- read.csv('PCA_Bio_Check2.csv', header=TRUE, sep=",")

#str of data
str(df2)

#54 rows and 24 columns in data set
dim(df2)

summary(df2)


#remove missingness
nomiss <- sum(is.na(df2check))

#remove categorical (season)
kelp_numeric <- df2check[1:54, 3:10]

PCA <- prcomp(kelp_numeric, scale. =TRUE) #true ensure data is standardised & avoid bias
summary(PCA)

#Elements of PCA object
names(PCA)

#standard deviations
PCA$sdev
#Eigenvectors -> provided by rotation -> contains loadings per variable per component
PCA$rotation
#std and mean of ORIGINAL variables
PCA$center
PCA$scale

#pca score
PCA$x

df2check <- read.csv('PCA_Bio_Check.csv', header=TRUE, sep=",")

#Hypothesis 1: Testing M/G relationship with seasonal variation
##PCA M/G vs seasonal variation 
df2check <- read.csv('PCA_Bio_Check.csv', header=TRUE, sep=",")
summary(df2check)

df2check <- read.csv('PCA_Bio_Check.csv')


#Check #levels in Season
table(df2check$Season)

#ANOVA - not significant
mg <-df2check$M.G
season.aov <- aov(mg~Season, data=df2check)
summary(season.aov)

# normality - normal
shapiro.test(season.aov$residuals)

# variance - equal variance
var.test(mg ~ Season, data = df2)

# t-test -> Two-sample t-test - not significant
t.test(mg~Season, data=df2, var.equal=TRUE)

#Scree plot of Variance
fviz_eig(PCA, addlabels=TRUE, ylim=c(0,60))
#>Biplot with default settings
fviz_pca_biplot(PCA)
#>Biplot with labeled variables only, not sample #s
fviz_pca_biplot(PCA,
                label="var")
#>Biplot with customised colour groups
fviz_pca_biplot(PCA,
                label="var",
                habillage=df2check$Development)

#Biplot with customised colour groups and variables - SEASON
fviz_pca_biplot(PCA,
                label="varPCA",
                habillage= df2check$Season,
                addEllipses = TRUE,
                ellipse.level = 0.95,
                col.var="black",
                palette = (values=c("royalblue4", "tomato4"))) +
  geom_point(aes(colour = df2check$Season,
                 shape  = df2check$Development), size = 2) +
  labs(x = "PC1 (44.7%)",
       y = "PC2 (19.7%)")+
  scale_shape_manual(name = "Development",
                     values = c("Mature"    = 16, 
                                "Senescent" = 17))+
  theme_classic()+
  theme(panel.border = element_rect(color="black",
                                    fill = NA,
                                    size=1))

y+geom_point(aes(colour = df2check$Season,
                 shape  = df2check$Development),
             size = 3) 


varsPCA<- df2check[c("Carbon", "Nitrogen", "Chlorophyll","Mannitol", "M.G", "Fucose", "Glucose", "X")]

+ geom_text_repel(
  data = loadings,
  aes(PC1, PC2, label = Variable),
  size = 4,
  box.padding = 0.4,
  point.padding = 0.3,
  segment.color = "grey50",
  max.overlaps = Inf
)

library(ggplot2)
library(ggrepel)

####PERMANOVA - Does season affect the multivariate structure of the data?
library(vegan)

#Create the matrix of variables used in the PCA
varsPCA<- df2check[c("Carbon", "Nitrogen", "Chlorophyll","Mannitol", "M.G", "Fucose", "Glucose", "Man.URA")]

#Run PERMANOVA
adonis2(varsPCA~Season, data = df2check, method ="euclidean", permutations = 999)

#Equal Dispersion check - check assumption with PERMDISP
dist_matrix <- dist(varsPCA)
disp_season <- betadisper(dist_matrix, df2check$Season)
anova(disp_season)

#Check #levels in Development
table(df2check$Development)
#ANOVA
mg <- df2check$M.G
dev.aov <-aov(mg~Development, data=df2check)
summary(dev.aov)
# normality - normal
shapiro.test(dev.aov$residuals)

#Optional: Plot normality - the data is normal
par(mfrow = c(1, 2))
hist(dev.aov$residuals)
library(car)
qqPlot(season.aov$residuals,
       id=FALSE)

# variance - non equal variance
var.test(mg ~ Development, data = df2)

# t-test -> Welch two sample t-test - significant
t.test(mg~Development, data=df2)

#For PERMANOVA
#PCA Biplot with customised colour groups and variables - DEVELOPMENT
fviz_pca_biplot(PCA,
                label="var",
                habillage= df2check$Development,
                addEllipses = TRUE,
                ellipse.level = 0.95,
                col.var="black",
                palette = (values=c("royalblue4", "tomato4"))) +
  labs(x = "PC1 (19.8%)",
       y = "PC2 (44.8%)")+
  theme_classic()+
  theme(panel.border = element_rect(color="black",
                                    fill = NA,
                                    size=1))


###PERMANOVA
library(vegan)
#Create the matrix of variables used in the PCA
varsPCA<- df2check[c("Carbon", "Nitrogen", "Chlorophyll","Mannitol", "M.G", "Fucose", "Glucose", "Man.URA")]
#Run PERMANOVA
adonis2(varsPCA~Development, data = df2check, method ="euclidean", permutations = 999)
#Equal Dispersion check - check assumption with PERMDISP
dist_matrix <- dist(varsPCA)
disp_dev <- betadisper(dist_matrix, df2check$Development)
anova(disp_dev)

##Mannuronic acid vs seasonal variation
##RB Thickness vs seasonal variation
df <- read.csv('Master_Kelp_data.csv')
#Blade thickness vs M:G
thick <- df$relative_thickness
mg1 <- df$Alginate.M.G.Ratio..area.based.
#ANOVA
aov.thick <-aov(mg~thick, data=df)
summary(aov.thick)

aggregate(mg~thick, data =df, FUN = mean)

confint(aov.thick)


#Relationship between season and thickness?
aov.tss <- aov(thick~season, data=df)
summary(aov.tss)
plot(aov.tss)

TukeyHSD(aov.ts)

#Relationship between development and thickness?
aov.tsd <- aov(thick~Development, data=df)
summary(aov.tsd)
plot(aov.tsd)

TukeyHSD(aov.tsd)

#Relationship btw development and mannuronic acid?
manu <- df$Percent_Mannuronic
aov.manu <- aov(manu~Development, data = df)
summary(aov.manu)

summary(lm(manu ~ age_d, data = df), correlation=TRUE)

#Relationship btw Season and MU?
aov.manu <- aov(manu~Season, data = df)
summary(aov.manu)


qnorm(0.025, mean=)
?qnorm

plot(mg~thick, data=df)

ggplot(df, aes(thick, mg1)) +
  geom_point(shape = 2, size=1)+
  geom_point(aes(colour = season)) +
  geom_smooth(method = "lm",color="darkgray",size=1) +
  scale_colour_manual(values = c(
    "Spring" = "royalblue4",
    "Summer" = "tomato"
  )) +
  xlab("Relative Blade Thickness") +
  ylab("M:G Alginate Ratio (Area-Based)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))

ggplot(df, aes(thick, mg1)) +
  geom_point(shape = 2, size=1)+
  geom_point(aes(colour = Development)) +
  geom_smooth(method = "lm",color="darkgray",size=1) +
  scale_colour_manual(values = c(
    "Mature" = "royalblue4",
    "Senescent" = "tomato"
  )) +
  xlab("Relative Blade Thickness") +
  ylab("M:G Alginate Ratio (Area-Based)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))

lm(formula=mg~thick, data=df)


#thickness vs mannitol anova


# Normality
shapiro.test(aov.thick$residuals)
#Variance
fligner.test(mg ~ thick, data = df)

library(readxl)

lmmgthick = lm(mg1~thick, data = df) #Create the linear regression
plot(lmmgthick)

summary(lmmgthick) #Review the results


##RB Thickness vs M/G
aov.manni <- aov(Mannitol~thick, data=df3)
summary(aov.manni)
plot(aov.manni)
plot(Mannitol~thick, data=df3)

ggplot(df2check, aes(thick, Mannitol)) +
  geom_point(shape = 2, size=1)+
  geom_point(aes(colour = Season)) +
  geom_smooth(method = "lm",color="darkgray",size=1) +
  scale_colour_manual(values = c(
    "Spring" = "royalblue4",
    "Summer" = "tomato"
  )) +
  xlab("Relative Blade Thickness") +
  ylab("Mannitol Concentration(%)") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1)) 


#Hypothesis 2: Testing M/G relationship with lifecycle development stage
##PCA M/G vs  lifecycle development stage
##RB Thickness vs  lifecycle development stage

#Hypothesis 3: Testing mannitol content with season and lifecycle development stage
## Mannitol vs seasonal variation and lifecycle development stage
f <- ggplot(df2check, aes(Development, Mannitol))

f+geom_boxplot(aes(color=Season))+
  scale_colour_manual(values = c(
    "Spring" = "royalblue4",
    "Summer" = "tomato"
  )) +
  labs(x="Kelp Development Stage", y="Mannitol Concentration (%)")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))
#ANOVA
aov.mannitol <- aov(Mannitol~Group, data=df3)
summary(aov.mannitol)
# normality - normal
shapiro.test(aov.mannitol$residuals)
#mean
aggregate(Percent_Mannitol~Group, data =df3, FUN = mean)
#tukey test
tukey.test <- TukeyHSD(aov.mannitol)
tukey.test
plot(tukey.test)

## Mannitol vs RB thickness
df <- read.csv('Master_Kelp_data.csv')

thick <- df$relative_thickness
plot(manu~thick, data=df)
lm(formula=mani~thick, data=df)

lm(formula=Percent_Mannuronic~thick, data=df)
plot(Percent_Mannuronic~thick, data=df)



