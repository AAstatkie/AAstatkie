id<-c(1, 2, 3, 4, 5)
gender<-c("Male", "Female", "Female", "Male", "Male")
age<-c(25, 22, 30, 29, 35)

id
gender
age

sampledf1<-data.frame(id, gender, age) # Creating a dataframe based on the three objects
  # created above.
sampledf1
sampledf2<-data.frame(id=c(1, 2, 3, 4, 5), gender=c("Male", "Female", "Female", 
  "Male", "Male"), age=c(25, 22, 30, 29, 35))
sampledf2
sampledf2$residence<-c("urban", "rural", "urban", "rural", "rural")
sampledf2

#===============Listing data sets available in R======================
# To list the data sets in all *available* packages
data(package = .packages(all.available = TRUE)) # NOTE: To make use of the daya set,
    # first you need to call (load) the package containing the data.

# To list only data sets from the base package
data()

#=============Importing data==========================================
# Importing CSV data
titanic<-read.csv("titanic.csv", header = TRUE)
titanic

#============Exporting data===========================================
# E.g., exporting to Stata
write_dta(titanic, "titanic.dta") # Exports 'titanic' data to Stata.

#==========Saving data in R==========================================
# E.g., to Save the 'titanic' data in R:
saveRDS(titanic, "titanic.rds") # This saves the 'titanic' data in R-format/
  # as an R object with the file name "titanic.rds" in the working directory.
readRDS("titanic.rds") # Reads or restores the RDS file.

#============Manipulating data======================
# Renaming variables.
# Read the dataset 'pima.tr2' from the 'MASS' package with object name 'dm', 
    # and rename the variable 'type' as 't2dm'.
library("MASS")
dm<-Pima.tr2
head(dm, 10)
colnames(dm) [8]="t2dm" #To change the variable "type" (8th column) to "t2dm".
head(dm, 10)

# Changing Variable types
# E.g., check the structure (type) of variables in the data frame 'dm' and change 
    # all 'integer' variable types into 'numeric' and 'character' variables into 
    # 'factor'.
str(dm)
#
dm$npreg<-as.numeric(dm$npreg)
dm$glu<-as.numeric(dm$glu)
dm$bp<-as.numeric(dm$bp)
dm$skin<-as.numeric(dm$skin)
dm$t2dm<-as.factor(dm$t2dm)
str(dm)

# Recoding numerical variables to factor
  # In the dm dataset, create BMI categories ("underweight", "normal", "overweight", 
  # "obese") and store it with a new variable name bmi_cat.
dm$bmi_cat<-cut(dm$bmi, breaks = c(-Inf, 18.49, 24.99, 29.9, Inf), 
      labels = c("Underweight", "normal", "Overweight", "Obese"),
        include.lowest=T, right = T)

# Recoding factor variables
# E.g., In the dataset dm, the variable t2dm is coded 'no=1' and 'yes=2'. 
    # Recode this variable so that 'no=0' and 'yes=1'.

# replace(x, list, values)
dm$t2dm<-as.numeric(dm$t2dm)
dm$t2dm <-replace(dm$t2dm, dm$t2dm == 1, 0) # replaced 1s by 0's
dm$t2dm <-replace(dm$t2dm, dm$t2dm == 2, 1) # replaced 2s by 1's
head(dm)

dm$t2dm<-as.factor(dm$t2dm)
levels(dm$t2dm)=c("No", "Yes")

# Changing the labels of factor variables
# E.g., in the dm dataset, change the label 'No' to "Absent" and 'Yes' to "Present". 

dm$t2dm<-recode_factor(dm$t2dm, No= "Absent", Yes = "Present") # recode works
    # equally well as recode_factor.

# Labeling categories of factor variables
# create the following data frame and label:
    # sex:  1=Male, 2=Female'
    # residence: 1=Urban, 2=Rural

exdf <- data.frame(age=c(22, 29, 30, 27, 35), 
                  sex=as.factor(c(1, 1, 2, 1, 2)),
                  residence=as.factor(c(1,2,1,2,2)))
levels(exdf$sex)=c("Male", "Female")
levels(exdf$residence)=c("Urban", "Rural")

# Computing new variables from existing ones
# E.g., In the data frame 'exdf' you created above create a new variable 'agemo' 
# which stands for age in months. Compute it as age*12.

exdf$agemo<-exdf$age*12

# Creating a new data frame from an existing one
# E.g., from the dm dataframe create a new dataframe  named 'dm_small' containing
    # only the variable 't2dm' (in 8th column) and 'bmi_cat' (in 9th column).

dm_small<-dm[, c(8,9)]
head(dm_small)
    
# Or, you can do it by directly specifying the column names as follows.

dm_small2<-dm[, c("t2dm","bmi_cat")]
head(dm_small2)

# The select function can be used to create a new dataset from an exiting one.
# E.g.,
dm_small3<-select(dm, "t2dm", "bmi_cat")
head(dm_small3)

#===============Filtering data===================
library(dplyr)
obese<-filter(dm, bmi_cat=="Obese") # Subsets only "obese" subjects from dm dataset.
obeseoverw<-filter(dm, bmi_cat=="Overweight" | bmi_cat=="Obese") ## Subsets those who
    # are either "overweight" or "obese" subjects from dm dataset. The subset
    # contains both categories (overweight and obese).
normalbmi<-filter(dm, bmi_cat=="normal") # Subsets only those with "normal" BMI.
obeset2dm<-filter(dm, t2dm=="Yes" & bmi_cat=="Obese") # Subsets those who have 
    # type-2 DM and are obese.

#============Sorting data========================
# The function sort (from the base package) sorts a vector into ascending or 
  # descending order. By default, it orders the data in ascending order. 
# E.g., sort the variable 'age' in the exdf dataset in ascending order.

sort(exdf$age)

# To sort the entire dataframe by 'age' use 'order' function as follows.

exdf[order(exdf$age),]

#===========Reshaping data=======================
SatisfactionData<-read.delim("Honeymoon Period.dat", head = TRUE)
saveRDS(SatisfactionData, "SatisfactionData.rds")

# E.g., the dataframe 'SatisfactionData.rds' a repeated-measures data in wide format.
# It has six columns: the person's ID; satisfaction measured at baseline, 
# 6 mo, 12 mo, and 18 mo; and gender.
# Reshape  this dataset to a long format.

SatisfactionData<-readRDS("SatisfactionData.rds")
SatisfactionStacked<-stack(SatisfactionData, 
    select = c("Satisfaction_Base", "Satisfaction_6_Months", "Satisfaction_12_Months",
    "Satisfaction_18_Months"))

# To unstack:
SatisfactionUnstacked<-unstack(SatisfactionStacked, values ~ ind)

# Reshaping data using 'melt()' and 'cast()'.
# Stack and unstack are not convenient for larger (and complex) datasets.
# Generally, better to use the melt and cast (from reshape2 package).
    # 'melt' changes wide data to long form.
    # 'cast' changes long data to wide form.
library(reshape2)
RestructuredData<-melt(SatisfactionData, id = c("Person", "Gender"), 
  measured = c("Satisfaction_Base", "Satisfaction_6_Months", "Satisfaction_12_
  Months", "Satisfaction_18_Months"))
    
    # NOTE: 1) id: This option specifies any variables in the dataframe that do not
    # vary over time.For these data we have two variables that don't vary over
    # time, the first is the person's identifier (Person), and the second is their
    # gender (Gender). We can specify these variables as id = c("Person", "Gender").
          # 2) measured: This option specifies the variables that do vary over time
    # or are repeated measures (i.e., scores within the same entity). In other words,
    # it specifies the names of variables currently in different columns that you
    # would like to be restructured so that they are in different rows. We have four
    # columns that we want to restructure (Satisfaction_Base, Satisfaction_6_Months,
    # Satisfaction_12_Months, Satisfaction_18_Months)

# To  change back to wide form:

WideData<-dcast(RestructuredData, Person + Gender ~ variable, value = "value")

#==========Exploring columns and rows============================
setwd("D:/Ayalew's_Folder/Wisdom_Barn/Biost_etc/Intro_to_Stata_and_R/
      Statistics_with_R")
dm<-readRDS("dm.rds")
nrow(dm) # number of rows (sample size)of the "dm" dataset.
ncol(dm) # number of columns (variables) of the "dm" dataset.
dm[1,4] # The observation in the first row, fourth column.
dm[,4] # All rows of the fourth column.
dm[2,] # All columns of the second row.

#==========Exploring missing observations=========================
is.na(dm) # To check for missing observation.
sum(is.na(dm)) # To check for missing observation in the entire dataset.
colSums(is.na(dmpima)) #A simpler way to get number of missing observations
# for each variable.
which(is.na(dm$bp)) # To identify the location of the missing value
which(is.na(dm$skin))
which(is.na(dm$bmi))

#===========Descriptive statistics=================================
summary(dm) # Simple summary statistics for all variables in the dm dataset.
summary(dm[,-c(8,9)]) # Simple summary statistics for all variables in the dm 
  # dataset excluding categorical variables.
lapply(dm, summary) # Simple summary statistics for all variables in the dm dataset.

# Calculating mean and standard deviation
mean(dm$glu) # Calculate the mean of the variable "glu" in the dm dataset.
sd(dm$glu) # Calculate the mean of the variable "glu" in the dm dataset.
cbind(Mean=mean(dm$glu), SD=sd(dm$glu)) # Combine the mean and SD output in rows
  # and columns
sapply(dm[,-c(8,9)], mean, na.rm=TRUE) # Calculates mean for all variables excluding the
  # variables in the 8th and 9th columns (categorical ones).
sapply(dm[,-c(8,9)], sd, na.rm=TRUE) # Calculates SD for all variables excluding the
  # variables in the 8th and 9th columns (categorical ones).
tapply(dm$glu, dm$bmi_cat, mean, na.rm=TRUE) # Calculates mean by levels of a 
  # categorical variable (mean of "glu" for different BMI categories).

# Descriptive statistics for categorical variables
table(dm$t2dm) # Frequency distribution for the variable "glu" (serum glucose level)
  # in the dataset dm.
table(dm$bmi_cat) # Frequency distribution for the variable "bmi_cat"
# in the dataset dm.
table(dm$bmi_cat)*100/(nrow(dm)) # Calculates the percentage (relative frequency)

# To obtain a frequency distribution table showing both absolute frequencies
  # and percentages use the package epiDisplay. 
epiDisplay::summ.factor(dm$t2dm)
epiDisplay::summ.factor(dm$t2dm)

#=====================Cross-tabulations===============================
table(dm$bmi_cat, dm$t2dm) # Cross-tab reporting cell counts
table(dm$bmi_cat, dm$t2dm)*100/nrow(dm) # To get percentage frequencies.

crosstab<- table(dm$bmi_cat, dm$t2dm)
prop.table(crosstab) # sum of all cells gives 1.
prop.table(crosstab, 1) # sum of rows cells gives 1.
prop.table(crosstab, 2) # Sum of column cells gives 1.

# "gmodels" has CrossTable function.
library(gmodels)

CrossTable(dm$bmi_cat, dm$t2dm, prop.chisq = TRUE)
CrossTable(dm$bmi_cat, dm$t2dm, digits=2, expected=TRUE, prop.chisq = TRUE, 
           chisq = TRUE,  fisher=TRUE, missing.include=FALSE, format=c("SPSS"))

CrossTable(dm$bmi_cat, dm$t2dm, digits=2, max.width = 5, expected=TRUE, 
           prop.r=TRUE, prop.c=TRUE, prop.t=TRUE, prop.chisq=TRUE, chisq = TRUE, 
           fisher=TRUE, mcnemar=FALSE,
           resid=FALSE, sresid=FALSE, asresid=FALSE,
           missing.include=FALSE,
           format=c("SPSS"), dnn = NULL)

chisq.test(dm$bmi_cat, dm$t2dm) # Chi-square test
fisher.test(dm$bmi_cat, dm$t2dm) # Fisher's exact test

epiDisplay::tabpct(dm$bmi_cat, dm$t2dm, decimal = 1, percent = "both")

# ==============Plotting with R==================================
# Bar plot using the based R packages.
barplot(table(dm$t2dm), col=c("RED", "BLUE"))
barplot(table(dm$t2dm), col=c("RED", "BLUE"), ylim = c(0,200), 
        main = "Bar plot of type 2 DM status", ylab = "Count", xlab = "Type 2 DM")
barplot(table(dm$t2dm)*100/nrow(dm), col=c("RED", "BLUE"), ylim = c(0,70),
        main = "Bar plot of type 2 DM", ylab = "Proportion", xlab = "Type 2 DM")

tab_bmicat<-table(dm$bmi_cat)*100/nrow(dm)
barplot(tab_bmicat, xlab="BMI category", ylab="Percent", 
        main = "BMI category among Pima native American Women",
        cex.axis = 1, cex.names = 1, cex.main=1.2,
        col=c("MAGENTA"))        

# Boxplot using the base R package

boxplot(dm$glu) # Without options

boxplot(dm$glu, range = 1.5, border = par("fg"), col = "burlywood1",
        xlab=NULL, ylab="Serum glucose level (mg/dl)",
        main="Blood glucose level of the participant",
        pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
        add = FALSE, at = NULL) # Embellished with options

# Box plot for sub-groups
boxplot(dm$glu ~ dm$t2dm, range = 1.5, border = par("fg"), col = "palegreen",
        xlab="Type 2 DM status",
        ylab="Serum glucose level (mg/dl)",
        main="Blood glucose level of the participant",
        pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
        add = FALSE, at = NULL) # Embellished with options

# Histogram of numeric variables. NOTE: Integers must be converted to numeric.
hist(dm$bmi, col = "PINK", xlab="BMI", main="Histogram of BMI", ylim = c(0,100))

hist(dm$bmi, col = "PINK", xlab="BMI", ylab = "Probability density", 
     main="Histogram of BMI", freq = FALSE) # Histogram based on probability 
      # densities

par(mfrow=c(1,2))
hist(dm$bmi, col = "PINK", xlab="BMI", main="Histogram of bmi")
hist(dm$glu, col = "PINK", xlab="Serum glucose", main="Histogram of serum glucose")

par(mfrow=c(1,1))
hist(dm$bmi, col = "PINK", xlab="BMI", main="Histogram of BMI", 
     xlim = c(10,60), freq=FALSE) # REPLACE FREQ WITH DENSITY
lines(density(dm$bmi, na.rm = TRUE, kernel = c("gaussian")), col = "BLUE", lwd=2) # To
# superimpose a density curve on the histogram.

# Scatter plot
plot(dm$glu~dm$bmi, col=c("BLUE"), xlab="BMI", ylab="Glucose", 
     pch=16)

# Scatter plot with a regression line
plot(dm$glu~dm$bmi, col=c("BLUE"), xlab="BMI", ylab="Glucose", 
     pch=16)
abline(lm(glu ~ bmi, data = dm), col = "magenta", lwd=3)

# Scatter plot with a regression line with 95% CI.
library(ggplot2)
scatter<-ggplot(dm, aes(bmi, glu))
scatter + geom_point(colour = "Blue", na.rm = TRUE) + 
  geom_smooth(method = "lm", colour = "Magenta", fill = "Magenta", na.rm = TRUE) + 
  labs(x = "Body Mass Index", y="Serum glucose level (mg/dl)")

# Scatter plot with a smoothing (lowess) curve with 95% CI.
scatter2<-ggplot(dm, aes(bmi, glu))
scatter2 + geom_point(colour = "Blue", na.rm = TRUE) + 
  geom_smooth(colour = "Magenta", fill = "Magenta", na.rm = TRUE) + 
  labs(x = "Body Mass Index", y="Serum glucose level (mg/dl)")

#==============Hypothesis testing=============================================
# t-test
t.test(dm$bp, mu=120, alternative = "two.sided")

# Normality check: graphical
par(mfrow=c(1,2))
hist(dm$bp, main = "Histogram of BP", freq = FALSE)
lines(density(dm$bp, na.rm = TRUE), col = "Blue", lwd = 2)
qqnorm(dm$bp, main = "Q-Q plot of BP")

# Statistical test of normality
shapiro.test(dm$bp) # Shapiro-Wilk test (often used for small sample sizes, n<50)
ks.test(dm$bp, pnorm) # Kolmogorov-Smirnov test (often for large sample size, n>/=50)

# Normality test for type-2 DM and non-type-2 DM patients
yes<-filter(dm, t2dm==1)
shapiro.test(yes$bp)

no<-filter(dm, t2dm=0)
shapiro.test(no$bp)

# Equality of variance test
var.test(dm$bp ~ dm$t2dm)

#Two-sample t-test
t.test(dm$bp ~ dm$t2dm, alternative = "two.sided", var.equal=TRUE)

#==============Paired-sample t-test============
library(datarium)

# Checking normality
par(mfrow=c(1,2))
hist(mice2$before, prob = TRUE, col = "magenta")
lines(density(mice2$before), col="Blue", lwd=2)
hist(mice2$after, prob = TRUE, col = "magenta")
lines(density(mice2$after), col="Blue", lwd=2)
shapiro.test(mice2$before)
shapiro.test(mice2$after)

# Equality of variance test
var.test(mice2$before, mice2$after) 

# Paired t-test        
t.test(mice2$before, mice2$after, alternative = "two.sided", paired=TRUE,
       na.rm=TRUE, var.equal=TRUE)

#===================One-way ANOVA==============================
# Example: A study is designed to test whether there is a difference in mean daily 
    # calcium intake in adults with normal bone density, 
    # adults with osteopenia (a low bone density which may lead to osteoporosis) 
    # and adults with osteoporosis. Adults 60 years of age with normal bone density, 
    # osteopenia and osteoporosis are selected at random from hospital records and 
    # invited to participate in the study. Each participant's daily calcium intake 
    # is measured based on reported food intake and supplements. 
    # The data are in the file "calcium_data.xlsx". 
# Is there a statistically significant difference in mean calcium intake in patients
    # with normal bone density as compared to patients with osteopenia 
    # and osteoporosis? (source: https://sphweb.bumc.bu.edu/otlt/MPH-Modules/
    # BS/BS704_HypothesisTesting-ANOVA/BS704_HypothesisTesting-Anova4.html)

# First, import the data into R.
library(readxl)
calcium_data <- read_excel("calcium_data.xlsx")

# The variable "Bone_density" is a character-data type. We should convert it to
    # factor-data type.
calcium_data$Bone_density<-as.factor(calcium_data$Bone_density)

# One-way ANOVA
ca_aov<-aov(Daily_ca_intake ~ Bone_density, data = calcium_data)
summary(ca_aov)

# Pair-wise comparison to see which groups have significant differences.This is done
    # only when the AOV shows a significant overall difference. In the current
    # there is no significant overall difference. We are doing pair-wise 
    # comparison here only to demonstrate how it is done.
TukeyHSD(ca_aov)

# To check the assumptions of ANOVA
par(mfrow=c(2,2))
plot(ca_aov)

library(car)
leveneTest(ca_aov)

# =======================Correlation====================================
# Correlation matrix of several variables
cor(dm[,-c(8,9)], use = "complete.obs")

# To limit the number of decimal places to only 3 digit:
round(cor(dm[,-c(8,9)], use = "complete.obs"), 3)

# To test whether the correlation coefficient (r) is significantly 
    # different from zero
cor.test(dm$bmi, dm$glu)

library(Hmisc)
rcorr(as.matrix(dm[,-c(8,9)]))

# =====================Regression=====================================
# ======Linear regression========

# Empty model
linear0<-lm(glu~1, data = dm, na.action = na.omit)
summary(linear0)

# Linear reg (Full model, with all IDVs).
linear1<-lm(glu~npreg+bp+bmi+ped+age, data = dm, na.action = na.omit)
summary(linear1)

# TO estimate 95% CI of the model (linear1).
confint(linear1, level=0.95)

# AIC & BIC to compare model goodness-of-fit.
AIC(linear0)
AIC(linear1)
#
BIC(linear0)
BIC(linear1)

# To obtain the VIF.
library(car)
vif(linear1)

# To obtain both tolerance and VIF.
library(olsrr)
ols_vif_tol(linear1)

# To plot model diagnostic models.
par(mfrow=c(2,2))
plot(linear1)

#Breusch-Pagan test test against heteroskedasticity.
library(lmtest)
bptest(linear1)

# ===================Logistic regression==============================
# First create a dataset without missing values. Missing values may create
    # different number of observations in different models which may make
    # some model dignostics such anova() not to work.
dm.complete<-na.omit(dm)

# Logistic regression (full model)
model1<-glm(t2dm~npreg+glu+bp+bmi+ped+age+skin, data = dm.complete,
            family = binomial, na.action = na.omit)
summary(model1)

#Logistic regression (condensed model)
model2<-glm(t2dm~npreg+glu+bp+bmi+ped+age, data = dm.complete, 
            family = binomial, na.action = na.omit)
summary(model2)

# Odds ratio with 95% CI.
or<-exp(coef(model2)) # Computes OR.
or
betaconfint<-confinterval<-confint(model2) # Computes 95% CI of coefficients.
betaconfint
aor<-exp(confinterval) # Computes 95% CI of OR.
aor

# CBIND: To bind together the AORs with 95% CIs.
exp(cbind(OR=coef(model2), confint(model2)))

# To get the output in the form of odds ratios with 95% CIs (for both CORs and
# AORs), the epiDisplay::logistic.display function provides a nice
# functionality.
library(epiDisplay)
epiDisplay::logistic.display(model2, alpha = 0.05, crude=TRUE, 
      crude.p.value = TRUE, decimal = 2, simplified = FALSE)

# Testing the goodness-of-fit of the logistic model
#=========LR test===============
# Subsequent logictic models can be compared using the likelihood ratio test.
# This can be done using the anova() function on subsequent nested
# logistic models, and asking R to yield a Chi-square test. Hence, now, to
# compare model1 with six IDVs with model2 containing five IDVs:
anova(model1, model2, test = "Chisq") 
#=========Pseudo-R square==========
library(DiscTools)
PseudoR2(model2, which = "Nagelkerke")

#=========Hosmer and Lemeshow test=====
library(ResourceSelection)
hoslem.test(model2$y, fitted(model2))

#=================== Saving R output================================
sink("Example.R.Output.txt")

# To unsink or stop sinking
sink()