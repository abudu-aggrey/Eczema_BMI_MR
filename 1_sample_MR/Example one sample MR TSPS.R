#Script for performing 1 sample MR with individual level data, where the outcome variable is binary
#MR method = two-stage predictor substitution. For more information see Bugress et al, 2015 (pmid:25773750)
#Modified from Tom Palmer

#********************************************************************************************

#install.packages("boot")
#install.packages("foreign")
#install.packages("tidyverse")

library(foreign)
library(boot)
library(tidyverse)

#read in data file (dataframe with phenotype data and genetic data (dosage) for each SNP)

#ideal format:
#ID Phenotype1 Phenotype2 SNP1 SNP2 SNP3   
#1
#2


fam_snps <- read.table("./fam_snps.txt", header=T)



#****************************************************************************************************************************************************************************************************

# function to return exp() IV ratio and CI using delta-method SE

delta <- function(x, y, u, g, PC_1, PC_2, PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10){
  bx = lm(x ~ g + u + PC1 +PC2 +PC3 +PC4 +PC5 
          +PC6 +PC7 +PC8 +PC9 +PC10, 
          subset=y==0)$coef[2]
  bxse = summary(lm(x ~ g + u + PC1 +PC2 +PC3 +PC4 +PC5 
                    +PC6 +PC7 +PC8 +PC9 +PC10, 
                    subset=y==0))$coef[2,2]
  by = glm(y ~ g + + u + PC1 +PC2 +PC3 +PC4 +PC5 
           +PC6 +PC7 +PC8 +PC9 +PC10, 
           family="binomial")$coef[2]
  byse = summary(glm(y ~ g + u + PC1 +PC2 +PC3 +PC4 +PC5 
                     +PC6 +PC7 +PC8 +PC9 +PC10, 
                     family="binomial"))$coef[2,2]
  theta <- 0
  se_ratio_approx = sqrt(byse^2/bx^2 + by^2*bxse^2/bx^4 - 2*theta*by/bx^3)
  ratio = by/bx
  or = exp(ratio)
  low = exp(ratio - 1.96*se_ratio_approx)
  upp = exp(ratio + 1.96*se_ratio_approx)
  se = se_ratio_approx
  return(list(or=or, low=low, upp=upp, ratio=ratio, se=se))
}

# read in data ----
#fam_snps = read_delim(file="./fam.txt", delim=" ")

snps = names(fam_snps[19:30]) #define columns for SNP names
n_snps = length(snps)

results.frame = data.frame(matrix(ncol=8, nrow=n_snps))
names(results.frame)=c("SNP", "Outcome", "Exposure","OR","OR_L95","OR_U95","Beta","SE")

#***********************************************************************************************************

# Loop over snps
for (i in 1:n_snps) {
  #  ----
  print(snps[i])
  results.frame[i,1] = snps[i]
  results.frame[i,2] = "Outcome" #outcome
  results.frame[i,3] = "Exposure" #exposure
  x = fam_snps$Exposure
  y <- fam_snps$Outcome # define outcome
  u <- fam_snps$GRSpheno.chip #if necessary for example control for genotyping chip
  
  PC1 = fam_snps$PC1
  PC2 =fam_snps$PC2
  PC3 =fam_snps$PC3
  PC4 =fam_snps$PC4
  PC5 =fam_snps$PC5
  PC6 =fam_snps$PC6
  PC7 =fam_snps$PC7
  PC8 =fam_snps$PC8
  PC9 =fam_snps$PC9
  PC10 =fam_snps$PC10
  
  #optional , if wish to perform MR with GRS
  #g <- fam_snps$SDscore

#to perform MR with individual SNPs
  g <- eval(parse(text=paste0("fam_snps$", snps[i]))) 
  res <- delta(x=x, y=y, u=u, g=g)
  
  # check estimate within CI limits
  if (((res$or < res$upp) & (res$or > res$low)) == FALSE) {
    print(paste("error for CI for", snps[i]))
  }
  
  results.frame[i,4] = res$or
  results.frame[i,5] = res$low
  results.frame[i,6] = res$upp
  results.frame[i,7] = res$ratio
  results.frame[i,8] = res$se
  res = NULL
}

#save results to excel file
library(xlsx)
write.xlsx(results.frame,file="./results.xlsx", sheetName="MR_results", row.names=F, append=T)


#write_csv(results.frame, path="./results.csv", col_names=TRUE)
