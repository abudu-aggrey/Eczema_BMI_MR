#Script for performing 1 sample MR with individual level data.
#MR method = two-stage least squares. 
#Modified from Tom Palmer

#********************************************************************************************

#install.packages("boot")
#install.packages("foreign")
#install.packages("tidyverse")

library(foreign)
library(boot)
library(tidyverse)

#read in data file

#ideal format:
#ID Phenotype1 Phenotype2 SNP1 SNP2 SNP3   
#1
#2

fam_snps <- read.table("./data.txt", header=T)

#****************************************************************************************************************************************************************************************************

# function to return exp() IV ratio and CI using delta-method SE

delta <- function(x, y, g, u, PC_1, PC_2, PC_3){
  bx = glm(x ~ g + u + PC1 +PC2 +PC3, 
           family=binomial(link="logit"))$coeff[2]
  
  bxse = summary(glm(x ~ g +u + PC1 +PC2 +PC3,
                     family=binomial(link="logit")))$coeff[2,2]
  
  by = lm(y ~ g + u + PC1 +PC2 +PC3)$coef[2]
  
  byse = summary(lm(y ~ g + u + PC1 +PC2 +PC3))$coef[2,2]
  theta <- 0
  se_ratio_approx = sqrt(byse^2/bx^2 + by^2*bxse^2/bx^4 - 2*theta*by/bx^3)
  ratio = by/bx
  #or = exp(ratio)
  low = ratio - 1.96*se_ratio_approx
  upp = ratio + 1.96*se_ratio_approx
  se = se_ratio_approx
  return(list(ratio=ratio, low=low, upp=upp, se=se))
}

# read in data ----
#fam_snps = read_delim(file="./data.txt", delim=" ")

snps = names(fam_snps[50:length(fam_snps)]) #define columns for SNP names
n_snps = length(snps)

results.frame = data.frame(matrix(ncol=7, nrow=n_snps))
names(results.frame)=c("SNP", "Outcome", "Exposure","Beta","L95","U95","SE")

#***********************************************************************************************************

# Loop over snps
for (i in 1:n_snps) {
  # ----
  print(snps[i])
  results.frame[i,1] = snps[i]
  results.frame[i,2] = "outcome" #outcome
  results.frame[i,3] = "exposure" #exposure
  x = fam_snps$Exposure
  y <- fam_snps$Outcome # define outcome
  u <- fam_snps$GRSpheno.chip #variable for genotyping chip
  
  PC1 = fam_snps$PC_1
  PC2 =fam_snps$PC_2
  PC3 =fam_snps$PC_3
  
  #optional , if wish to perform MR with GRS
  #g <- fam_snps$SDscore
  
  #to perform MR with individual SNPs
    g <- eval(parse(text=paste0("fam_snps$", snps[i]))) 
  res <- delta(x=x, y=y, u=u, g=g)
  
  # check estimate within CI limits
  if (((res$ratio < res$upp) & (res$ratio > res$low)) == FALSE) {
    print(paste("error for CI for", snps[i]))
  }
  
  results.frame[i,4] = res$ratio
  results.frame[i,5] = res$low
  results.frame[i,6] = res$upp
  results.frame[i,7] = res$se
  res = NULL
}

#save results to excel file
library(xlsx)
write.xlsx(results.frame,file="./TSLS_MR_analysis.xlsx", sheetName="MR_results", row.names=F, append=F)



