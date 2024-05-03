library(madbito)

#go to the folder
setwd('~/research/crea-za_2024_geno_pheno_arabidopsis/data/DaysToFlowering1/')


# DATA LOAD ---------------------------------------------------------------

#read files
SNP_test = read_SNP_table('X_test_703.csv.zip')
SNP_train = read_SNP_table('X_train_703.csv.zip')

#merge
SNP = rbind(SNP_test, SNP_train)
rm(SNP_test, SNP_train)

# COMPUTE KINSHIP ---------------------------------------------------------
library(madbito)

writeLines('computing kinship')
k = kinship.AB(data = SNP, ploidy = 1, verbose = TRUE)

writeLines('saving')
outfile = '../kinship_703.csv.gz'
fp = gzfile(description = outfile, open = 'w')
write.csv(k, file = fp)
close(fp)

