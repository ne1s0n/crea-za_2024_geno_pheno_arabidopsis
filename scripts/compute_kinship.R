library(madbito)

#go to the folder
setwd('~/research/crea-za_2024_geno_pheno_arabidopsis/data/Cauline Leaf number (CL)')


# DATA LOAD ---------------------------------------------------------------
read_SNP_table = function(infile){
  infile_decompressed = gsub(pattern = '.zip', replacement = '', x = infile)
  unzip(infile)
  
  #reading the header
  header = scan(file = infile_decompressed, what = character(), nlines = 1)
  header = strsplit(header, split = ',')[[1]]
  
  #reading the data
  res = matrix(data = scan(file = infile_decompressed, what = integer(), skip = 1, sep = ','), 
              ncol = length(header), byrow = TRUE)
  colnames(res) = header
  
  #rownames from accession id
  rownames(res) = paste(sep='', 'accession_', res[,1])
  res = res[,-1]
  
  #cleanup
  file.remove(infile_decompressed)
  
  return(res)
}

#read files
SNP_test = read_SNP_table('X_test_705.csv.zip')
SNP_train = read_SNP_table('X_train_705.csv.zip')

#merge
SNP = rbind(SNP_test, SNP_train)
rm(SNP_test, SNP_train)

# COMPUTE KINSHIP ---------------------------------------------------------
library(madbito)

writeLines('computing kinship')
k = kinship.AB(data = SNP, ploidy = 1, verbose = TRUE)

writeLines('saving')
outfile = '../kinship_705.csv.gz'
fp = gzfile(description = outfile, open = 'w')
write.csv(k, file = fp)
close(fp)

