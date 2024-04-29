library(madbito)

#loads (and, if necessary, computes) kinship matrix for the passed trait
load_kinship = function(infolder){
  res = list()
  
  #target files
  infile_train = file.path(infolder, 'kinship_train_gwas.csv.gz')
  infile_test = file.path(infolder, 'kinship_test_gwas.csv.gz')
  
  #if either don't exist we compute the kinship
  if (!file.exists(infile_train) | !file.exists(infile_test)){
    #no kinship available, let's compute it. 
    
    #load train genotypes
    candidates = list.files(infolder, pattern = 'X_train_gwas', full.names = TRUE, recursive = FALSE)
    stopifnot(length(candidates) == 1)
    SNP_train = read.csv(candidates, stringsAsFactors = FALSE)
    
    #load test genotypes
    candidates = list.files(infolder, pattern = 'X_test_gwas', full.names = TRUE, recursive = FALSE)
    stopifnot(length(candidates) == 1)
    SNP_test = read.csv(candidates, stringsAsFactors = FALSE)
    
    #sanity
    stopifnot(all(colnames(SNP_train) == colnames(SNP_test)))
    SNP = rbind(SNP_train, SNP_test)
    accessions = paste(sep='', 'accession_', SNP$accession_id)
    SNP$accession_id = NULL
    
    #computing kinship
    kinship = kinship.AB(SNP)
    
    #splitting back in train and test
    train_selector = rep(TRUE, nrow(kinship))
    train_selector[(nrow(SNP_train)+1):nrow(kinship)] = FALSE
    kinship_train = data.frame(cbind(accession_id=SNP_train$accession_id, kinship[train_selector,]))
    kinship_test  = data.frame(cbind(accession_id=SNP_test$accession_id,  kinship[!train_selector,]))
    
    #names
    colnames(kinship_train) = c('accession_id', accessions)
    colnames(kinship_test) = c('accession_id', accessions)
    rownames(kinship_train) = accessions[train_selector]
    rownames(kinship_test)  = accessions[!train_selector]
    
    #saving
    fp = gzfile(description = infile_train, open = 'w')
    write.csv(x = kinship_train, file = fp, row.names = TRUE)
    close(fp)
    fp = gzfile(description = infile_test, open = 'w')
    write.csv(x = kinship_test, file = fp, row.names = TRUE)
    close(fp)
  }
  
  #storing the data
  res$train = read.csv(file = infile_train, stringsAsFactors = FALSE, row.names = 1)  
  res$test = read.csv(file = infile_test, stringsAsFactors = FALSE, row.names = 1)  
  
  return(res)
}

#load train and test phenotypes from the passed data folder
load_phenotypes = function(infolder){
  res = list()
  
  #train
  candidates = list.files(infolder, pattern = 'y_train', full.names = TRUE, recursive = FALSE)
  stopifnot(length(candidates) == 1)
  res$train = read.csv(candidates, stringsAsFactors = FALSE)
  
  #test
  candidates = list.files(infolder, pattern = 'y_test', full.names = TRUE, recursive = FALSE)
  stopifnot(length(candidates) == 1)
  res$test = read.csv(candidates, stringsAsFactors = FALSE)
  
  #sanity
  tmp = intersect(res$y_train$accession_id, res$y_test$accession_id)
  stopifnot(tmp == 0)
  
  return(res)    
}

#load train and test genotypes from the passed data folder
load_genotypes = function(infolder){
  res = list()
  
  #train
  candidates = list.files(infolder, pattern = 'X_train_gwas', full.names = TRUE, recursive = FALSE)
  stopifnot(length(candidates) == 1)
  res$train = read.csv(candidates, stringsAsFactors = FALSE)
  
  #test
  candidates = list.files(infolder, pattern = 'X_test_gwas', full.names = TRUE, recursive = FALSE)
  stopifnot(length(candidates) == 1)
  res$test = read.csv(candidates, stringsAsFactors = FALSE)
  
  #sanity
  tmp = intersect(res$y_train$accession_id, res$y_test$accession_id)
  stopifnot(tmp == 0)
  
  return(res)
}

#load genos and phenos for the specified trait   
load_all = function(
    basefolder = '/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/data', 
    trait = c('Cauline Leaf number (CL)', 'DaysToFlowering1', 'DaysToFlowering2', 'DaysToFlowering3', 'Rosette Leaf number (RL)')
){
  #check for valid traits
  trait = match.arg(trait)
  
  #room for results
  res = list()
  res$trait = trait
  
  #load
  genos = load_genotypes(file.path(basefolder, trait))
  phenos = load_phenotypes(file.path(basefolder, trait))
  kinship = load_kinship(file.path(basefolder, trait)) 
  
  #sanity for individuals order
  stopifnot(all(genos$train$accession_id == phenos$train$accession_id))
  stopifnot(all(genos$test$accession_id == phenos$test$accession_id))
  stopifnot(all(kinship$train$accession_id == phenos$train$accession_id))
  stopifnot(all(kinship$train$accession_id == phenos$train$accession_id))
  
  #taking notes of individuals
  res$train_accession_id = genos$train$accession_id
  res$test_accession_id = genos$test$accession_id
  genos$train$accession_id = NULL
  genos$test$accession_id = NULL
  phenos$train$accession_id = NULL
  phenos$test$accession_id = NULL
  kinship$train$accession_id = NULL
  kinship$test$accession_id = NULL
  
  #final data storage
  res$genos_train = genos$train
  res$genos_test = genos$test
  res$phenos_train = phenos$train
  res$phenos_test = phenos$test
  res$kinship_train = kinship$train
  res$kinship_test = kinship$test
  
  #done
  return(res)
}

load_kinship_705 = function(infolder){
  infile = file.path(infolder, 'kinship_705.csv.gz')
  return(read.csv(infile, row.names = 1))
}

load_phenotypes_705 = function(infolder){
  tmp = load_phenotypes(infolder)
  tmp = rbind(tmp$train, tmp$test)
  tmp$accession_id = paste(sep='', 'accession_', tmp$accession_id)
  return(tmp)
}

load_splits_705 = function(infolder){
  #stub function
  return(NULL)
}

load_all_705 = function(
  basefolder = '/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/data', 
  trait = c('Cauline Leaf number (CL)', 'DaysToFlowering1', 'DaysToFlowering2', 'DaysToFlowering3', 'Rosette Leaf number (RL)')
){
  #check for valid traits
  trait = match.arg(trait)
  
  #room for results
  res = list()
  res$trait = trait
  
  #load data
  res$phenos = load_phenotypes_705(file.path(basefolder, trait))
  res$kinship = load_kinship_705(basefolder) 
  res$splits = load_splits_705(file.path(basefolder, trait))
  
  #sanity
  stopifnot(all(res$phenos$accession_id %in% rownames(res$kinship)))
  stopifnot(all(rownames(res$kinship) %in% res$phenos$accession_id))
  
  #same order
  res$kinship = res$kinship[res$phenos$accession_id, res$phenos$accession_id]
  
  #done
  return(res)
}