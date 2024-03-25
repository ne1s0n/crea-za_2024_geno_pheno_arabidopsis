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
  
  #sanity for individuals order
  stopifnot(all(genos$train$accession_id == phenos$train$accession_id))
  stopifnot(all(genos$test$accession_id == phenos$test$accession_id))
  
  #taking notes of individuals
  res$train_accession_id = genos$train$accession_id
  res$test_accession_id = genos$test$accession_id
  genos$train$accession_id = NULL
  genos$test$accession_id = NULL
  phenos$train$accession_id = NULL
  phenos$test$accession_id = NULL
  
  #final data storage
  res$genos_train = genos$train
  res$genos_test = genos$test
  res$phenos_train = phenos$train
  res$phenos_test = phenos$test
  
  #done
  return(res)
}
