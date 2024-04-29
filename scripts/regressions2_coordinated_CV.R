library(GROAN)

#target traits
trait_list = c('Cauline Leaf number (CL)', 'DaysToFlowering1', 'DaysToFlowering2', 'DaysToFlowering3', 'Rosette Leaf number (RL)')

#setup, data management library
basefolder = '/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/'
source(file.path(basefolder, 'scripts', 'data_manager.R'))

#room for results
outfolder = file.path(basefolder, 'results')
outfile = 'accuracy_regression2_coordinated_CV.csv'
dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

#Bayesian regressors
nIter = 15000

# SUPPORT FUNCTIONS -------------------------------------------------------
get_split_data = function(all_data, split_current){
  res = list()
  
  #extracting the accession IDs for the current split
  IDs_train = subset(all_data$splits, split == split_current & train)$accession_id
  IDs_test  = subset(all_data$splits, split == split_current & !train)$accession_id
  
  #extracting the data
  res$kinship_train = all_data$kinship[IDs_train,]
  res$kinship_test = all_data$kinship[IDs_test,]
  res$phenos_train = subset(all_data$phenos, accession_id %in% IDs_train)
  res$phenos_test = subset(all_data$phenos, accession_id %in% IDs_test)
  
  #sanity
  stopifnot(all(rownames(res$kinship_train) %in% res$phenos_train$accession_id))
  stopifnot(all(rownames(res$kinship_test) %in% res$phenos_test$accession_id))
  
  #same order
  res$kinship_train = res$kinship_train[res$phenos_train$accession_id,]
  res$kinship_test = res$kinship_test[res$phenos_test$accession_id,]
  
  #done
  return(res)
}

# ACTUAL SCRIPT -----------------------------------------------------------
for (trait in trait_list){
  #loading data from the full SNP dataset
  all_data = load_all_703(trait = trait)
  
  #looping over all available splits
  for (split in unique(all_data$splits$split)){
    #we are doing by-hand special fold split 
    split_dataset = get_split_data(all_data, split)
    
    #kinship based, full set
    wb = createWorkbench(outfolder = outfolder, outfile.name = outfile, folds = NULL, regressor.name = 'KIN GBLUP')
    wb = addRegressor(wb, regressor.name = 'KIN RKHS', regressor = phenoRegressor.BGLR, type='RKHS', nIter = nIter)
    nds_train = createNoisyDataset(name = paste(sep='', trait, ' split=', split, ' train'), covariance = split_dataset$kinship_train, phenotypes = split_dataset$phenos_train$phenotype_value)
    nds_test = createNoisyDataset(name = paste(sep='', trait, ' split=', split, ' test'), covariance = split_dataset$kinship_test, phenotypes = split_dataset$phenos_test$phenotype_value)
    tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
  }
  
}


