library(GROAN)

#target traits
trait_list = c('Cauline Leaf number (CL)', 'DaysToFlowering1', 'DaysToFlowering2', 'DaysToFlowering3', 'Rosette Leaf number (RL)')

#setup, data management library
basefolder = '/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/'
source(file.path(basefolder, 'scripts', 'data_manager.R'))

#room for results
outfolder = file.path(basefolder, 'results')
dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

#Bayesian regressors
nIter = 15000

#pick the desired computations
DO_SNP_BASED_REGRESSIONS = FALSE
DO_KINSHIP_BASED_REGRESSIONS = TRUE

# SNP BASED REGRESSIONS ---------------------------------------------------
if (DO_SNP_BASED_REGRESSIONS){
  for (trait in trait_list){
    #load data
    all_data = load_all(basefolder=file.path(basefolder, 'data'), trait=trait)
    
    #SNP based, crossvalidation
    wb = createWorkbench(outfolder = outfolder, regressor.name = 'rrBLUP')
    wb = addRegressor(wb, regressor.name = 'BL', regressor = phenoRegressor.BGLR, type='BL', nIter = nIter)
    wb = addRegressor(wb, regressor.name = 'BayesC', regressor = phenoRegressor.BGLR, type='BayesC', nIter = nIter)
    nds_train = createNoisyDataset(name = paste(sep='', trait, 'CVtrain'), genotypes = all_data$genos_train, phenotypes = all_data$phenos_train$phenotype_value)
    nds_test = createNoisyDataset(name = paste(sep='', trait, 'CVtest'), genotypes = all_data$genos_test, phenotypes = all_data$phenos_test$phenotype_value)
    tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
    
    #SNP based, full set
    wb = createWorkbench(outfolder = outfolder, folds = NULL, regressor.name = 'rrBLUP')
    wb = addRegressor(wb, regressor.name = 'BL', regressor = phenoRegressor.BGLR, type='BL', nIter = nIter)
    wb = addRegressor(wb, regressor.name = 'BayesC', regressor = phenoRegressor.BGLR, type='BayesC', nIter = nIter)
    nds_train = createNoisyDataset(name = paste(sep='', trait, 'FULLtrain'), genotypes = all_data$genos_train, phenotypes = all_data$phenos_train$phenotype_value)
    nds_test = createNoisyDataset(name = paste(sep='', trait, 'FULLtest'), genotypes = all_data$genos_test, phenotypes = all_data$phenos_test$phenotype_value)
    tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
  }
}

# KINSHIP BASED REGRESSIONS -----------------------------------------------
if(DO_KINSHIP_BASED_REGRESSIONS){
  for (trait in trait_list){
    #load data
    all_data = load_all(basefolder=file.path(basefolder, 'data'), trait=trait)
    
    #kinship based, crossvalidation
    wb = createWorkbench(outfolder = outfolder, regressor.name = 'KIN GBLUP')
    wb = addRegressor(wb, regressor.name = 'KIN RKHS', regressor = phenoRegressor.BGLR, type='RKHS', nIter = nIter)
    nds_train = createNoisyDataset(name = paste(sep='', trait, 'CVtrain'),   covariance = all_data$kinship_train, phenotypes = all_data$phenos_train$phenotype_value)
    nds_test = createNoisyDataset(name = paste(sep='', trait, 'CVtest'), covariance = all_data$kinship_test, phenotypes = all_data$phenos_test$phenotype_value)
    tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
    
    #kinship based, full set
    wb = createWorkbench(outfolder = outfolder, folds = NULL, regressor.name = 'KIN GBLUP')
    wb = addRegressor(wb, regressor.name = 'KIN RKHS', regressor = phenoRegressor.BGLR, type='RKHS', nIter = nIter)
    nds_train = createNoisyDataset(name = paste(sep='', trait, 'FULLtrain'), covariance = all_data$kinship_train, phenotypes = all_data$phenos_train$phenotype_value)
    nds_test = createNoisyDataset(name = paste(sep='', trait, 'FULLtest'), covariance = all_data$kinship_test, phenotypes = all_data$phenos_test$phenotype_value)
    tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
  }
  
}
