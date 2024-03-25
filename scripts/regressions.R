library(GROAN)

#target traits
trait_list = c('Cauline Leaf number (CL)', 'DaysToFlowering1', 'DaysToFlowering2', 'DaysToFlowering3', 'Rosette Leaf number (RL)')

#setup, data management library
basefolder = '/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/'
source(file.path(basefolder, 'scripts', 'data_manager.R'))

#room for results
outfolder = file.path(basefolder, 'results')
dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

# REGRESSIONS -------------------------------------------------------------
for (trait in trait_list){
  #load data
  all_data = load_all(basefolder=file.path(basefolder, 'data'), trait=trait)
  
  #SNP based, crossvalidation
  wb = createWorkbench(outfolder = outfolder)
  nds_train = createNoisyDataset(name = paste(sep='', trait, 'CVtrain'), genotypes = all_data$genos_train, phenotypes = all_data$phenos_train$phenotype_value)
  nds_test = createNoisyDataset(name = paste(sep='', trait, 'CVtest'), genotypes = all_data$genos_test, phenotypes = all_data$phenos_test$phenotype_value)
  tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
  
  #SNP based, full set
  wb = createWorkbench(outfolder = outfolder, folds = NULL)
  nds_train = createNoisyDataset(name = paste(sep='', trait, 'FULLtrain'), genotypes = all_data$genos_train, phenotypes = all_data$phenos_train$phenotype_value)
  nds_test = createNoisyDataset(name = paste(sep='', trait, 'FULLtest'), genotypes = all_data$genos_test, phenotypes = all_data$phenos_test$phenotype_value)
  tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
}

