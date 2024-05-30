library(GROAN)

#target traits
trait_list = c('Cauline Leaf number (CL)', 'DaysToFlowering1', 'DaysToFlowering2', 'DaysToFlowering3', 'Rosette Leaf number (RL)')

#setup, data management library
basefolder = '/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/'
source(file.path(basefolder, 'scripts', 'data_manager.R'))

#room for results
outfolder = file.path(basefolder, 'results', 'test_regression_results')
outfile = 'accuracy_regression2_coordinated_CV_SNP.csv'
dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

#Bayesian regressors
nIter = 15000

# DATA LOAD ---------------------------------------------------------------
for (trait in trait_list){
  #loading data from the full SNP dataset
  all_data = load_all_703(trait = trait, load_SNPs = TRUE)
  
  #looping over all available splits
  for (split_current in unique(all_data$splits$split)){
    
    #samples for train and test set
    samples_train = subset(all_data$splits, split == split_current & train)$accession_id
    samples_test = subset(all_data$splits, split == split_current & !train)$accession_id
    
    #kinship based, full set
    wb = createWorkbench(outfolder = outfolder, outfile.name = outfile, folds = NULL, regressor.name = 'rrBLUP', regressor = phenoRegressor.rrBLUP, saveExtraData = TRUE, reps = 1)
    nds_train = createNoisyDataset(name = paste(sep='', trait, ' split=', split_current, ' train'), genotypes = all_data$SNP[samples_train,], phenotypes = all_data$phenos[samples_train, 'phenotype_value'])
    nds_test = createNoisyDataset(name = paste(sep='', trait, ' split=', split_current, ' test'),   genotypes = all_data$SNP[samples_test,],  phenotypes = all_data$phenos[samples_test, 'phenotype_value'])
    tmp = GROAN.run(nds=nds_train, nds.test=nds_test, wb=wb)
  }
}