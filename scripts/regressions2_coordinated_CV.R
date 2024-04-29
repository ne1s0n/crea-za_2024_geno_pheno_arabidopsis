library(GROAN)

#target traits
trait_list = c('Cauline Leaf number (CL)', 'DaysToFlowering1', 'DaysToFlowering2', 'DaysToFlowering3', 'Rosette Leaf number (RL)')

#setup, data management library
basefolder = '/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/'
source(file.path(basefolder, 'scripts', 'data_manager.R'))

#room for results
outfolder = file.path(basefolder, 'results2_coordinated_CV')
dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

for (trait in trait_list){
  #loading data from the full SNP dataset
  all_data = load_all_705(trait = trait)
  
  #we are doing by-hand special fold split 
  stop('here')
}


