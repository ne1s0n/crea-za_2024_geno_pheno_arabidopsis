setwd('/home/nelson/research/crea-za_2024_geno_pheno_arabidopsis/results/test_regression_results/extra_data/')

# SUPPORT FUNCTIONS -------------------------------------------------------
get_trait_name = function(x){
  i = grepl(pattern = 'split', x = x)
  res = x[i]
  trait = gsub(pattern = ' split.*', replacement = '', res)
  trait = gsub(pattern = 'GROAN_phenotypes_predicted_', replacement = '', trait)
  return(c(res, trait))
}

# ACTUAL SCRIPT -----------------------------------------------------------
res = data.frame()
for (f in list.files()){
  #load the file
  load(f)
  
  #extract the trait name
  tmp = get_trait_name(names(extradata))
  field = tmp[1]
  trait = tmp[2]
  
  res[names(extradata[[field]]), trait] = extradata[[field]]
}

write.csv(x = res, file = '../predictions.csv')