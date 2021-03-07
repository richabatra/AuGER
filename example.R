# start fresh
rm(list = ls())
setwd('path_to_your_working_dir')
source('auger.R')
load('testdata100.rds')
# For finding associations between molecular profiles and clinical profiles with multiple imputations of clinical variables
test_mult <- auger_multi_impute(
  molprofiles=moldata, # dataframe with samples in rows and molecular variable in columns
  list_clinprofiles=list_clindata, # list of dataframes with samples in rows and clinical variable in columns
  clin_datatype=clindatatype # dataframe with name and datatype of clinical variables
)
# For finding associations between molecular profiles and clinical profiles for both cases and controls
test_cc <- auger_casecontrol(
  case_molprofiles=moldata, # dataframe with case samples in rows and molecular variable in columns
  case_clinprofiles=clindata, # dataframe with case samples in rows and clinical variable in columns
  ctrl_molprofiles=ctrl_moldata, # dataframe with control samples in rows and molecular variable in columns
  ctrl_clinprofiles=clindata, # dataframe with control samples in rows and clinical variable in columns
  clin_datatype=clindatatype # dataframe with name and datatype of clinical variables
  )
# For finding associations between molecular profiles and clinical profiles. It is called by auger_casecontrol and auger_multi_impute functions.
test <- auger (
  molprofiles=moldata, # dataframe with samples in rows and molecular variable in columns
  clinprofiles=clindata, # dataframe with samples in rows and clinical variable in columns
  clin_datatype=clindatatype, # dataframe with name and datatype of clinical variables
  save_each_model=T # boolean to save each model?
)
# To plot the associations in the network format using the association dataframes.
get_network (asso_df=test$asso_df, cytoscape=F, network_title = 'test') 
  