# computes associations between genes/proteins/variables and clinical/patient attributes
# author: Richa Batra rib4003@med.cornell.edu

library(glmnet)
library(tidyverse)
library(magrittr)
library(RCy3)
library(visNetwork)
library(igraph)

# For computing coefficients between molecular profiles and clinical profiles. It is called by auger function.
# @param molprofiles: a dataframe with samples in rows and molecular variable in columns
# @param clinprofiles: a dataframe with samples in rows and clinical variable in columns
# Rows of the above two dataframes should represent the same sample
# Columns names should be set to identify the molecular and clinical variables
# @param clin_datatype: a dataframe with clinical variable in rows and two columns 'name' - name of clinial variable which should name the column name of clinprofiles
# 'datatype' - 'ordinal' / 'categorical', 'numeric', 'binary' representing datatype of the clinical variable
# @param clin_cov: boolean - 'T' or 'F' to correct for other clinical variables in glmnet or not, default is T
# @param alpha_value: level of regularization to be used in glmnet default is set to 0.5
# @param cv_folds: cross validation folds to be used in glmnet default is set to 10
# @param opt_lambda: optimized model to select - options are '1se' or 'min, default is set to '1se'
# @param save_each_model: boolean - 'T' or 'F' to save each glmnet model or not, default is F
# @param seed_value: seed to set before glmnet, default is set to 809
# @returns: a list for each clinical variable with (a) models (if opted for) (b) coefficients 
get_coefficients <- function(molprofiles, # dataframe with samples in rows and molecular variable in columns
                             clinprofiles, # dataframe with samples in rows and clinical variable in columns
                             clin_datatype, # dataframe with name and datatype of clinical variables
                             clin_cov=T, # boolean to correct for other clinical variables in glmnet or not
                             alpha_value=0.5, # regularization parameter to be used in glmnet
                             cv_folds=10, # crossvalidation fold to be used in glmnet
                             opt_lambda='1se', # optimized model to be used '1se' or 'min'
                             save_each_model=F, # boolean to save each model?
                             seed_value=809) # seed to set before glmnet
                             {
  
  # is alpha in acceptable range?
  if(alpha_value>1 | alpha_value<0){
    warning('Alpha out of bounds, using alpha 0.5')
    alpha_value <- 0.5
  }

# loop over clinical variables
all_coefs_list <- lapply(1:nrow(clin_datatype), FUN=function(i){
  
  # value of clinical variable to be used as response
  this_y <- clinprofiles %>% select(clin_datatype$name[i])
  # value of clinical variables to be used as covariates
  this_cov <- clinprofiles %>% select(-clin_datatype$name[i])
  # type of clinical variable
  this_type <- clin_datatype$datatype[i]
  # family type based on the datatype of clinical variable
  if(this_type%in%"numeric"){
    fam_type <- "gaussian"
    this_y <- as.numeric(as.matrix(this_y))
  } else if(this_type%in%"binary"){
    fam_type <- "binomial"
    this_y <- as.factor(as.matrix(this_y))
  } else if(this_type%in%"ordinal" | this_type%in%"categorical"){
    fam_type <- "multinomial"
    this_y <- as.factor(as.matrix(this_y))
  }
  
  # correct for rest of the clinical variables?
  if(clin_cov) {
    this_x<- data.frame(cbind(this_cov, molprofiles))
  } else{
    this_x<- molprofiles
  }
  
  # set seed
  set.seed(seed_value)
  # model matrix
  mod_mat <- model.matrix(~0+.,this_x)
  this_fit <- this_coefs <- NULL
  # try in case it fails for one variable
  try({
    # glmnet fit 
    this_fit <- cv.glmnet(x=mod_mat, y=this_y, alpha=alpha_value, family=fam_type, nfolds=cv_folds)
    
    # which opt model?
    if(opt_lambda=='1se'){
      this_coefs <- coef(this_fit, this_fit$lambda.1se)
    } else if(opt_lambda=='min'){
      this_coefs <- coef(this_fit, this_fit$lambda.min)
    } else {
      warning('Unrecognized opt_lambda, returning values of  "1se"') 
      this_coefs <- coef(this_fit, this_fit$lambda.1se)
    } 
  })
  # what all to return?
  if(save_each_model){
    res <- list(model=this_fit, coefs=this_coefs)
  } else{
    res <- list(model=NULL, coefs=this_coefs)
  }
  return(res)
})
# name the list after the clinical variables
names(all_coefs_list) <- clin_datatype$name
return(all_coefs_list)
}

# Extract coefficients from the output of get_coefficients function. It is called by get_assodf function.
# @param all_coefs_list: a list of the type get_coefficients outputs
# @returns coefs_list: a list for each clinical variable in the input list with only nonzero variables
extract_nonzero_coefs <- function(all_coefs_list) {
  # loop over the clinical variables
  coefs_list <- lapply(all_coefs_list, FUN=function(x) {
    # dataframe of coefficients
    if(mode(x$coefs)=='list'){
      y <- do.call(cbind, x$coefs) %>% data.matrix() %>% data.frame()
      names(y) <- paste0('coef', names(x$coefs))
    } else{
      y <- x$coefs %>% data.matrix() %>% data.frame()
      # name the columns
      names(y) <- 'coef'
    }
    # add sum of coefficients as a column and # add rownames as columns
    y %<>% mutate(sumcoefs = rowSums(abs(.)), variables=rownames(y)) %>%
          # filter out non zero variables
          filter(sumcoefs > 0)
  })
  return(coefs_list)
}
# Generate association data frame from the output of get_coefficients function. It is called by auger and auger_casecontrol functions
# @param all_coefs_list: a list of the type get_coefficients outputs
# @returns asso_df: a dataframe of associations between molecules and clinical variables
get_assodf <- function(all_coefs_list) {
  # extract non zero coefs
  coefs_list <- extract_nonzero_coefs (all_coefs_list) 
  # loop over the clinical variables in coefs_list
  edge_list <- lapply(1:length(coefs_list), FUN=function(i) {
        # remove the sumcoefs column 
    y <- coefs_list[[i]] %>% select(-sumcoefs) %>% 
      # reformat the data to get long format
      reshape2::melt(id='variables') %>% 
      # add the name of clinical variable to the variable levels
      mutate(res_var=paste0(names(coefs_list)[i], '_', variable),
        # remove 'coef' from the name
        res_var= sub('coef', '', res_var), res_var = sub('_$', '', res_var)) %>% 
      filter(!variables%in% '(Intercept)') %>%
      # filter non zero and select relevant columns
      filter(value!=0) %>% select(res_var, variables, value) %>% 
      # rename columns
      dplyr::rename(pred_var=variables, coef=value)
    return(y)
  })
  # combine the lists in a dataframe
  asso_df <- do.call(rbind, edge_list)
  return(asso_df)
}

# For finding associations between molecular profiles and clinical profiles. It is called by auger_casecontrol and auger_multi_impute functions.
# @param molprofiles: a dataframe with samples in rows and molecular variable in columns
# @param clinprofiles: a dataframe with samples in rows and clinical variable in columns
# Rows of the above two dataframes should represent the same sample
# Columns names should be set to identify the molecular and clinical variables
# @param clin_datatype: a dataframe with clinical variable in rows and two columns 'name' - name of clinial variable which should name the column name of clinprofiles
# 'datatype' - 'ordinal' / 'categorical', 'numeric', 'binary' representing datatype of the clinical variable
# @param clin_cov: boolean - 'T' or 'F' to correct for other clinical variables in glmnet or not, default is T
# @param alpha_value: level of regularization to be used in glmnet default is set to 0.5
# @param cv_folds: cross validation folds to be used in glmnet default is set to 10
# @param opt_lambda: optimized model to select - options are '1se' or 'min, default is set to '1se'
# @param save_each_model: boolean - 'T' or 'F' to save each glmnet model or not, default is F
# @param seed_value: seed to set before glmnet, default is set to 809
# @returns: a list for each clinical variable with (a) models (if opted for) (b) dataframe of associations between molecules and clinical variables
auger <- function(
                  molprofiles, # dataframe with samples in rows and molecular variable in columns
                  clinprofiles, # dataframe with samples in rows and clinical variable in columns
                  clin_datatype, # dataframe with name and datatype of clinical variables
                  clin_cov=T, # boolean to correct for other clinical variables in glmnet or not
                  alpha_value=0.5, # regularization parameter to be used in glmnet
                  cv_folds=10, # crossvalidation fold to be used in glmnet
                  opt_lambda='1se', # optimized model to be used '1se' or 'min'
                  save_each_model=F, # boolean to save each model?
                  seed_value=809 # seed to set before glmnet
                  ) {
  
  # get coeffiecients
  all_coefs_list <- get_coefficients (molprofiles=molprofiles, 
                 clinprofiles=clinprofiles, 
                 clin_datatype=clin_datatype, 
                 clin_cov=clin_cov, 
                 alpha_value=alpha_value, 
                 cv_folds=cv_folds, 
                 opt_lambda=opt_lambda, 
                 save_each_model=save_each_model, seed_value = seed_value)
  # extract associations from only non zero coefficients
  asso_df <- get_assodf (all_coefs_list)
  # extract models
  mod_list <- lapply(all_coefs_list, FUN=function(x) x$model)
  
  # what all to output?
  if(save_each_model){
    out <- list(mod_list=mod_list, asso_df=asso_df)
  } else{
    out <- asso_df
  }
  return(out)
}

# Compare the associations in cases and controls for the same set of 
# clinical variables using the output from auger function. It is called by auger_casecontrol function.
# @param case_asso: output of auger function only the asso_df for cases
# @param ctrl_asso: output of auger function only the asso_df for controls
# @return cc_asso: combined dataframe of the above two with information on the associations
compare_case_control <- function(case_asso, ctrl_asso) {
  # create unique assocation identifier and add sample type information
  case_asso <- case_asso %>% mutate(uid=paste0(res_var, '(asso with)', pred_var), sample_type='cases')
  ctrl_asso <- ctrl_asso %>% mutate(uid=paste0(res_var, '(asso with)', pred_var), sample_type='controls')
  # intersection of associations between the two sample types
  commons <- intersect(case_asso %>% pull(uid), ctrl_asso %>% pull(uid))
  # add assocation information in the dataframes
  case_asso %<>% mutate(asso_info = ifelse(uid%in% commons, 'also associate in control', 'unique to cases'))
  ctrl_asso %<>% mutate(asso_info = ifelse(uid%in% commons, 'also associate in cases', 'unique to controls'))
  # combine the two data frames
  cc_asso <- bind_rows(case_asso, ctrl_asso) %>% .[order(.$uid), ]
  return(cc_asso)
}
# For finding associations between molecular profiles and clinical profiles for both cases and controls
# @param case_molprofiles: a dataframe with case samples in rows and molecular variable in columns
# @param case_clinprofiles: a dataframe with case samples in rows and clinical variable in columns
# Rows of the above two dataframes should represent the same sample
# Columns names should be set to identify the molecular and clinical variables
# @param ctrl_molprofiles: optional - a dataframe with control samples in rows and molecular variable in columns
# @param ctrl_clinprofiles: optional - a dataframe with control samples in rows and clinical variable in columns
# Rows of the above two dataframes should represent the same sample
# Columns names should be set to identify the molecular and clinical variables
# @param clin_datatype: a dataframe with clinical variable in rows and two columns 'name' - name of clinial variable which should name the column name of clinprofiles
# 'datatype' - 'ordinal' / 'categorical', 'numeric', 'binary' representing datatype of the clinical variable
# @param clin_cov: boolean - 'T' or 'F' to correct for other clinical variables in glmnet or not, default is T
# @param alpha_value: level of regularization to be used in glmnet default is set to 0.5
# @param cv_folds: cross validation folds to be used in glmnet default is set to 10
# @param opt_lambda: optimized model to select - options are '1se' or 'min, default is set to '1se'
# @param seed_value: seed to set before glmnet, default is set to 809
# @returns: a dataframe of associations from both cases and controls
auger_casecontrol <- function(
  case_molprofiles, # dataframe with case samples in rows and molecular variable in columns
  case_clinprofiles, # dataframe with case samples in rows and clinical variable in columns
  ctrl_molprofiles, # dataframe with control samples in rows and molecular variable in columns
  ctrl_clinprofiles, # dataframe with control samples in rows and clinical variable in columns
  clin_datatype, # dataframe with name and datatype of clinical variables
  clin_cov=T, # boolean to correct for other clinical variables in glmnet or not
  alpha_value=0.5, # regularization parameter to be used in glmnet
  cv_folds=10, # crossvalidation fold to be used in glmnet
  opt_lambda='1se', # optimized model to be used '1se' or 'min'
  save_each_model=F, # boolean to save each model?
  seed_value=809 # seed to set before glmnet
) {
  
# case associations
 case_assodf <- auger (
    molprofiles=case_molprofiles, # dataframe with samples in rows and molecular variable in columns
    clinprofiles=case_clinprofiles, # dataframe with samples in rows and clinical variable in columns
    clin_datatype=clin_datatype, # dataframe with name and datatype of clinical variables
    clin_cov=clin_cov, # boolean to correct for other clinical variables in glmnet or not
    alpha_value=alpha_value, # regularization parameter to be used in glmnet
    cv_folds=cv_folds, # crossvalidation fold to be used in glmnet
    opt_lambda=opt_lambda, # optimized model to be used '1se' or 'min'
    save_each_model=F, # boolean to save each model?
    seed_value=seed_value)
 # control associations
 ctrl_assodf <- auger (
   molprofiles=ctrl_molprofiles, # dataframe with samples in rows and molecular variable in columns
   clinprofiles=ctrl_clinprofiles, # dataframe with samples in rows and clinical variable in columns
   clin_datatype=clin_datatype, # dataframe with name and datatype of clinical variables
   clin_cov=clin_cov, # boolean to correct for other clinical variables in glmnet or not
   alpha_value=alpha_value, # regularization parameter to be used in glmnet
   cv_folds=cv_folds, # crossvalidation fold to be used in glmnet
   opt_lambda=opt_lambda, # optimized model to be used '1se' or 'min'
   save_each_model=F, # boolean to save each model?
   seed_value=seed_value)
 # compare the assocations
 cc_assodf <- compare_case_control(case_asso=case_assodf, ctrl_asso=ctrl_assodf)
 return(cc_assodf)
}

# To plot the associations in the network format using the association dataframes.
# Note: can be slow or crash depending on the size of the network
# @param asso_df: a association dataframe similar to the output from auger
# @param cytoscape: boolean - can be set to T if you have cytoscape running in the machine default is F
# @returns a network: either in cytoscape and/or outputs a html using visnetwork.
get_network <- function(asso_df, cytoscape=F, network_title='net') {
  net <- asso_df %>%
    # rename the columns
    dplyr::rename(source = res_var, target = pred_var, edge_wt =coef) %>% 
    # generate the column edge type based on the sign of coef
    mutate(edge_type = case_when(edge_wt > 0 ~ 'Positive', TRUE ~ 'Negative'), 
           dashes = case_when(edge_wt > 0 ~ TRUE, TRUE ~ FALSE))
  # node attributes based on net and a prior information about node types in case of a hybrid network
  node_attributes <- data.frame(node_name=unique(c(net$source, net$target))) %>% 
    mutate(group = case_when(node_name%in%net$target ~ "Predictor", 
                             node_name%in%net$source ~ "Response"), label = node_name)
  # create numeric node ids
  node_attributes %<>% mutate(id=1:nrow(node_attributes))
  # add the same numeric node ids to net
  net <- node_attributes %>% dplyr::select(id, node_name) %>% 
    left_join(net,.,by=c("source"= "node_name")) %>% dplyr::rename(from=id)
  net <- node_attributes %>% dplyr::select(id, node_name) %>% 
    left_join(net,.,by=c("target"= "node_name")) %>% dplyr::rename(to=id)
  # edge legend for visnetwork
  ledges <- data.frame(color = c("#17202A"), label = c("Negative", "Positive"), dashes =c(TRUE, FALSE))
  # visNetwork object
  visNetwork(nodes=node_attributes, edges = net, main='Auger based associations') %>%
    visIgraphLayout(layout='layout.davidson.harel', physics = F, smooth = F) %>%
    visPhysics(stabilization = FALSE) %>% 
    # node shapes based on node types defined in node_attributes
    visGroups(groupname = "Response", color = "#FADBD8", shape = "diamond") %>% 
    visGroups(groupname = "Predictor", color = "#D6EAF8", shape = "dot") %>% 
    # edge legend as defined above
    visLegend(addEdges = ledges, useGroups = T) %>%
    visNodes(borderWidth = 3) %>% 
    visEdges(smooth = FALSE, shadow = TRUE, color = "#ABB2B9") %>%
    visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group")%>%
    visInteraction(navigationButtons = T) %>%
    visConfigure(enabled=T) %>%
    visSave(file=paste0(network_title, '.html'), selfcontained = TRUE, background = "#FDFEFE")
  # print in cytoscape?
  if(cytoscape){

    # convert source and target to numeric values based on ids of nodes
    edge_list <- net %>% dplyr::rename(tmp1=source, tmp2=target, 
                                       source=from, target=to) %>% 
                          dplyr::rename(from=tmp1, to=tmp2)
    # data frame specification needed for cytoscape
    edge_list <- data.frame(edge_list, stringsAsFactors=FALSE)
    node_attributes <- data.frame(node_attributes, row.names = node_attributes$id, stringsAsFactors = F)
    # name styles of cytoscape networks based on the outcome_name
    this_style <- this_net <- network_title
    collection_name <- 'Auger based associations'

    # new graph object
    g <- graph_from_data_frame(edge_list[,c('source', 'target')], directed = F) %>% 
      igraph::simplify() # remove self loops
    # load the network to cytoscape
    createNetworkFromIgraph(igraph = g, title=this_net, collection=collection_name)
    # workaround for a bug in cytoscape
    edge_list %<>% mutate(key=paste(source, "(interacts with)", target))
    # get SUID for edges and match with edge attributes
    cy_edges <- getTableColumns(table = 'edge')
    cy_edges <- cy_edges[order(cy_edges$name), ]
    edge_list <- edge_list[order(edge_list$key), ]
    edge_list$cpSUID <- cy_edges$SUID
    
    # add edge attributes
    loadTableData(subset(edge_list, select = -c(source, target)), 
                  data.key.column = 'cpSUID', table.key.column = 'SUID', table = 'edge')
    # add node attributes
    loadTableData(node_attributes, data.key.column = 'id', table = 'node')
    #then prepare style variables
    style_name <- this_style
    nodeLabels <- mapVisualProperty('node label','node_name','p')
    nodeShapes <- mapVisualProperty('node shape','group','d',
                                    c('Response',  'Predictor'), c ('DIAMOND', 'ELLIPSE'))
    nodeFills <- mapVisualProperty('node fill color', 'group', 'd',
                                   c('Response',  'Predictor'), c ('#FADBD8', '#D6EAF8'))
    edgeStyles <- mapVisualProperty('edge line type', 'edge_type', 'd', 
                                    c("Negative", "Positive"), c("LONG_DASH", "SOLID"))
    #and then create the style
    createVisualStyle(style.name=style_name, base.url = 'http://localhost:1234/v1',
                      mappings = list(nodeLabels,nodeShapes,edgeStyles, nodeFills)) 
    #finish by applying the style
    setVisualStyle(style.name=style_name, base.url = 'http://localhost:1234/v1',
                   network = this_net)
    setEdgeColorDefault('#ABB2B9', style.name = this_style,base.url =  'http://localhost:1234/v1')
  }
}
# For finding associations between molecular profiles and clinical profiles with multiple imputations of clinical variables
# @param molprofiles: a dataframe with samples in rows and molecular variable in columns
# @param list_clinprofiles: a list of dataframes with samples in rows and clinical variable in columns
# Rows of the above two dataframes should represent the same sample
# Columns names should be set to identify the molecular and clinical variables
# @param clin_datatype: a dataframe with clinical variable in rows and two columns 'name' - name of clinial variable which should name the column name of clinprofiles
# 'datatype' - 'ordinal' / 'categorical', 'numeric', 'binary' representing datatype of the clinical variable
# @param clin_cov: boolean - 'T' or 'F' to correct for other clinical variables in glmnet or not, default is T
# @param alpha_value: level of regularization to be used in glmnet default is set to 0.5
# @param cv_folds: cross validation folds to be used in glmnet default is set to 10
# @param opt_lambda: optimized model to select - options are '1se' or 'min, default is set to '1se'
# @param seed_value: seed to set before glmnet, default is set to 809
# @returns: a dataframe of associations ranked by the frequency of dicovery in the imputations
auger_multi_impute <- function(
  molprofiles, # dataframe with samples in rows and molecular variable in columns
  list_clinprofiles, # list of dataframes with samples in rows and clinical variable in columns
  clin_datatype, # dataframe with name and datatype of clinical variables
  clin_cov=T, # boolean to correct for other clinical variables in glmnet or not
  alpha_value=0.5, # regularization parameter to be used in glmnet
  cv_folds=10, # crossvalidation fold to be used in glmnet
  opt_lambda='1se', # optimized model to be used '1se' or 'min'
  save_each_model=F, # boolean to save each model?
  seed_value=809 # seed to set before glmnet
) {
  # loop over list of clinprofiles
  assodf_list <- lapply(list_clinprofiles, FUN=function(clinprofiles){
    # associations
    asso_df <- auger (
      molprofiles=molprofiles, # dataframe with samples in rows and molecular variable in columns
      clinprofiles=clinprofiles, # dataframe with samples in rows and clinical variable in columns
      clin_datatype=clin_datatype, # dataframe with name and datatype of clinical variables
      clin_cov=clin_cov, # boolean to correct for other clinical variables in glmnet or not
      alpha_value=alpha_value, # regularization parameter to be used in glmnet
      cv_folds=cv_folds, # crossvalidation fold to be used in glmnet
      opt_lambda=opt_lambda, # optimized model to be used '1se' or 'min'
      save_each_model=F, # boolean to save each model?
      seed_value=seed_value)
    })
  # bind the replicates
  all_asso_df <- do.call(rbind, assodf_list) %>% 
    # unique assocation id
    mutate(uid = paste0(res_var, '(asso with)', pred_var))
  # frequency of discovery of association
  asso_freq <- all_asso_df %>% group_by(uid) %>% count()
  # coefs string per association
  coefs <- do.call(rbind, all_asso_df %>% group_by(uid) %>% group_rows()) %>% 
   apply(., 1, FUN=function(x) toString(all_asso_df$coef[unique(x)])) %>% data.frame(coefs=.)
  # bind asso_freq and coefs
  tmp <- bind_cols(asso_freq, coefs)
  # remove the coef column and duplicates
  out <- all_asso_df %>% left_join (tmp, by='uid') %>% select(-coef) %>% unique()
  return(out)
}
