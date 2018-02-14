
#' Title
#'
#' This function is used to covert the dataset to a contingincy table.
#'
#' @param data the raw dataset
#'
#'
#'
#'
#'



convert_ct = function(data)
{
  ct = data
  ct[ct!=0]=1
  return(ct)
}

prelim = function(data) ##where data is a data frame which as least has 2 column
{
  name = colnames(data)
  num =length(name)
  for(i in 2:num)
  {
    data = cbind(data,ifelse(data[name[i]]!=0,1,0))
    colnames(data)[ncol(data)]=paste("I",name[i],sep = "_")
  }
  name = colnames(data)


  ######
  #Ind = data[,(num+1):ncol(data)]
  #Ind = cbind(Ind,1)
  #colnames(Ind)[ncol(Ind)]="num"
  #Ind = aggregate(num~.,data=Ind,FUN = "sum")
  #m = Ind$num
  ####
  #p = rowSums(Ind[1:(ncol(Ind)-1)])+1
  #reg = lm(as.formula(paste(name[1],"~.",sep = "")),data = data)
  #res = cbind(residuals(reg)^2,data[(num+1):ncol(data)])
  #colnames(res)[1]="res"
  #resi = aggregate(res~.,data = res,FUN = "sum")
  #sigma2 = resi$res/(m-p)
  #sigma2[m<=p]=sigma2[m>p]%*%m[m>p]/sum(m)
  #################
  ## with interaction

  m = nrow(data)
  reg = lm(as.formula(paste(paste(name[1],"~",paste(name[2:num],collapse = "+")),"+",paste(name[(num+1):length(name)],collapse = "*"))),data=data)

  #combination = combn(length(name)-num,2)
  #f1 = paste(name[1],"~",paste(name[2:num],collapse = "+"))
  #name_ind = name[(num+1):length(name)]
  #f2 = paste(apply(combination,2,FUN = function(x){paste(name_ind[x[1]],":",name_ind[x[2]],sep = "")}),collapse = "+")
  #reg = lm(as.formula(paste(f1,"+",f2,sep = "")),data = data)

  p = length(reg$coefficients)
  sigma2 = sigma(reg)^2





  ####
  #sigma2 = rep(0,length(m))
  #p = rep(0,length(m))
  #for(i in 1:length(m))
  #{
  #  indicator = Ind[i,1:(ncol(Ind)-1)]
  #  p[i] = 1 + sum(indicator==1)
  #  d_matrix = as.data.frame(selectdata_prelim(data,as.numeric(indicator),m[i]))
  #
  #  linear_reg = lm(as.formula(paste(name[1],"~.",sep = "")),d_matrix[c(T,indicator==1)])
  #  sigma2[i] = sigma(linear_reg)^2
  #  if(nrow(d_matrix)<=(sum(indicator)+1))
  #  {
  #    sigma2[i] = 1
  #  }
  #}

  return(list(m = m, sigma2 = sigma2,p=p))
}



local_score = function(x,ct=NULL,data=NULL,num_p,type = c("discrete","gaussian","zero_inflated"))
{
  if(is.null(ct))
  {
    name = setdiff(colnames(data),x)
  }
  else{name = setdiff(colnames(ct),x)}

  n = rep(0,num_p+1)
  index = NULL

  for( j in 0:num_p)
  {
    te = apply(combn(length(name),j),2,function(x){t = ifelse(1:length(name) %in% x, 1,0);paste(t,collapse = "")})
    n[j+1] = length(te)
    index = c(index,te)
  }

  result = matrix(0,nrow=sum(n),ncol = 2)
  result = data.frame(result)
  colnames(result) = c("parents","score")

  result[,1] = index

  if(type =="discrete")
  {
    for(i in 1:sum(n))
    {

      indicator = (unlist(strsplit(result[i,1],""))==1)
      variables = name[indicator]
      logit = glm(as.formula(paste(x,"~.")),family = binomial,data = ct[c(x,variables)])
      prob = fitted(logit)
      fit_prob = ifelse(ct[x]==1,prob,1-prob)
      result[i,2] = sum(log(fit_prob))-(length(variables)+1)/2*log(length(prob))
    }
  }

  else if(type =="gaussian")
  {
    for( i in 1:sum(n))
    {
      indicator = (unlist(strsplit(result[i,1],""))=="1")
      variables = name[indicator]
      l = lm(as.formula(paste(x,"~.")),data = data[c(x,variables)])
      sigma2 = sigma(l)^2
      m = nrow(data)
      result[i,2] = -m/2*log(2*pi*sigma2)-1/2*(m-length(variables)-1)-(length(variables)+1)/2*log(m)
    }
  }
  else if(type=="zero_inflated")
  {
    for(i in 1:sum(n))
    {
      indicator = (unlist(strsplit(result[i,1],""))==1)
      variables = name[indicator]
      logit = glm(as.formula(paste(x,"~.")),family = binomial,data = ct[c(x,variables)],model = F,x=F)
      prob = fitted(logit)
      fit_prob = ifelse(ct[x]==1,prob,1-prob)
      result[i,2] = sum(log(fit_prob))-(length(variables)+1)/2*log(length(prob))


      design = data[c(x,variables)][data[x]!=0,]
      if(length(variables)==0&&length(design)>0)
      {
        u = mean(design)
        variance = var(design)
        num_non_zero = length(design)
        result[i,2] = result[i,2] -num_non_zero/2*log(2*pi*variance)-0.5*(num_non_zero-1) -1/2*log(num_non_zero)
      }

      else if(length(variables)>0&&nrow(design)>0)
      {
        param = prelim(design)
        num_non_zero = nrow(design)
        result[i,2] = result[i,2]- 0.5*num_non_zero*log(2*pi)-0.5*param$m%*%log(param$sigma2)-0.5*(sum(param$m)-sum(param$p))
        result[i,2] = result[i,2]- (param$p)%*%log(param$m)/2
      }
    }
  }
  return(list(table = result,parent_set =name,obj = x ))
}


get_all_local_score_parallel = function(ct=NULL,data = NULL,num_p,type=c("discrete","gaussian","zero_inflated"),clustersize)
{
  if(type!="discrete"&&type!="gaussian"&&type!="zero_inflated")
  {
    return("type is wrong")
  }

  if(is.null(data))
  {
    name = colnames(ct)
  }
  else
  {
    name = colnames(data)
  }
  ncols = length(name)
  cluster = makeCluster(clustersize,type = "SOCK")
  clusterEvalQ(cluster,Rcpp::sourceCpp('cpp_code3.cpp'))
  clusterExport(cluster,list = list("prelim"))
  result = clusterApply(cluster,name[1:ncols],local_score,ct,data,num_p,type)
  return(result)
  stopCluster(cluster)
}







best_parents = function(x,all_names,num_p) ## score_all is need
{
  variables = setdiff(all_names,x)
  for(i in 1:length(all_names))
  {
    if(score_all[[i]]$obj==x)
    {
      score_file = score_all[[i]]$table
    }
  }
  result = comparison_bps(score_file,num_p)
  return(list(table = result,parent_set = variables,obj = x))
}

best_parents_all = function(all_names,num_p,clustersize) ## score_all is need
{
  cluster = makeCluster(clustersize,type="SOCK")
  clusterExport(cluster,list("score_all"))
  clusterEvalQ(cluster,Rcpp::sourceCpp('cpp_code3.cpp'))
  result = clusterApply(cluster,all_names[1:length(all_names)],best_parents,all_names,num_p)
  return(result)
  stopCluster(cluster)
}



get_sinks = function(W,all_names) ## sinks_all is needed
{

  indicator = as.character(ifelse(all_names%in%W,1,0))
  indicator = paste(indicator,collapse = "")

  int = Tode(indicator)

  score = sinks_all[int,3]
  sink = sinks_all[int,2]
  return(list(sink = sink,score = score))
}

get_ordering = function(all_names)
{
  names = all_names
  ord = rep(0,length(all_names))
  for(i in length(all_names):1)
  {
    ord[i] = get_sinks(names,all_names)$sink
    names = setdiff(names,ord[i])
  }
  return(ord)
}

#' Title
#'
#' This function is used to find the best bayesian network for dicrete, gaussian and zero-inflated data.
#'
#' @param ct the contingincy table converted from data
#' @param data dataset
#' @param num_p the maximum number of parents
#' @param type a character representing the type of your data. Choosing from "discrete","gausssian","zero_inflated"
#' @param clustersize the number of parallel
#' @details
#' This function is used to find the best bayesian network for dicrete, gaussian and zero-inflated data by using dynamic programming.
#'
#' @return
#' @param DAG represents the best bayesian network
#' @param score represents the score of the bayesian network
#'
#' @author Yejiong Zhu
#'
#'
#'
dag = function(ct = NULL,data = NULL,num_p,type=c("discrete","gaussian","zero_inflated"),clustersize)
{
  if(is.null(ct))
  {
    all_names = colnames(data)
  }
  else
  {
    all_names = colnames(ct)
  }
  score_all <<- get_all_local_score_parallel(ct,data,num_p = num_p,type = type,clustersize = clustersize)
  bps_all <<- best_parents_all(all_names = all_names,num_p = num_p,clustersize = clustersize)
  sinks_all <<- sinks(all_names = all_names)
  result = list()
  W = all_names
  for(i in 1:length(all_names))
  {
    sink = get_sinks(W,all_names = all_names)$sink
    W = setdiff(W,sink)
    k = which(all_names%in%sink)
    result[[k]] = get_best_parents(sink,W,all_names)$bps
  }
  score = get_sinks(all_names,all_names)$score
  return(list(DAG = result,score = score))
}


simulation_gaussian = function(sample_size,num)
{
  result = rep(0,num)
  for(i in 1:sample_size)
  {
    x = rep(0,num)
    x[1] = rnorm(1)
    for(j in 2:num)
    {
      x[j] = rnorm(1,x[j-1],j)
    }
    result = rbind(result,x)
  }
  return(data.frame(result[2:(sample_size+1),]))
}
