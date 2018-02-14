#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

std::string toBinary(int n,int k)
{
  std::string r;
  while(k!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;k=k-1;}
  return r;
}

CharacterVector split_string( std::string s)
{
  CharacterVector  result;
  for(int i=0;i<s.size();i++)
  {
    result.push_back(s[i]);
  }
  return(result);
}

String convert(CharacterVector s)
{
  String result;
  int n = s.size();
  for(int i=0;i<n;i++)
  {
    result+=s[i];
  }
  return(result);
}

CharacterVector all_subset(CharacterVector s)
{
  CharacterVector result;
  CharacterVector t;

  for(int i =0;i<s.size();i++)
  {
    if(s[i]=="1")
    {
      t = clone(s);
      t[i] = "0";
      result.push_back(convert(t));

    }
  }
  return(result);
}

//[[Rcpp::export]]
int Tode(String s)
{
  int result=0;
  CharacterVector t = split_string(s);
  int k = t.size()-1;
  for(int i =0;i<t.size();i++)
  {
    if(t[i]=="1")
    {
      result+= std::pow(2,k);
    }
    k=k-1;
  }
  return(result);
}

NumericVector Todecimal(CharacterVector s)
{
  NumericVector result(s.size());
  for(int i=0;i<s.size();i++)
  {
    result[i] = Tode(s[i]);
  }
  return(result);
}

//[[Rcpp::export]]
DataFrame selectdata_prelim(DataFrame df, NumericVector indicator,int m)
{
  StringVector name = df.names();
  int ncol = name.size();
  int num = ncol-indicator.size();
  NumericMatrix t(df.nrows(),ncol);
  for(int i=0;i<ncol;i++)
  {
    t(_,i)= NumericVector(df[i]);
  }
  NumericMatrix res(m,num);
  bool bo;
  int k=0;

  for(int j =0;j<df.nrows();j++)
  {
    bo=true;
    for(int i=0;i<indicator.size();i++)
    {
      if(t(j,num+i)!=indicator[i])
      {
        bo=false;
        break;
      }
    }
    if(bo)
    {
      for(int h=0;h<num;h++)
      {
        res(k,h) = t(j,h);
      }
      k++;
    }
  }
  DataFrame result;
  String s;
  for(int i=0;i<num;i++)
  {
    s = name[i];
    result.push_back(res(_,i),s);
  }

  return(result);
}






//[[Rcpp::export]]
DataFrame comparison_bps (DataFrame df,int num_p)
{
  CharacterVector parents=df[0];
  NumericVector score=df[1];

  CharacterVector t = split_string(Rcpp::as<std::string>(parents[0]));
  int num_variables = t.size();
  CharacterVector p(pow(2,num_variables));
  CharacterVector bps(p.size());
  NumericVector bps_score(p.size());
  for(int i=0;i<p.size();i++)
  {
    p[i] = toBinary(i,num_variables);
  }


  bps[0] = p[0];
  bps_score[0] = score[0];

  CharacterVector allset;
  NumericVector index;
  double sub_score;
  String sub_parents;

  for(int i = 1;i<p.size();i++)
  {
    sub_score = R_NegInf;
    if(bps_score[i]==0)
    {

      allset = all_subset(split_string(Rcpp::as<std::string>(p[i])));

      index = Todecimal(allset);

      for(int j=0;j<index.size();j++)
      {
        if(sub_score<bps_score[index[j]])
        {
          sub_score = bps_score[index[j]];
          sub_parents = bps[index[j]];
        }
      }

      t = split_string(Rcpp::as<std::string>(p[i]));
      int s = 0;
      for(int j = 0;j<t.size();j++)
      {
        if(t[j]=="1")
        {
          s=s+1;
        }
      }

      double sc=0;
      if(s<=num_p)
      {
        for(int j = 0;j<parents.size();j++)
        {
          if(parents[j]==p[i])
          {
            sc = score[j];
            break;
          }
        }

        if(sub_score<sc)
        {

          sub_score = sc;
          sub_parents = p[i];
        }
      }


      bps_score[i] = sub_score;
      bps[i] = sub_parents;
    }
  }
  DataFrame result = DataFrame::create(Named("parents")=p,Named("bps")=bps,Named("score") = bps_score,_["stringsAsFactors"] = false);
  return(result);
}



bool Inset(String x,CharacterVector set)
{
  bool result =FALSE;
  for(int i=0;i<set.size();i++)
  {
    if(x==set[i])
    {
      result= TRUE;
    }
  }
  return(result);
}

//[[Rcpp::export]]
List get_best_parents ( String x, CharacterVector C, CharacterVector all_names)
{
  Environment env = Environment::global_env();
  List bps_all = env["bps_all"];
  DataFrame bps;
  List t;
  List result;
  double score;
  CharacterVector p;
  String s;
  for(int i=0;i<bps_all.size();i++)
  {
    t = bps_all[i];
    s = Rcpp::as<String>(t[2]);
    if(s==x)
    {
      bps = t[0];
      break;
    }
  }

  NumericVector scores = bps["score"];
  CharacterVector best_p = bps["bps"];
  if(C.size()==0)
  {
    return(List::create(Named("bps")=R_NaN,Named("score")=scores[0]));
  }

  CharacterVector variables = all_names;
  for(int i=0;i<variables.size();i++)
  {
    if(variables[i]==x)
      variables.erase(i);
  }

  CharacterVector index(variables.size());
  for(int i=0;i<variables.size();i++)
  {
    index[i] = (Inset(variables[i],C)?"1":"0");
  }

  int k = Tode(convert(index));

  score = scores[k];
  String index_p = best_p[k];

  CharacterVector index_pa = split_string(index_p);

  for(int j =0 ;j<index_pa.size();j++)
  {
    if(Rcpp::as<std::string>(index_pa[j])=="1")
    {
      p.push_back(variables[j]);
    }
  }

  if(p.size()==0)
  {
    p=R_NaN;
  }

  result = List::create(Named("bps")=p,Named("score")=score);

  return(result);
}


//[[Rcpp::export]]
DataFrame sinks(CharacterVector all_names)
{
  Environment env = Environment::global_env();
  List bps_all = env["bps_all"];
  CharacterVector subset(pow(2,all_names.size())-1);
  CharacterVector Sinks(subset.size());
  NumericVector score(subset.size());
  CharacterVector t;
  NumericVector index;
  CharacterVector sub;
  NumericVector sub_index;
  double s,obj;
  int num=0;
  CharacterVector variable2(1);
  for(int i =0;i<subset.size();i++)
  {
    subset[i] = toBinary(i+1,all_names.size());
    t = split_string(Rcpp::as<std::string>(subset[i]));
    for(int j = 0;j<t.size();j++)
    {
      if(t[j]=="1")
      {
        num=num+1;
        index.push_back(j);
      }
    }
    if(num==1)
    {
      Sinks[i] = all_names[index[0]];
      score[i] = get_best_parents(Sinks[i],CharacterVector::create(),all_names)[1];
    }
    else
    {
      sub = all_subset(t);
      sub_index = Todecimal(sub);
      s= R_NegInf;
      CharacterVector variable(sub_index.size());
      String sink;
      for(int j=0;j<sub_index.size();j++)
      {
        variable[j] = all_names[index[j]];
      }

      for(int k =0;k<sub_index.size();k++)
      {
        variable2[0] = variable[k];
        obj= get_best_parents(variable[k],setdiff(variable,variable2),all_names)[1];
        obj= obj+score[sub_index[k]-1];
        if(s<obj)
        {
          s = obj;
          sink = variable[k];
        }
      }

      Sinks[i] = sink;
      score[i] = s;
    }
    index.erase(0,index.size());
  }
  return(DataFrame::create(Named("subset")= subset,Named("sink")=Sinks,Named("score")=score,_["stringsAsFactors"] = false));
}










