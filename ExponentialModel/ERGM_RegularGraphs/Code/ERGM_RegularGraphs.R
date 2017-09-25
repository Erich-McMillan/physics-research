library(Rcpp)
library(igraph)

cppFunction('
List genERGM(NumericVector ds)
        {
            int N=ds.size();
            int D = ds[1];
            double lg = log((N-1)/(double)D - 1.0);
            double pr = 1.0/(1.0+exp(lg));

            std::vector<int> startlist;
            std::vector<int> endlist;
            
            int i,j;
            
            for(i=0;i<N;i++)
            {
              for(j=i+1;j<N;j++)
              {
                if(R::runif(0,1)<pr)
                {
                  startlist.push_back(i);
                  endlist.push_back(j);
                }
              }
            }
              return List::create(startlist,endlist);
        }
')

zERGM=function(ds)
{
  el_list=genERGM(ds)
  el=cbind(el_list[[1]],el_list[[2]])+1
  #print(el)
  g=graph_from_edgelist(el,directed=F)
  return(g)
}

N = 1000
D = 20
Sq=0
SqTotal=10000
count = 0

deg = 0
ds=rep(D,N)

for(Sq in seq(SqTotal)) 
{
  
  g=zERGM(ds)
  degr = degree(g)
  #print(mean(degr))
  if(is.nan(mean(degr)))
  {
    print(degree(g))
    print()
  }
  deg = deg + mean(degree(g))
  
}

print("mean degree:")
print(deg/SqTotal)
