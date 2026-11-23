#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
int main(int argc,char**argv){
 const int genes=33694,cells=111668;
 if(argc!=4)return 1;
 const int groups=std::atoi(argv[1]);
 FILE* labels=std::fopen(argv[2],"r");
 std::vector<int> label(cells); for(auto &v:label) if(std::fscanf(labels,"%d",&v)!=1)return 2; std::fclose(labels);
 std::vector<double> sums((size_t) genes*groups);
 std::vector<int> ix(genes),val(genes);
 int g,c,n,current=1,k=0,total=0;
 long long seen=0;
 char buffer[128];
 while(std::fgets(buffer,sizeof buffer,stdin) && buffer[0]=='%'){}
 long long expected=0;int ng,nc;
 if(std::sscanf(buffer,"%d%d%lld",&ng,&nc,&expected)!=3 || ng!=genes || nc!=cells)return 3;
 auto flush=[&](){
   if(label[current-1]>=0){double scale=10000.0/(total?total:1);
     for(int j=0;j<k;j++)sums[(size_t)label[current-1]*genes+ix[j]-1]+=std::log1p(scale*val[j]);}
   total=0;k=0;
 };
 while(std::fgets(buffer,sizeof buffer,stdin)){
   char *p=buffer;g=std::strtol(p,&p,10);c=std::strtol(p,&p,10);n=std::strtol(p,&p,10);
   if(c!=current){if(c<current || c>cells)return 4;flush();current=c;}
   if(g<1||g>genes||n<0||k>=genes)return 5;
   ix[k]=g;val[k++]=n;total+=n;seen++;
 }
 flush();
 if(seen!=expected)return 6;
 FILE* out=std::fopen(argv[3],"wb");if(!out)return 7;
 const auto written=std::fwrite(sums.data(),sizeof(double),sums.size(),out);
 std::fclose(out);return written==sums.size()?0:8;
}
