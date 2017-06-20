

.onLoad = function(libname, pkgname) {

code_swapC =
"srand(time(NULL));
int ok,vc1,vc2,vr1,vr2,a1,a2,npos;
/*int pos [(*nr)-1];*/
int pos [(*nc)-1];
int K=(*Kb);
for(int S=0;S<(*nsim);S++){
             for(int j=0;j<((*nr)*(*nc));j++){
                     mat[j+(S*(*nr)*(*nc))]=mat0[j];}
     int k=0;
     if((S==1)&&(*recursive==1)){
                                K=(*Kt);}
     while(k<(K+1)){
                     ok=0;             
                     vr1=(rand()%(*nr)); 
                     vr2=vr1;                    
                     while((vr2==vr1)||(strata[vr1]!=strata[vr2])){
                                     vr2=(rand()%(*nr));}
                     vc1=(rand()%(*nc));
                     if(mat[vr1+vc1*(*nr)+(S*(*nr)*(*nc))]!=mat[vr2+vc1*(*nr)+(S*(*nr)*(*nc))]){
                                                                npos=0;
                                                                for(int i=0;i<(*nc);i++){
                                                                        if((mat[vr1+i*(*nr)+(S*(*nr)*(*nc))]==mat[vr2+vc1*(*nr)+(S*(*nr)*(*nc))])&&(mat[vr2+i*(*nr)+(S*(*nr)*(*nc))]==mat[vr1+vc1*(*nr)+(S*(*nr)*(*nc))])&&(i!=vc1)){    
                                                                                                                                                                         pos[npos]=i;
                                                                                                                                                                         npos+=1;}}
                                                                if(npos>0){
                                                                           vc2=pos[rand()%npos];
                                                                           ok=1;}
                                                                if(ok==1){
                                                                          a1=mat[vr1+vc1*(*nr)+(S*(*nr)*(*nc))]==mat[vr2+vc2*(*nr)+(S*(*nr)*(*nc))];
                                                                          a2=mat[vr2+vc1*(*nr)+(S*(*nr)*(*nc))]==mat[vr1+vc2*(*nr)+(S*(*nr)*(*nc))];
                                                                          if(a1&&a2){
                                                                                        mat[vr1+vc1*(*nr)+(S*(*nr)*(*nc))]=mat[vr2+vc1*(*nr)+(S*(*nr)*(*nc))];
                                                                                        mat[vr2+vc1*(*nr)+(S*(*nr)*(*nc))]=mat[vr2+vc2*(*nr)+(S*(*nr)*(*nc))];
                                                                                        mat[vr1+vc2*(*nr)+(S*(*nr)*(*nc))]=mat[vr2+vc2*(*nr)+(S*(*nr)*(*nc))];
                                                                                        mat[vr2+vc2*(*nr)+(S*(*nr)*(*nc))]=mat[vr1+vc1*(*nr)+(S*(*nr)*(*nc))];
                                                                                        k+=1;
                                                                                        perms[vr1+vc1*(*nr)+(S*(*nr)*(*nc))]=perms[vr1+vc1*(*nr)+(S*(*nr)*(*nc))]+1;
                                                                                        perms[vr1+vc2*(*nr)+(S*(*nr)*(*nc))]=perms[vr1+vc2*(*nr)+(S*(*nr)*(*nc))]+1;
                                                                                        perms[vr2+vc1*(*nr)+(S*(*nr)*(*nc))]=perms[vr2+vc1*(*nr)+(S*(*nr)*(*nc))]+1;
                                                                                        perms[vr2+vc2*(*nr)+(S*(*nr)*(*nc))]=perms[vr2+vc2*(*nr)+(S*(*nr)*(*nc))]+1;
                                                                                        }
                                                                          }
                     }
     }
if(*recursive==1){
                 for(int j=0;j<((*nr)*(*nc));j++){
                         mat0[j]=mat[j+(S*(*nr)*(*nc))];}}
}"

sign_swapC = signature(
mat0="numeric", 
mat="numeric",
nsim = "integer",
recursive = "integer",
nr = "integer",
nc = "integer",
Kb = "integer",
Kt = "integer",
strata = "integer",
perms = "numeric"
)

code_NBP = 
"int a,b;     
for(int i=0;i<(*nCombSP);i++){
        for(int j=0;j<(*nCombCL);j++){
                a=0;
                b=0;
                for(int k=0;k<(*nli);k++){
                        
                        if((mat[k+dsgSP[2*i]*(*nli)]==dsgCL[2*j])&&(mat[k+dsgSP[2*i+1]*(*nli)]==dsgCL[2*j+1])){
                                                                                               a+=1;}
                        if((mat[k+dsgSP[2*i]*(*nli)]==dsgCL[2*j+1])&&(mat[k+dsgSP[2*i+1]*(*nli)]==dsgCL[2*j])){
                                                                                               b+=1;}
                }
                res[j+(*nCombCL)*i]=min(a,b);
                        
        }        
}"

sign_NBP = signature(
mat="numeric",
dsgSP = "integer",
dsgCL="numeric",
nCombSP="integer",
nCombCL="integer",
nli="integer",
res="numeric")
setCMethod(c("swapClass_swapC","swapClass_NBP"),list(sign_swapC,sign_NBP),body=list(code_swapC,code_NBP),convention=".C",language="C",includes=c("#include <stdlib.h>","#include <time.h>","#define min(a,b) (a<=b?a:b)"),
	where=baseenv())
}


nullModel=function(mat,nsim=100,recursive=TRUE,burnin=NULL,thin=NULL,strata=NULL){
f1=function(x) t(apply(x,1,list))
f2=function(x) x[[1]]
f3=function(x,nc) matrix(x,ncol=nc,byrow=FALSE)
if(!is.matrix(mat) & !is.data.frame(mat)) stop("mat has to be a matrix or a data frame") 
if(!is.matrix(mat)) mat=as.matrix(mat)
if(!is.numeric(mat)) stop("mat has to be numerical")
if(is.null(strata)) strata=rep(1,nrow(mat))
if(!is.vector(strata)) stop("strata has to be a vector")
if(length(strata)!=nrow(mat)) stop("strata has to be a vector of length egual to number of rows in mat")
if(is.null(burnin)) burnin=prod(dim(mat))*10
if(is.null(thin)) thin=max(c(prod(dim(mat)),1000))
if((nrow(mat)<2) | (ncol(mat)<2)) stop("mat has to have more than 1 row and 1 column")
nbp=nbPerm(mat)
if(nbp==0) stop("mat is not swappable")
if(nbp<0.5) warning("be carefull: nbPerm index of mat is lower than 0.5")
strata=as.numeric(strata)
nc=dim(mat)[2];nr=dim(mat)[1]
mat0=as.vector(matrix(mat,1)[1,])
mat=rep(0,length(mat0)*nsim)
perms=rep(0,length(mat))
res=swapClass_swapC(mat0=as.double(mat0),mat=as.double(mat),nsim=as.integer(nsim),recursive=as.integer(recursive*1),nr=as.integer(nr),nc=as.integer(nc),Kb=as.integer(burnin),Kt=as.integer(thin),strata=as.integer(strata),perms=as.double(perms))
mat=lapply(f1(matrix(res$mat,nrow=nsim,byrow=T)),f2)
perms=lapply(f1(matrix(res$perms,nrow=nsim,byrow=T)),f2)
return(list(sim=lapply(mat,f3,nc=nc),perms=lapply(perms,f3,nc=nc)))}

nbPerm=function(mat){
n=prod(dim(mat))
dsgSP=combn(dim(mat)[2],2)-1
nCombSP=dim(dsgSP)[2]
dsgSP=matrix(dsgSP,ncol=1)[,1]
nli=dim(mat)[1]
mat=matrix(mat,ncol=1)[,1]
abm=as.numeric(levels(as.factor(mat)))
if(length(abm)>1){
dsgCL0=combn(abm,2)
nCombCL=dim(dsgCL0)[2]
dsgCL=matrix(dsgCL0,ncol=1)[,1]
res=rep(0,nCombSP*nCombCL)
res=swapClass_NBP(mat=as.double(mat),dsgSP=as.integer(dsgSP),dsgCL=as.double(dsgCL),nCombSP=as.integer(nCombSP),nCombCL=as.integer(nCombCL),nli=as.integer(nli),res=as.double(res))$res
res=colSums(matrix(res,ncol=nCombCL,byrow=T))/2
res=sum(res/n)}
else{
res=0}
return(res)}



