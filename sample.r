s<-read.delim('./sample.input')
ptm <- proc.time()
s<-subset(s,s$cell==23)
Unique<-unique(s$sp)
n.obs<-length(Unique)
iter<-500
dbks.out<-rep(NA,n.obs)
pvalue.out<-rep(NA,n.obs)
#The Loop
for (k in 1:n.obs){
  s.k <- s[s$sp == Unique[k], ]
  if (nrow(s.k)>30 & sum(s.k$seed)!=0 & sum(s.k$tree != 0)){
    x<-s.k$tree
    y<-s.k$seed
    n<-length(x)
    obs<-matrix(NA,n,4)
    exp<-matrix(NA,n,4)
#True Dbks
    for(i in 1:n){
      x0<-x[i]
      y0<-y[i]
      q1<-table(factor(x>x0,levels=c(F,T)),factor(y>=y0,levels=c(F,T)))["TRUE","TRUE"]/n
      q2<-table(factor(x<=x0,levels=c(F,T)),factor(y>=y0,levels=c(F,T)))["TRUE","TRUE"]/n
      q3<-table(factor(x<=x0,levels=c(F,T)),factor(y<y0,levels=c(F,T)))["TRUE","TRUE"]/n
      q4<-table(factor(x>x0,levels=c(F,T)),factor(y<y0,levels=c(F,T)))["TRUE","TRUE"]/n
      obs[i,]<-cbind(q1,q2,q3,q4)
      eq1<-table(factor(x>x0,levels=c(F,T)))["TRUE"]/n*table(factor(y>=y0,levels=c(F,T)))["TRUE"]/n
      eq2<-table(factor(x<=x0,levels=c(F,T)))["TRUE"]/n*table(factor(y>=y0,levels=c(F,T)))["TRUE"]/n
      eq3<-table(factor(x<=x0,levels=c(F,T)))["TRUE"]/n*table(factor(y<y0,levels=c(F,T)))["TRUE"]/n
      eq4<-table(factor(x>x0,levels=c(F,T)))["TRUE"]/n*table(factor(y<y0,levels=c(F,T)))["TRUE"]/n
      exp[i,]<-cbind(eq1,eq2,eq3,eq4)
    }
 
    dbks.out[k]<-max(abs(obs-exp))

#Data Randomization (Distribution of Dbks)
    obsr<-matrix(NA,n,4)
    expr<-matrix(NA,n,4)
    dbksr<-rep(NA,iter)
    for(j in 1:iter){
      yr<-sample(y,n,replace=F)
      for(l in 1:n){
        x0<-x[l]
        y0<-yr[l]
        q1<-table(factor(x>x0,levels=c(F,T)),factor(yr>=y0,levels=c(F,T)))["TRUE","TRUE"]/n
        q2<-table(factor(x<=x0,levels=c(F,T)),factor(yr>=y0,levels=c(F,T)))["TRUE","TRUE"]/n
        q3<-table(factor(x<=x0,levels=c(F,T)),factor(yr<y0,levels=c(F,T)))["TRUE","TRUE"]/n
        q4<-table(factor(x>x0,levels=c(F,T)),factor(yr<y0,levels=c(F,T)))["TRUE","TRUE"]/n
        obsr[l,]<-cbind(q1,q2,q3,q4)
        eq1<-table(factor(x>x0,levels=c(F,T)))["TRUE"]/n*table(factor(yr>=y0,levels=c(F,T)))["TRUE"]/n
        eq2<-table(factor(x<=x0,levels=c(F,T)))["TRUE"]/n*table(factor(yr>=y0,levels=c(F,T)))["TRUE"]/n
        eq3<-table(factor(x<=x0,levels=c(F,T)))["TRUE"]/n*table(factor(yr<y0,levels=c(F,T)))["TRUE"]/n
        eq4<-table(factor(x>x0,levels=c(F,T)))["TRUE"]/n*table(factor(yr<y0,levels=c(F,T)))["TRUE"]/n
        expr[l,]<-cbind(eq1,eq2,eq3,eq4)
      }
      dbksr[j]<-max(abs(obsr-expr))
    }
    pvalue.out[k]<-table(factor(dbksr>dbks.out[k],levels=c(T,F)))["TRUE"]/iter
  }
}
label<-as.character(Unique)
twodksout<-cbind(label,dbks.out,pvalue.out)
proc.time() - ptm
write.table(twodksout,"2DKSCell23.txt",sep="\t",row.names=F)
