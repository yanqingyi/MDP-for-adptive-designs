##### https://github.com/yanqingyi/MDP-for-adptive-designs.git

##### The R-code for binary responses: notation and inputs #################
r # runs of simulation.

d # the number of subjects alloation by adaptive randomization; d=n-2, where n the totoal number of subjects in the trial.

pa # success probability on A.
pb # success probability on B.

h # gamma: the minimum allocatiom proportion.

eta<-0.1 # tuning parameter for MDP.

zeta<-0.05 # tuning parameter for MDP.

cutZ<-1.645 # critical value for type I error rate/statistical power.

##### output from the mdpb function: Statistical power can be read from 'powerZ'; the distribution of N_A can be obtained from 'probNA' thus the distribution of N_A/n.

######################## the mdpb function #################################
mdpb<-function(pa,pb,eta,zeta,d,r,cutZ,h){
  naCount<-rep(0,d+1)
  xNA<-1:(d+1)
  sCount<-rep(0,d+3)
  xS<-0:(d+2)
  
  zCount<-0
  
  zDistr<-rep(0,r)
  
  for (n in 1:r){
    sa<-0
    sb<-0
    na<-1 
    nb<-1
    
    u1<-runif(2,0,1)
    if (u1[1]<pa) {sa<-sa+1}
    if (u1[2]<pb) {sb<-sb+1}
    
    paHat<-(sa+1)/(na+2)
    pbHat<-(sb+1)/(nb+2)
    
    v<-rep(0,3)
    vm<-rep(0,3)
    ub<-pbHat+eta
    lb<-pbHat-eta
    p<-1/2
    
    ee<-0
    
    if (paHat>ub){p<-1-h
    v[1]<-v[1]+1
    vm[1]<- 2*(p*paHat+(1-p)*pbHat)
    ee<-1}
    if (paHat<lb){p<-h
    v[2]<-v[2]+1
    vm[2]<-2*(p*paHat+(1-p)*pbHat)
    ee<-1}
    if (ee==0){p<-1/2
    v[3]<-v[2]+1
    vm[3]<-2*(p*paHat+(1-p)*pbHat)}
    
    
    pas<-(sa+1+1)/(na+1+2) 
    ae<-0
    
    if (pas>ub) {as<- vm[1]
    ae<-1}
    if (pas<lb) {as<- vm[2]
    ae<-1}
    if (ae==0) {as<- vm[3]}
    
    paf<-(sa+1)/(na+1+2) 
    ae<-0
    if (paf>ub) {af<-vm[1]
    ae<-1}
    if (paf<lb) {af<-vm[2]
    ae<-1}
    if (ae==0) {af<-vm[3]}
    
    av<-paHat+paHat*as+(1-paHat)*af
    
    pbs<-(sb+1+1)/(nb+1+2) 
    ubt<-pbs+eta
    lbt<-pbs-eta
    be<-0
    
    if (paHat>ubt) {bs<- vm[1]
    be<-1}
    if (paHat<lbt) {bs<- vm[2]
    be<-1}
    if (be==0) {bs<- vm[3]}
    
    pbf<-(sb+1)/(nb+1+2) 
    ubt<-pbf+eta
    lbt<-pbf-eta
    be<-0
    if (paHat>ubt) {bf<- vm[1]
    be<-1}
    if (paHat<lbt) {bf<- vm[2]
    be<-1}
    if (be==0) {bf<- vm[3]}
    
    bv<-pbHat+pbHat*bs+(1-pbHat)*bf
    
    ev<-0
    if(av>bv+zeta){p<-1-h
    ev<-1}
    if (av<bv-zeta){p<-h
    ev<-1}
    if (ev==0){p<-1/2}
    
    for (i in 1:d){
      u<-runif(3,0,1)
      e<-0
      if (u[1]<p){na<-na+1
      if (u[2]<pa){sa<-sa+1}
      e<-1}
      if (e==0){nb<-nb+1
      if (u[3]<pb){sb<-sb+1}}
      
      paHat<-(sa+1)/(na+2)
      pbHat<-(sb+1)/(nb+2)
      
      ub<-pbHat+eta # for optimal value vm
      lb<-pbHat-eta
      ee<-0
      
      if (paHat>ub){p<-1-h
      v[1]<-v[1]+1
      vm[1]<-(i+2)*(p*paHat+(1-p)*pbHat)
      ee<-1}
      if (paHat<lb){p<-h
      v[2]<-v[2]+1
      vm[2]<-(i+2)*(p*paHat+(1-p)*pbHat)
      ee<-1}
      if (ee==0){p<-1/2
      v[3]<-v[2]+1
      vm[3]<-(i+2)*(p*paHat+(1-p)*pbHat)}
      
      
      pas<-(sa+1+1)/(na+1+2) 
      ae<-0
      
      if (pas>ub) {as<- vm[1]
      ae<-1}
      if (pas<lb) {as<- vm[2]
      ae<-1}
      if (ae==0) {as<- vm[3]}
      
      paf<-(sa+1)/(na+1+2) 
      ae<-0
      if (paf>ub) {af<-vm[1]
      ae<-1}
      if (paf<lb) {af<-vm[2]
      ae<-1}
      if (ae==0) {af<-vm[3]}
      
      av<-paHat+paHat*as+(1-paHat)*af
      
      
      pbs<-(sb+1+1)/(nb+1+2) 
      ubt<-pbs+eta
      lbt<-pbs-eta
      be<-0
      
      if (paHat>ubt) {bs<- vm[1]
      be<-1}
      if (paHat<lbt) {bs<- vm[2]
      be<-1}
      if (be==0) {bs<- vm[3]}
      
      pbf<-(sb+1)/(nb+1+2) 
      ubt<-pbf+eta
      lbt<-pbf-eta
      be<-0
      if (paHat>ubt) {bf<- vm[1]
      be<-1}
      if (paHat<lbt) {bf<- vm[2]
      be<-1}
      if (be==0) {bf<- vm[3]}
      
      bv<-pbHat+pbHat*bs+(1-pbHat)*bf
      
      ev<-0
      if(av>bv+zeta){p<-1-h
      ev<-1}
      if (av<bv-zeta){p<-h
      ev<-1}
      if (ev==0){p<-1/2}
      
    }
    s<-sa+sb
    naCount[na]<-naCount[na]+1
    sCount[s+1]<-sCount[s+1]+1
    
    paHat<-(sa+1)/(na+2)
    pbHat<-(sb+1)/(nb+2)
    
    
    Def<-paHat-pbHat
    varDef<-sqrt(paHat*(1-paHat)/na+pbHat*(1-pbHat)/nb)
    
    z<-Def/varDef
    zDistr[n]<-z
    
    c2<-0
    if (z<cutZ){c2<-1} 
    if (c2==0){zCount<-zCount+1}
    
  }
  probNA<-naCount/r
  probS<-sCount/r
  
  powerZ<-zCount/r
  
  return(list(zDistr=zDistr,xNA=xNA,probNA=probNA,probS=probS,
              powerZ=powerZ,pa=pa,pb=pb,cutZ=cutZ, h=h, d=d, rep=r, eta=eta, zeta=zeta))}


##### example: the codes for Figure 3
r <- 1000000
eta <- 0.1 # tuning parameter of MDP.
zeta <- 0.05 # tuning parameter of MDP.
n <- 200
d <- 200-2

h <- 0.25

pa <- 0.9
pb <- 0.8
cutZ<-1.645

aa <- mdpb(pa,pb,eta,zeta,d,r,cutZ,h)

pNA <- aa$xNA/n
plot(pNA, aa$probNA)
matplot(pNA, aa$probNA, xlab=expression(N[A]/n),ylab=" ", xlim=c(0,1),type="l", lty=1)

##################################################################################
########## R-code for exponential distributions: notation and input ##############
r # runs of simulation

eta<-0.1 # tuning parameter of MDP.
zeta<-0.05 # tuning parameter of MDP.

h # gamma: the minimum allocation proportion. 

r1<-1000 # the number of repeatition used to calculate the intergral for general distributions

d1 # number of subjects assigned to each treatment by complete randomization
d # number of subjectes assigned using adaptive randomization; d=n-2*d1, where n is the total number of subjects in the trial.

ma # mean of the distribution for treatment A
mb # mean of the distribution for treatment B

cutZ<-1.82 # critical value for type I error rate/statisitcal power

##### output from the mdpc2 function: Statistical power can be read from 'powerZ'; the distribution of N_A can be obtained from 'probNA' thus the distribution of N_A/n.

############## the mdpc2 fucntion #######################################

mdpc2<-function(ma,mb,eta,zeta,d,d1,r,r1,cutZ,h){
  naCount<-rep(0,d+d1)
  xNA<-1:(d+d1)
  zCount<-0
  
  zDistr<-rep(0,r)
  
  for (n in 1:r){
    na<-d1 
    nb<-d1
    
    u1<-rexp(d1,1/ma)
    u2<-rexp(d1,1/mb)
    
    maHat<-mean(u1)
    mbHat<-mean(u2)
    
    v<-rep(0,3)
    vm<-rep(0,3)
    ub<-mbHat+eta
    lb<-mbHat-eta
    p<-1/2
    
    ee<-0
    d2<-2*d1
    if (maHat>ub){p<-1-h
    v[1]<-v[1]+1
    vm[1]<-d2*(p*maHat+(1-p)*mbHat)
    ee<-1}
    if (maHat<lb){p<-h
    v[2]<-v[2]+1
    vm[2]<-d2*(p*maHat+(1-p)*mbHat)
    ee<-1}
    if (ee==0){p<-1/2
    v[3]<-v[2]+1
    vm[3]<-d2*(p*maHat+(1-p)*mbHat)}
    
    at<-rep(0,r1)
    ra<-rexp(r1,1/maHat)
    mat<-(na*maHat+ra)/(na+1)
    ind1<-which(mat>ub)
    if(length(ind1)>0){at[ind1]<-vm[1]}
    
    ind2<-which(mat<lb)
    if(length(ind2)>0){at[ind2]<-vm[2]}
    
    att1<-which(mat<=ub)
    att2<-which(mat>=lb)
    ind3<-intersect(att1,att2)
    if(length(ind3)>0){at[ind3]<-vm[3]}
    
    av<-maHat+sum(at)
    
    bt<-rep(0,r1)
    rb<-rexp(r1,1/mbHat)
    mbt<-(nb*mbHat+rb)/(nb+1)
    ubt<-mbt+eta
    lbt<-mbt-eta
    indb1<-which(maHat>ubt)
    if(length(indb1)>0){bt[indb1]<-vm[1]}
    
    indb2<-which(maHat<lbt)
    if(length(indb2)>0){bt[indb2]<-vm[2]}
    
    btt1<-which(maHat<=ubt)
    btt2<-which(maHat>=lbt)
    indb3<-intersect(btt1,btt2)
    if(length(indb3)>0){bt[indb3]<-vm[3]}
    
    bv<-mbHat++sum(bt)
    
    ev<-0
    if(av>bv+zeta){p<-1-h
    ev<-1}
    if (av<bv-zeta){p<-h
    ev<-1}
    if (ev==0){p<-1/2}
    
    for (i in 1:d){	u<-runif(1,0,1)
    e<-0
    if (u<p){na<-na+1
    raa<-rexp(1,1/ma)
    maHat<-((na-1)*maHat+raa)/na
    e<-1}
    if (e==0){nb<-nb+1
    rbb<-rexp(1,1/mb)
    mbHat<-((nb-1)*mbHat+rbb)/nb}
    
    ub<-mbHat+eta # for optimal value vm
    lb<-mbHat-eta
    ee<-0
    
    vm<-rep(0,3)
    
    if (maHat>ub){p<-1-h
    v[1]<-v[1]+1
    vm[1]<-(i+d2)*(p*maHat+(1-p)*mbHat)
    ee<-1}
    if (maHat<lb){p<-h
    v[2]<-v[2]+1
    vm[2]<-(i+d2)*(p*maHat+(1-p)*mbHat)
    ee<-1}
    if (ee==0){p<-1/2
    v[3]<-v[2]+1
    vm[3]<-(i+d2)*(p*maHat+(1-p)*mbHat)}
    
    at<-rep(0,r1)
    ra<-rexp(r1,1/maHat)
    mat<-(na*maHat+ra)/(na+1)
    ind1<-which(mat>ub)
    if(length(ind1)>0){at[ind1]<-vm[1]}
    
    ind2<-which(mat<lb)
    if(length(ind2)>0){at[ind2]<-vm[2]}
    
    att1<-which(mat<=ub)
    att2<-which(mat>=lb)
    ind3<-intersect(att1,att2)
    if(length(ind3)>0){at[ind3]<-vm[3]}
    
    av<-maHat+sum(at)
    
    bt<-rep(0,r1)
    rb<-rexp(r1,1/mbHat)
    mbt<-(nb*mbHat+rb)/(nb+1)
    ubt<-mbt+eta
    lbt<-mbt-eta
    indb1<-which(maHat>ubt)
    if(length(indb1)>0){bt[indb1]<-vm[1]}
    
    indb2<-which(maHat<lbt)
    if(length(indb2)>0){bt[indb2]<-vm[2]}
    
    btt1<-which(maHat<=ubt)
    btt2<-which(maHat>=lbt)
    indb3<-intersect(btt1,btt2)
    if(length(indb3)>0){bt[indb3]<-vm[3]}
    
    bv<-mbHat+sum(bt)
    
    ev<-0
    if(av>bv+zeta){p<-1-h
    ev<-1}
    if (av<bv-zeta){p<-h
    ev<-1}
    if (ev==0){p<-1/2}
    
    }
    naCount[na]<-naCount[na]+1
    
    Def<-maHat-mbHat
    sd<-sqrt(maHat^2/na+mbHat^2/nb)
    
    z<-Def/sd
    zDistr[n]<-z
    
    c2<-0
    if (z<cutZ){c2<-1} 
    if (c2==0){zCount<-zCount+1}
    
  }
  probNA<-naCount/r
  
  powerZ<-zCount/r
  
  return(list(zDistr=zDistr,xNA=xNA,probNA=probNA,powerZ=powerZ,ma=ma,mb=mb,cutZ=cutZ, h=h, d=d, rep=r, eta=eta, zeta=zeta))}

####### Example: the codes for Figure 2
r<-1000000

cutZ<-1.82 

eta<-0.1
zeta<-0.05
h<-0.25

r1<-1000 # to calculate the intergral for general distributions

d1<-10
d<-80

ma<-7
mb<-5
c75d10<-mdpc2(ma,mb,eta,zeta,d,d1,r,r1,cutZ,h)

pNA<-c75d10$xNA/100
matplot(pNA,c75d10$probNA,  xlab=expression(N[A]/n), ylab=' ', xlim=c(0,1), ylim=c(0,0.08), type='l',lty=1)

