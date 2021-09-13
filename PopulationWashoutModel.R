
##General functions
#Base set of parameters; note that per capita reproduction R corresponds to r in the paper; gammma is emergence rate in insect models
#Sigma denotes the proportion of individuals invulnerable to advective dispersal; fixed to 0 throughout.
parms=c(d=0.08,mu=0.3,rho=0.25,sigma=0,R=0.25,L=100,m=1,n=2,gamma=0)

#Simple function to update parameter vector. Note that parmsNew must be a vector with named entries
setParms=function(parmsOld,parmsNew){ parmsOld[names(parmsNew)]=parmsNew; return(parmsOld); }

#Function to manipulate the off-diagonal entries of (square) projection matrices x.
#k is either (a) a single index of diagonal to extract, or 
#(b) indeces of diagonals to replace with vec. k<0 refers to upper diagonals.
#vec is the replacement for when modifying off-diagonal elements
offDiag=function(x,k,vec=NULL){
	if(nrow(x)!=ncol(x)){ print("offDiag: input matrix not square!!"); return(NA); }
	#if vec not given, return specified off-diagonal:
	if(is.null(vec)){
		y=x
		x[row(x)==col(x)-k]=NA
		return(y[is.na(x)])  
	}
	#if no off-diagonals given (k empty), return original matrix:
	if(length(k)==0) return(x)
	#if vec AND multiple off-diagonals given, replace all:
	if(length(k)>1){ for(i in k) x[row(x)==col(x)-i]=vec; return(x); }
	#if vec and single k given, return replacement:
	x[row(x)==col(x)-k]=vec; return(x);
}

UpstreamDispersBias=DispersBias=0 #DispersalBias is number of patches by which all dispersers move upstream (values <0) or downstream (values >0) of home patch
SpatialMatr=function(psi,npatch){
	D=diag(npatch); dists=-npatch:npatch; if(psi==0) return(D);
	distsPrs=pnorm(dists+0.5, mean=DispersBias, sd=psi*(pi/2)^0.5) - pnorm(dists-0.5, mean=DispersBias, sd=psi*(pi/2)^0.5)
	for(i in 1:npatch) D[,i]=distsPrs[(npatch+2-i):(2*npatch+1-i)]
	D=D*(!lower.tri(D)) + UpstreamDispersBias*t(D*lower.tri(D)) + (1-UpstreamDispersBias)*D*lower.tri(D)
	return(D)
}



#Function to plot heatmaps
require(fields)
HeatmapPlot=function(spacetime,Xax=1:dim(spacetime)[2],Yax=1:dim(spacetime)[1],XaxN="X",YaxN="Y",figtitle="Title",Zlim=c(0,max(spacetime)),cont=0,cexAx=1,SciPen=-10,ContCex=1.5,ContPlot=t(spacetime),Nlevels=5){  
	image.plot(x=Xax,y=Yax,z=t(spacetime), zlim=Zlim,xlab=XaxN,ylab=YaxN,cex.axis=cexAx,cex.lab=cexAx,main=figtitle, 
	col=ColPallette);  box();  
	options(scipen=SciPen)
	contour(x=Xax,y=Yax,z=ContPlot,add=T,col=cont,labcex=ContCex,nlevels=Nlevels,lwd=2.5)
	options(scipen=0)
}
ColPallette=rev(rainbow(1e3,start=0,end=0.7))

##Function to build life history matrices
#Throughout, propInvalid denotes the proportion of patches affected by permanent habitat loss (parameter E in paper),
#and HabSeed is used to run multiple treatments on the same set of habitat loss configurations
Lmatr=function(parms,propInvalid=0,HabSeed=NULL,LHfun="Mussel"){
	Npatch=as.numeric(parms["L"])
	#Number of patches affected by permanent habitat loss:
	NpatchInvalid=round(Npatch*propInvalid,0); set.seed(HabSeed);
	#Assign patches with habitat loss
	E=sample(rep(c(0,1), c(NpatchInvalid,Npatch-NpatchInvalid)))
	#Matrix describing advective dispersal:
	D=offDiag((1-parms["rho"]+parms["rho"]*parms["sigma"])*diag(parms["L"]), -(0:parms["m"])[-1], (1/parms["m"])*(1-parms["mu"])*parms["rho"]*(1-parms["sigma"]))
	#Matrix describing active (recruitment) dispersal:
	F=SpatialMatr(parms["n"],parms["L"])

	#Return projection matrix:
	if(LHfun=="Mussel") return(D%*%((1-parms["d"])*E*diag(Npatch) + parms["R"]*F))
	if(LHfun=="Insect") return(D%*%((1-parms["d"])*E*(1-parms["gamma"])*diag(Npatch) + parms["R"]*parms["gamma"]*F))
	if(LHfun=="Fish"){
		#This matrix ensures adults reproduce only in suitable habitats:
		M=diag((rep(1,Npatch)%*%(F*E%*%t(E)+1e-6))[1,]^-1)
		return(((1-parms["d"])*F*E%*%t(E)%*%M + parms["R"]*D))
	}
}






##Functions to simulate different scenarios of habitat loss:
#Function to calculate stochastic growth rate under temporary patch extinctions:
StochLambda=function(x,trim=2){ Ns=rowMeans(x[,-(1:trim)]); return(mean(Ns[-1]/Ns[-length(Ns)], na.rm=T)); }

#Function to run stochastic simulations (if applicable) and calculate population growth rate:
#propInvalid is proportion of patches affected by permanent habitat loss (parameter E in text)
#ExtProp is proportion of vulnerable patches affected by temporary habitat loss (exinction) each year
#CVext is variation in the long-term disturbance frequency of patches;
#we only explore extreme case with CVext=0 (patches either vulnerable to disturbance or not)
tries=function(parms,out=1,runl=1e2, propInvalid=0, ExtProp=0, CVext=0, HabSeed=NULL, LHfun="Mussel"){
	#Build projection matrix for given life history
	L=Lmatr(parms,propInvalid,HabSeed,LHfun)

	#For testing effects of temporary habitat loss in the form of patch extinctions:
	if(ExtProp>0){
		N=matrix(0,runl,parms["L"]) #Matrix to store simulation results
		N[1,]=rep(1,parms["L"]) #Initially 1 individual in every patch (this is arbitrary)
		#Mean disturbance frequencies for every patch in system:
		Es=matrix(rnorm(parms["L"],(1-ExtProp),CVext*ExtProp), 1,parms["L"]); Es[Es<0 | Es>1]=ExtProp;
		#Generate disturbace history throughout landscape for entire run (1=patch survives year)
		E=apply(Es, 2, function(x) rbinom(runl,1,prob=x))
		#Run stochastic simulation and return long-term stochastic growth rate:
		for(t in 2:runl) N[t,]=L%*%N[t-1,]*E[t,]
		return(StochLambda(N))
	}

	#If instead simulating permanent habitat loss, return dominant eigenvalue of projection matrix:
	return(max(Re(eigen(L)$values)))
}



###Figure 3: effects of permanent habitat loss on population growth:
parms=c(d=0.1,mu=0.3,rho=0.37,sigma=0,R=0.35,L=1e2,m=1,n=2,gamma=0) #Base mussel parameters from Table A1

#propInvalid is proportion of patches affected by permanent habitat loss (parameter E in text)
tries=function(parms,out=1,propInvalid=0, HabSeed=NULL, LHfun="Mussel"){
	#Build projection matrix for given life history
	L=Lmatr(parms,propInvalid,HabSeed,LHfun)
	#Return dominant eigenvalue of projection matrix:
	return(max(Re(eigen(L)$values)))
}

dev.new(); par(mfrow=c(1,2));

##Small-scale habitat loss:
ntests=25 #Number of parameter sets to test
Rhos=seq(0,0.5,length=ntests) #set of advective dispersal levels
propsInvalids=seq(0,0.99,length=ntests) #set of habitat loss levels
nreps=50; HabSeeds=rnorm(nreps); #Number of times to replicate each analysis (to aggregate across different spatial configurations of habitat loss)
Toplot3=array(0,dim=c(ntests,length(Rhos),nreps)) #Object to to hold results
for(i in 1:ntests) for(j in 1:length(Rhos)) for(k in 1:nreps) Toplot3[i,j,k]=tries(setParms(parms,c(rho=Rhos[j])),propInvalid=propsInvalids[i],HabSeed=HabSeeds[k])
#Plot mean poulation growth rates:
LocalGrowthRelatives=t(apply(apply(Toplot3,c(1,2),mean),1,c))/max(Toplot3)
RunAvg=function(x,n=168) filter(x,rep(1/n,n),sides=1) #For clarity, smooth out contours using a running average
LocalGrowthRelativesSmoo=apply(LocalGrowthRelatives,2,RunAvg,3); LocalGrowthRelativesSmoo[1:2,]=LocalGrowthRelatives[1:2,];
HeatmapPlot(LocalGrowthRelatives,Zlim=c(0.39,1),Yax=Rhos, cont=1,SciPen=-1,ContCex=1.5,ContPlot=LocalGrowthRelativesSmoo,Nlevels=6)



##Large-scale habitat loss:
ntests=25 #Number of parameter sets to test
Rhos=seq(0,0.5,length=ntests) #set of advective dispersal levels
Ls=1:15 #Sizes of remnant habitat domains
#nreps=50 #Number of times to replicate each analysis (to aggregate across different spatial configurations of habitat loss)
Toplot=matrix(0,length(Ls),length(Rhos)) #Object to to hold results
for(i in 1:length(Ls)) for(j in 1:length(Rhos)) Toplot[i,j]=tries(setParms(parms,c(rho=Rhos[j],L=Ls[i])))
LargeGrowthRelatives=t(Toplot[15:1,])/max(Toplot)
HeatmapPlot(LargeGrowthRelatives,Zlim=c(0.375,1),Yax=Rhos, cont=1,SciPen=-1,ContCex=1.5,Nlevels=6)



###Figure 3.1: effects of temporary habitat loss (i.e., local extinctions from disturbances) on population growth:
parms=c(d=0.1,mu=0.3,rho=0.37,sigma=0,R=0.4,L=1e2,m=1,n=2,gamma=0) #Base mussel parameters from Table A1

#Function to run stochastic simulations and calculate population growth rate
#propInvalid is proportion of patches affected by permanent habitat loss (parameter E in text)
#runl is length of stochastic simulations; much longer values can lead to rounding errors at high disturbance levels
#ExtProp is proportion of vulnerable patches affected by temporary habitat loss (exinction) each year
#permfreq is frequency of patches invulnerable to disturbance
#CVext is variation in the long-term disturbance frequency of vulnerable patches;
#we only explore extreme case with CVext=0 (patches either vulnerable to disturbance or not)
triesStoch=function(parms,out=1,runl=1e2, propInvalid=0, ExtProp=0, CVext=0, permfreq=0, RhoLevs){
	#Mean disturbance frequencies for every patch in system:
	#When some patches invlunerable, adjust ExtProp to keep system-wide disturbance frequency identical to permfreq=0
	Es=matrix(rnorm(parms["L"],(1-ExtProp/(1-permfreq)),CVext*ExtProp), 1,parms["L"]); Es[Es<0 | Es>1]=ExtProp;
	#Generate disturbace history throughout landscape for entire run (1=patch survives year)
	E=apply(Es, 2, function(x) rbinom(runl,1,prob=x))
	#Invunerable patches never disturbed:
	E[,rbinom(parms["L"],1,permfreq)*(1:parms["L"])]=1
	
	#Matrix to store simulation results
	N=matrix(0,runl,parms["L"])
	#Initially 1 individual in every patch (this is arbitrary)	
	N[1,]=rep(1,parms["L"])
	#Vector to hold stochastic growth rates across all dispersal levels
	StochLambdas=double(length(RhoLevs))
	for(i in 1:length(RhoLevs)){
		L=Lmatr(setParms(parms,c(rho=RhoLevs[i])),propInvalid);
		for(t in 2:runl) N[t,]=L%*%N[t-1,]*E[t,]
		StochLambdas[i]=StochLambda(N)
	}
	return(StochLambdas)
}

dev.new(); par(mfrow=c(1,2));

ntests=25 #Number of parameter sets to test
nreps=50 #Number of times to replicate each analysis (to aggregate across different histories of disturbance)
EPs=seq(0.01,0.4,length=ntests) #set of disturbance frequencies loss levels
Rhos=seq(0,.5,length=ntests) #set of advective dispersal levels
RndmStoch3=array(0,dim=c(ntests,length(Rhos),nreps)) #Object to to hold results
for(i in 1:ntests) for(k in 1:nreps) RndmStoch3[i,,k]=triesStoch(parms,ExtProp=EPs[i],RhoLevs=Rhos)
nreps=4*nreps; StochWithPerms3=array(0,dim=c(ntests,length(Rhos),nreps));
for(i in 1:ntests) for(k in 1:nreps) StochWithPerms3[i,,k]=triesStoch(parms,ExtProp=EPs[i],permfreq=0.3,RhoLevs=Rhos)

StochWithPerms3[StochWithPerms3>1.25]=1.256 #rounding errors in a few cases yeild inaccurate estimates of population growth
AllStoch=apply(apply(RndmStoch3,c(1,2),mean),1,c)
StochWithPerms=apply(apply(StochWithPerms3,c(1,2),mean),1,c)
HeatmapPlot(AllStoch/max(AllStoch),Zlim=c(0.519,1), Yax=Rhos,Xax=EPs,cont=1,SciPen=-1,ContCex=1.5,Nlevels=5)
HeatmapPlot(StochWithPerms/StochWithPerms[1,1],Zlim=c(0.67,1), Yax=Rhos,Xax=EPs,cont=1,SciPen=-1,ContCex=1.5,Nlevels=4)



###Figure 1: washout effects under different life histories
#Function to compare effects of each stressor on population growth for a given life history:
relatives=function(parms,LHfun,HablossLevel=0.85,nreps=20){
	#Growth rate under ideal conditions
	lamBase=tries(setParms(parms,c(rho=0,L=1e2)),propInvalid=0,LHfun=LHfun)
	#Growth rate under advective dispersal only
	lamAD=tries(setParms(parms,c(L=1e2)),propInvalid=0,LHfun=LHfun)

	lamHL=lamADHL=double(nreps); for(i in 1:nreps){
		#Growth rate under habitat loss only
		lamHL[i]=tries(setParms(parms,c(rho=0,L=1e2)),propInvalid=HablossLevel,HabSeed=i,LHfun=LHfun)
		#Growth rate under habitat loss and advective dispersal
		lamADHL[i]=tries(setParms(parms,c(L=1e2)),propInvalid=HablossLevel,HabSeed=i,LHfun=LHfun)
	}
	#return percent declines in lambda in each case, compared to ideal conditions (lamBase)
	out=1-c(AD=lamAD,HL=mean(lamHL),ADHL=mean(lamADHL))/lamBase 
	return(c(BaseLambda=lamBase,out))
}

dev.new(); HL=0.85;
Insect=relatives(parms=c(d=0.5,mu=0.3,rho=0.5,sigma=0,R=1.43,L=100,m=2,n=0.35,gamma=0.8),LHfun="Insect",HablossLevel=HL)
Mussel=relatives(parms=c(d=0.10,mu=0.3,rho=0.37,sigma=0,R=0.35,L=100,m=1,n=2,gamma=0),LHfun="Mussel",HablossLevel=HL)
MusselLong=relatives(parms=c(d=0.025,mu=0.3,rho=0.37,sigma=0,R=0.125,L=100,m=1,n=2,gamma=0),LHfun="Mussel",HablossLevel=HL)
Fish=relatives(parms=c(d=0.25,mu=0.3,rho=0.5,sigma=0,R=0.5,L=100,m=4,n=5,gamma=0),LHfun="Fish",HablossLevel=HL)
FishLong=relatives(parms=c(d=0.15,mu=0.3,rho=0.5,sigma=0,R=0.26,L=100,m=4,n=5,gamma=0),LHfun="Fish",HablossLevel=HL)
Bplot2=rbind(Insect,Mussel,MusselLong,Fish,FishLong); colnames(Bplot)=c("BaseLam","Advective Dispersal","Habitat Loss","Advective Dispersal and Habitat Loss")
barplot(100*t(Bplot[,-1]),beside=T,ylab="Percent decline in population growth",legend=T, args.legend=list(box.col=0))
#This part plots sensitivity results and requires all remaining code to be run first:
Lams=100*cbind(InsectGSA[,12:14],MusselGSA[,12:14],MusselLongGSA[,12:14],FishGSA[,12:14],FishLongGSA[,12:14])
segments(x0=0.5+c(1:19)[c(T,T,T,F)],y0=apply(Lams,2,quantile,0.25),y1=apply(Lams,2,quantile,0.75),lwd=4)




#Appendix 5: effects of bias in reproductive dispersal on washout effects and population growth:
UpstreamDispersBias=0
DispersBias=-0.35; InsectBF=relatives(parms=c(d=0.5,mu=0.3,rho=0.5,sigma=0,R=1.43,L=1e2,m=2,n=0.35,gamma=0.8),LHfun="Insect",HablossLevel=HL)
DispersBias=-2; MusselBF=relatives(parms=c(d=0.10,mu=0.3,rho=0.37,sigma=0,R=0.35,L=1e2,m=1,n=2,gamma=0),LHfun="Mussel",HablossLevel=HL)
DispersBias=-2; MusselLongBF=relatives(parms=c(d=0.025,mu=0.3,rho=0.37,sigma=0,R=0.125,L=1e2,m=1,n=2,gamma=0),LHfun="Mussel",HablossLevel=HL)
DispersBias=-5; FishBF=relatives(parms=c(d=0.25,mu=0.3,rho=0.5,sigma=0,R=0.5,L=1e2,m=4,n=5,gamma=0),LHfun="Fish",HablossLevel=HL)
DispersBias=-5; FishLongBF=relatives(parms=c(d=0.15,mu=0.3,rho=0.5,sigma=0,R=0.26,L=1e2,m=4,n=5,gamma=0),LHfun="Fish",HablossLevel=HL)
BplotBF=rbind(InsectBF,MusselBF,MusselLongBF,FishBF,FishLongBF); colnames(BplotBF)=c("BaseLam","Advective Dispersal","Habitat Loss","Advective Dispersal and Habitat Loss")
barplot(100*t(BplotBF[,-1]),ylim=c(0,60),beside=T,ylab="Percent decline in population growth",legend=T, args.legend=list(box.col=0))
Xval=0.5+(1:20)[c(T,T,T,F)]; points(Xval,100*t(Bplot[,-1]),cex=2,pch=16,col=2); points(Xval,100*t(Bplot[,-1]),cex=2,lwd=2,col=1);
DispersBias=0 #re-set global parameter





###Sensitivity analysis across life histories for Appendix C:
require(doParallel); require(foreach); require(randomForest);
parmsBound=function(parms){
	#ensure that L, m, and n are integers
	parms=round(parms,c(rep(4,5),rep(0,2),rep(4,3)))
	#ensure tested parameters are reasonable (e.g., <100% mortality and >0 reproduction)
	lbounds=c(d=0.01,mu=0.,rho=0,sigma=0,R=0.01,L=1e2,m=1,n=0.1,gamma=0.01,E=0)
	ubounds=c(d=0.98,mu=0.9,rho=1,sigma=0,R=6,L=1e2,m=10,n=10,gamma=1,E=0.98)
	parms[(parms-lbounds)<0]=lbounds[(parms-lbounds)<0]
	parms[(ubounds-parms)<0]=ubounds[(ubounds-parms)<0]
	return(parms)
}

relativesGSA=function(parms,LHfun,HablossLevel=0.85,nreps=1,nsims=4e3){
	GSAparms0=apply(cbind(0.45*c(parms,E=HablossLevel),1.55*c(parms,E=HablossLevel)), 1, function(x) runif(n=nsims,min=x[1],max=x[2]))
	GSAparms=t(apply(GSAparms0, 1, parmsBound))
	library(doParallel); library(foreach); registerDoParallel(cores=3);
	ParmsLambdas=foreach(n=seq(1:dim(GSAparms)[1]), .combine=rbind, .export=c("Lmatr","tries","SpatialMatr","relatives","setParms","offDiag")) %dopar%{ 
		c(GSAparms[n,],round(relatives(GSAparms[n,1:9],LHfun=LHfun,HablossLevel=GSAparms[n,10]),4)) }
	return(ParmsLambdas)
}

HL=0.85; Nsims=4e3;
InsectGSA=relativesGSA(parms=c(d=0.5,mu=0.3,rho=0.5,sigma=0,R=1.43,L=100,m=2,n=0.5,gamma=0.8),LHfun="Insect",HablossLevel=HL,nsims=Nsims)
MusselGSA=relativesGSA(parms=c(d=0.10,mu=0.3,rho=0.37,sigma=0,R=0.35,L=100,m=1,n=2,gamma=0),LHfun="Mussel",HablossLevel=HL,nsims=Nsims)
MusselLongGSA=relativesGSA(parms=c(d=0.025,mu=0.3,rho=0.37,sigma=0,R=0.125,L=100,m=1,n=2,gamma=0),LHfun="Mussel",HablossLevel=HL,nsims=Nsims)
FishGSA=relativesGSA(parms=c(d=0.25,mu=0.3,rho=0.5,sigma=0,R=0.5,L=100,m=4,n=5,gamma=0),LHfun="Fish",HablossLevel=HL,nsims=Nsims)
FishLongGSA=relativesGSA(parms=c(d=0.15,mu=0.3,rho=0.5,sigma=0,R=0.26,L=100,m=4,n=5,gamma=0),LHfun="Fish",HablossLevel=HL,nsims=Nsims)


#Analyze and plot GSA results:
GSAsummary=function(ParmsLambdas,VarsInclude=1:4,parmsExclude=c(4,6),plotgive=NULL){
	testGSA=as.data.frame(ParmsLambdas)[,(1:10)[-parmsExclude]]; Vars=tail(colnames(ParmsLambdas),4);
	ImportancesAll=matrix(0,nrow=4,ncol=(dim(testGSA)[2]+3)); colnames(ImportancesAll)=c(names(testGSA),"Q1","MeanEffect","Q2"); rownames(ImportancesAll)=Vars;
	for(i in VarsInclude){
		Response=ParmsLambdas[,Vars[i]]; EffectRng=round(c(quantile(Response,0.25),mean(Response),quantile(Response,0.75)),2);
		RFresult=randomForest(Response ~ ., data=testGSA, importance=TRUE)
		Importances=pmax(importance(RFresult,scale=FALSE)[,1],0); ImportancesAll[i,]=c(Importances/sum(Importances),EffectRng);
		if(plotgive==1) ImpOrd=(order(Importances,decreasing=F)) #In this case, parameters plotted in order of importance
		if(!is.null(plotgive)) barplot(Importances[ImpOrd],names=(parmnames[-parmsExclude])[ImpOrd],horiz=T,main=paste0("||",Vars[i],"||=(",EffectRng[1],", ",EffectRng[2],", ",EffectRng[3],")"), las=1,col="Indianred2",xlab="Relative importance value")
	}
	return(ImportancesAll)
}
parmnames=c("Mortality d",expression(paste("AD mortality ",mu)),expression(paste("AD frequency ",rho)),expression(paste("AD refugia ",sigma)),
   "Recruitment rate R","Domain size L","AD distance m","Recruitment distance n",expression(paste("Maturation rate ",gamma)),"Habitat loss E")
linch=max(strwidth(parmnames,"inch"),na.rm=T); par(mai=c(1.02,linch,0.82,0.42),mfrow=c(3,2))
dev.new(); par(mfrow=c(3,2)); 
GSAsummary(InsectGSA,plotgive=1,VarsInclude=4);
GSAsummary(MusselGSA,plotgive=1,VarsInclude=4,parmsExclude=c(4,6,9))
GSAsummary(MusselLongGSA,plotgive=1,VarsInclude=4,parmsExclude=c(4,6,9))
GSAsummary(FishGSA,plotgive=1,VarsInclude=4,parmsExclude=c(4,6,9))
GSAsummary(FishLongGSA,plotgive=1,VarsInclude=4,parmsExclude=c(4,6,9))






