### Initialization of parameters

library(rARPACK); #install.packages('rARPACK');
parms = c(d = 0.1, mu = 0.3, rho = 0.37, R = 0.25, L=100, S=2, m=1, n=2, g=0.5); #initialization of model parameters

#function for changing parameter value.
setParms = function(d = parms["d"], mu = parms["mu"], rho = parms["rho"], R = parms["R"], L = parms["L"], S = parms["S"], m = parms["m"], n = parms["n"], g = parms["g"]){
	parms["d"] = d; #natural mortality
	parms["mu"] = mu; #dispersal mortality
	parms["rho"] = rho; #advective dispersal rate
	parms["R"] = R; #recruitment rate
	parms["L"] = L; #number of patches
	parms["S"] = S; #number of stages
	parms["m"] = m; #distance of advective dispersal 
	parms["n"] = n; #distance of fish dispersal
	parms["g"] = g; #rate of stage transition
	return(parms);
}

#initialization of functional parameters
StageStruct = 0; # 1 for structrued model and 0 for unstructrued model
NoWashout = 0; #dispersal on a torus: mussels being washed out from patchs furthest downstream are transported to the furthest upstream patches in a loop
propInvalid = 0; #proportion of permanent patch loss due to human disturbance
propExt = 0; #proportion of temporary patch loss due to natural disturbance
propPerm = 0; #proportion of patches invulnarable to natural disturbance
runl = 100; #time period of simulation
runt = 50; #time period of simulation to eliminate transient dynamics



###Construction of projection matrix L

#calculate vec-permutation matrix
Pmatr = function(parms){
	#initialize matrices
	p = matrix (0, parms["L"] * parms["S"], parms["L"] * parms["S"]); 
	E = matrix(0, parms["S"], parms["L"]); 
	#calculate according to Caswell paper
	for(i in 1 : parms["S"]){
		for (j in 1 : parms["L"]){
			E[i,j] = 1;
			p = p + kronecker(E, t(E));
			E[,] = 0;
		}
	}
	return (p)
}

p = Pmatr(parms); #initialize vec-permutation matrix

#Create projection matrix
Lmatr = function(parms, propInvalid = 0, StageStruct = 0, NoWashout = 0){
	if (length(p) != parms["L"] * parms["S"] * parms["L"] * parms["S"]) p = Pmatr(parms); #Recalculate corresponding vec-permutation matrix for changed L
	
	#define stage-related parameters and matrices
	if(StageStruct == 0){ #Unstructured model
		R1 = parms["R"]; R2 = parms["R"]; #R1 = R2 = R
		U = matrix(c((1 - parms["d"]), 0, 0, (1-parms["d"])), ncol=2, byrow=TRUE); #g = 0
		F = matrix(c((1 - parms["d"]) * R1, (1 - parms["d"]) * R2, 0, 0), ncol=2, byrow=TRUE);	
	} 
	if(StageStruct == 1){ #2-stage structured
		R1 = 0; R2 = parms["R"]; 
		U = matrix(c((1 - parms["d"]) * (1 - parms["g"]), 0, (1 - parms["d"]) * parms["g"], (1 - parms["d"])), ncol=2, byrow=TRUE);
		F = matrix(c((1 - parms["d"]) * R1, (1 - parms["d"]) * R2, 0, 0), ncol = 2, byrow = TRUE);	
	}
	if(StageStruct == 2){ #3-stage structured
		parms["S"] = 3;
		R1 = 0; R2 = 0.5 * parms["R"]; R3 = parms["R"]; 
		U = matrix(0, parms["S"], parms["S"]);
		U[1, 1] = U[2, 2] = (1 - parms["d"]) * (1 - parms["g"]); U[2,1] = U[3,2] = (1 - parms["d"]) * parms["g"]; U[3, 3] = (1 - parms["d"]);
		F = matrix(0, parms["S"], parms["S"]);
		F[1, ] = c((1 - parms["d"]) * R1, (1 - parms["d"]) * R2, (1 - parms["d"]) * R3);
	}
	
	#create spatial transition matrices
	#advective dispersal
	D = parms["rho"] * (1 - parms["mu"]) / parms["m"];
	Du = diag(parms["L"]) * (1 - parms["rho"]);
	if(parms["m"] > 0 && parms["L"] > 1) {
		for(i in 1 : (parms["L"] - 1)){
			Du[(i + 1) : min(parms["L"], i+parms["m"]), i] = D;
		};
		if (NoWashout == 1) { #turn off effects of adults being washed out from last patches by connecting the patches in a loop
			for (j in (parms["L"] - parms["m"] + 1) : parms["L"]) Du[1 : (parms["m"] - (parms["L"] - j)), j] = D;
		}
	}
	#fish dispersal
	Df = diag(parms["L"]) * 1 / (2 * parms["n"] + 1);
	if(parms["n"] > 0 && parms["L"] > 1) {	
		for(i in 1 : parms["L"]){
			Df[max(1, i-parms["n"]): min(parms["L"], i+parms["n"]), i] = 1 / (2 * parms["n"] + 1);
		};
	}

	#For testing effects of making patches permanently unviable (but invalid patches still trap advective dispersers)
	if(propInvalid > 0){  
		Npatch = as.numeric(parms["L"]); 
		Delpatches = sample( #generate a sequence of patches ramdomly chosen to become invalid
			c(rep(1, round(Npatch*propInvalid,0)), rep(0,(Npatch-round(Npatch*propInvalid,0))))
		) * (1:Npatch);
		Du[, Delpatches] = 0; Df[, Delpatches] = 0; #all individuals transported to selected patches are dead
	};
	
	return (t(p) %*% (kronecker(diag(parms["S"]), Du)) %*% p %*% (kronecker(diag(parms["L"]), U))
			+ t(p) %*% (kronecker(diag(parms["S"]),Df)) %*% p %*% (kronecker(diag(parms["L"]), F)))
}




###Calculate the dominant eigenvalue of projection matrix

tries = function(parms, out = 1, propInvalid = 0, StageStruct = 0, NoWashout = 0){
		L=Lmatr(parms, propInvalid, StageStruct, NoWashout)
		if(out == 0){ return(Re(eigs(L,1,which = "LR")$values)) } 
		if(out == 1){ return(max(Re(eigen(L)$values))) }
		if(out == 2){ return(Re(eigen(L)$vectors[,1])) } 
		if(out == 3){ return(which.max(nrm(Re(eigen(L)$vectors[,1])))) }; 
		if(out == 4){ for(t in 2:runl) N[t,] = L%*%N[t-1,]; return(mean(N[runl,]/N[runl-1,])); }; 
}



###Calculate population growth rate from numerical simulations for stochastic model

#generate a matrix with 0 corresponding entries for disturbed patches each year
StochE = function(parms, propExt = 0, propPerm = 0){
	Es = matrix(1 - propExt, 1, parms["L"]); 
	E = apply(Es, 2, function(x) rbinom(runl + runt, 1, prob=x)); #proportion of peroExt patches are choosen to be disturbed each year (by row)
	E[,rbinom(parms["L"], 1, propPerm) * (1 : parms["L"])] = 1; #choose proportion of propPerm patches to be invulnarable to disturbance
	return(E);
}

#Calculate population growth rate																																													on growth rate from numerical simulations for stochastic model
StochLambda = function(parms, propInvalid = 0, StageStruct = 0, E){ 
	L = Lmatr(parms, propInvalid, StageStruct); #create projection matrix
	
	#initialization: each patch has the same number of adults and no juveniles
	N = matrix(0, parms["L"] * parms["S"], runl + runt); 
	N[, 1] = rep(c(0, 1 / (parms["L"])), parms["L"]);
	
	lambda = rep(0, runl); #initialize a vector to store growth rate of each year
	
	#simulate for runt years to eliminate transient dynamics
	for(t in 2 : runt) {
		N[, t] = L %*% N[, t-1] * rep(E[t, ], each = 2);
		tempLambda = sum(N[,t]); #calculate population growth rate of year t
		N[, t] = N[, t] / tempLambda; #normalization for next year
	}
	
	#simulate for runl years
	for(t in (runt + 1) : (runl + runt)) {
		N[, t] = L %*% N[, t-1] * rep(E[t,], each = 2);
		lambda[t - runt] = sum(N[, t]); #calculate population growth rate of year t
		N[, t] = N[, t] / lambda[t - runt]; #normalization for next run
	}
	return(exp(mean(log(lambda)))); #calculate long-term population growth rate by calculating the geometric mean of growth rate for each year
}






#########PART I: EFFECT OF FLOW DISPERSAL ON POP GROWTH#########

#initialization
parms = setParms(d = 0.1,  mu = 0.3, rho = 0.37, R = 0.25, L = 100, S = 2, m = 1, n = 2, g = 0.5)
L = matrix(0, parms["L"] * parms["S"], parms["L"] * parms["S"]);
mus = c(0, 0.3, 0.6, 1); 
ntests = 4; #number of rhos tested
Rhos = seq(0, 0.5, length = ntests);

#calculate difference between lambda and 1
LL2 = function(Rtry, rhotry, mutry, StageStruct) {
	parms = setParms(mu = mutry, rho = rhotry, R = Rtry);
	abs(1 - tries(parms, StageStruct = StageStruct))
}

#calculate required R for lambda = 1
fitR2 = function(rhotry, mutry, StageStruct = 0) 
	optim(par = 0.25, fn = LL2, method = "L-BFGS-B", lower = .01, upper = 4, control = list(maxit = 100), 
		rhotry = rhotry, mutry = mutry, StageStruct = StageStruct) $par;

#plot required R vs rho for mu = 0, 0.3, 0.6, 1
plotF1 = function(StageStruct = 0){
	Toplot=matrix(0, ntests, length(mus));
	for(i in 1 : ntests){ 
		for(j in 1 : length(mus)){ 
			Toplot[i, j] = fitR2(Rhos[i], mus[j], StageStruct = StageStruct);
		}	
	}
	return (Toplot);
}

###Unstructured model
dev.new()
par(mfrow = c(1,2));
Toplot1 = plotF1();
matplot(Rhos, Toplot1, type = "b", lty = 1, lwd = 2, col = c(3, 4, 2, 1), pch = 1, ylim = c(0, 1.2),  
	main = "Unstructured model",  ylab = "Per capita fecundity R required for persistence", xlab = "Proportion of individuals annually dispersed by flow,  Rho");
legend("topleft", lty = 1, col = c(1, 2, 4, 3), lwd = 2, seg.len = 3, box.col = 0, 
	c("Dispersal mortality 100%",  "Dispersal mortality 60%",  "Dispersal mortality 30%",  "Dispersal mortality 0%"))

###Structured model
Toplot2 = plotF1(StageStruct = 1);
matplot(Rhos, Toplot2, type = "b", lty = 1, lwd = 2, col = c(3, 4, 2, 1), pch = 1, ylim = c(0, 4),  
	main = "Two-Stage Structured",  ylab = "Per capita fecundity R required for persistence", xlab = "Proportion of individuals annually dispersed by flow,  Rho");
legend("topleft", lty = 1, col = c(1, 2, 4, 3), lwd = 2, seg.len = 3, box.col = 0, 
	c("Dispersal mortality 100%",  "Dispersal mortality 60%",  "Dispersal mortality 30%",  "Dispersal mortality 0%"))

###Plot both for compare
dev.new()
#unstructured
matplot(Rhos, Toplot1, type = "b", lty = 2, lwd = 2, col = "grey", pch = c(1, 17, 22, 4), ylim = c(0, 4),  
	ylab = "Per capita fecundity R required for persistence", xlab = "Proportion of individuals annually dispersed by flow,  Rho");
legend("topleft", lty = 1, col = 1, pch = c(4, 22, 17, 1), lwd = 2, seg.len = 3, box.col = 0,  
	c("Dispersal mortality 100%",  "Dispersal mortality 60%",  "Dispersal mortality 30%",  "Dispersal mortality 0%"))
#structured
par(new = TRUE)
matplot(Rhos, Toplot2, type = "b", lty = 1, lwd = 2, col = 1, pch = c(4, 22, 17, 1), ylim = c(0, 4),  
	ylab = "Per capita fecundity R required for persistence", xlab = "Proportion of individuals annually dispersed by flow,  Rho");
legend(-0.025, 3.3, lty = c(1, 2),  col = c(1, "grey"), lwd = 2, seg.len = 3, box.col = 0,  
	c("Structured with Transition rate 0.5",  "Unstructured"))

	
	
	
	
	
	
	
#########PART II:  EFFECT OF PATCH LOSS AND SMALL DOMAINS ON POP GROWTH#########

#initialization
parms = setParms(d = 0.1,  mu = 0.3, rho = 0.37, R = 0.25, L = 100, S = 2, m = 1, n = 2, g = 0.5)
L = matrix(0, parms["L"] * parms["S"], parms["L"] * parms["S"]);
Rhos = c(0., 0.25, .5); 
Mu = parms["mu"];

#plot lambda vs proportion of patch loss
ntests = 15; #number of tests for increasing proportion of invalid patches
nreps = 15; #number of repeats for the same level of propInvalid
propsInvalids = seq(0, 0.9, length = ntests); 
plotF2a = function (StageStruct = 0, NoWashout = 0, muZero = 0, compare = 0) {
	Toplot3 = array(0, dim = c(ntests, length(Rhos), nreps)); #plot for patch loss
	if (parms["L"] != 100) {
		parms = setParms(L = 100);
		L = matrix(0, parms["L"] * parms["S"], parms["L"] * parms["S"]);
	}
	if (compare == 1) Rhos = 0.37; #when plotting for comparison, set rho to default value
	if (muZero == 1) Mu = 0; #turn off effects of washout
	if (muZero == 1 && NoWashout == 1) Rhos = 0; #turn off dispersal
	for(i in 1 : ntests){ 
		for(j in 1 : length(Rhos)){ 
			for(k in 1 : nreps){ #calculate lambda
				parms = setParms(rho = Rhos[j], mu = Mu);
				Toplot3[i, j, k] = tries(parms, propInvalid = propsInvalids[i], StageStruct = StageStruct, NoWashout = NoWashout); 
			} 
		}
	}
	return(Toplot3);
}


#plot lambda vs patch number in a small domain
maxPatches = 15; 
Ls = 1 : maxPatches; #number of patches for test
plotF2b = function(StageStruct = 0, NoWashout = 0, muZero = 0, compare = 0){
	Toplot = matrix(0, length(Ls), length(Rhos)); #plot for small domains
	if (compare == 1) Rhos = 0.37;
	if (muZero == 1) Mu = 0; #turn off effects of washout
	if (muZero == 1 && NoWashout == 1) Rhos = 0; #turn off dispersal
	for(i in 1 : length(Ls)){ 
		for(j in 1 : length(Rhos)){ #calculate lambda
			parms = setParms(L = Ls[i], rho = Rhos[j], mu = Mu); #change patch number and rho
			Toplot[i, j] = tries(parms, propInvalid = 0, StageStruct = StageStruct, NoWashout = NoWashout);
		} 
	}
	return(Toplot);
}


#Unstructured model & two-stage structured model
dev.new(); 
par(mfrow = c(1, 2));
Toplota1 = plotF2a(); #unstructured, increasing patch loss
Toplota2 = plotF2a(StageStruct = 1); #structured, increasing patch loss
matplot(propsInvalids, apply(Toplota1, c(1, 2), mean), type = "l", lty = 1, lwd = 3, pch = 1, ylim = c(0.5, 1.2), yaxt = "n", 
	ylab = "Population growth rate", xlab = "Proportion of patches lost"); axis(2, at = seq(0.5, 1.4, by = 0.2));
par(new = TRUE);
matplot(propsInvalids, apply(Toplota2, c(1, 2), mean), type = "l", lty = 2, lwd = 3, pch = 1, ylim = c(0.5, 1.2), yaxt = "n", 
	ylab = "Population growth rate", xlab = "Proportion of patches lost"); axis(2, at = seq(0.5, 1.4, by = 0.2));
legend("topleft", lty = 1, col = c(1, 2, 3), lwd = 2, seg.len = 3, box.col = 0, 
	c("Dispersal rate 0%", "Dispersal rate 25%", "Dispersal rate 50%"))
legend("bottomleft", lty = c(1, 2), lwd = 2, seg.len = 3, box.col = 0, 
	c("Unstructured model", "Two-stage structured model"))	

Toplotb1 = plotF2b(); #unstructured, increasing patch number in a small domain
Toplotb2 = plotF2b(StageStruct = 1); #structured, increasing patch number in a small domain
matplot(Ls, Toplotb1, type = "b", lty = 1, lwd = 2, pch = 1, ylim = c(0.5, 1.2), yaxt = "n", xlim = c(0, maxPatches), 
	ylab = "Population growth rate", xlab = "Number of patches in population, L"); axis(2, at = seq(0.5, 1.3, by = 0.2));
par(new = TRUE);
matplot(Ls, Toplotb2, type = "b", lty = 1, lwd = 2, pch = 17, ylim = c(0.5, 1.2), yaxt = "n", xlim = c(0, maxPatches), 
	ylab = "Population growth rate", xlab = "Number of patches in population, L"); axis(2, at = seq(0.5, 1.3, by = 0.2));
legend("topleft", lty = 1, col = c(1, 2, 3), lwd = 2, seg.len = 3, box.col = 0, 
	c("Dispersal rate 50%", "Dispersal rate 25%", "Dispersal rate 0%"))
legend("bottomright", lty = 1, pch = c(1, 17), lwd = 2, seg.len = 3, box.col = 0, 
	c("Unstructured model", "Two-stage structured model"))	




#compare effects of mu and washout on unstructured model
Ls = (parms["m"] + 1) : maxPatches; #cannot turn off washout if patch number < dispersal distance
dev.new();
par(mfrow = c(1, 2));
#patch loss
ToplotA1 = array(0, dim = c(4, ntests, 3, nreps));
ToplotA1[1, , , ] = plotF2a (compare = 1); #all effects
ToplotA1[3, , , ] = plotF2a (muZero = 1, compare = 1); #only washout
ToplotA1[4, , , ] = plotF2a (NoWashout = 1, muZero = 1, compare = 1); #no dispersal
ToplotA1[2, , , ] = ToplotA1[4, , , ] - (ToplotA1[3, , , ] - ToplotA1[1, , , ]); #only mu (also turn off washout from patches in middle to invalid patches)
for (i in 1 : 4) {
	if (i > 1) par(new = TRUE);
	matplot(propsInvalids, apply(ToplotA1[i, , 1, ], 1, mean), type = "l", lty = 1, col = i, lwd = 3, pch = 1, ylim = c(0.7, 1.3), yaxt = "n", 
	ylab = "Population growth rate", xlab = "Proportion of patches lost"); axis(2, at = seq(0.7, 1.3, by = 0.2));
}
legend("topleft", lty = 1, col = c(1, 2, 3, 4), lwd = 2, seg.len = 3, box.col = 0, 
	c("All effects of dispersal", "Only effect of dispersal mortality", "Only effect of washout", "No dispersal"))
	
#small domains
maxPatches = 40;
Ls = 1 : maxPatches; #number of patches for test
ToplotB1 = array(0, dim = c(4, length(Ls), 3));
ToplotB1[1, , ] = plotF2b (compare = 1); #all effects
ToplotB1[2, , ] = plotF2b (NoWashout = 1, compare = 1); #only mu
ToplotB1[3, , ] = plotF2b (muZero = 1, compare = 1); #only washout
ToplotB1[4, , ] = plotF2b (NoWashout = 1, muZero = 1, compare = 1); #no dispersal
for (i in 1 : 4) {
	if (i > 1) par(new = TRUE);
	matplot(Ls, ToplotB1[i, , 1], type = "l", lty = 1, col = i, lwd = 2, pch = 1, ylim = c(0.5, 1.2), yaxt = "n", xlim = c(0, maxPatches), 
	ylab = "Population growth rate", xlab = "Number of patches in population, L"); axis(2, at = seq(0.5, 1.3, by = 0.2));
}
legend("bottomright", lty = 1, col = c(1, 2, 3, 4), lwd = 2, seg.len = 3, box.col = 0, 
	c("All effects of dispersal", "Only effect of dispersal mortality", "Only effect of washout", "No dispersal"))
	

#compare difference between structured and unstructured model
Ls = (parms["m"] + 1) : maxPatches; #cannot turn off washout if patch number < dispersal distance
dev.new();
par(mfrow = c(1, 2));
#patch loss
ToplotA2 = array(0, dim = c(4, ntests, 3, nreps));
ToplotA2[1, , , ] = plotF2a (StageStruct = 1, compare = 1); #all effects
ToplotA2[3, , , ] = plotF2a (StageStruct = 1, muZero = 1, compare = 1); #only washout
ToplotA2[4, , , ] = plotF2a (StageStruct = 1, NoWashout = 1, muZero = 1, compare = 1); #no dispersal
ToplotA2[2, , , ] = ToplotA2[4, , , ] - (ToplotA2[3, , , ] - ToplotA2[1, , , ]); #only mu (also turn off washout from patches in middle to invalid patches)
A = ToplotA1 - ToplotA2; #plot difference between structured and unstructured model
for (i in 1 : 4) {
	if (i > 1) par(new = TRUE);
	matplot(propsInvalids, apply(A[i, , 1, ], 1, mean), type = "l", lty = 1, col = i, lwd = 3, pch = 1, ylim = c(0, 0.15), yaxt = "n", 
	ylab = "Difference in population growth rate", xlab = "Proportion of patches lost"); axis(2, at = seq(0, 0.15, by = 0.05));
}

#small domains
maxPatches = 40;
Ls = 1 : maxPatches; #number of patches for test
ToplotB2 = array(0, dim = c(4, length(Ls), 3));
ToplotB2[1, , ] = plotF2b (StageStruct = 1, compare = 1); #all effects
ToplotB2[2, , ] = plotF2b (StageStruct = 1, NoWashout = 1, compare = 1); #only mu
ToplotB2[3, , ] = plotF2b (StageStruct = 1, muZero = 1, compare = 1); #only washout
ToplotB2[4, , ] = plotF2b (StageStruct = 1, NoWashout = 1, muZero = 1, compare = 1); #no dispersal
B = ToplotB1 - ToplotB2; #plot difference between structured and unstructured model
for (i in 1 : 4) {
	if (i > 1) par(new = TRUE);
	matplot(Ls, B[i, , 1], type = "l", lty = 1, col = i, lwd = 2, pch = 1, ylim = c(0, 0.15), yaxt = "n", xlim = c(0, maxPatches), 
	ylab = "Difference in opulation growth rate", xlab = "Number of patches in population, L"); axis(2, at = seq(0, 0.15, by = 0.05));
}
legend("bottomright", lty = 1, col = c(1, 2, 3, 4), lwd = 2, seg.len = 3, box.col = 0, 
	c("All effects of dispersal", "Only effect of dispersal mortality", "Only effect of washout", "No dispersal"));
	
#print the difference for largest patch number
print(ToplotB1[1, length(Ls), 1] - ToplotB2[1, length(Ls), 1]);
print(ToplotB1[2, length(Ls), 1] - ToplotB2[2, length(Ls), 1]);
print(ToplotB1[3, length(Ls), 1] - ToplotB2[3, length(Ls), 1]);
print(ToplotB1[4, length(Ls), 1] - ToplotB2[4, length(Ls), 1]);








#########PART III: METAPOP DYNS

#initialization
parms = setParms(d = 0.1, mu = 0.3, rho = 0.37, R = 0.25, L=100, S=2, m=1, n=2, g=0.5);
L = matrix(0, parms["L"]*parms["S"], parms["L"]*parms["S"]);
ntests = 10; #number of tests for increasing disturbance frequency
nreps = 10; #number of repeats for each level of disturbance frequency
EPs = seq(0.01, 0.5, length = ntests); #levels of disturbance frequencies for test
Rhos = c(0., 0.25, .5); #levels of rhos for test
RndmStoch3 = StochWithPerms3 = array(0, dim = c(ntests, length(Rhos), nreps)); 
runl = 1e2; #simulation time length

#plot lambda vs disturbance frequency
F3 = function (StageStruct = 0){
	#plot for random disturbance 
	for(i in 1:ntests){ 
		for(k in 1:nreps){ 
			E1 = StochE (parms, propExt = EPs[i], propPerm = 0); #calculate vectors of disturbance distribution
			for (j in 1:length(Rhos)){
				parms = setParms(rho = Rhos[j]); #change level of rho
				RndmStoch3[i, j, k] = StochLambda(parms, StageStruct = StageStruct, E = E1); #calculate lambda from simulation
			}
		} 
	}
	#plot for cases where some patches are invulneralble
	for(i in 1:ntests){ 
		for(k in 1:nreps){ 
			E2 = StochE (parms, propExt = EPs[i], propPerm = 0.3); #calculate vectors of disturbance distribution
			for (j in 1:length(Rhos)){
				parms = setParms(rho = Rhos[j]); #change level of rho
				StochWithPerms3[i, j, k] = StochLambda(parms, StageStruct = StageStruct, E = E2); #calculate lambda from simulation
			}
		} 
	}
	dev.new(); 
	par(mfrow = c(1, 2)); 
	include = 1:(ntests-2);
	matplot(EPs[include], apply(RndmStoch3, c(1, 2), mean)[include, ], type = "l", lty = 1, lwd = 3, pch = 1, ylim = c(0.5, 1.15), yaxt = "n", 
		ylab = "Population growth rate", xlab = "Proportion of patches lost annually"); axis(2, at = seq(0.5, 1.15, by = 0.2));
	matlines(EPs[include], apply(RndmStoch3, c(1, 2), function(x) mean(x)-sd(x))[include, ], type = "l", lty = 2, lwd = 1); 
	matlines(EPs[include], apply(RndmStoch3, c(1, 2), function(x) mean(x)+sd(x))[include, ], type = "l", lty = 2, lwd = 1);
	legend("topleft", lty = 1, col = c(1, 2, 3), lwd = 2, seg.len = 3, box.col = 0, 
		c("Dispersal rate 0%", "Dispersal rate 25%", "Dispersal rate 50%"))
	matplot(EPs[include], apply(StochWithPerms3, c(1, 2), mean)[include, ], type = "l", lty = 1, lwd = 3, pch = 1, ylim = c(0.5, 1.15), yaxt = "n", 
		ylab = "Population growth rate", xlab = "Proportion of patches lost annually"); axis(2, at = seq(0.5, 1.15, by = 0.2));
	matlines(EPs[include], apply(StochWithPerms3, c(1, 2), function(x) mean(x)-sd(x))[include, ], type = "l", lty = 2, lwd = 1); 
	matlines(EPs[include], apply(StochWithPerms3, c(1, 2), function(x) mean(x)+sd(x))[include, ], type = "l", lty = 2, lwd = 1);
}

###Unstructured
F3();
title("Unstructured model");

###structured
F3(StageStruct = 1);
title("Two-stage structured model");




