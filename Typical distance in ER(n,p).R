# Studying the Typical Distance in Erdős–Rényi Random Graph through Simulations
# a project by Aditya Ghosh and Sayak Chatterjee, done under the supervision of Prof. Antar Bandyopadhyay

library(igraph)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# 				   Sparse but super-critical Regime (p = c/n, c > 1) 
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-----------------------------------------------------------------------------------
# 								I. Data Generation
#-----------------------------------------------------------------------------------

size = c(20, 60, 100, 150, 250, 400, 675, 1000, 2000)
a = seq(1.1, 2.5, by = 0.2)
rpt = 1000 		
PlayER<-function(n, p, N = 1000){
	d = rep(0,N)
	for(i in 1:N){
		G = sample_gnp(n, p)
		s = sample(1:n, 2)
		d[i] = distances(G, s[1], s[2])
	}
	return(d)
}
# PlayER generates realizations of H_n for ER(n, p) for N times and returns the vector of the observed typical distances

set.seed(140)
num = matrix(nrow = length(size), ncol = length(a), dimnames=list(size,a)) 
# num stores how many times the observed distance is finite
D = array(dim = c(rpt, length(size), length(a)), dimnames = list(1:rpt,size,a)) 
# D is the 3D array for storing the observed typical distances
for(j in 1:length(a)){
	for(i in 1:length(size)){
	n = size[i]
	p = a[j]/n
	d = PlayER(n, p, N = rpt)
	D[,i,j] = d
	num[i,j] = length(subset(d, d<Inf))
	}
}
# Next we calculate the sample mean and s.d. for each pair of n and c = np
mcbyn = matrix(nrow=length(size),ncol=length(a))
for(j in 1:length(a)) mcbyn[,j] = apply(D[,,j], 2, function(dat) mean(dat[dat<Inf]))
dimnames(mcbyn) = list(size, a)
scbyn = matrix(nrow=length(size),ncol=length(a))
for(j in 1:length(a)) scbyn[,j] = apply(D[,,j], 2, function(dat) sd(dat[dat<Inf]))
dimnames(scbyn) = list(size, a)
# Saving the data for future use:
for(j in 1:length(a))
write.csv(D[ , , j],paste("yourpath/cbyn",a[j],".csv",sep=""), row.names = FALSE)

#-----------------------------------------------------------------------------------
# 								    II. Histograms
#-----------------------------------------------------------------------------------
	
for(j in 1:length(a)){
	jpeg(file=paste("hist",a[j],".jpeg",sep=""), width=1200, height=800, pointsize = 20, quality=100)
	layout(matrix(1:9,nrow=3,byrow=T))
	for(i in 1:length(size)){
		d = D[,i,j][!is.infinite(D[,i,j])]
		hist(d, breaks=seq(min(d-0.5),max(d+0.5),by=1), col="peachpuff", 
		border="black", prob = TRUE, xlab = "typical distance", 
		main = paste("for n =",size[i]," & p =",a[j],"/n"))
		x = seq(min(d-0.5), max(d+0.5), length=1000)
		lines(x,dnorm(x,mean=mcbyn[i,j],sd=scbyn[i,j]), lwd=2, col="dark blue")
		lines(x,dnorm(x,mean=log(size[i])/log(a[j]),sd=scbyn[i,j]), lwd=2, col="dark red")
	}
	dev.off()
}
# Standardized Histograms 
for(j in 1:length(a)){
	jpeg(file=paste("sdhist",a[j],".jpeg",sep=""), width=1200, height=800, pointsize = 20, quality=100)
	layout(matrix(1:9,nrow=3,byrow=T))
	for(i in 1:length(size)){
		d = D[,i,j][!is.infinite(D[,i,j])]
		d = (d - mean(d))/sd(d)
		h = hist(d, breaks=round(log(size[i]))*2-1, col="peachpuff", 
		border="black", prob = TRUE, xlab = "typical distance", 
		main = paste("for n =",size[i]," & p =",a[j],"/n"))
		xfit = seq(min(h$breaks),max(h$breaks),length=40)
		yfit = dnorm(xfit, mean=0, sd=1)
		lines(xfit, yfit, col="dark blue", lwd=2)
	}
	dev.off()
}

#-----------------------------------------------------------------------------------
# 								 III. Growth Plots
#-----------------------------------------------------------------------------------

error = mcbyn
for(i in 1:length(size)) 
	for(j in 1:length(a)) 
		error[i,j] = mcbyn[i,j]*log(a[j])/log(size[i]) - 1 
# error is the matrix that stores the o(1) term (observed)

quartz.options(width=10, height=7)
par(family="System Font")
k = length(a)
#-------------- Plotting the o(1) term = (sample mean*log(c)/log(n) - 1)
matplot(error, ylim=c(min(error),max(error)), type = "o", pch=20, lty=1, lwd=3, col = rainbow(k+1)[-5], xlab="Graph size n", main="The o(1) term in the theoretical mean of the typical distance in ER(n, c/n)", ylab="sample mean*log(c)/log(n) - 1", xaxt = "n")
axis(side=1,at=1:length(size),labels=size)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright",inset=c(-0.18,0), legend=paste("c = ",rev(a),sep=""),
		col=rev(rainbow(k+1)[-5]),lwd=3)

#-------------- Plotting the sample s.d. (unscaled)
matplot(scbyn, ylim=c(min(scbyn),max(scbyn)), type = "o", pch=20, lty=1, lwd=3, col = (rainbow(k+1)[-5]), xlab="Graph size n", main="Sample s.d. for the typical distance in ER(n, c/n)", ylab="Sample s.d.", xaxt = "n")
axis(side=1,at=1:length(size),labels=size)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright",inset=c(-0.18,0), legend=paste("c = ",a,sep=""),
		col=(rainbow(k+1)[-5]),lwd=3)

#-------------- Plotting the proportion of times the typical distance is finite
matplot(num/rpt, type="o", pch=20, lty=1, lwd=3, col = rainbow(k+1)[-5], xlab="Graph size n", main="Proportion of times the typical distance is finite, for p = c/n (where c > 1)", ylab="Proportion of times the typical distance is finite", xaxt = "n", yaxt = "n")
axis(side=1, at=1:length(size), labels=size)
axis(side=2, at=seq(0.1, max(num/rpt), by=0.1))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright", inset=c(-0.18,0), legend=paste("c =",rev(a),sep=" "), col=rev(rainbow(k+1)[-5]), lty=1, lwd=3)

#-------------- Column mean and s.d. of the '#times d is finite' matrix
par(mfrow=t(1:2), mai = c(1, 0.5, 0.4, 0.4))
plot(apply(num, 2, mean)~a, main="column means", xlab="value of c=np", col="Dark blue", pch=20, cex=1.5, xaxt="n")
axis(side=1, at=a)
plot(apply(num, 2, sd)~a, main="column sd", xlab="value of c=np", col="Maroon", pch=20, cex=1.5, xaxt="n")
axis(side=1, at=a)

#-----------------------------------------------------------------------------------
#								IV. Testing Normality
#-----------------------------------------------------------------------------------

require(nortest)
mychisq.test<-function(dat){
	sdat = (dat - mean(dat))/sd(dat)
	f.os = table(cut(sdat, breaks=c(-Inf,-2:2,Inf)))
	f.ex = (pnorm(c(-2:2,Inf)) - pnorm(c(-Inf,-2:2)))*length(sdat)
	chi.stat = sum((f.os-f.ex)^2/f.ex)
	return(pchisq(chi.stat, df=length(f.os)-1, lower.tail=FALSE))
}
myks.test<-function(dat){
	sdat = (dat - mean(dat))/sd(dat) 
	jdat = sdat + rnorm(length(sdat), mean=0, sd=0.001) # jittered to remove ties
	y = rnorm(length(sdat))
	pval = ks.test(jdat, y)$p.value
	return(pval)
}
norm.testing<-function(dat){
	sw = shapiro.test(sdat)$p.value
	sf = sf.test(sdat)$p.value
	ad = ad.test(sdat)$p.value
	cvm = cvm.test(sdat)$p.value
	ks = myks.test(sdat)
	psn = mychisq.test(sdat)
	return(c(psn, ks, sw, cvm, sf, ad))
}
testlist = c("Pearson chi-square", "Kolmogorov-Smirnov", "Shapiro-Wilk", "Cramer-von Mises", "Shapiro-Francia", "Anderson-Darling")
P = array(dim = c(6, length(size), length(a)), dimnames = list(testlist,size,a)) 
for(j in 1:length(a)){
	for(i in 1:length(size)){
	dat = D[, i, j][!is.infinite(D[, i, j])] 
	sdat = (dat - mean(dat))/sd(dat)
	P[,i,j] = norm.testing(sdat)
	}
}   # P is a 3D array for storing the p-values of normality testing
for(i in 1:6){
	cat("\n Name of the test : ", testlist[i],"\n\n")
	print(round(P[i, ,], 5))
	cat("\n")
}

#-----------------------------------------------------------------------------------
# We repeated above for dist(1,v) where v is randomly chosen. The results were very similar to above.
#-----------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# 							   Connectivity Regime								   
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-----------------------------------------------------------------------------------
# 							    I. Data Generation
#-----------------------------------------------------------------------------------

size = c(20,60,100,150,250,400,675,1000,2000)	
a = seq(1.1, 2.5, by = 0.2)	
rpt = 500

numc = matrix(nrow = length(size), ncol = length(a), dimnames=list(size,a)) 
# numc stores how many times the observed distance is finite
Dc = array(dim = c(rpt, length(size), length(a)), dimnames = list(1:rpt,size,a)) 
# Dc is a 3D array for storing the observed typical distances
for(j in 1:length(a)){
	for(i in 1:length(size)){
	n = size[i]
	p = a[j]*log(n)/n
	d = PlayER(n, p, N = rpt)
	Dc[,i,j] = d
	numc[i,j] = length(subset(d, d<Inf))
	}
}
# Next we calculate the sample mean and s.d. for each pair of n and c = np
mconn = matrix(nrow=length(size),ncol=length(a))
for(j in 1:length(a)) mconn[,j] = apply(Dc[,,j], 2, function(dat) mean(dat[dat<Inf]))
dimnames(mconn) = list(size, a)
sconn = matrix(nrow=length(size),ncol=length(a))
for(j in 1:length(a)) sconn[,j] = apply(Dc[,,j], 2, function(dat) sd(dat[dat<Inf]))
dimnames(sconn) = list(size, a)

for(j in 1:length(a))
write.csv(Dc[,,j],paste("yourpath/conn",a[j],".csv",sep=""), row.names = FALSE)

#-----------------------------------------------------------------------------------
# 								  II. Histograms
#-----------------------------------------------------------------------------------
	
for(j in 1:length(a)){
	jpeg(file=paste("histc",a[j],".jpeg",sep=""), width=1200, height=800, pointsize = 20, quality=100)
	layout(matrix(1:9,nrow=3,byrow=T))
	for(i in 1:length(size)){
		d = Dc[,i,j][!is.infinite(Dc[,i,j])]
		hist(d, breaks=seq(min(d-0.5),max(d+0.5),by=1), col="peachpuff", 
		border="black", prob = TRUE, xlab = "typical distance", 
		main = paste("for n =",size[i]," & p =",a[j],"log(n)/n"))
		x = seq(min(d-0.5), max(d+0.5), length=1000)
		lines(x,dnorm(x,mean=mconn[i,j],sd=sconn[i,j]), lwd=2, col="dark blue")
		lines(x,dnorm(x,mean=log(size[i])/log(a[j]*log(size[i])),sd=sconn[i,j]), lwd=2, col="dark red")
	}
	dev.off()
}
# Standardized Histograms 
for(j in 1:length(a)){
	jpeg(file=paste("sdhistc",a[j],".jpeg",sep=""), width=1200, height=800, pointsize = 20, quality=100)
	layout(matrix(1:9,nrow=3,byrow=T))
	for(i in 1:length(size)){
		d = Dc[,i,j][!is.infinite(Dc[,i,j])]
		d = (d - mean(d))/sd(d)
		h = hist(d, breaks=round(log(size[i])), col="peachpuff", 
		border="black", prob = TRUE, xlab = "typical distance", 
		main = paste("for n =",size[i]," & p =",a[j],"log(n)/n"))
		xfit = seq(min(h$breaks),max(h$breaks),length=40)
		yfit = dnorm(xfit, mean=0, sd=1)
		lines(xfit, yfit, col="dark blue", lwd=2)
	}
	dev.off()
}

#-----------------------------------------------------------------------------------
# 								III. Growth Plots
#-----------------------------------------------------------------------------------

error = mconn
for(i in 1:length(size)) 
	for(j in 1:length(a)) 
		error[i,j] = mconn[i,j]*log(a[j]*log(size[i]))/log(size[i]) - 1
		
k = length(a)
#-------------- Plotting the o(1) term = (sample mean - theo. mean)*log(np)/log(n)
matplot(error, ylim=c(min(error),max(error)), type = "o", pch=20, lty=1, lwd=3, col = rainbow(k+1)[-5], xlab="Graph size n", main="The o(1) term in the theoretical mean of the typical distance in ER(n, c*log(n)/n)", ylab="sample mean*log(c*log n)/log n - 1", xaxt = "n")
axis(side=1,at=1:length(size),labels=size)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright",inset=c(-0.18,0), legend=paste("c = ",rev(a),sep=""),
		col=rev(rainbow(k+1)[-5]),lwd=3)

#-------------- Plotting the sample s.d. (unscaled)
matplot(sconn, ylim=c(min(sconn),max(sconn)), type = "o", pch=20, lty=1, lwd=3, col = (rainbow(k+1)[-5]), xlab="Graph size n", main="Sample s.d. for the typical distance in ER(n, c*log(n)/n)", ylab="Sample s.d.", xaxt = "n")
axis(side=1,at=1:length(size),labels=size)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright",inset=c(-0.18,0), legend=paste("c = ",a,sep=""),
		col=(rainbow(k+1)[-5]),lwd=3)

#-------------- Plotting the proportion of times the typical distance is finite
matplot(numc/rpt, type="o", pch=20, lty=1, lwd=3, col = rainbow(k+1)[-5], xlab="Graph size n", main="Proportion of times the typical distance is finite, for p = c*log(n)/n (where c > 1)", ylab="Proportion of times the typical distance is finite", xaxt = "n", yaxt = "n")
axis(side=1, at=1:length(size), labels=size)
axis(side=2, at=seq(min(numc/rpt), max(numc/rpt), by=0.01))
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright", inset=c(-0.18,0), legend=paste("c =",rev(a),sep=" "), col=rev(rainbow(k+1)[-5]), lty=1, lwd=3)

#-------------- Column mean and s.d. of the '#times d is finite' matrix
par(mfrow=t(1:2), mai = c(1, 0.5, 0.4, 0.4))
plot(apply(num, 2, mean)~a, main="column means", xlab="value of c=np", col="Dark blue", pch=20, cex=1.5, xaxt="n")
axis(side=1, at=a)
plot(apply(num, 2, sd)~a, main="column sd", xlab="value of c=np", col="Maroon", pch=20, cex=1.5, xaxt="n")
axis(side=1, at=a)

#-----------------------------------------------------------------------------------
#								IV. Testing normality
#-----------------------------------------------------------------------------------

testlist = c("Pearson chi-square", "Kolmogorov-Smirnov", "Shapiro-Wilk", "Cramer-von Mises", "Shapiro-Francia", "Anderson-Darling")
P = array(dim = c(6, length(size), length(a)), dimnames = list(testlist,size,a)) 
for(j in 1:length(a)){
	for(i in 1:length(size)){
	dat = Dc[, i, j][!is.infinite(Dc[, i, j])] 
	sdat = (dat - mean(dat))/sd(dat)
	P[ , i, j] = norm.testing(sdat)
	}
}	# P is a 3D array for storing the testing p-values
for(i in 1:6){
	cat("\n Name of the test : ", testlist[i],"\n\n")
	print(round(P[i, ,], 5))
	cat("\n")
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# 							   Constant probability
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-----------------------------------------------------------------------------------
# 							     I. Data Generation
#-----------------------------------------------------------------------------------

set.seed(125)
size = c(20,60,100,150,250,400,675,1000,2000)
cp = seq(0.1,0.9,by=0.2)
rpt = 200 		

numct = matrix(nrow = length(size), ncol = length(cp), dimnames=list(size,cp)) 
Dct = array(dim = c(rpt, length(size), length(cp)), dimnames = list(1:rpt,size,cp))
for(j in 1:length(cp)){
	for(i in 1:length(size)){
	n = size[i]
	d = PlayER(n, cp[j], N = rpt)
	Dct[,i,j] = d
	numct[i,j] = length(subset(d, d<Inf))
	}
}
# Next we calculate the sample mean and s.d. for each pair of n and p
k = length(cp)
mcnst = matrix(nrow=length(size), ncol=k)
for(j in 1:k) mcnst[,j] = apply(Dct[,,j], 2, function(dat) mean(dat[dat<Inf]))
dimnames(mcnst) = list(size, cp)
scnst = matrix(nrow=length(size), ncol=k)
for(j in 1:k) scnst[,j] = apply(Dct[,,j], 2, function(dat) sd(dat[dat<Inf]))
dimnames(scnst) = list(size, cp)

for(j in 1:length(cp))
write.csv(Dct[ , , j],paste("yourpath/const",cp[j],".csv",sep=""), row.names = FALSE)

#-----------------------------------------------------------------------------------
# 								 II. Histograms
#-----------------------------------------------------------------------------------
	
for(j in 1:k){
	jpeg(file=paste("histc",cp[j],".jpeg",sep=""), width=1200, height=800, pointsize = 20, quality=100)
	layout(matrix(1:9, nrow=3, byrow=T))
	for(i in 1:length(size)){
		d = Dct[,i,j][!is.infinite(Dct[,i,j])]
		hist(d, breaks=seq(min(d-0.5),max(d+0.5),by=1), col="peachpuff", 
		border="black", prob = TRUE, xlab = "typical distance", 
		main = paste("for n =",size[i]," & p =", cp[j]))
		x = seq(min(d-0.5), max(d+0.5), length=1000)
		lines(x,dnorm(x,mean=mcnst[i,j],sd=scnst[i,j]), lwd=2, col="dark blue")
		lines(x,dnorm(x,mean=log(n)/log(cp[j]*n),sd=scnst[i,j]), lwd=2, col="dark red")
	}
	dev.off()
}
# Standardized Histograms 
for(j in 1:k){
	jpeg(file=paste("sdhistc",cp[j],".jpeg",sep=""), width=1200, height=800, pointsize = 20, quality=100)
	layout(matrix(1:9, nrow=3, byrow=T))
	for(i in 1:length(size)){
		d = Dct[,i,j][!is.infinite(Dct[,i,j])]
		d = (d - mean(d))/sd(d)
		h = hist(d, breaks=round(log(size[i]))*2, col="peachpuff", 
		border="black", prob = TRUE, xlab = "typical distance", 
		main = paste("for n =",size[i]," & p =",cp[j]))
		xfit = seq(min(h$breaks),max(h$breaks),length=40)
		yfit = dnorm(xfit, mean=0, sd=1)
		lines(xfit, yfit, col="dark blue", lwd=2)
	}
	dev.off()
}

#-----------------------------------------------------------------------------------
#								  III. Growth Plots
#-----------------------------------------------------------------------------------

error = mcnst
for(i in 1:length(size)) 
	for(j in 1:length(cp)) 
		error[i,j] = mcnst[i,j]*log(cp[j]*size[i])/log(size[i]) - 1

cols = c("red","orange","green","blue","purple")
#-------------- Plotting the o(1) term = (sample mean - theo. mean)*log(np)/log(n)
matplot(error, ylim=c(min(error),max(error)), type = "o", pch=20, lty=1, lwd=3, col = cols, xlab="Graph size n", main="The o(1) term in the theoretical mean of typical distance in ER(n, p) (p constant)", ylab="sample mean*log(np)/log(n) - 1", xaxt = "n")
axis(side=1,at=1:length(size),labels=size)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright",inset=c(-0.18,0), legend=paste("p = ",cp,sep=""),
		col=cols,lwd=3)

#-------------- Plotting the sample s.d. (unscaled)
matplot(scnst, ylim=c(min(scnst),max(scnst)), type = "o", pch=20, lty=1, lwd=3, col = (cols), xlab="Graph size n", main="Sample s.d. for typical distance of ER(n, p) (p constant)", ylab="Sample s.d.", xaxt = "n")
axis(side=1,at=1:length(size),labels=size)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright",inset=c(-0.18,0), legend=paste("p = ",cp,sep=""),
		col=(cols),lwd=3)

#-----------------------------------------------------------------------------------
#								IV. Testing Normality
#-----------------------------------------------------------------------------------

testlist = c("Pearson chi-square", "Kolmogorov-Smirnov", "Shapiro-Wilk", "Cramer-von Mises", "Shapiro-Francia", "Anderson-Darling")
P = array(dim = c(6, length(size), length(cp)), dimnames = list(testlist,size,cp)) # 3D array for storing the testing p-values
for(j in 1:k){
	for(i in 1:length(size)){
	dat = Dct[,i,j][!is.infinite(Dct[,i,j])] 
	sdat = (dat - mean(dat))/sd(dat)
	P[,i,j] = norm.testing(sdat)
	}
}
for(i in 1:6){
	cat("\n Name of the test : ", testlist[i],"\n\n")
	print(round(P[i, ,], 6))
	cat("\n")
}

#-----------------------------------------------------------------------------------
#									Testing Symmetry
#-----------------------------------------------------------------------------------

f<-function(a,b,c) sign(a+b-2*c)+sign(b+c-2*a)+sign(c+a-2*b)
B1<-function(x,t){
	n=length(x)
	g=0
	if(t==1){
		for(j in (t+1):(n-1))
			for(k in (j+1):n)
				g=g+f(x[t],x[j],x[k])
	}
	else if(t==n){
		for(i in 1:(t-2))
			for(j in (i+1):(t-1))
				g=g+f(x[i],x[j],x[t])
	}
	else if(t==2){
		for(j in (t+1):(n-1))
			for(k in (j+1):n)
				g=g+f(x[t],x[j],x[k])
		for(i in 1:(t-1))
			for(k in (t+1):n)
				g=g+f(x[i],x[t],x[k])
	}
	else if(t==(n-1)){
		for(i in 1:(t-1))
			for(k in (t+1):n)
				g=g+f(x[i],x[t],x[k])

		for(i in 1:(t-2))
			for(j in (i+1):(t-1))
				g=g+f(x[i],x[j],x[t])
	}
	else{
		for(j in (t+1):(n-1))
			for(k in (j+1):n)
				g=g+f(x[t],x[j],x[k])
		for(i in 1:(t-1))
			for(k in (t+1):n)
				g=g+f(x[i],x[t],x[k])

		for(i in 1:(t-2))
			for(j in (i+1):(t-1))
				g=g+f(x[i],x[j],x[t])
	}
	return(g)
}
B2<-function(x,s,t){
	n=length(x)
	g=0
	if((t-s)==1){
		if(s==1){
			for(k in 3:n)
				g=g+f(x[s],x[t],x[k])
		}
		else if(t==n){
			for(i in 1:(n-2))
				g=g+f(x[i],x[s],x[t])
		}
		else{
			for(i in 1:(s-1))
				g=g+f(x[i],x[s],x[t])
			for(k in (t+1):n)
				g=g+f(x[s],x[t],x[k])
		}
	}
	else{
		if(s==1 & t<n){
			for(j in (s+1):(t-1))
				g=g+f(x[s],x[j],x[t])
			for(k in (t+1):n)
				g=g+f(x[s],x[t],x[k])
		}
		else if(s>1 & t==n){
			for(i in 1:(s-1))
				g=g+f(x[i],x[s],x[t])
			for(j in (s+1):(t-1))
				g=g+f(x[s],x[j],x[t])
		}
		else if(s==1 & t==n){
			for(j in (s+1):(t-1))
				g=g+f(x[s],x[j],x[t])
		}
		else{
			for(i in 1:(s-1))
				g=g+f(x[i],x[s],x[t])
			for(j in (s+1):(t-1))
				g=g+f(x[s],x[j],x[t])
			for(k in (t+1):n)
				g=g+f(x[s],x[t],x[k])
		}
	}
	return(g)
}
test_sym_large<-function(x){
	n=length(x)
	b1=0
	b=0
	for(t in 1:n){
		k=B1(x,t)
		b1=b1+k^2
		b=b+k
	}
	T=b/3
	b2=0
	for(s in 1:(n-1))
		for(t in (s+1):n)
			b2=b2+(B2(x,s,t))^2
	s_sq=((n-3)*(n-4)*b1)/((n-1)*(n-2))+((n-3)*b2)/(n-4)+(n*(n-1)*(n-2))/6-(1-((n-3)*(n-4)*(n-5))/(n*(n-1)*(n-2)))*(T^2)
	Z=T/sqrt(s_sq)
	p_val=2*(1-pnorm(abs(Z),0,1))
	ret.list=list(Z,p_val)
	names(ret.list)=c("obs.stat","p.val")
	return(ret.list)
}
P1 = array(dim = c(length(size), length(a)), dimnames = list(size, a)) 
for(j in 1:length(a)){
	for(i in 1:length(size)){
	dat = D[,i,j][!is.infinite(D[,i,j])] 
	sdat = (dat - mean(dat))/sd(dat)
	P1[i,j] = test_sym_large(sdat)$p.val
	}
}
P2 = array(dim = c(length(size), length(a)), dimnames = list(size, a)) 
for(j in 1:length(a)){
	for(i in 1:length(size)){
	dat = Dc[,i,j][!is.infinite(Dc[,i,j])] 
	sdat = (dat - mean(dat))/sd(dat)
	P2[i,j] = test_sym_large(sdat)$p.val
	}
}
P3 = array(dim = c(length(size), length(cp)), dimnames = list(size,cp)) 
for(j in 1:length(cp)){
	for(i in 1:length(size)){
	dat = Dct[,i,j][!is.infinite(Dct[,i,j])] 
	sdat = (dat - mean(dat))/sd(dat)
	P3[i,j] = test_sym_large(sdat)$p.val
	}
}
# P-values of the test of Symmetry
print(round(P1, 4)) # for the p = c/n case
print(round(P2, 4)) # for the p = c*log(n)/n case
print(round(P3, 4)) # for the P = constant case 
