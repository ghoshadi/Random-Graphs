library(igraph)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#                Lattice Random Graph
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-----------------------------------------------------------------------------------
#                              I. Data Generation
#-----------------------------------------------------------------------------------

vert<-function(i,j,n)
{
	return((i-1)*(2*n+1)+j)
}

Lattice<-function(n,p)
{
	g=make_empty_graph(directed=F)+vertices(1:((2*n+1)^2),color="black")
	for(i in 1:(2*n))
	{
		for(j in 1:(2*n))
		{
			if(rbinom(1,1,p)==1)(g=g+edges(vert(i,j,n),vert(i+1,j,n)))
			if(rbinom(1,1,p)==1)(g=g+edges(vert(i,j,n),vert(i,j+1,n)))
		}
		if(rbinom(1,1,p)==1)(g=g+edges(vert(i,2*n+1,n),vert(i+1,2*n+1,n)))
	}
	for(j in 1:(2*n))
	{
		if(rbinom(1,1,p)==1)(g=g+edges(vert(2*n+1,j,n),vert(2*n+1,j+1,n)))
	}
	plot(g,layout=layout_on_grid(g),vertex.size=2,vertex.col="black")
	return(g)
}

PlayLattice<-function(n,p,N=1000){
	d=rep(0,N)
	for(i in 1:N){
		g=Lattice(n, p)
		s=sample(1:((2*n+1)^2),2,replace=T)
		d[i]=distances(g,s[1],s[2])
	}
	return(d)
}

# PlayLattice generates realizations of H_n for Lat(n, p) for N times and returns the vector of the observed typical distances

n=seq(5,40,by=5)
p=seq(0.55,0.95,by=0.05)
N=1000
H = array(0,c(length(p),length(n),N))
t1<-proc.time()
for(i in 1:length(p))
{
	for(j in 1:length(n))
	{
		H[i,j,]=PlayLattice(n[j],p[i],1000)
		write.table("your.path/n_jp_i",header=T)
	}
}
t2<-proc.time()

# Generating and storing the data in a folder

#-----------------------------------------------------------------------------------
#                               II. Histograms
#-----------------------------------------------------------------------------------

for(j in 1:length(p)
{
	layout(matrix(1:8,nrow=2,byrow=T))
	for(i in 1:length(n))
	{
		a=unlist(read.table("your.path/n_jp_i",header=T))
		d=subset(a,a<Inf)
		hist(d,breaks=10,col="peachpuff",border="black",prob=TRUE,xlab="Typical distance",main = paste("For n =",n[i],"& p =",p[j]))
		x=seq(min(d-0.5),max(d+0.5),length=1000)
		lines(x,dnorm(x,mean=Mean_mat[i,j],sd=Sd_mat[i,j]),lwd=2,col="dark blue")
	}
}

#-----------------------------------------------------------------------------------
#                               III. Growth Plots
#-----------------------------------------------------------------------------------

Mean_mat=matrix(nrow=length(n),ncol=length(p))
Sd_mat=matrix(nrow=length(n),ncol=length(p))
Prop_mat=matrix(nrow=length(n),ncol=length(p))
dimnames(Mean_mat)=list(n,p)
dimnames(Sd_mat)=list(n,p)
dimnames(Prop_mat)=list(n,p)

for(i in 1:length(n))
{
	for(j in 1:length(p))
	{
		a=unlist(read.table("your.path/n_jp_i",header=T))
		b=subset(a,a<Inf)
		Mean_mat[i,j]=mean(b)
		Sd_mat[i,j]=sd(b)
		Prop_mat[i,j]=length(b)/length(a)
	}
}

k=length(p)

#-------------- Plotting the sample mean when the typical distance is finite

matplot(Mean_mat, ylim=c(min(Mean_mat),max(Mean_mat)),type="o",pch=20,lty=1,lwd=1.9,col=rainbow(k+1)[-5],xlab="Graph dimension (n)",main="The mean typical distance in lattice when its finite",ylab="Mean typical distance",xaxt="n")
axis(side=1,at=1:length(n),labels=n)
par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
legend("topright",inset=c(-0.11,0), legend=paste("p=",p,sep=""),col=rainbow(k+1)[-5],lwd=3)

#-------------- Plotting the sample s.d. when the typical distance is finite

matplot(Sd_mat, ylim=c(min(Sd_mat),max(Sd_mat)),type="o",pch=20,lty=1,lwd=1.9,col=rainbow(k+1)[-5],xlab="Graph dimension (n)",main="The sd of typical distance in lattice when its finite",ylab="Sd of typical distance",xaxt="n")
axis(side=1,at=1:length(n),labels=n)
par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
legend("topright",inset=c(-0.11,0), legend=paste("p=",p,sep=""),col=rainbow(k+1)[-5],lwd=3)

#-------------- Plotting the proportion of times the typical distance is finite

matplot(Prop_mat, ylim=c(min(Prop_mat),max(Prop_mat)),type="o",pch=20,lty=1,lwd=1.9,col=rainbow(k+1)[-5],xlab="Graph dimension (n)",main="The proportion of times when typical distance in lattice is finite",ylab="Proportion when typical distance is finite",xaxt="n")
axis(side=1,at=1:length(n),labels=n)
par(mar=c(5.1,4.1,4.1,8.1),xpd=T)
legend("topright",inset=c(-0.11,0), legend=paste("p=",rev(p),sep=""),col=rev(rainbow(k+1)[-5]),lwd=3)

#-----------------------------------------------------------------------------------
#                              IV. Testing Normality
#-----------------------------------------------------------------------------------

for(j in 1:length(p)
{
	layout(matrix(1:8,nrow=2,byrow=T))
	for(i in 1:length(n))
	{
		a=unlist(read.table("your.path/n_jp_i",header=T))
		b=subset(a,a<Inf)
		qqnorm((b-mean(b))/sd(b))
		qqline((b-mean(b))/sd(b))
	}
}

# QQ plots

mychisq.test<-function(dat)
{
	sdat = (dat - mean(dat))/sd(dat)
	f.os = table(cut(sdat, breaks=c(-Inf,-2:2,Inf)))
	f.ex = (pnorm(c(-2:2,Inf)) - pnorm(c(-Inf,-2:2)))*length(sdat)
	chi.stat = sum((f.os-f.ex)^2/f.ex)
	return(pchisq(chi.stat, df=length(f.os)-1, lower.tail=FALSE))
}

myks.test<-function(dat)
{
	sdat = (dat - mean(dat))/sd(dat) 
	jdat = sdat + rnorm(length(sdat), mean=0, sd=0.001) # jittered to remove ties
	y = rnorm(length(sdat))
	pval = ks.test(jdat, y)$p.value
	return(pval)
}

norm.testing<-function(sdat)
{
	sw = shapiro.test(sdat)$p.value
	ks = myks.test(sdat)
	psn = mychisq.test(sdat)
	return(c(psn,ks,sw))
}

testlist = c("Pearson chi-square", "Kolmogorov-Smirnov", "Shapiro-Wilk")
P = array(dim = c(3, length(n), length(p)), dimnames = list(testlist,n,p))

for(j in 1:length(p))
{
	for(i in 1:length(n))
	{
		a=unlist(read.table("your.path/n_jp_i",header=T))
		dat=subset(a,a<Inf)
		sdat = (dat - mean(dat))/sd(dat)
		P[,i,j] = norm.testing(sdat)
	}
}

for(i in 1:3)
{
	cat("\n Name of the test : ", testlist[i],"\n\n")
	print(round(P[i, ,], 5))
	cat("\n")
}

#-----------------------------------------------------------------------------------
#                                 Testing Symmetry
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

P1 = array(dim = c(length(n), length(p)), dimnames = list(n,p)) 
for(j in 1:length(p))
{
	for(i in 1:length(n))
	{
		a=unlist(read.table("your.path/n_jp_i",header=T))
		dat=subset(a,a<Inf)
		sdat = (dat - mean(dat))/sd(dat)
		P1[i,j] = test_sym_large(sdat)$p.val
	}
}

# P-values of the test of Symmetry
print(round(P1, 5))