
.onAttach <- 
		function(lib, pkg, ...)
{
	pkgDescription <- packageDescription(pkg)
	pkgVersion <- pkgDescription$Version
	pkgDate <- pkgDescription$Date
	pkgName <- pkgDescription$Package
	pkgTitle <- pkgDescription$Title
	pkgAuthor <- pkgDescription$Author
	pkgMaintainer <- pkgDescription$Maintainer
	packageStartupMessage(paste("\n", pkgName, ": ", pkgTitle, sep = ""))
	packageStartupMessage(paste("Version ", pkgVersion, " (", pkgDate, ") installed", sep = ""))
	packageStartupMessage(paste("Author: ", pkgAuthor, sep = ""))
	packageStartupMessage(paste("Maintainer: ", pkgMaintainer, "\n", sep = ""))
	packageStartupMessage('Use citation("bounds") to know how to cite this work.\n')
}

bounds.arvid <- function(formula, X, Z, data, monotonicity = FALSE, weights, R, cc){
	if(missing(weights)){
		weights=NULL
	}
	est=bounds.inner(formula,X,Z,data,monotonicity,weights,cc)
	out=matrix(NA,nrow=4,ncol=5) 
	rownames(out)=c("lower","upper","p.lower","p.upper") 
	colnames(out)=c("Estimate","Standard error","2.5th percentile ","97.5th percentile", "p.value")
	out[1:length(est),1]=est
	if(!missing(R)){
		est.boot=matrix(rep(NA,length(est)*R),nrow=length(est),ncol=R)
		n=nrow(data)
		r=1
		while(r<=R){
			data.r=data[sample(1:n,n,replace=TRUE),]
			temp=bounds.inner(formula,X,Z,data.r,monotonicity,weights,cc)
			if(is.numeric(temp)){
				est.boot[,r]=temp
				r=r+1
			}
		}
		out[1:length(est),2]=apply(est.boot,1,sd)
		out[1:length(est),3]=apply(est.boot,1,quantile,probs=0.025)
		out[1:length(est),4]=apply(est.boot,1,quantile,probs=0.975)
		out[1:length(est),5] <- pchisq((out[1:length(est),1]/out[1:length(est),2])**2, 1, lower.tail = FALSE)
	}
	out  
}

bounds.inner <- function(formula,X,Z,data,monotonicity,weights,cc) {
	
	#check that all levels of X and X are represented in the data set,
	#which is not necessarily the case when you bootstrap. If not, then
	#redo the bootstrap sample
	if(any(is.na(match(levels(data$X),unique(data$X)))) |
			any(is.na(match(levels(data$Z),unique(data$Z)))))
		return("redo")
	
	#check that both exposures are factors 
	if(!is.factor(data[,X]) | !is.factor(data[,Z])){
		return("Error: both exposures must be factors")  
	}
	#check that the outcome is numeric 0/1
	Y=all.vars(formula)[attr(terms(formula),"response")]
	if(!is.numeric(data[,Y])){
		return("Error: the outcome must be numeric")
	} 
	else{
		if(sum(data[,Y]!=0 & data[,Y]!=1,na.rm=TRUE)>0){
			return("Error: the outcome must be coded as 0/1")
		}
	}
	
	#fit logistic regression model for the outcome
	if(is.null(weights)){
		fit=glm(formula=formula,family=binomial,data=data)
	}
	else{
		w=data[,weights]
		fit=glm(formula=formula,family=quasibinomial,data=data,weights=w) 
	}
	
	#use the fitted model to obtain estimates of p.xz and put these into matrix P
	levX = levels(data[,X]) 
	KX=length(levX)
	levZ = levels(data[,Z]) 
	KZ=length(levZ)
	P=matrix(rep(NA,KX*KZ),nrow=KX,ncol=KZ)
	for(i in 1:KX){
		for(j in 1:KZ){
			newdata=data
			newdata[,X]=factor(levX[i])
			newdata[,Z]=factor(levZ[j])
			if (!cc) {
				P[i,j]=mean(predict.glm(object=fit,newdata=newdata,type="respons"),
						na.rm=TRUE)
			} else {
				P[i,j]=mean(predict.glm(object=fit,newdata=newdata[newdata[,Y] == 0,],type="respons"),
					    na.rm=TRUE)  
			}
		}
	}
	
	#calculate bounds
	lower=0
	upper=1
	
	if(monotonicity){
		#lower bound
		for(x in 1:KX){
			xf=1
			while(xf<=x){
				xff=1
				while(xff<=xf){
					for(z in 1:KZ){
						zf=1
						while(zf<=z){
							zff=1
							while(zff<=zf){
								pxz=P[x,z]
								pxzf=P[x,zf]
								pxzff=P[x,zff]
								pxfz=P[xf,z]
								pxfzf=P[xf,zf]
								pxfzff=P[xf,zff]
								pxffz=P[xff,z]
								pxffzf=P[xff,zf]
								pxffzff=P[xff,zff]
								psi.00=pxz-pxfz-pxzf+pxfzf
								psi.01=pxzf-pxfzf-pxzff+pxfzff
								psi.10=pxfz-pxffz-pxfzf+pxffzf
								psi.11=pxfzf-pxffzf-pxfzff+pxffzff
								lower=max(lower,
										psi.00,
										-psi.00,
										(psi.00+psi.01),
										-(psi.00+psi.01),
										(psi.00+psi.10),
										-(psi.00+psi.10),
										psi.00-psi.01,
										psi.00-psi.10,
										(psi.00+psi.11),
										-(psi.00+psi.11),
										(psi.00-psi.11),
										-(psi.00-psi.11),
										psi.10-psi.01-psi.11,
										psi.01-psi.10-psi.11,
										(psi.00+psi.10+psi.11),
										-(psi.00+psi.10+psi.11),
										(psi.00+psi.01+psi.11),
										-(psi.00+psi.01+psi.11),
										psi.00+psi.10-psi.01,
										psi.00-psi.10+psi.01,
										psi.00+psi.10+psi.01-psi.11,
										psi.00+psi.10-psi.01-psi.11,
										psi.00-psi.10+psi.01-psi.11,
										psi.00-psi.10-psi.01-psi.11)
								zff=zff+1
							}
							zf=zf+1
						}
					}
					xff=xff+1
				}
				xf=xf+1
			}
		}
		#upper bound
		eX=P[KX,KZ]-P[1,1]
		power.set.X=powerSet(1:(KX-1))
		#note: need not consider the first element in power.set.X,
		#since this is the empty set, which is already taken care of 
		#by setting eX=P[KX,KZ]-P[1,1]
		for(i in 2:length(power.set.X)){
			wX=power.set.X[[i]]
			eX=min(eX,P[KX,KZ]-P[1,1]+sum(P[wX,KZ]-P[wX+1,1]))
		}
		eZ=P[KX,KZ]-P[1,1]
		power.set.Z=powerSet(1:(KZ-1))
		#note: need not consider the first element in power.set.Z,
		#since this is the empty set, which is already taken care of 
		#by setting eZ=P[KX,KZ]-P[1,1]
		for(i in 2:length(power.set.Z)){
			wZ=power.set.Z[[i]]
			eZ=min(eZ,P[KX,KZ]-P[1,1]+sum(P[KX,wZ]-P[1,wZ+1]))
		}
		upper=min(P[KX,KZ]-P[1,1]+(P[2,1]-P[1,2]),
				P[KX,KZ]-P[1,1]-(P[2,1]-P[1,2]),
				eX,
				eZ)
		#proportion of total effect due to causal interaction 
		p.lower=lower/(P[KX,KZ]-P[1,1])
		p.upper=upper/(P[KX,KZ]-P[1,1])
		out=c(lower,upper,p.lower,p.upper)
	}
	else{
		#lower bound
		for(x in 1:KX){
			xf=1
			while(xf<=x){
				for(z in 1:KZ){
					zf=1
					while(zf<=z){
						pxz=P[x,z]
						pxzf=P[x,zf]
						pxfz=P[xf,z]
						pxfzf=P[xf,zf]
						psi=pxz-pxfz-pxzf+pxfzf
						lower=max(lower,
								0.5*psi,
								-0.5*psi,
								psi-pxfzf,
								-psi-pxzf,
								-psi-pxfz,
								psi-pxz,
								-psi+pxfzf-1,
								psi+pxzf-1,
								psi+pxfz-1,
								-psi+pxz-1)
						zf=zf+1
					}
				}
				xf=xf+1
			}
		}
		#upper bound
		eX=1
		power.set.X=powerSet(1:KX)
		for(i in 1:length(power.set.X)){
			wX=power.set.X[[i]]
			v=rep(0,KX)
			v[wX]=1
			I=matrix(v,nrow=KX,ncol=KZ)
			Ic=1-I
			eX=min(eX,sum(I*(1-P)+Ic*P))
		}
		eZ=1
		power.set.Z=powerSet(1:KZ)
		for(i in 1:length(power.set.Z)){
			wZ=power.set.Z[[i]]
			v=rep(0,KZ)
			v[wZ]=1
			I=matrix(v,nrow=KX,ncol=KZ,byrow=TRUE)
			Ic=1-I
			eZ=min(eZ,sum(I*(1-P)+Ic*P))
		}
		upper=min(eX,eZ)
		out=c(lower,upper)
	}
	out   
}


bounds.suff.arvid=function(formula,X,x,Z,z,data,weak=FALSE,monotonicity=FALSE,weights,R,cc){
	if(missing(weights)){
		weights=NULL
	}
	est=bounds.suff.inner(formula,X,x,Z,z,data,weak,monotonicity,weights,cc)
	out=matrix(NA,nrow=4,ncol=5) 
	rownames(out)=c("lower","upper","p.lower","p.upper") 
	colnames(out)=c("Estimate","Standard error","2.5th percentile ","97.5th percentile","p.value")
	out[1:length(est),1]=est
	r=0
	if(!missing(R)){
		est.boot=matrix(rep(NA,4*R),nrow=4,ncol=R)
		n=nrow(data)
		while(r<=R){
			data.r=data[sample(1:n,n,replace=TRUE),]
			temp=bounds.suff.inner(formula,X,x,Z,z,data.r,weak,monotonicity,weights,cc)
			if(is.numeric(temp)){
				est.boot[,r]=temp
				r=r+1
			}
		}
		out[1:length(est),2]=apply(est.boot,1,sd)
		out[1:length(est),3]=apply(est.boot,1,quantile,probs=0.025)
		out[1:length(est),4]=apply(est.boot,1,quantile,probs=0.975)
		out[1:length(est),5] <- pchisq((out[1:length(est),1]/out[1:length(est),2])**2, 1, lower.tail = FALSE)
		#bias=matrix(rowMeans(est.boot)-est,nrow=4,ncol=R)
		#out[1:length(est),3]=apply(est.boot-bias,1,quantile,probs=0.025)
		#out[1:length(est),4]=apply(est.boot-bias,1,quantile,probs=0.975)
	}
	out  
}

bounds.suff.inner=function(formula,X,x,Z,z,data,weak,monotonicity,weights,cc){
	
	#check that all levels of X and X are represented in the data set,
	#which is not necessarily the case when you bootstrap. If not, then
	#redo the bootstrap sample
	if(any(is.na(match(levels(data$X),unique(data$X)))) |
			any(is.na(match(levels(data$Z),unique(data$Z)))))
		return("redo")
	
	#check that both exposures are factors 
	if(!is.factor(data[,X]) | !is.factor(data[,Z])){
		return("Error: both exposures must be factors")  
	}
	#check that x and z are levels of X and Z, respectively
	if(is.na(match(x,levels(data[,X]))) | is.na(match(z,levels(data[,Z])))){
		return("Error: x and z must be levels of X and Z, respectively")  
	}
	#check that the outcome is numeric 0/1
	Y=all.vars(formula)[attr(terms(formula),"response")]
	if(!is.numeric(data[,Y])){
		return("Error: the outcome must be numeric")
	} 
	else{
		if(sum(data[,Y]!=0 & data[,Y]!=1,na.rm=TRUE)>0){
			return("Error: the outcome must be coded as 0/1")
		}
	}
	
	#fit logistic regression model for the outcome
	if(is.null(weights)){
		fit=glm(formula=formula,family=binomial,data=data)
	}
	else{
		w=data[,weights]
		fit=glm(formula=formula,family=quasibinomial,data=data,weights=w)
	}
	
	#use the fitted model to obtain estimates of p.xz and put these into matrix P
	levX = levels(data[,X]) 
	KX=length(levX)
	levZ = levels(data[,Z]) 
	KZ=length(levZ)
	P=matrix(rep(NA,KX*KZ),nrow=KX,ncol=KZ)
	for(i in 1:KX){
		for(j in 1:KZ){
			newdata=data
			newdata[,X]=factor(levX[i])
			newdata[,Z]=factor(levZ[j])
			if (!cc) {
				P[i,j]=mean(predict.glm(object=fit,newdata=newdata,type="respons"),
						na.rm=TRUE)
			} else {
				P[i,j]=mean(predict.glm(object=fit,newdata=newdata[newdata[,Y] == 0,],type="respons"),
						na.rm=TRUE)  
			} 
		}
	}
	
	x=match(x,levels(data[,X]))
	z=match(z,levels(data[,Z]))
	
	#calculate bounds
	lower=0
	upper=1
	
	if(monotonicity){
		if(weak){
			#lower bound
			xf=1
			while(xf<=x){
				xff=1
				while(xff<=xf){
					for(z in 1:KZ){
						zf=1
						while(zf<=z){
							zff=1
							while(zff<=zf){
								pxz=P[x,z]
								pxzf=P[x,zf]
								pxzff=P[x,zff]
								pxfz=P[xf,z]
								pxfzf=P[xf,zf]
								pxfzff=P[xf,zff]
								pxffz=P[xff,z]
								pxffzf=P[xff,zf]
								pxffzff=P[xff,zff]
								psi.a=pxz-pxfz-pxzf+pxfzf
								psi.b=pxz-pxffz-pxzff+pxffzff
								psi.c=pxfzf-pxffzf-pxfzff+pxffzff
								lower=max(lower,psi.a,psi.b-psi.c)
								zff=zff+1
							}
							zf=zf+1
						}
					}
					xff=xff+1
				}
				xf=xf+1
			}
			#upper bound
			upper=min(upper,P[x,z]-P[1,z],P[x,z]-P[x,1])
		}
		else{
			if(x!=KX | z!=KZ){
				upper=0
			}
			else{
				#lower
				lower=max(lower,P[KX,KZ]-P[KX-1,KZ]-P[KX,KZ-1]+P[KX-1,KZ-1])
				#upper
				upper=min(upper,P[KX,KZ]-P[KX-1,KZ],P[KX,KZ]-P[KX,KZ-1])
			}
		}
		#proportion of total effect due to sufficient-cause interaction 
		#p.lower=lower/(P[KX,KZ]-P[1,1])
		#p.upper=upper/(P[KX,KZ]-P[1,1])
		#out=c(lower,upper,p.lower,p.upper)
	}
	else{
		if(weak){
			#lower bound
			for(xf in (1:KX)[-x]){
				for(zf in (1:KZ)[-z]){
					lower=max(lower,P[x,z]-P[xf,z]-P[x,zf])
				}
			}
			#upper bound
			upper=min(upper,P[x,z],sum(1-P[(1:KX)[-x],z]),sum(1-P[x,(1:KZ)[-z]]))
		}
		else{
			#lower bound
			lower=max(lower,P[x,z]-sum(P[(1:KX)[-x],z])-sum(P[x,(1:KZ)[-z]]))
			#upper bound
			upper=min(upper,P[x,z])
			for(xf in (1:KX)[-x]){
				for(zf in (1:KZ)[-z]){
					upper=min(upper,1-P[xf,z],1-P[x,zf])
				}
			}
		} 
		#out=c(lower,upper)
	} 
	#proportion of P_{xz} due to sufficient-cause interaction 
	p.lower=lower/P[x,z]
	p.upper=upper/P[x,z]
	
	out=c(lower,upper,p.lower,p.upper)
	out
} 