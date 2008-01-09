####
# PROBLEMS:
# Consider optionally making use of the coda or mcmc package facilities
# Change output variable names to something better
# Needs input checking

library(survival)
library(MASS)

##############################################################
# Generate simulated data
##############################################################

generaterandom<-function(n,type,params)
{
    if(!(type%in%c("fixed","weibull","gamma","normal","lognormal","normmix","lognormmix"))) stop("Invalid distribution type")
    if(type=="fixed"){
        if(!("value"%in%names(params))) stop("Parameter value not specified for type fixed")
        return(rep(params$value,n))
    }
    if(type=="weibull"){
        if(!all(c("lambda0","gamweib")%in%names(params))) stop("Parameters lambda0, gamweib not specified for type weibull")
        lambda0<-params$lambda0
        gamweib<-params$gamweib
        return(rweibull(n,shape=gamweib,scale=lambda0^(-1/gamweib)))
    }
    if(type=="gamma"){
        if(!all(c("mu","sigma2")%in%names(params))) stop("Parameters mu, sigma2 not specified for type gamma")
        mu<-params$mu
        sigma2<-params$sigma2
        return(rgamma(n,shape=mu^2/sigma2,scale=sigma2/mu))
    }
    if(type=="normal"){
        if(!all(c("mu","sigma2")%in%names(params))) stop("Parameters mu, sigma2 not specified for type normal")
        mu<-params$mu
        sigma2<-params$sigma2
        return(rnorm(n,mean=mu,sd=sqrt(sigma2)))
    }
    if(type=="lognormal"){
        if(!all(c("mu","sigma2")%in%names(params))) stop("Parameters mu, sigma2 not specified for type lognormal")
        mu<-params$mu
        sigma2<-params$sigma2
        sigma2prime<-log(1+sigma2/mu^2)        
        muprime<-log(mu)-1/2*sigma2prime
        return(rlnorm(n,meanlog=muprime,sdlog=sqrt(sigma2prime)))
    }
    if(type=="normmix"){
        if(!all(c("mu","sigma2","w")%in%names(params))) stop("Parameters mu, sigma2, w not specified for type normmix")
        mu<-params$mu
        sigma2<-params$sigma2
        w<-params$w/sum(params$w)
        if(length(w)==1) w=rep(w,length(mu))
        if(length(mu)!=length(sigma2) | length(mu)!=length(w) | length(mu)<2) stop("Bad parameter lengths for type normmix")
        out<-mvrnorm(n,mu,mdiag(sigma2)) 
        return(t(out)[findInterval(runif(n),cumsum(w))+1+0:(n-1)*length(w)])
    }
    if(type=="lognormmix"){
        if(!all(c("mu","sigma2","w")%in%names(params))) stop("Parameters mu, sigma2, w not specified for type lognormmix")
        mu<-params$mu
        sigma2<-params$sigma2
        w<-params$w/sum(params$w)
        if(length(w)==1) w=rep(w,length(mu))
        if(length(mu)!=length(sigma2) | length(mu)!=length(w) | length(mu)<2) stop("Bad parameter lengths for type lognormmix")
        sigma2prime<-log(1+sigma2/mu^2)        
        muprime<-log(mu)-1/2*sigma2prime
        out<-mvrnorm(n,muprime,mdiag(sigma2prime)) 
        return(exp(t(out)[findInterval(runif(n),cumsum(w))+1+0:(n-1)*length(w)]))   
    }
}

generateevents<-function(m,Ji,beta,Ui,Zij,type,params)
{
    Tij<-rep(0,sum(Ji))
    Uij<-rep(Ui,Ji)
    if(type=="weibull"){
        if(!all(c("lambda0","gamweib")%in%names(params))) stop("Parameters lambda0, gamweib, w not specified for type weibull")
        lambda0<-params$lambda0
        gamweib<-params$gamweib
        for(ind in 1:sum(Ji)){
            Tij[ind]<-((-exp(-beta*Zij[ind])*log(1-runif(1))/(lambda0*Uij[ind]))^(1/gamweib))
        }
    }
    if(type=="stepfunction"){
        breaks<-params$breaks
        haz<-params$haz
        if(length(haz)!=length(breaks)+1) stop("Step function params: haz should be one longer than breaks")
        for(ind in 1:sum(Ji)){
            accept<-FALSE
            maxhaz<-max(haz)*Uij[ind]*exp(beta*Zij[ind])
            while(!accept){
                Tijprop<- -1/maxhaz*log(runif(1))
                thishaz<-haz[findInterval(Tijprop,breaks)+1]*Uij[ind]*exp(beta*Zij[ind])
                if(maxhaz*runif(1)<thishaz){
                    Tij[ind]<-Tijprop
                    accept<-TRUE
                }
            }
        }
    }
    if(type=="bspline"){
        b<-params$b
        w<-params$w
        rbound<-attr(b,"Boundary.knots")[2]
        survfn<-bs.survfn(b,w,1000,rbound)
        if(survfn>.25) warning(paste("Baseline survival over B-spline support is high:",format(survfn,digits=3)))
        for(ind in 1:sum(Ji)){
            accept<-FALSE
            Tijprop<-rbound+1
            while(!accept){
                maxhaz<-max(w)*Uij[ind]*exp(beta*Zij[ind])
                while(Tijprop>rbound) Tijprop<- -1/maxhaz*log(runif(1))
                thishaz<-predict(b,Tijprop)%*%w*Uij[ind]*exp(beta*Zij[ind])
                if(maxhaz*runif(1)<thishaz){
                    Tij[ind]<-Tijprop
                    accept<-TRUE
                }
            }
        }
    }
    return(Tij)
}

bs.survfn<-function(b,w,n,t=NULL)
{
    rbound<-attr(b,"Boundary.knots")[2]
    x<-seq(from=0,to=rbound,length=n)
    h<-(predict(b,x)%*%w)/n
    if(is.null(t)){
        H<-cumsum(h)
        S<-exp(-H)
        return(S)
    }else{
        Ht<-sum(h[x<t])
        St<-exp(-Ht)
        return(St)
    }
}

makeagdata<-function(m,Ji,Tij,deltaij,Zij)
{
    agdata<-data.frame(
        i=rep(1:m,Ji),
        j=c(sapply(Ji,function(x) return(1:x)),recursive=TRUE),
        time=Tij,
        delta=deltaij
    )
    agdata<-cbind(agdata,Zij)
    return(agdata)  
}

sim.sample<-function(m=10,Ji=rep(5,10),params=NULL)
{
    if(length(Ji)==1) Ji<-rep(Ji,m)
    if(length(Ji)!=m) stop("Length of Ji does not match m")
    params.in<-params
    params.default<-list(
        beta=1,
        haz.type="weibull",
        haz.params=list(lambda0=1,gamweib=1.8),
        frail.type="lognormal",
        frail.params=list(mu=1,sigma2=.25),
        Z.type="normal",
        Z.params=list(mu=0,sigma2=1),
        C.type="weibull",
        C.params=list(lambda0=1,gamweib=1.8)
    )
    params<-params.default
    if(!is.null(params.in)){
        for(n in names(params.in)) eval(parse(text=paste("params$",n,"<-params.in$",n,sep="")))
        if(params.in$haz.type=="bspline" & is.null(params.in$haz.params)){
            tmax=3;
            N=4
            b<-bs(0,knots=seq(from=0,to=tmax,length=N+1)[-(N+1)],Boundary.knots=c(0,2),degree=3)
            w<-c(.3,.2,.4,.2,.6,.8,1)*4
            params$haz.params<-list(b=b,w=w)
        }
    }
    #attach(params)
    Zij<-generaterandom(sum(Ji),params$Z.type,params$Z.params)
    Cij<-generaterandom(sum(Ji),params$C.type,params$C.params)
    Ui<-generaterandom(m,params$frail.type,params$frail.params)
    Tij<-generateevents(m,Ji,params$beta,Ui,Zij,params$haz.type,params$haz.params)
    deltaij<-as.numeric(Tij<Cij)
    Tij[deltaij==0]<-Cij[deltaij==0]
    #detach(params)
    
    agdata<-makeagdata(m,Ji,Tij,deltaij,data.frame(Z=Zij))
    return(list(agdata=agdata,Ui=Ui,params=params))
}

##############################################################
# B-spline integrals and other matrices
##############################################################

ki<-function(knots,j) return(match(j,attr(knots,"i")))

evalBinte<-function(knots,ord)
{
    K<-sum(knots>attr(knots,"b")[1] & knots<attr(knots,"b")[2])
    Binte<-rep(0,K+ord)
    for(j in 1:(K+ord)){
        Binte[j]<-(knots[ki(knots,j)]-knots[ki(knots,j-ord)])/ord
    }
    return(Binte)
}

evalCinte<-function(knots,ord,obs,Binte)
{
    K<-sum(knots>attr(knots,"b")[1] & knots<attr(knots,"b")[2])
    Cinte<-matrix(0,length(obs),K+ord)
    knots2<-c(attr(knots,"b")[1],knots,attr(knots,"b")[2])
    attr(knots2,"i")<-c(min(attr(knots,"i"))-1,attr(knots,"i"),max(attr(knots,"i"))+1)
    attr(knots2,"b")<-attr(knots,"b")
    Bordp1<-splineDesign(knots2, x=obs, outer.ok = TRUE, ord=ord+1)
    for(i in 1:length(obs)){
        for(j in 1:(K+ord)){
            if(obs[i]>=knots[ki(knots,j)]) Cinte[i,j]<-Binte[j]
            if(obs[i]<knots[ki(knots,j)] & obs[i]>=knots[ki(knots,j-ord)]) Cinte[i,j]<-Binte[j]*sum(Bordp1[i,(j+1):(K+ord+1)])
        }
    }
    return(Cinte)
}

evalEinte<-function(knots,ord)
{
    K<-sum(knots>attr(knots,"b")[1] & knots<attr(knots,"b")[2])
    Einte<-rep(0,K+ord)
    for(j in 1:(K+ord)){
        Einte[j]<-nBsmom(1,ord,j,knots)
    }
    return(1-Einte)
}

nBsmom<-function(N,ord,j,knots)
{
    lknot<-knots[ki(knots,j-ord)]
    rknot<-knots[ki(knots,j)]
    if(lknot==rknot) return(rknot^N)
    return(ord/(rknot-lknot)/(N+1)*(-nBsmom(N+1,ord-1,j-1,knots)+nBsmom(N+1,ord-1,j,knots)))
}

splineconv<-function(k,n1,j1,n2,j2,knots)
{
    if(j1-n1>=j2 | j2-n2>=j1) return(0)
    if(n1==1 & n2==1){
        out<-1/(k+1)*(knots[ki(knots,j1)]^(k+1)-knots[ki(knots,j1-1)]^(k+1))
        return(out)
    }
    if(n2>n1){
        n3<-n1; n1<-n2; n2<-n3
        j3<-j1; j1<-j2; j2<-j3
    }
    out<-0
    denom1<-knots[ki(knots,j1-1)]-knots[ki(knots,j1-n1)]
    denom2<-knots[ki(knots,j1)]-knots[ki(knots,j1-n1+1)]
    if(denom1>0){
        out<-out+1/denom1*splineconv(k+1,n1-1,j1-1,n2,j2,knots)
        out<-out-knots[ki(knots,j1-n1)]/denom1*splineconv(k,n1-1,j1-1,n2,j2,knots)
    }
    if(denom2>0){
        out<-out+knots[ki(knots,j1)]/denom2*splineconv(k,n1-1,j1,n2,j2,knots)
        out<-out-1/denom2*splineconv(k+1,n1-1,j1,n2,j2,knots)
    }
    return(out)
}

splinederivint<-function(l1,n1,j1,l2,n2,j2,knots)
{
    if(j1-n1>=j2 | j2-n2>=j1) return(0)
    if(l1==0 & l2==0) return(splineconv(0,n2,j1,n2,j2,knots))
    if(l2>l1){
        l3<-l1; l1<-l2; l2<-l3
        n3<-n1; n1<-n2; n2<-n3
        j3<-j1; j1<-j2; j2<-j3
    }
    out<-0
    denom1<-knots[ki(knots,j1-1)]-knots[ki(knots,j1-n1)]
    denom2<-knots[ki(knots,j1)]-knots[ki(knots,j1-n1+1)]
    if(denom1>0) out<-out+(n1-1)/denom1*splinederivint(l1-1,n1-1,j1-1,l2,n2,j2,knots)
    if(denom2>0) out<-out-(n1-1)/denom2*splinederivint(l1-1,n1-1,j1,l2,n2,j2,knots)
    return(out)
}


mysplineDesign<-function (knots, x, ord = 4) 
{
    nk<-length(knots)
    x <- as.numeric(x)
    derivs<-integer(1)
    temp <- .Call("spline_basis", knots, ord, x, derivs, PACKAGE = "splines")
    ncoef <- nk - ord
    design <- rep(0,ncoef)
    jj <- 1:ord + attr(temp, "Offsets")
    design[jj] <- temp
    dim(design)<-c(1,ncoef)
    design
}


makePenalty.2diff<-function(K){
    D<-matrix(0,K-2,K)
    for(i in 1:(K-2)){
        D[i,i]<-1
        D[i,i+1]<--2
        D[i,i+2]<-1
    }
    P<-t(D)%*%D
    return(P)
}

makePenalty.2deriv<-function(ord,knots)
{
    nspline<-sum(knots>attr(knots,"b")[1])
    out<-matrix(0,nspline,nspline)
    for(j1 in 1:nspline){
        for(j2 in j1:nspline){
            out[j1,j2]<-splinederivint(2,ord,j1,2,ord,j2,knots)
        }
    }
    for(j1 in 2:nspline){
        for(j2 in 1:(j1-1)) out[j1,j2]<-out[j2,j1]
    }
    return(out)
}


plotspline<-function(knots,theta,ord,npoints=1000,plotknots=TRUE,plotmean=FALSE,plotsplines=FALSE,norm=FALSE,xlim=NULL,ylim=NULL,...)
{
    if(is.null(xlim)) xlim=attr(knots,"b")
    x=seq(from=xlim[1],to=xlim[2],length=npoints)
    dx=diff(xlim)/npoints
    spl<-splineDesign(knots, x=x, ord=ord)
    if(norm){
        Bint<-evalBinte(knots,ord)
        for(i in 1:dim(spl)[1]) spl[i,]<-spl[i,]/Bint
    }
    splsum<-spl%*%exp(theta) 
    if(norm) splsum<-splsum/sum(exp(theta))
    #if(plotsplines) ymax<-max(max(spl),max(splsum)) else 
    ymax<-max(splsum)
    if(is.null(ylim)) ylim<-c(0,ymax)
    if(plotsplines){
        matplot(x,spl/max(spl)*ylim[2]/2,type="l",lty=2,xlim=xlim,ylim=ylim,...)
        lines(x,splsum,col="red",lwd=2,lty=1)
    }else{
        plot(x,splsum,col="red",type="l",lwd=2,lty=1,xlim=xlim,ylim=ylim,...)
    }
    if(plotknots) abline(v=knots,col="grey")
    if(plotmean){
        Bint<-colSums(spl)/dx
        spl.norm<-spl%*%mdiag(1/Bint)
        Eint<-rep(0,dim(spl)[2])
        for(i in 1:dim(spl)[2]) Eint[i]<-sum(spl.norm[,i]*x)/dx
        E<-Eint%*%exp(theta)/sum(exp(theta))
        abline(v=E,col="grey",lwd=2)
    }
}

plotspline.meta<-function(fit,which="base",...)
{
    plotbase<-FALSE;plotfrail<-FALSE;plotdens<-FALSE
    if(which=="base") plotbase<-TRUE
    if(which=="frail") plotfrail<-TRUE
    if(which=="coef") plotdens<-TRUE
    if(which=="all") {
        par(mfrow=c(3,1))
        plotbase<-TRUE;plotfrail<-TRUE;plotdens<-TRUE
    }
    if(plotbase){
        knots<-fit$spline$knots.haz
        ord<-fit$spline$ord.haz
        theta<-fit$posterior.mean$splinepar.haz
        plotspline(knots,theta,ord,...)
    }
    if(plotfrail){
        knots<-fit$spline$knots.frail
        ord<-fit$spline$ord.frail
        theta<-fit$posterior.mean$splinepar.frail
        plotspline(knots,theta,ord,norm=TRUE,...)
    }
    if(plotdens){
        burnin<-fit$control$burnin
        if(is.null(burnin)) burnin<-fit$call$control$burnin
        if(is.null(burnin)) burnin<-dim(fit$history$coefficients)[2]
        for(i in 1:dim(fit$history$coefficients)[2]){
            betai<-fit$history$coefficients[burnin:(dim(fit$history$coefficients)[1]),i]
            plot(density(betai),...)
        }
    }
}

accrate.predict.s<-function(x,s) predict(s,x)$y-.25
accrate.predict.s2<-function(x,s) predict(s,x)$y-.25
accrate.predict.lm<-function(x,m) (predict(m,data.frame(x=x))-.25)^2
accrate.predict.lm2<-function(x,m) (predict(m,data.frame(x=x))-.25)^2

##############################################################
# Likelihoods, gradients and hessians
##############################################################

mdiag<-function(x) if(is.vector(x) && length(x)==1 && x<1) return(matrix(x,1,1)) else return(diag(x))

################# For frailties

mklik.Ui<-function(u,whichi,i.all,delta,beta,Z,C.l,theta.l,B.u,theta.u)
{
    iind<-i.all==whichi
    delta<-delta[iind]
    Z<-Z[iind,]; dim(Z)<-c(sum(iind),length(beta))
    C.l<-C.l[iind,]; dim(C.l)<-c(sum(iind),length(theta.l))
    B.u<-B.u[whichi,]; dim(B.u)<-c(1,length(theta.u))
    #Z<-matrix(Z[iind,],sum(iind),length(beta))
    #C.l<-matrix(C.l[iind,],sum(iind),length(theta.l))
    #B.u<-matrix(B.u[whichi,],1,length(theta.u))
    
    lik.Ui<-0
    lik.Ui<-lik.Ui+sum(delta*log(u))
    lik.Ui<-lik.Ui-sum(u*(C.l%*%exp(theta.l))*exp(Z%*%beta))
    lik.Ui<-lik.Ui+log(B.u%*%exp(theta.u))
    
    lik.Ui<-as.numeric(lik.Ui)
    return(lik.Ui)
}

################# For beta

mklik.b<-function(beta,delta,Z,U.u,C.l,theta.l,sigma2.b)
{
    lik.b<-0
    lik.b<-lik.b+t(delta)%*%Z%*%beta
    lik.b<-lik.b-U.u%*%((C.l%*%exp(theta.l))*exp(Z%*%beta))
    lik.b<-lik.b-t(beta)%*%beta/(2*sigma2.b)
    lik.b<-as.numeric(lik.b)
    return(lik.b)
}

mkgrad.b<-function(beta,delta,Z,U.u,C.l,theta.l,sigma2.b)
{
    grad.b<-rep(0,length(beta))
    grad.b<-grad.b+t(Z)%*%delta
    grad.b<-grad.b-t(Z)%*%((C.l%*%exp(theta.l))*exp(Z%*%beta)*U.u)
    grad.b<-grad.b-beta/sigma2.b
    grad.b<-as.numeric(grad.b)
    return(grad.b)
}

mkhess.b<-function(beta,delta,Z,U.u,C.l,theta.l,sigma2.b)
{
    hess.b<-matrix(0,length(beta),length(beta))
    hess.b<-hess.b-t(Z)%*%mdiag(as.vector((C.l%*%exp(theta.l))*(exp(Z%*%beta))*U.u))%*%Z
    hess.b<-hess.b-mdiag(length(beta))/sigma2.b
    return(hess.b)
}


################# For theta_lambda

mklik.l<-function(theta.l,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty)
{   
    lik.l<-0
    lik.l<-lik.l+t(log(B.l%*%exp(theta.l)))%*%delta
    lik.l<-lik.l-t(U.u*C.l%*%exp(theta.l))%*%exp(Z%*%beta)
    lik.l<-lik.l-smoothpen.l(theta.l,P.l,penalty,0)/(2*sigma2.l)
    lik.l<-as.numeric(lik.l)
    return(lik.l)
}

mkgrad.l<-function(theta.l,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty)
{
    grad.l<-rep(0,length(theta.l))
    grad.l<-grad.l+exp(theta.l)*t(B.l)%*%(1/(B.l%*%exp(theta.l))*delta)
    grad.l<-grad.l-exp(theta.l)*t(C.l)%*%(exp(Z%*%beta)*U.u)
    grad.l<-grad.l-smoothpen.l(theta.l,P.l,penalty,1)/(2*sigma2.l)
    grad.l<-as.numeric(grad.l)
    return(grad.l)
}

mkhess.l<-function(theta.l,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty)
{
    hess.l<-matrix(0,length(theta.l),length(theta.l))
    hess.l<-hess.l+mdiag(exp(theta.l)*as.vector(t(B.l)%*%(1/(B.l%*%exp(theta.l))*delta)))
    hess.l<-hess.l-mdiag(exp(theta.l))%*%t(B.l)%*%mdiag(delta/as.vector(B.l%*%exp(theta.l))^2)%*%B.l%*%mdiag(exp(theta.l))
    hess.l<-hess.l-mdiag(exp(theta.l)*as.vector(t(C.l)%*%(exp(Z%*%beta)*U.u)))
    hess.l<-hess.l-smoothpen.l(theta.l,P.l,penalty,2)/(2*sigma2.l)
    return(hess.l)
}

################# For theta_u

mklik.u<-function(theta.u,B.u,E.u,sigma2.u,P.u,M,penalty,fix)
{
    theta.u<-theta.u.repair(theta.u,fix)
    m<-dim(B.u)[1]
    lik.u<-0
    lik.u<-lik.u+sum(log(B.u%*%exp(theta.u)))
    lik.u<-lik.u-m*log(sum(exp(theta.u)))
    lik.u<-lik.u-smoothpen.u(theta.u,P.u,penalty,0)/(2*sigma2.u)
    lik.u<-lik.u-M*( (t(E.u)%*%exp(theta.u))/sum(exp(theta.u)) )^2
    lik.u<-as.numeric(lik.u)
    return(lik.u)
}

mkgrad.u<-function(theta.u,B.u,E.u,sigma2.u,P.u,M,penalty,fix)
{
    theta.u<-theta.u.repair(theta.u,fix)
    m<-dim(B.u)[1]
    grad.u<-rep(0,length(theta.u))
    grad.u<-grad.u+exp(theta.u)*(t(B.u)%*%(1/B.u%*%exp(theta.u)))
    grad.u<-grad.u-m/sum(exp(theta.u))*exp(theta.u)
    grad.u<-grad.u-smoothpen.u(theta.u,P.u,penalty,1)/(2*sigma2.u)
    grad.u<-grad.u-2*M*t(E.u)%*%exp(theta.u)/sum(exp(theta.u))^2*E.u*exp(theta.u)
    grad.u<-grad.u+2*M*(t(E.u)%*%exp(theta.u))^2/sum(exp(theta.u))^3*exp(theta.u)
    grad.u<-as.numeric(grad.u)
    return(grad.u[-fix])
}

mkhess.u<-function(theta.u,B.u,E.u,sigma2.u,P.u,M,penalty,fix)
{
    theta.u<-theta.u.repair(theta.u,fix)
    m<-dim(B.u)[1]
    hess.u<-matrix(0,length(theta.u),length(theta.u))
    hess.u<-hess.u+mdiag(exp(theta.u))%*%mdiag(as.vector(t(B.u)%*%mdiag(as.vector(1/B.u%*%exp(theta.u)))%*%rep(1,m)))
    hess.u<-hess.u-mdiag(exp(theta.u))%*%t(B.u)%*%mdiag(1/as.vector(B.u%*%exp(theta.u))^2)%*%B.u%*%mdiag(exp(theta.u))
    hess.u<-hess.u-m/sum(exp(theta.u))^2*(sum(exp(theta.u))*mdiag(exp(theta.u))-exp(theta.u)%*%t(exp(theta.u)))
    hess.u<-hess.u-smoothpen.u(theta.u,P.u,penalty,2)/(2*sigma2.u)
    hess.u<-hess.u-2*M/sum(exp(theta.u))^3 * ( sum(exp(theta.u))*((exp(theta.u)*E.u)%*%t((exp(theta.u)*E.u))+as.numeric(t(E.u)%*%exp(theta.u))*mdiag(E.u*exp(theta.u))) 
                                                -2*as.numeric(t(E.u)%*%exp(theta.u))*mdiag(E.u)%*%exp(theta.u)%*%t(exp(theta.u)) )
    hess.u<-hess.u+2*M/sum(exp(theta.u))^4 * ( sum(exp(theta.u))*( 2*as.numeric(t(E.u)%*%exp(theta.u))*exp(theta.u)%*%t(exp(theta.u))%*%mdiag(E.u)+as.numeric(t(E.u)%*%exp(theta.u))^2*mdiag(exp(theta.u)))
                                                -3*as.numeric(t(E.u)%*%exp(theta.u))^2*exp(theta.u)%*%t(exp(theta.u)) )
    hess.u<-hess.u[-fix,-fix]
    return(hess.u)
}

theta.u.repair<-function(theta.u,fix)
{
    if(fix==1) return(c(1,theta.u))
    if(fix==length(theta.u)) return(c(theta.u,1))
    return(c(theta.u[1:(fix-1)],1,theta.u[fix:length(theta.u)]))
}

################# Smoothness penalty functions and gradients/hessians

smoothpen.l<-function(theta.l,P.l,penalty,der)
{
    if(penalty=="2diff"){
        if(der==0) pen<-t(theta.l)%*%P.l%*%theta.l
        if(der==1) pen<-2*P.l%*%theta.l
        if(der==2) pen<-2*P.l
    }
    if(penalty=="2deriv"){
        if(der==0) pen<-t(exp(theta.l))%*%P.l%*%exp(theta.l)
        if(der==1) pen<-2*mdiag(exp(theta.l))%*%P.l%*%exp(theta.l)
        if(der==2) {
            pen<-2*mdiag(exp(theta.l))%*%P.l%*%mdiag(exp(theta.l))
            pen<-pen+2*mdiag(as.vector(P.l%*%exp(theta.l)))%*%mdiag(exp(theta.l))
        }
    }
    return(pen)       
}

smoothpen.u<-function(theta.u,P.u,penalty,der)
{
    if(penalty=="2diff"){
        if(der==0) pen<-t(theta.u)%*%P.u%*%theta.u
        if(der==1) pen<-2*P.u%*%theta.u
        if(der==2) pen<-2*P.u
    }
    if(penalty=="2deriv"){
        st<-sum(exp(theta.u))
        tpt<-as.numeric(t(exp(theta.u))%*%P.u%*%exp(theta.u))
        if(der==0) pen<-tpt/st^2
        if(der==1) {
            dpt<-mdiag(exp(theta.u))%*%P.u%*%exp(theta.u)
            pen<-2/st^2 * ( dpt - tpt/st*exp(theta.u) )
        }
        if(der==2) { 
            tt<-exp(theta.u)%*%t(exp(theta.u)) 
            dp<-mdiag(exp(theta.u))%*%P.u             
            dptt<-dp%*%tt
            dpd<-dp%*%mdiag(exp(theta.u))
            pen<-2/st^2 * (dpd + mdiag(as.vector(P.u%*%exp(theta.u)))%*%mdiag(exp(theta.u)) )
            pen<-pen - 2*tpt/st^3 *(mdiag(exp(theta.u))-3/st*tt)
            pen<-pen - 4/st^3 *(dptt+t(dptt))
        }
    }
    return(pen)
}

##############################################################
# Metropolis-Hastings routines
##############################################################

mh.U<-function(Ui,i.all,delta,beta,Z,C.l,theta.l,B.u,theta.u,nu2.U,gamma.U,knots.u,ord.u,Bint.u)
{
    acc<-NULL
    for(i in 1:length(Ui)){
        u<-Ui[i]
        baselik<-mklik.Ui(u,i,i.all,delta,beta,Z,C.l,theta.l,B.u,theta.u)
        #cand<-abs(rnorm(1,u,gamma.U*nu2.U))
        v<-gamma.U*nu2.U
        cand<-rgamma(1,shape=u^2/v,scale=v/u)
        B.c<-B.u
        if(is.nan(cand) || cand>attr(knots.u,"b")[2] | cand<attr(knots.u,"b")[1]) {
            candlik<- -Inf
        }else{
            B.c[i,]<-mysplineDesign(knots.u,x=cand,ord=ord.u)/Bint.u
            candlik<-mklik.Ui(cand,i,i.all,delta,beta,Z,C.l,theta.l,B.c,theta.u)
        }
        if(is.nan(candlik)) candlik<--Inf
        puc<-suppressWarnings(dgamma(cand,shape=u^2/v,scale=v/u)) # transition u->cand
        pcu<-suppressWarnings(dgamma(u,shape=cand^2/v,scale=v/cand)) # transition cand->u
        r<-exp(candlik-baselik)*pcu/puc
        p.accept<-min(1,r)
        if(is.nan(p.accept)) p.accept<-0
        if(runif(1)<p.accept){
            Ui[i]<-cand
            acc<-c(acc,1)
        }else{
            acc<-c(acc,0)
        }
        if(Ui[i]>attr(knots.u,"b")[2] | Ui[i]<attr(knots.u,"b")[1]) browser()
    }
    return(list(Ui=Ui,acc=mean(acc)))
}

mh.b<-function(beta,delta,Z,U.u,C.l,theta.l,sigma2.b,Sigma.b,gamma.b)
{
    baselik<-mklik.b(beta,delta,Z,U.u,C.l,theta.l,sigma2.b)
    cand<-mvrnorm(1,beta,Sigma.b*gamma.b)
    candlik<-mklik.b(cand,delta,Z,U.u,C.l,theta.l,sigma2.b)
    if(is.nan(candlik)) candlik<--Inf
    r<-exp(candlik-baselik)
    p.accept<-min(1,r)
    if(runif(1)<p.accept){
        beta<-cand
        acc<-1
    }else{
        acc<-0
    }
    return(list(beta=beta,acc=acc))
}

mh.l<-function(theta.l,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty,Sigma.l,gamma.l)
{
    baselik<-mklik.l(theta.l,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty)
    cand<-mvrnorm(1,theta.l,Sigma.l*gamma.l)
    candlik<-mklik.l(cand,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty)
    if(is.nan(candlik)) candlik<--Inf
    r<-exp(candlik-baselik)
    p.accept<-min(1,r)
    if(runif(1)<p.accept){
        theta.l<-cand
        acc<-1
    }else{
        acc<-0
    }
    return(list(theta.l=theta.l,acc=acc))
}

mh.u<-function(theta.u,B.u,E.u,sigma2.u,P.u,M,penalty,fix,Sigma.u,gamma.u)
{
    K.u<-length(theta.u)
    baselik<-mklik.u(theta.u[-fix],B.u,E.u,sigma2.u,P.u,M,penalty,fix)
    cand<-mvrnorm(1,theta.u[-fix],Sigma.u*gamma.u)
    candlik<-mklik.u(cand,B.u,E.u,sigma2.u,P.u,M,penalty,fix)
    if(is.nan(candlik)) candlik<--Inf
    r<-exp(candlik-baselik)
    p.accept<-min(1,r)
    if(runif(1)<p.accept){
        theta.u<-theta.u.repair(cand,fix)
        acc<-1
    }else{
        acc<-0
    }
    return(list(theta.u=theta.u,acc=acc))
}

##############################################################
# S3 Methods for fitting, printing, summary, etc
##############################################################

################# Methods for fitting

splinesurv<-function(x,...)
{
    UseMethod("splinesurv")
}

splinesurv.data.frame<-function(x,...)
{
    call<-match.call()
    if(dim(x)[2]>4 && all(colnames(x)[1:4]==c("i","j","time","delta"))) {
        x<-x[order(x$i),]
        out<-splinesurv.agdata(x,...)
        out$call<-call
        return(out)
    } else {
        stop("input data frame needs to have colnames i,j,time,delta")
    }
}

splinesurv.formula<-function(formula,data=parent.frame(),...)
{
    # in part based on coxph function
    call<-match.call()
    m<-match.call(expand.dots=FALSE)
    if(is.matrix(eval(m$data,sys.parent()))) m$data<-as.data.frame(data)
    m$...<-NULL
    m[[1]]<-as.name("model.frame")
    special<-"cluster"
    Terms <- if (missing(data)) terms(formula, special) else terms(formula, special, data = data)    
    m$formula<-Terms
    m<-eval(m,sys.parent())
    n<-nrow(m)
    resp <- model.extract(m, "response")
    if (!is.Surv(resp)) stop("model response must be a Surv object")
    if(attr(resp,"type")!="right") stop("right-censored survival data only")
    time<-resp[,"time"]
    delta<-resp[,"status"]
    clusterind<-attr(Terms,"specials")$cluster
    dropx<-NULL
    clusternames<-NULL
    if(length(clusterind)>0){
        cluster<-m[,clusterind]
        if(is.factor(cluster)) clusternames<-levels(cluster)
        if(is.numeric(cluster)) clusternames<-as.character(unique(cluster))
        i<-as.numeric(as.factor(cluster))
        tempc <- untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1)) stop("Cluster can not be used in an interaction")
        dropx <- c(dropx,tempc$terms)
    }else{
        i<-rep(1,n)
    }
    Ji<-table(i)
    j<-as.vector(sapply(Ji,function(x) 1:x))
    newTerms <- if(length(dropx))  Terms[-dropx] else Terms
    X <- model.matrix(newTerms, m)
    X<-X[,-1,drop=FALSE]
    agdata<-as.data.frame(cbind(i,j,time,delta,X))
    agdata[,-2]<-agdata[order(agdata$i),-2]
    class(agdata)<-c("agdata","data.frame")
    fit<-splinesurv.agdata(agdata,...)
    fit$call<-call
    colnames(fit$history$frailty)<-clusternames
    if(!is.null(fit$posterior.mean)) names(fit$posterior.mean$frailty)<-clusternames
    return(fit)
}

################# Methods for printing

print.splinesurv<-function(x,...)
{
    cat("Regression parameter posterior means:\n")
    coef<-as.matrix(x$posterior.mean$coef,ncol=1)
    colnames(coef)<-"coef"
    print(coef,...)
    invisible(x)
}

summary.splinesurv<-function(object,quantiles=c(.025,.975),...)
{
    x<-object
    out<-NULL
    out$call<-x$call
    out$coef<-as.matrix(x$posterior.mean$coefficients,ncol=1)
    colnames(out$coef)<-"mean"
    out$iter<-x$control$iter
    out$burnin<-x$control$burnin
    out$nknots.haz<-length(x$spline$knots.haz)
    out$nknots.frail<-length(x$spline$knots.frail)
    out$ord.haz<-x$spline$ord.haz
    out$ord.frail<-x$spline$ord.frail    
    out$knotspacing.haz<-x$spline$knotspacing.haz
    out$knotspacing.frail<-x$spline$knotspacing.frail
    out$knotrange.haz<-attr(x$spline$knots.haz,"boundary")
    out$knotrange.frail<-attr(x$spline$knots.frail,"boundary")
    out$quantiles<-NULL
    goodcoef<-x$history$coefficients[(out$burnin+2):(out$iter+1),,drop=FALSE]
    if(length(quantiles)){
        for(q in quantiles){
            out$quantiles<-cbind(out$quantiles,apply(goodcoef,2,function(x) quantile(x,q)))
        }
        colnames(out$quantiles)<-paste(quantiles*100,"%",sep="")
    }
    out$dots<-as.list(substitute(list(...)))[-1]
    class(out)<-"summary.splinesurv"   
    return(out)
}

print.summary.splinesurv<-function(x,...)   
{ 
    cat("Call: \n")
    print(x$call)
    cat("\nIterations: ",x$iter," (",x$burnin," discarded as burn-in)\n",sep="")
    cat("\nRegression parameter posterior:\n")
    printpars<-paste(names(x$dots),unlist(x$dots),sep="=",collapse=",")
    if(nchar(printpars)) printpars<-paste(",",printpars)
    eval(parse(text=paste("print(cbind(x$coef,x$quantiles,...)",printpars,")")))
    cat("\nBaseline hazard spline: Order ",x$ord.haz," B-spline\n\twith ", x$nknots.haz-2*x$ord.haz+2," interior knots distributed ", if(x$knotspacing.haz=="equal") "evenly" else "by quantiles",
        " on (",paste(format(x$knotrange.haz,digits=2),collapse=","),")",sep="")
    cat("\nFrailty density spline: Order ",x$ord.frail," B-spline\n\twith ", x$nknots.frail-2*x$ord.frail+2," interior knots distributed ", if(x$knotspacing.frail=="equal") "evenly" else "by quantiles",
        " on (",paste(format(x$knotrange.frail,digits=2),collapse=","),")\n",sep="")
    
    invisible(x)
}

################# Methods for plotting

plot.splinesurv<-function(x,which=c("base","frail","coef","all"),...)
{
    if(length(which)>=1 && which[1]%in%c("base","frail","coef","all")){
        plotspline.meta(x,which[1],...)
    }
}


##############################################################
# MAIN FUNCTION
##############################################################


splinesurv.agdata<-function(x,verbose=3,initial=NULL,control=NULL,spline=NULL,coda=FALSE,...)
{
        
    if(verbose>=1) cat("Initializing...\n")
    
    agdata<-x
    
    call<-match.call()
    m<-length(unique(agdata$i))
    Ji<-table(agdata$i)
    
    if(verbose>=2) cat("\tSetting initial parameters...\n")
    
    # Parse input (control)
    control.in<-control
    control.default<-list(
        burnin=500, # Length of the burn-in period
        maxiter=1000, # Max number of iterations
        cal.gamma=TRUE, # Auto-calibrate tuning parameters
        calint=1e2, # Interval for calibration of the acceptance rate
        alpha=rep(.01,6), # Hyperparameters
        M=1e3,  # Penalty parameter for density mean
        gamma=rep(1,4) # Tuning parameters
    )
    control<-control.default
    controlnames<-c("burnin","burnin.max","burnin.min","terminate","maxiter","cal.gamma","calint","alpha","M","gamma")
    innames<-names(control.in)
    validnames<-intersect(innames,controlnames)
    if(!is.null(control.in)){
        for(n in validnames) eval(parse(text=paste("control$",n,"<-control.in$",n,sep="")))
    }
    burnin<-control$burnin; maxiter<-control$maxiter; cal.gamma<-control$cal.gamma
    calint<-control$calint; alpha<-control$alpha; M<-control$M;gamma<-control$gamma
        
    # Parse input (spline)
    spline.in<-spline
    spline.default<-list(
        ord.haz=4, # Order of the baseline hazard spline
        ord.frail=4, # Order of the frailty density spline
        knotspacing.haz="quantile", # Spacing of the knots (if automatically chosen)
        knotspacing.frail="equal", # Spacing of the knots
        nknots.haz=NULL, # Number of knots (hazard) 
        nknots.frail=NULL, # Number of knots (frailty)
        knots.haz=NULL, # knot positions (hazard)
        knots.frail=NULL, # knot positions (frailty)
        penalty.haz="2deriv",
        penalty.frail="2diff"
    )
    spline<-spline.default
    if(!is.null(spline.in)){
        for(n in names(spline.in)) eval(parse(text=paste("spline$",n,"<-spline.in$",n,sep="")))
    }
    #for(n in c("ord.haz","ord.frail","knotspacing.haz","knotspacing.frail","nknots.haz","nknots.frail","knots.haz","knots.frail")) eval(parse(text=paste(n,"<-spline$",n,sep="")))
    ord.l<-spline$ord.haz; knotspacing.l<-spline$knotspacing.haz; nknots.l<-spline$nknots.haz; knots.l<-spline$knots.haz; penalty.l<-spline$penalty.haz
    ord.u<-spline$ord.frail; knotspacing.u<-spline$knotspacing.frail; nknots.u<-spline$nknots.frail; knots.u<-spline$knots.frail; penalty.u<-spline$penalty.frail
        
    ### Hyperparameters
    if(length(alpha)!=6) stop("alpha needs to have length 6")
    alpha.b<-alpha[1:2]
    alpha.l<-alpha[3:4]
    alpha.u<-alpha[5:6]
    
    # sigma^2 initial values
    sigma2.b<-.1
    sigma2.l<-.1
    sigma2.u<-.1
    
    # Tuning parameters default values
    gamma.U<-gamma[1]
    gamma.b<-gamma[2]
    gamma.l<-gamma[3]
    gamma.u<-gamma[4]
    
    # Automatic number of knots
    if(is.null(nknots.l)) nknots.l<-min(round(sum(Ji)/4),35)
    if(is.null(nknots.u)) nknots.u<-min(round(m/4),35)
    
    if(verbose>=2) cat("\tFitting Cox survival models...\n")

    
    # Cox fit with gamma frailties for initial values of Ui and beta
    varnames<-colnames(agdata)[-(1:4)]
    qvarnames<-paste("`",varnames,"`",sep="")
    if(m>1){
        coxfit<-coxph(as.formula(paste("Surv(time,delta)~",paste(qvarnames,collapse="+"),"+frailty(i)")),data=agdata)
        Ui<-exp(coxfit$frail)
        if(var(Ui)<1e-5) {
            Ui<-2*runif(length(Ui))
            Ui<-1+(Ui-mean(Ui))
            Ui[Ui<0]<-1
        }
    }else{
        coxfit<-coxph(as.formula(paste("Surv(time,delta)~",paste(qvarnames,collapse="+"))),data=agdata)
        Ui<-1
    }        
    beta<-coxfit$coef

    if(verbose>=2) cat("\tInitializing spline bases...\n")

    # Vector lenghts
    K.u<-nknots.u+ord.u
    K.l<-nknots.l+ord.l
    p<-length(beta)
    
    # Knots for the splines
    bounds.l<-c(min(agdata$time),max(agdata$time))
    if(is.null(knots.l)){
        if(knotspacing.l=="quantile"){
            ibounds.l<-c(min(agdata$time[agdata$delta==1]),max(agdata$time[agdata$delta==1]))
            nintknots<-nknots.l; lrep<-ord.l; rrep<-ord.l
            if(ibounds.l[1]==bounds.l[1]) {nintknots<-nintknots+1; lrep<-ord.l-1}
            if(ibounds.l[2]==bounds.l[2]) {nintknots<-nintknots+1; rrep<-ord.l-1}
            knots.l<-quantile(agdata$time[agdata$delta==1],seq(from=0,to=1,length=nintknots))
            knots.l<-c(rep(bounds.l[1],lrep),knots.l,rep(bounds.l[2],rrep))
        }
        if(knotspacing.l=="equal"){
            knots.l<-seq(from=bounds.l[1],to=bounds.l[2],length=nknots.l+2)
            knots.l<-c(rep(bounds.l[1],ord.l-1),knots.l,rep(bounds.l[2],ord.l-1))
        }
    }
    attr(knots.l,"boundary")<-bounds.l
    attr(knots.l,"index")<-seq(from=-(ord.l-1),length=length(knots.l),by=1)

    #browser()

    bounds.u<-c(0,2*max(Ui))
    if(bounds.u[1]<0) bounds.u[1]<-0
    if(is.null(knots.u)){
        if(knotspacing.u=="quantile"){
            stop("Quantile knot spacing not allowed for frailty density")
        }
        if(knotspacing.u=="equal"){    
            knots.u<-seq(from=bounds.u[1],to=bounds.u[2],length=nknots.u+2)
            knots.u<-c(rep(bounds.u[1],ord.u-1),knots.u,rep(bounds.u[2],ord.u-1))
        }
    }
    attr(knots.u,"boundary")<-bounds.u
    attr(knots.u,"index")<-seq(from=-(ord.u-1),length=length(knots.u),by=1)
    
    # Evaluate the splines and integrals
    B.u<-splineDesign(knots.u,x=Ui,ord=ord.u)
    Bint.u<-evalBinte(knots.u,ord.u)
    for(i in 1:dim(B.u)[1]) B.u[i,]<-B.u[i,]/Bint.u
    E.u<-evalEinte(knots.u,ord.u)
    
    B.l<-splineDesign(knots.l,x=agdata$time,ord=ord.l)
    Bint.l<-evalBinte(knots.l,ord.l)
    C.l<-evalCinte(knots.l,ord.l,agdata$time,Bint.l)
    
    if(verbose>=2) cat("\tInitializing penalty matrices...\n")

    # Penalty matrices
    if(penalty.l=="2diff") P.l<-makePenalty.2diff(K.l)
    if(penalty.u=="2diff") P.u<-makePenalty.2diff(K.u)
    if(penalty.l=="2deriv") P.l<-makePenalty.2deriv(ord.l,knots.l)
    if(penalty.u=="2deriv") P.u<-makePenalty.2deriv(ord.u,knots.u)/(Bint.u%*%t(Bint.u))        

    # Other matrices and vectors
    delta<-agdata$delta
    Z<-agdata[,-(1:4)]; Z<-as.matrix(Z,sum(Ji),p)
    U.u<-rep(Ui,Ji)
    i.all<-agdata$i
    
    if(verbose>=2) cat("\tObtaining initial values for spline parameters...\n")

    # Initial values for the theta vectors
    theta.l<-rep(1,K.l)
    theta.u<-rep(1,K.u)
    #fixthetau<-round(K.u/2)
    fixthetau<-max(1,which.min((knots.u-mean(Ui))^2)-ord.u)
    #if(penalty.l=="2deriv") sigma2.l<-1/rgamma(1,K.l/2+alpha.l[1],1/(as.numeric(smoothpen.l(theta.l,P.l,penalty.l,0))/2+alpha.l[2]))
    #if(penalty.u=="2deriv") sigma2.u<-1/rgamma(1,K.u/2+alpha.u[1],1/(as.numeric(smoothpen.u(theta.u,P.u,penalty.u,0))/2+alpha.u[2]))                
    opt.theta.l<-optim(theta.l,fn=mklik.l,gr=mkgrad.l,method="BFGS",control=list(fnscale=-1),delta=delta,beta=beta,Z=Z,U.u=U.u,B.l=B.l,C.l=C.l,sigma2.l=sigma2.l,P.l=P.l,penalty=penalty.l,hessian=FALSE)
    theta.l<-opt.theta.l$par
    opt.theta.u<-optim(theta.u[-fixthetau],fn=mklik.u,gr=mkgrad.u,method="BFGS",control=list(fnscale=-1),B.u=B.u,E.u=E.u,sigma2.u=sigma2.u,P.u=P.u,M=M,penalty=penalty.u,hessian=FALSE,fix=fixthetau)    
    theta.u<-theta.u.repair(opt.theta.u$par,fixthetau)
    
    #browser()
    
    # Evaluate variances and hessians for candidate generation
    nu2.U<-diff(range(Ui))^2/6
    hess.b<-mkhess.b(beta,delta,Z,U.u,C.l,theta.l,sigma2.b)
    Sigma.b<-solve(-hess.b)
    hess.l<-mkhess.l(theta.l,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty.l)
    Sigma.l<-try(solve(-hess.l),silent=TRUE)
    d<-10
    while(inherits(Sigma.l,"try-error")){ Sigma.l<-try(solve(-(hess.u-10^(-d)*mdiag(K.u-1))),silent=TRUE); d<-d-1}
    d<-10
    while(!all(eigen(Sigma.l)$val>0)){ Sigma.l<-solve(-(hess.l-10^(-d)*mdiag(K.l))) ;d<-d-1  }
    #if(det(Sigma.l<0)) stop("Bad Hessian for baseline spline likelihood")
    hess.u<-mkhess.u(theta.u[-fixthetau],B.u,E.u,sigma2.u,P.u,M,penalty.u,fixthetau)
    Sigma.u<-try(solve(-hess.u),silent=TRUE)
    d<-10
    while(inherits(Sigma.u,"try-error")){ Sigma.u<-try(solve(-(hess.u-10^(-d)*mdiag(K.u-1))),silent=TRUE); d<-d-1}
    d<-10
    while(!all(eigen(Sigma.u)$val>0)){ Sigma.u<-solve(-(hess.u-10^(-d)*mdiag(K.u-1))); d<-d-1}

    #if(det(Sigma.u<0)) stop("Bad Hessian for frailty density spline likelihood")
    
    # Store initial values in parameter history
    paramhist<-list(Ui=Ui,beta=beta,theta.l=theta.l,theta.u=theta.u,sigma2.b=sigma2.b,sigma2.l=sigma2.l,sigma2.u=sigma2.u,acc.U=0,acc.b=0,acc.l=0,acc.u=0)
    avg.gammahist<-NULL
    avg.accepthist<-NULL
    
    ######################
    # main<-function() Main Loop
        
    #browser()
    if(verbose>=1) cat("Starting MCMC...\n")
    
    iter<-0
    acc.U<-0;acc.b<-0;acc.l<-0;acc.u<-0
    while(iter<maxiter)
    {
        iter<-iter+1
        
        # MH update of frailties
        mhout.U<-mh.U(Ui,i.all,delta,beta,Z,C.l,theta.l,B.u,theta.u,nu2.U,gamma.U,knots.u,ord.u,Bint.u)
        Ui<-mhout.U$Ui
        acc.U<-mhout.U$acc
        U.u<-rep(Ui,Ji)
        B.u<-splineDesign(knots.u,x=Ui,ord=ord.u)
        for(i in 1:dim(B.u)[1]) B.u[i,]<-B.u[i,]/Bint.u
        
        # MH update of regression parameters
        mhout.b<-mh.b(beta,delta,Z,U.u,C.l,theta.l,sigma2.b,Sigma.b,gamma.b)
        beta<-mhout.b$beta
        acc.b<-mhout.b$acc
        
        # MH update of baseline parameters
        mhout.l<-mh.l(theta.l,delta,beta,Z,U.u,B.l,C.l,sigma2.l,P.l,penalty.l,Sigma.l,gamma.l)
        theta.l<-mhout.l$theta.l
        acc.l<-mhout.l$acc
        
        # MH update of frailty density parameters
        mhout.u<-mh.u(theta.u,B.u,E.u,sigma2.u,P.u,M,penalty.u,fixthetau,Sigma.u,gamma.u)
        theta.u<-mhout.u$theta.u
        acc.u<-mhout.u$acc
        
        # Update of the sigmas
        sigma2.b<-1/rgamma(1,length(beta)/2+alpha.b[1],rate=as.numeric(t(beta)%*%beta)/2+alpha.b[2])
        sigma2.l<-1/rgamma(1,K.l/2+alpha.l[1],rate=smoothpen.l(theta.l,P.l,penalty.l,0)/2+alpha.l[2])
        sigma2.u<-1/rgamma(1,K.u/2+alpha.u[1],rate=smoothpen.u(theta.u,P.u,penalty.u,0)/2+alpha.u[2])
        #if(penalty.l=="2diff") sigma2.l<-1/rgamma(1,K.l/2+alpha.l[1],1/(as.numeric(t(theta.l)%*%P.l%*%theta.l)/2+alpha.l[2]))
        #if(penalty.u=="2diff") sigma2.u<-1/rgamma(1,K.u/2+alpha.u[1],1/(as.numeric(t(theta.u)%*%P.u%*%theta.u)/2+alpha.u[2]))
        #if(penalty.l=="2deriv") sigma2.l<-1/rgamma(1,K.l/2+alpha.l[1],1/(as.numeric(t(exp(theta.l))%*%P.l%*%exp(theta.l))/2+alpha.l[2]))
        #if(penalty.l=="2deriv") sigma2.u<-1/rgamma(1,K.u/2+alpha.u[1],1/(as.numeric(t(exp(theta.u))%*%P.u%*%exp(theta.u))/2+alpha.u[2]))        
        
        # Update parameter history
        paramhist$Ui<-rbind(paramhist$Ui,Ui)
        paramhist$beta<-rbind(paramhist$beta,beta)
        paramhist$theta.l<-rbind(paramhist$theta.l,theta.l)
        paramhist$theta.u<-rbind(paramhist$theta.u,theta.u)
        paramhist$acc.U<-c(paramhist$acc.U,acc.U)
        paramhist$acc.b<-c(paramhist$acc.b,acc.b)
        paramhist$acc.l<-c(paramhist$acc.l,acc.l)
        paramhist$acc.u<-c(paramhist$acc.u,acc.u)
        paramhist$sigma2.b<-c(paramhist$sigma2.b,sigma2.b)
        paramhist$sigma2.l<-c(paramhist$sigma2.l,sigma2.l)
        paramhist$sigma2.u<-c(paramhist$sigma2.u,sigma2.u)
        
        
       # if(iter%%1000==0) browser()
        
        # Periodic calibration check
        if(iter%%calint==0 & iter<maxiter){
            
            if(verbose==1 | verbose==2) cat(iter," ")
            
            #if(iter>=burnin.max) burnin<-iter
        
            # Calibration of the tuning parameters for acceptance rate
            if(cal.gamma & iter<=burnin){
               if(verbose>=3) cat("\n Calibration ...\n")
              # if(iter==500) browser()
               gamma<-c(gamma.U,gamma.b,gamma.l,gamma.u)
                avg.gammahist<-rbind(avg.gammahist,gamma)
                avg.accepthist<-rbind(avg.accepthist,c(mean(paramhist$acc.U[(iter-calint+1):iter]),mean(paramhist$acc.b[(iter-calint+1):iter]),
                                                       mean(paramhist$acc.l[(iter-calint+1):iter]),mean(paramhist$acc.u[(iter-calint+1):iter])))
                for(g in 1:4){
                    if(all(avg.accepthist[,g]>.25)) gamma[g]<-gamma[g]*2
                    if(all(avg.accepthist[,g]<.25)) gamma[g]<-gamma[g]/2
                    if(any(avg.accepthist[,g]>.25) & any(avg.accepthist[,g]<.25)){
                        fit<-lm(y~x,data.frame(x=avg.gammahist[,g],y=avg.accepthist[,g]))
                        gamma[g]<-max(.01,nlm(accrate.predict.lm,gamma[g],m=fit)$est)
                    }
                }
                gamma.U<-gamma[1];gamma.b<-gamma[2];gamma.l<-gamma[3];gamma.u<-gamma[4]
                if(verbose>=4){
                    cat("Acceptance rates (U,b,l,u): ",avg.accepthist[dim(avg.accepthist)[1],],"\n")
                    cat("Tuning parameters (U,b,l,u): ",gamma,"\n")
                }
            }
            # Check whether burnin is done
            if(verbose>=4){
                thisUi<-matrix(paramhist$Ui[(iter-calint+1):iter,],calint,m)
                thisbeta<-matrix(paramhist$beta[(iter-calint+1):iter,],calint,length(beta))
                thistheta.l<-matrix(paramhist$theta.l[(iter-calint+1):iter,],calint,K.l)
                thistheta.u<-matrix(paramhist$theta.u[(iter-calint+1):iter,],calint,K.u)
                avgUi<-apply(thisUi,2,mean)
                avgbeta<-apply(thisbeta,2,mean)
                avgtheta.l<-apply(thistheta.l,2,mean)
                avgtheta.u<-apply(thistheta.u,2,mean)
                avgsigma2.b<-mean(paramhist$sigma2.b[(iter-calint+1):iter])
                avgsigma2.u<-mean(paramhist$sigma2.u[(iter-calint+1):iter])
                avgsigma2.l<-mean(paramhist$sigma2.l[(iter-calint+1):iter])
                if(verbose>=5){
                    cat("Ui: ",avgUi,"\n")
                    cat("beta: ",avgbeta,"\n")
                    cat("theta.l: ",avgtheta.l,"\n")
                    cat("theta.u: ",avgtheta.u,"\n")
                    cat("sigma2: ",avgsigma2.b,avgsigma2.l,avgsigma2.u,"\n")
                }
            }
        }
        
        if(verbose>=3) cat(iter," ")

    }
    if(verbose>0) cat("Done!\n")
    outspline=list(knots.u=knots.u,knots.l=knots.l,ord.u=ord.u,ord.l=ord.l,knotspacing.l=knotspacing.l,knotspacing.u=knotspacing.u)
    if(burnin<iter){
        sub<-(1:(iter+1))>(burnin+1)
        npost<-sum(sub)
        postmean.Ui<-apply(matrix(paramhist$Ui[sub,],npost,m),2,mean)
        postmean.beta<-apply(matrix(paramhist$beta[sub,],npost,length(beta)),2,mean)
        postmean.l<-apply(matrix(paramhist$theta.l[sub,],npost,K.l),2,mean)
        postmean.u<-apply(matrix(paramhist$theta.u[sub,],npost,K.u),2,mean)
        postmean<-list(Ui=postmean.Ui,beta=postmean.beta,theta.l=postmean.l,theta.u=postmean.u)
        names(postmean$beta)<-varnames
    }else{
        postmean=NULL
    }
    posterior.mean<-list(frailty=postmean$Ui,coefficients=postmean$beta,splinepar.haz=postmean$theta.l,splinepar.frail=postmean$theta.u)
    outspline2<-list(knots.haz=knots.l,knots.frail=knots.u,ord.haz=ord.l,ord.frail=ord.u,knotspacing.haz=knotspacing.l,knotspacing.frail=knotspacing.u)
    history<-list(frailty=paramhist$Ui,coefficients=paramhist$beta,splinepar.haz=paramhist$theta.l,splinepar.frail=paramhist$theta.u,
        sigma2=cbind(paramhist$sigma2.b,paramhist$sigma2.l,paramhist$sigma2.u),acc=cbind(paramhist$acc.U,paramhist$acc.b,paramhist$acc.l,paramhist$acc.u))
    if(coda){
        library(coda)
        history$frailty<-mcmc(history$frailty); history$coefficients<-mcmc(history$coefficients)
        history$splinepar.haz<-mcmc(history$splinepar.haz); history$splinepar.frail<-as.mcmc(history$splinepar.frail)
        history$sigma2<-mcmc(history$sigma2); history$acc<-mcmc(history$acc)
    }
    colnames(history$sigma2)<-c("sigma2.coef","sigma2.par.haz","sigma2.par.frail")
    colnames(history$acc)<-c("acc.frail","acc.coef","acc.par.haz","acc.par.frail")
    rownames(history$frailty)<-rownames(history$coefficients)<-rownames(history$splinepar.haz)<-rownames(history$splinepar.frail)<-NULL
    control$iter<-iter
    out<-list(call=call,history=history,posterior.mean=posterior.mean,spline=outspline2,control=control)
    class(out)<-"splinesurv"
    return(out)
}