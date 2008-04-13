# Todo: (high) try normalizing the penalties for the frailty estimates
# Todo: (med) allow plotting spline and parametric component separately
# Todo: (med) allow plotting credible bands for the splines
# Todo: (low) Needs input checking
# Todo: (low) allow knots to birth/death/move
# Todo: (low) update history to also include tuning parameter data
# Todo: (low) allow using tkrplot to explore estimates by iteration

#{{{ #Header
#library(survival)
#library(MASS)
#dyn.load("C/init.so")
#}}}

{{{ #Simulation
##############################################################
# \section{Simulation} Generate simulated data
##############################################################

# my MYmvrnorm for testing (to make sure normal variates in C code are the same)
# Todo: Remove this function eventually
MYmvrnorm<- function (n = 1, mu, Sigma) {
    l<-length(mu)
    candnorm<-matrix(rnorm(n*l),nrow=n)
    out<-candnorm%*%chol(Sigma, pivot=TRUE)
    out<-out + rep(mu,rep(n,l)) 
    if(n==1) out<-as.numeric(out)
    out
}

generaterandom<-function(n,type,params)
{
    if(!(type%in%c("fixed","weibull","gamma","normal","lognormal","normmix","lognormmix","unifmix"))) stop("Invalid distribution type")
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
        out<-MYmvrnorm(n,mu,mdiag(sigma2)) 
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
        out<-MYmvrnorm(n,muprime,mdiag(sigma2prime)) 
        return(exp(t(out)[findInterval(runif(n),cumsum(w))+1+0:(n-1)*length(w)]))   
    }
    if(type=="unifmix"){
        if(!all(c("w","bounds")%in%names(params))) stop("Parameters w, bounds not specified for type unifmix")
        w<-params$w/sum(params$w)
        bounds<-matrix(params$bounds,ncol=2,byrow=TRUE)
        out<-rep(0,n)
        for(i in 1:n){
            which<-sum(runif(1)>cumsum(w))+1
            out[i]<-runif(1,bounds[which,1],bounds[which,2])
        }
        return(out)
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
            Tijprop<-0
            maxhaz<-max(haz)*Uij[ind]*exp(beta*Zij[ind])
            while(!accept){
                Tijprop<- Tijprop -1/maxhaz*log(runif(1))
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
            Tijprop<-0
            while(!accept){
                maxhaz<-max(w)*Uij[ind]*exp(beta*Zij[ind])
                Tijprop<- min(Tijprop -1/maxhaz*log(runif(1)),rbound-1e-5)
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
    h<-(predict(b,x)%*%w)
    if(is.null(t)){
        H<-cumsum(h*rbound/n)
        S<-exp(-H)
        return(S)
    }else{
        Ht<-sum(h[x<t]*rbound/n)
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

}}}

{{{ #Utility
##############################################################
# \section{Utility} Utility and Extraction functions
##############################################################

hasspline<-function(curve) return(curve$type=="spline" | curve$type=="both")
haspar<-function(curve) return(curve$type=="parametric" | curve$type=="both")
mdiag<-function(x) if(is.vector(x) && length(x)==1 && x<1) return(matrix(x,1,1)) else return(diag(x))

repairfrailtypar<-function(par,ind)
{
    if(ind==1) return(c(0,par))
    if(ind==length(par)) return(c(par,0))
    return(c(par[1:(ind-1)],0,par[ind:length(par)]))
}

inverthessian<-function(hess){
    K<-dim(hess)[1]
    Sigma<-try(solve(-hess),silent=TRUE)
    d<-10
    while(inherits(Sigma,"try-error")){ Sigma<-try(solve(-(hess-10^(-d)*mdiag(K))),silent=TRUE); d<-d-1}
    while(!all(eigen(Sigma)$val>0)){ Sigma<-solve(-(hess-10^(-d)*mdiag(K))) ;d<-d-1  }
    return(Sigma)
}

numHess.par<-function(param.par,fun,eps=1e-5,...)
{
    numDer.par<-function(param.par,fun,eps=1e-5,...)
    {        
        lik1<-fun(param.par,...)
        nd<-rep(0,length(param.par))
        for(i in 1:length(nd))
        {
            param.par2<-param.par
            param.par2[i]<-param.par2[i]+eps
            lik2<-fun(param.par2,...)
            nd[i]<-(lik2-lik1)/eps
        }
        return(nd)
    }
    nh<-matrix(0,length(param.par),length(param.par))
    gr1<-numDer.par(param.par,fun,eps,...)
    for(i in 1:length(param.par))
    {
        param.par2<-param.par
        param.par2[i]<-param.par2[i]+eps
        gr2<-numDer.par(param.par2,fun,eps,...)
        nh[i,]<-(gr2-gr1)/eps
    }
    return(nh)
}

inithistory<-function(hazard,frailty,regression,control)
{
    history<-NULL
    maxiter<-control$maxiter
    history$frailty<-matrix(0,maxiter,length(frailty$x))
    history$coefficients<-matrix(0,maxiter,length(regression$coefficients))
    if(hazard$hasspline) history$hazard.spline.par<-matrix(0,maxiter,length(hazard$spline.par))
    if(frailty$hasspline) history$frailty.spline.par<-matrix(0,maxiter,length(frailty$spline.par))
    if(hazard$haspar) history$hazard.param.par<-matrix(0,maxiter,length(hazard$param.par))
    if(frailty$haspar) history$frailty.param.par<-matrix(0,maxiter,length(frailty$param.par))
    if(hazard$hasspline & hazard$haspar) history$hazard.weight<-matrix(0,maxiter,1)
    if(frailty$hasspline & frailty$haspar) history$frailty.weight<-matrix(0,maxiter,1)
    history$priorvar<-matrix(0,maxiter,7); 
    colnames(history$priorvar)<-c("coefficients","hazard.spline","frailty.spline",
            "hazard.param","frailty.param","hazard.weight","frailty.weight")
    history$accept<-matrix(0,maxiter,8)
    colnames(history$accept)<-c(colnames(history$priorvar),"frailty")
    history<-updatehistory(history,1,hazard,frailty,regression)
    return(history)
}

updatehistory<-function(history,i,hazard,frailty,regression)
{
    history$frailty[i,]<-frailty$x
    history$coefficients[i,]<-regression$coefficients
    if(hazard$hasspline) history$hazard.spline.par[i,]<-hazard$spline.par
    if(frailty$hasspline) history$frailty.spline.par[i,]<-frailty$spline.par
    if(hazard$haspar) history$hazard.param.par[i,]<-hazard$param.par
    if(frailty$haspar) history$frailty.param.par[i,]<-frailty$param.par
    if(hazard$hasspline & hazard$haspar) history$hazard.weight[i]<-hazard$weight
    if(frailty$hasspline & frailty$haspar) history$frailty.weight[i]<-frailty$weight
    history$priorvar[i,]<-c(regression$priorvar,hazard$spline.priorvar,frailty$spline.priorvar,
        hazard$param.priorvar,frailty$param.priorvar,hazard$weight.priorvar,frailty$weight.priorvar)
    history$accept[i,]<-c(regression$accept,hazard$spline.accept,frailty$spline.accept,
        hazard$param.accept,frailty$param.accept,hazard$weight.accept,frailty$weight.accept,frailty$accept)

    return(history)
}

rinvgamma<-function(n,shape,scale=1)
{

    return(1/rgamma(n,shape,rate=scale))
} 

logit<-function(p) log(p/(1-p))
invlogit<-function(x) 1/(1+exp(-x))

accrate.predict.lm<-function(x,m) (predict(m,data.frame(x=x))-.25)^2

submean<-function(x,subset,f=mean) {
    if(is.null(x)) return(NULL)
    if(!is.null(dim(x))) return(apply(x[subset,,drop=F],2,f))
    else return(f(x[subset]))
}

makeoutputcurve<-function(curve)
{
    outcurve<-list(name=curve$name,
                    type=curve$type,
                    spline.knots=curve$spline.knots,
                    spline.nknots=curve$spline.nknots,
                    spline.knotspacing=curve$spline.knotspacing,
                    spline.ord=curve$spline.ord,
                    spline.norm=curve$spline.norm,
                    spline.penalty=curve$spline.penalty,
                    param.dist=curve$param.dist
                )
    return(outcurve)   
}
}}}

{{{ #Initialize
##############################################################
# \section{Initialize} Initialize knots, penalties, etc for spline components
##############################################################

makeknots<-function(curve,x,bounds=NULL)
{
    if(!curve$hasspline) return(curve)
    knots<-curve$spline.knots; nknots<-curve$spline.nknots;
    knotspacing<-curve$spline.knotspacing; ord<-curve$spline.ord
    K<-ord+nknots
    if(is.null(bounds)) bounds<-c(min(x),max(x))
    if(is.null(knots)){
        if(knotspacing=="quantile"){
            ibounds<-c(min(x),max(x))
            nintknots<-nknots; lrep<-ord; rrep<-ord
            if(ibounds[1]==bounds[1]) {nintknots<-nintknots+1; lrep<-ord-1}
            if(ibounds[2]==bounds[2]) {nintknots<-nintknots+1; rrep<-ord-1}
            knots<-quantile(x,seq(from=0,to=1,length=nintknots))
            knots<-c(rep(bounds[1],lrep),knots,rep(bounds[2],rrep))
        }
        if(knotspacing=="equal"){
            knots<-seq(from=bounds[1],to=bounds[2],length=nknots+2)
            knots<-c(rep(bounds[1],ord-1),knots,rep(bounds[2],ord-1))
        }
        if(knotspacing=="mindiff"){
            ibounds<-c(min(x),max(x))
            nintknots<-nknots; lrep<-ord; rrep<-ord
            if(ibounds[1]==bounds[1]) {nintknots<-nintknots+1; lrep<-ord-1}
            if(ibounds[2]==bounds[2]) {nintknots<-nintknots+1; rrep<-ord-1}
            p<-seq(from=0,to=1,length=nintknots)
            q<-quantile(x,p)
            while(min(diff(q))<diff(ibounds/100))
            {
                mind<-which.min(diff(q))
                p[mind+1]<-(p[mind+1]+p[mind+2])/2
                q<-quantile(x,p)
            }
            knots<-c(rep(bounds[1],lrep),q,rep(bounds[2],rrep))
        }
    }
    attr(knots,"boundary")<-bounds
    attr(knots,"index")<-seq(from=-(ord-1),length=length(knots),by=1)
    attr(knots,"order")<-ord
    curve$spline.knots<-knots
    return(curve)
}

makesplinebasis<-function(curve,quick=FALSE,usec=TRUE)
{
    if(!curve$hasspline) return(curve)
    knots<-curve$spline.knots; ord<-curve$spline.ord; x<-curve$x
    if(usec) B<-csplinedesign(knots,x=x,ord=ord) else B<-splineDesign(knots,x=x,ord=ord)
    if(usec) Bint<-cevalBinte(knots,ord) else Bint<-evalBinte(knots,ord)
    if(curve$spline.norm) for(i in 1:dim(B)[1]) B[i,]<-B[i,]/Bint
    if(!quick) { 
        if(curve$name=="hazard")  {
            if (usec) C<-cevalCinte(knots,ord,x,Bint)
            else C<-evalCinte(knots,ord,x,Bint)
        }
        else C<-NULL
        if(curve$name=="frailty") {
            if(usec) E<-cevalEinte(knots,ord)
            else E<-evalEinte(knots,ord)
        }
        else E<-NULL
        curve$spline.basiscum<-C
        curve$spline.basisexp<-E
    }
    curve$spline.basisint<-Bint
    curve$spline.basis<-B
    return(curve)
}

ki<-function(knots,j) return(j+attr(knots,"o"))

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

makepenalty<-function(curve,usec=TRUE)
{
    if(!curve$hasspline) return(curve)
    penalty<-curve$spline.penalty
    ord<-curve$spline.ord; nknots<-curve$spline.nknots; knots<-curve$spline.knots
    if(penalty=="2diff") P<-makePenalty.2diff(ord+nknots)
    if(penalty=="2deriv" | penalty=="log2deriv") {
        if(usec) P<-cmakePenalty.2deriv(ord,knots)
        else  P<-makePenalty.2deriv(ord,knots)
        if(curve$spline.norm){
            Bint<-curve$spline.basisint
            P<-P/(Bint%*%t(Bint))
        }
    }
    if(penalty=="none") P<-0
    curve$spline.penaltymatrix<-P
    return(curve)
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

smoothpen<-function(curve,der=0)
{
    type<-curve$spline.penalty
    name<-curve$name
    theta<-curve$spline.par
    P<-curve$spline.penaltymatrix
    sigma2<-curve$spline.priorvar
    if(der>=2) stop("second derivative not implemented")
    if(type=="2diff"){
        if(der==0) return(max( t(theta)%*%P%*%theta / (2*sigma2),0))
        if(der==1) return( P%*%theta /sigma2)
        #if(der==2) return( P / sigma2)
    }
    if(type=="2deriv"){
        et<-exp(theta)
        if(der==0) return(max( t(et)%*%P%*%et / (2*sigma2),0))
        if(der==1) return( mdiag(et)%*%P%*%et / sigma2 )
        #if(der==2) {
        #    pen<-mdiag(et)%*%P%*%mdiag(et)
        #    pen<-pen+mdiag(as.vector(P%*%et))%*%mdiag(et)
        #    return(pen / sigma2 )
        #}
    }
    if(type=="log2deriv"){
        et<-exp(theta)
        ePe<- as.numeric(t(et)%*%P%*%et) 
        if(der==0) return(max(log(ePe+1)/ (2*sigma2),0))
        if(der==1) return( mdiag(et)%*%P%*%et / sigma2 /(ePe+1))
    }
    if(type=="none"){
        if(der==0) return(0)
        if(der==1) return(rep(0,length(theta)))
        #if(der==2) return(matrix(0,length(theta),length(theta)))
    }
}


}}}

{{{ #CurveUpdate
##############################################################
# \section{CurveUpdate} Curve updating routines
##############################################################

fitparametric<-function(curve,x)
{
    name<-curve$name
    dist<-curve$param.dist
    if(dist=="none") return(curve)
    if(name=="frailty")
    {
        Ui<-x
        if(dist=="gamma")
            par<-log(var(Ui))
        if(dist=="lognormal")
        {
            varu<-var(Ui)
            par<-log(log(varu+1))
        }
        curve$param.par<-par
        curve$x<-Ui
    }
    if(name=="hazard")
    {
        agdata<-x
        varnames<-colnames(agdata)[-(1:4)]
        qvarnames<-paste("`",varnames,"`",sep="")
        fit<-survreg(as.formula(paste("Surv(time,delta)~",paste(qvarnames,collapse="+"))),data=agdata,dist=dist)
        if(dist=="exponential"){
            par<-log(fit$icoef[1])
        }
        if(dist=="weibull"){
            lambda<-exp(-fit$icoef[1])
            gamma<-1/exp(fit$icoef[2])
            par<-c(log(lambda),log(gamma))
        }
        if(dist=="lognormal"){
            par<-c(fit$icoef[1],fit$icoef[2])
        }
        names(par)<-NULL
        curve$param.par<-par
        curve$x<-agdata$time
    }
    curve<-evalparametric(curve)
    return(curve)
}

evalparametric<-function(curve,i=0)
{
    if(!curve$haspar) return(curve)
    if(i==0) ind<-1:length(curve$x) else ind<-i
    name<-curve$name
    dist<-curve$param.dist
    if(dist=="none") return(curve)
    par<-curve$param.par
    x<-curve$x[ind]
    if(name=="hazard"){
        if(dist=="exponential"){
            lambda<-exp(par)
            y<-rep(lambda,length(x))
            ycum<-x*lambda
        }
        if(dist=="weibull"){
            lambda<-exp(par[1])
            alpha<-exp(par[2])
             y<-alpha*lambda*x^(alpha-1)
             ycum<-lambda*x^alpha
        }
        if(dist=="lognormal")
            stop("lognormal distribution currently not fully supported")
    }
    if(name=="frailty"){
        ycum<-NULL
        if(dist=="gamma"){
            alpha<-exp(-par)
            y<-dgamma(x,shape=alpha,rate=alpha)
        }
        if(dist=="lognormal"){
            alpha<-exp(par)
            y<-exp(-(log(x)+alpha/2)^2/(2*alpha))/(x*sqrt(2*pi*alpha))
        }
    }
    curve$param.y[ind]<-y 
    curve$param.ycum[ind]<-ycum
    if(curve$hasspline) {
         curve<-weightcurve(curve,i)
    }else{
        curve$y[ind]<-curve$param.y[ind]
        curve$ycum[ind]<-curve$param.ycum[ind]
    }  
    return(curve)
}

evalspline<-function(curve,i=0,quick=FALSE)
{
    if(!curve$hasspline) return(curve)
    if(i==0) {
		ind<-1:length(curve$x)
		curve$spline.y<-NULL
	}else{ind<-i}
    spline.par<-curve$spline.par
    if(curve$spline.norm) 
        curve$spline.y[ind]<-as.vector(curve$spline.basis[ind,,drop=FALSE]%*%exp(spline.par)/sum(exp(spline.par)))
    else
        curve$spline.y[ind]<-as.vector(curve$spline.basis[ind,,drop=FALSE]%*%exp(spline.par))
    if(curve$name=="hazard" & !quick)
        curve$spline.ycum[ind]<-as.vector(curve$spline.basiscum[ind,,drop=FALSE]%*%exp(spline.par))
    else
        curve$spline.ycum<-NULL
    if(curve$haspar) {
        curve<-weightcurve(curve,i)
    }else{
        curve$y[ind]<-curve$spline.y[ind]
        curve$ycum[ind]<-curve$spline.ycum[ind]
    }
    return(curve)
}

updatespline<-function(curve,spline.par)
{
    if(!curve$hasspline) return(curve)
    if(curve$name=="frailty") spline.par<-repairfrailtypar(spline.par,curve$spline.fixedind)
    curve$spline.par<-spline.par
    curve<-evalspline(curve)
    return(curve)
}

updatecurvex<-function(curve,i)
{
    if(curve$name!="frailty") stop("Only frailty bases can be updated")
    if(curve$hasspline){
        knots<-curve$spline.knots; ord<-curve$spline.ord; x<-curve$x[i]
        curve$spline.basis[i,]<-mysplineDesign(knots,x,ord)/curve$spline.basisint 
        curve<-evalspline(curve,i)
    }
    curve<-evalparametric(curve,i)
    return(curve)   
}

updateparametric<-function(curve,param.par)
{
    if(!curve$haspar) return(curve)
    curve$param.par<-param.par
    curve<-evalparametric(curve)
    return(curve)
}

weightcurve<-function(curve,i=0)
{
    if(i==0) ind<-1:length(curve$x) else ind<-i
    curve$y[ind]<-curve$weight*curve$spline.y[ind]+(1-curve$weight)*curve$param.y[ind]
    if(curve$name=="hazard")
        curve$ycum[ind]<-curve$weight*curve$spline.ycum[ind]+(1-curve$weight)*curve$param.ycum[ind]
    return(curve)
}

updateregression<-function(regression,coef)
{
    regression$coefficients<-coef
    regression$lp<-as.vector(regression$covariates%*%regression$coefficients)
    return(regression)
}

}}}

{{{ #Likelihood
##############################################################
# \section{Likelihood} Likelihoods, gradients and hessians
##############################################################

mklik.coef<-function(coef,hazard,frailty,regression)
{
    status<-regression$status
    regression<-updateregression(regression,coef)
    lp<-regression$lp
    frailrep<-rep(frailty$x,regression$Ji)
    lik<-status%*%lp
    lik<-lik-sum(frailrep*hazard$ycum*exp(lp))
    lik<-lik-sum(coef^2)/(2*regression$priorvar)
    return(as.numeric(lik))
}

mkhess.coef<-function(coef,hazard,frailty,regression)
{
    status<-regression$status
    regression<-updateregression(regression,coef)
    lp<-regression$lp
    frailrep<-rep(frailty$x,regression$Ji)
    Z<-regression$covariates
    hess<--t(Z)%*%(rep(frailrep*exp(lp)*hazard$ycum,dim(Z)[2])*Z)
    hess<-hess-mdiag(rep(1,length(coef)))/regression$priorvar
    return(hess)
}

mklik.frail<-function(i,hazard,frailty,regression)
{
    ind<-which(regression$cluster==i)
    Ui<-frailty$x[i]
    status<-regression$status[ind]
    lp<-regression$lp[ind]
    cumhaz<-hazard$ycum[ind]
    lik<-sum(status*log(Ui)-Ui*cumhaz*exp(lp))+log(frailty$y[i])
    return(lik)   
}

mklik.spline.haz<-function(spline.par,hazard,frailty,regression)
{
    if(!hazard$hasspline) return(0)
    status<-regression$status
    lp<-regression$lp
    hazard<-updatespline(hazard,spline.par)
    frailrep<-rep(frailty$x,regression$Ji)
    lik<-sum(status*log(hazard$y) - frailrep*hazard$ycum*exp(lp))
    lik<-lik-hazard$spline.penaltyfactor*smoothpen(hazard,0)
    lik<-lik-sum(ifelse(spline.par< hazard$spline.min,(spline.par-hazard$spline.min)^2,0))    
    lik<-as.numeric(lik)
    return(lik)
}

mkgr.spline.haz<-function(spline.par,hazard,frailty,regression)
{
    if(!hazard$hasspline) return(rep(0,length(spline.par)))
    status<-regression$status
    lp<-regression$lp
    hazard<-updatespline(hazard,spline.par)
    frailrep<-rep(frailty$x,regression$Ji)
    Det<-diag(exp(hazard$spline.par))
    gr<-Det%*%t(hazard$spline.basis)%*%(status/hazard$y)
    gr<-gr-Det%*%t(hazard$spline.basiscum)%*%(frailrep*exp(lp))
    gr<-hazard$weight*gr
    gr<-gr-hazard$spline.penaltyfactor*smoothpen(hazard,1)
    gr<-gr-ifelse(spline.par< hazard$spline.min,2*(spline.par-hazard$spline.min),0)    
    gr<-as.numeric(gr)
    return(gr)
}


mklik.spline.frail<-function(spline.par,hazard,frailty,regression)
{
    if(!frailty$hasspline) return(0)
    frailty<-updatespline(frailty,spline.par)
    M<-frailty$spline.meanpenalty
    lik<-sum(log(frailty$y))
    lik<-lik-frailty$spline.penaltyfactor*smoothpen(frailty,0)
    lik<-lik-M*(frailty$spline.basisexp%*%exp(frailty$spline.par))^2
    lik<-lik-sum(ifelse(spline.par< frailty$spline.min,(spline.par-frailty$spline.min)^2,0))    
    lik<-as.numeric(lik)
    return(lik)
}

mklik.param.haz<-function(param.par,hazard,frailty,regression)
{
    if(!hazard$haspar) return(0)
    hazard<-updateparametric(hazard,param.par)
    status<-regression$status
    lp<-regression$lp
    frailrep<-rep(frailty$x,regression$Ji)
    lik<-sum(status*log(hazard$y) - frailrep*hazard$ycum*exp(lp))
    lik<-lik-sum(param.par^2)/(2*hazard$param.priorvar)
    return(lik)
}

mklik.param.frail<-function(param.par,hazard,frailty,regression)
{
    if(!frailty$haspar) return(0)
    frailty<-updateparametric(frailty,param.par)
    lik<-sum(log(frailty$y))-sum(param.par^2)/(2*frailty$param.priorvar)
    return(lik)
}

mklik.weight.haz<-function(weight,hazard,frailty,regression)
{
    hazard$weight<-weight
    hazard<-weightcurve(hazard)
    status<-regression$status
    lp<-regression$lp
    frailrep<-rep(frailty$x,regression$Ji)
    lik<-sum(status*log(hazard$y)-frailrep*hazard$ycum*exp(lp))
    #lik<-lik-hazard$weight^2/(2*hazard$weight.priorvar)
    lik <- lik + (hazard$weight.hyper[1]-1)*log(hazard$weight) + (hazard$weight.hyper[2]-1)*log(1-hazard$weight)
    return(lik) 
}

mklik.weight.frail<-function(weight,hazard,frailty,regression)
{
    frailty$weight<-weight
    frailty<-weightcurve(frailty)
    lik<-sum(log(frailty$y))
    #lik<-lik-frailty$weight^2/(2*frailty$weight.priorvar)
    lik <- lik + (frailty$weight.hyper[1]-1)*log(frailty$weight) + (frailty$weight.hyper[2]-1)*log(1-frailty$weight)
    return(lik) 
}

}}}

{{{ #Metropolis
##############################################################
# \section{Metropolis} Metropolis-Hastings
##############################################################

acceptreject<-function(baselik,candlik,ratio=1)
{
    if(is.nan(candlik)) candlik<--Inf
    r<-exp(candlik-baselik)*ratio
    p.accept<-min(1,r)
    if(is.nan(p.accept)) p.accept<-0
    if(runif(1)<p.accept){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

mh<-function(par,fun,candcov,tun,...)
{
    baselik<-fun(par,...)
    cand<-MYmvrnorm(1,par,candcov*tun)
    candlik<-fun(cand,...)
    acc<-acceptreject(baselik,candlik)
    if(acc) out<-cand else out<-par
    return(list(par=out,acc=acc))    
}


mh.frail<-function(hazard,frailty,regression)
{
    acc<-rep(0,regression$m)
    for(i in 1:regression$m){
        baselik<-mklik.frail(i,hazard,frailty,regression)
        u<-frailty$x[i]
        v<-frailty$tun
        cand<-rgamma(1,shape=u^2/v,scale=v/u)
        if(is.nan(cand) || (frailty$hasspline && 
          (cand>attr(frailty$spline.knots,"b")[2] | cand<attr(frailty$spline.knots,"b")[1]))) {
            candlik<- -Inf
        }else{
            temp<-frailty
            temp$x[i]<-cand
            temp<-updatecurvex(temp,i)
            candlik<-mklik.frail(i,hazard,temp,regression)
        }
        puc<-suppressWarnings(dgamma(cand,shape=u^2/v,scale=v/u)) # transition u->cand
        pcu<-suppressWarnings(dgamma(u,shape=cand^2/v,scale=v/cand)) # transition cand->u
        acc[i]<-acceptreject(baselik,candlik,pcu/puc)
        if(acc[i]) frailty<-temp
    }
    frailty$accept<-mean(acc)
    return(frailty)
}

mh.frailty.spline<-function(hazard,frailty,regression)
{
    if(!frailty$hasspline) return(frailty)
    mhout<-mh(frailty$spline.par[-frailty$spline.fixedind],mklik.spline.frail,frailty$spline.candcov,frailty$spline.tun,
        hazard=hazard,frailty=frailty,regression=regression)
    if(mhout$acc)  frailty<-updatespline(frailty,mhout$par)
    frailty$spline.accept<-mhout$acc
    return(frailty)
}

mh.hazard.spline<-function(hazard,frailty,regression)
{
    if(!hazard$hasspline) return(hazard)
    mhout<-mh(hazard$spline.par,mklik.spline.haz,hazard$spline.candcov,hazard$spline.tun,
        hazard=hazard,frailty=frailty,regression=regression)
    if(mhout$acc)  hazard<-updatespline(hazard,mhout$par)
    hazard$spline.accept<-mhout$acc
    return(hazard)
}

mh.frailty.param<-function(hazard,frailty,regression)
{
    if(!frailty$haspar) return(frailty)
    mhout<-mh(frailty$param.par,mklik.param.frail,frailty$param.candcov,frailty$param.tun,
        hazard=hazard,frailty=frailty,regression=regression)
    if(mhout$acc) frailty<-updateparametric(frailty,mhout$par)
    frailty$param.accept<-mhout$acc
    return(frailty)
}

mh.hazard.param<-function(hazard,frailty,regression)
{
    if(!hazard$haspar) return(hazard)
    mhout<-mh(hazard$param.par,mklik.param.haz,hazard$param.candcov,hazard$param.tun,
        hazard=hazard,frailty=frailty,regression=regression)
    if(mhout$acc) hazard<-updateparametric(hazard,mhout$par)
    hazard$param.accept<-mhout$acc
    return(hazard)
}

mh.coef<-function(hazard,frailty,regression)
{
    mhout<-mh(regression$coefficients,mklik.coef,regression$candcov,regression$tun,
        hazard=hazard,frailty=frailty,regression=regression)
    if(mhout$acc) regression<-updateregression(regression,mhout$par)
    regression$accept<-mhout$acc
    return(regression)
}

mh.weight<-function(which,hazard,frailty,regression)
{
    #browser()
    which<-match.arg(which,c("hazard","frailty"))
    if(which=="frailty"){
    	curve<-frailty
	fun<-mklik.weight.frail
    }
    if(which=="hazard"){
    	curve<-hazard
	fun<-mklik.weight.haz
    }
    if(!curve$haspar | !curve$hasspline) return(curve)
    w<-min(max(curve$weight,.01),.99)
    v<-curve$weight.tun
    alpha<-w*(w*(1-w)/v-1)
    beta<-(1-w)/w*alpha
    cand<-rbeta(1,alpha,beta)
    if(is.nan(cand)){
    	curve$weight.accept<-FALSE;
        return(curve)
    }
    alphac<-cand*(cand*(1-cand)/v-1)
    betac<-(1-cand)/cand*alphac
    baselik<-fun(w,hazard,frailty,regression)
    candlik<-fun(cand,hazard,frailty,regression)
    puc<-suppressWarnings(dbeta(cand,alpha,beta))
    pcu<-suppressWarnings(dbeta(w,alphac,betac))
    acc<-acceptreject(baselik,candlik,pcu/puc)
    if(acc){
	curve$weight<-cand
	curve<-weightcurve(curve)
    }
    curve$weight.accept<-acc
    return(curve)
}


updatepostvar.curve<-function(curve)
{
    if(curve$hasspline) curve$spline.priorvar<-rinvgamma(1,length(curve$spline.par)/2+curve$spline.hyper[1],
        scale=curve$spline.penaltyfactor*smoothpen(curve)*curve$spline.priorvar+curve$spline.hyper[2])
    if(curve$haspar) curve$param.priorvar<-rinvgamma(1,length(curve$param.par)/2+curve$param.hyper[1],
        scale=sum(curve$param.par^2)/2+curve$param.hyper[2])
#    if(curve$hasspline & curve$haspar) curve$weight.priorvar<-rinvgamma(1,1/2+curve$weight.hyper[1],
#        scale=curve$weight^2/2+curve$weight.hyper[2])
    return(curve)
}

updatepostvar.coef<-function(regression)
{
    regression$priorvar<-rinvgamma(1,length(regression$coefficients)/2+regression$hyper[1],
        scale=sum(regression$coefficients^2)/2+regression$hyper[2])
    return(regression)
}

}}}

{{{ #C-wrappers
    
csplinedesign<-function(knots,x,ord)
{
    K<-length(knots)-2*ord
    design<-matrix(0,length(x),K+ord)
    out<-.C("csplinedesign",
            des=as.double(design),
            x=as.double(x),
            nx=as.integer(length(x)),
            knots=as.double(knots),
            ord=as.integer(ord),
            K=as.integer(K)
        )
    des<-matrix(out$des,length(x),K+ord)
    return(des)
}

cevalEinte<-function(knots,ord,N=1)
{
    K<-length(knots)-2*ord
    einte<-rep(0,K+ord);
    out<-.C("cevalEinte",
            einte=as.double(einte),
            knots=as.double(knots),
            ord=as.integer(ord),
            K=as.integer(K),
            N=as.integer(N)
           )
    einte<-out$einte
    return(einte)
}

cevalBinte<-function(knots,ord)
{
    K<-length(knots)-2*ord;
    binte<-rep(0,K+ord);
    out<-.C("cevalBinte",
            binte=as.double(binte),
            knots=as.double(knots),
            ord=as.integer(ord),
            K=as.integer(K)
          )
    binte<-out$binte
    return(binte)
}

cevalCinte<-function(knots,ord,obs,Binte)
{
    K<-length(knots)-2*ord;
    cinte<-matrix(0,length(obs),length(Binte))
    out<-.C("cevalCinte",
            cinte=as.double(cinte),
            x=as.double(obs),
            nx=as.integer(length(obs)),
            knots=as.double(knots),
            ord=as.integer(ord),
            K=as.integer(K),
            binte=as.double(Binte)
        )
    cinte<-matrix(out$cinte,length(obs),K+ord)
    return(cinte)
}

cmakePenalty.2deriv<-function(ord,knots){
    K<-length(knots)-2*ord;
    P<-matrix(0,K+ord,K+ord)
    out<-.C("cMakePenalty2diff",
        P=as.double(P),
        knots=as.double(knots),
        ord=as.integer(ord),
        K=as.integer(K)
    )
    P<-matrix(out$P,K+ord,K+ord)
    return(P)
}

rmklik.spline.haz<-function(spline.par,status,lp,frailrep,hazParY,hazParYcum,weight,B,C,P,penaltyType,sigma2)
{
    hazSplineY <- B%*%exp(spline.par)
    hazY <- weight * hazSplineY + (1-weight) * hazParY
    hazSplineYcum <- C%*%exp(spline.par)
    hazYcum <- weight * hazSplineYcum + (1-weight) * hazParYcum
    lik<-sum(status*log(hazY) - frailrep*hazYcum*exp(lp))
    lik<-as.numeric(lik)
    if(penaltyType==1) lik <- lik-t((spline.par))%*%P%*%(spline.par)/(2*sigma2)
    if(penaltyType==2) lik <- lik-t(exp(spline.par))%*%P%*%exp(spline.par)/(2*sigma2)
    return(lik)
}

cmklik.spline.haz<-function(par,status,lp,frailrep,hazParY,hazParYcum,weight,B,C,P,penaltyType,sigma2,min)
{
    lik<-as.double(rep(0,1))
    out<-.C("cInitLikHazSpline",
            lik=lik,par=par,status=status,lp=lp,frailrep=frailrep,hazParY=hazParY,
            hazParYcum=hazParYcum,weight=weight,B=B,C=C,P=P,
            penaltyType=penaltyType,sigma2=sigma2,ny=as.integer(length(lp)),
            nj=as.integer(length(par)),DUP=FALSE)
    lik<-out$lik
    lik<-lik-sum(ifelse(par< min,(par-min)^2,0))    
    return(lik)
}

rmkgr.spline.haz<-function(spline.par,status,lp,frailrep,hazParY,hazParYcum,weight,B,C,P,penaltyType,sigma2)
{
    B<-matrix(B,length(lp),length(spline.par))
    C<-matrix(C,length(lp),length(spline.par))
    hazSplineY <- B%*%exp(spline.par)
    hazY <- weight * hazSplineY + (1-weight) * hazParY
    hazSplineYcum <- C%*%exp(spline.par)
    hazYcum <- weight * hazSplineYcum + (1-weight) * hazParYcum
    gr<-rep(0,length(spline.par))
    status=status/hazY
    lp=exp(lp)*frailrep
    gr<-gr+t(B)%*%status
    gr<-gr-t(C)%*%lp
    gr<-gr*exp(spline.par)
    return(gr)
}

cmkgr.spline.haz<-function(par,status,lp,frailrep,hazParY,hazParYcum,weight,B,C,P,penaltyType,sigma2,min)
{
    gr<-as.double(rep(0,length(par)))
    out<-.C("cInitGrHazSpline",
            gr=gr,par=par,status=status,lp=lp,frailrep=frailrep,hazParY=hazParY,
            hazParYcum=hazParYcum,weight=weight,B=B,C=C,P=P,
            penaltyType=penaltyType,sigma2=sigma2,ny=as.integer(length(lp)),
            nj=as.integer(length(par)),DUP=FALSE)
    gr<-out$gr
    gr<-gr-ifelse(par< min,2*(par-min),0)    
    return(gr)
}

cmklik.spline.frail<-function(par,fixedind, frailParY,weight,B, E, M, P, penaltyType, sigma2, min)
{
    par<-as.double(repairfrailtypar(par,fixedind))
    lik<-as.double(0);
    out<-.C("cInitLikFrailSpline", lik=lik, par=par, frailParY=frailParY,weight=weight, B=B, E=E, M=M, P=P, penaltyType=penaltyType, sigma2=sigma2, ny=as.integer(length(frailParY)), nj=as.integer(length(par)),DUP=FALSE)
    lik<-out$lik
    lik<-lik-sum(ifelse(par< min,(par-min)^2,0))    
    return(lik)
}
}}}

{{{ #S3Methods
##############################################################
# \section{S3Methods} S3 Methods for fitting, printing, summary, etc
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
    gcout<-gc()
    fit$call<-call
    colnames(fit$history$frailty)<-clusternames
    if(!is.null(fit$posterior.mean)) names(fit$posterior.mean$frailty)<-clusternames
    fit$terms<-newTerms
    return(fit)
}

################# Methods for printing

print.splinesurv<-function(x,...)
{
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
    out$hazard<-x$hazard
    out$frailty<-x$frailty
    out$frailty$spline.fvar<-spline.fvar(x)
    out$frailty$param.fvar<-if(haspar(x$frailty)) exp(x$posterior.mean$frailty.param.par) else NULL
    out$posterior.mean<-x$posterior.mean
    out$quantiles<-NULL
    if(out$iter<out$burnin) out$burnin<-0
    goodcoef<-x$history$coefficients[(out$burnin+1):(out$iter),,drop=FALSE]
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

spline.fvar<-function(x)
{
    if(!hasspline(x$frailty)) return(NULL)
    theta<-exp(x$posterior.mean$frailty.spline.par)
    moment2<-(1-cevalEinte(x$frailty$spline.knots,x$frailty$spline.ord,N=2))%*%theta/sum(theta)
    return(moment2-1)
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
    cat("\nBaseline hazard:")
    printcurvesummary(x$hazard,x$posterior.mean$hazard.weight,x$posterior.mean$hazard.param.par)
    cat("\nFrailty density:")
    printcurvesummary(x$frailty,x$posterior.mean$frailty.weight,x$posterior.mean$frailty.param.par)
    invisible(x)
}

printcurvesummary<-function(curve,w=NULL,param.par=NULL)
{
    haspar<-curve$type=="both" || curve$type=="parametric"
    hasspline<-curve$type=="both" || curve$type=="spline"
    if(hasspline) {
        cat("\n\tSpline")
        if(haspar) cat(" ( weight =", format(w,digits=3),"):") else cat(":")
        cat("\n\t\tOrder:",curve$spline.ord)
        cat("\n\t\tInterior knots:",curve$spline.nknots-2*curve$spline.ord+2)
        cat("\n\t\tKnot distribution:", curve$spline.knotspacing)
        cat("\n\t\tKnot boundaries: ", paste(format(attr(curve$spline.knots,"b"),digits=2),collapse=","))
        if(curve$name=="frailty") cat("\n\t\tVariance: ", format(curve$spline.fvar,digits=2))
    }
    if(haspar){
        cat("\n\tParametric")
        if(hasspline) cat(" ( weight =",format(1-w,digits=3),"):") else cat(":")
        dist<-curve$param.dist
        cat("\n\t\tDistribution:",curve$param.dist)
        if(curve$name=="hazard" & dist=="exponential"){
            cat("\n\t\tRate:",format(exp(param.par),digits=3))
        }
        if(curve$name=="hazard" & dist=="weibull"){
            cat("\n\t\tRate:",format(exp(param.par[1]),digits=3))
            cat("\n\t\tScale:",format(exp(param.par[2]),digits=3)) 
        }
        if(curve$name=="frailty" & (dist=="gamma" | dist=="lognormal")){
            cat("\n\t\tVariance:",format(exp(param.par),digits=3))
        }
        
    }
}


################# Methods for plotting


plot.splinesurv<-function(x,which=c("hazard","survival","frailty","coef","all"),newdata=NULL,iter=NULL,plotknots=TRUE,npoints=100,legend=NULL,
    lty=1,col=2,lwd=2,lty.knots=1,col.knots=8,lwd.knots=1,xlab=NULL,ylab=NULL,main=NULL,xlim=NULL,ylim=NULL,tk=FALSE,...)
{
    if(tk) splinesurvtkplot(x,newdata=NULL,plotknots=TRUE,npoints=100,legend=NULL,
    lty=1,col=2,lwd=2,lty.knots=1,col.knots=8,lwd.knots=1,xlab=NULL,ylab=NULL,main=NULL,xlim=NULL,ylim=NULL,tk=FALSE,...)
    oldask<-par("ask")
    which<-match.arg(which)
    if(which=="all") par(ask=TRUE)
    if(!is.null(iter) && iter<=0) iter<-NULL
    if(which=="hazard" | which=="all"){
        knots<-x$hazard$spline.knots
        if(is.null(knots)) knots<-range(x$data$time)
        if(is.null(xlim)) xlim1=range(x$data$time) else xlim1<-xlim
        times=seq(from=max(xlim1[1],min(knots)),to=min(xlim1[2],max(knots)),length=npoints)
        if(is.null(xlab)) xlab1<-"Time" else xlab1<-xlab
        if(is.null(ylab)) ylab1<-"Hazard" else ylab1<-ylab
        if(is.null(newdata)){
            haz<-predict(x,x=times,iter=iter)
            if(is.null(main)) main1<-"Baseline hazard" else main1<-main
            plot(haz,type="l",lty=lty,col=col,lwd=lwd,main=main1,xlab=xlab1,ylab=ylab1,xlim=xlim1,ylim=ylim,...)
            if(plotknots) {
                abline(v=knots,col=col.knots,lty=lty.knots,lwd=lwd.knots,...)
                lines(haz,lty=lty,col=col,lwd=lwd)
            }
        }else{
            if(is.null(main)) main1<-"Hazard" else main1<-main
            haz<-predict(x,type="risk",x=times,newdata=newdata,iter=iter)
            if(length(col)==1 & length(lty)==1 & length(lwd)==1) col=1:dim(newdata)[1]
            matplot(haz[1],haz[-1],type="l",col=col,lwd=lwd,lty=lty,main=main1,xlab=xlab1,ylab=ylab1,xlim=xlim1,ylim=ylim,...)
            if(is.null(legend)) legend<-rownames(newdata)
            if(plotknots) abline(v=knots,col=col.knots,lty=lty.knots,lwd=lwd.knots,...)
            legend("topleft",legend=legend,col=col,lty=lty,lwd=lwd)
        }
    }
    if(which=="survival" | which=="all"){
        knots<-x$hazard$spline.knots
        if(is.null(knots)) knots<-range(x$data$time)
        if(is.null(xlim)) xlim1=range(x$data$time) else xlim1<-xlim
        times=seq(from=max(xlim1[1],min(knots)),to=min(xlim1[2],max(knots)),length=npoints)
        if(is.null(xlab)) xlab1<-"Time" else xlab1<-xlab
        if(is.null(ylab)) ylab1<-"Survival" else ylab1<-ylab
        if(is.null(newdata)){
            haz<-predict(x,x=times,iter=iter)
            dx<-mean(diff(times))
            survival<-exp(-cumsum(haz$haz*dx))
            if(is.null(main)) main1<-"Survivor function" else main1<-main
            plot(times,survival,type="l",lty=lty,col=col,lwd=lwd,main=main1,xlab=xlab1,ylab=ylab1,xlim=xlim1,ylim=ylim,...)
            if(plotknots){
                abline(v=knots,col=col.knots,lty=lty.knots,lwd=lwd.knots,...)
                lines(times,survival,type="l",lty=lty,col=col,lwd=lwd)
            }
        }else{
            if(is.null(main)) main1<-"Survivor function" else main1<-main
            haz<-predict(x,type="risk",x=times,newdata=newdata,iter=iter)[,-1,drop=FALSE]
            survival<-exp(-apply(haz*diff(times),2,cumsum))
            if(length(col)==1 & length(lty)==1 & length(lwd)==1) col=1:dim(newdata)[1]
            matplot(times,survival,type="l",col=col,lwd=lwd,lty=lty,main=main1,xlab1=xlab,ylab1=ylab,xlim1=xlim,ylim=ylim,...)
            if(is.null(legend)) legend<-rownames(newdata)
            if(plotknots) abline(v=knots,col=col.knots,lty=lty.knots,lwd=lwd.knots,...)
            legend("topright",legend=legend,col=col,lty=lty,lwd=lwd)
        }
    }
    if(which=="frailty" | which=="all"){
        knots<-x$frailty$spline.knots
        if(is.null(knots)) knots<-range(x$posterior.mean$frailty)
        if(is.null(xlim)) {
               if(hasspline(x$frailty)) xlim1=attr(knots,"b") 
               else xlim1<-range(x$posterior.mean$frailty)
       }else{xlim1<-xlim}
        Ui=seq(from=max(xlim1[1],min(knots)),to=min(xlim1[2],max(knots)),length=npoints)
        if(is.null(xlab)) xlab1<-"x" else xlab1<-xlab
        if(is.null(ylab)) ylab1<-"Density" else ylab1<-ylab
        if(is.null(main)) main1<-"Frailty density" else main1<-main
        dens<-predict(x,type="frailty",x=Ui,iter=iter)
        plot(dens,type="l",lty=lty,col=col,lwd=lwd,main=main1,xlab=xlab1,ylab=ylab1,xlim=xlim1,ylim=ylim,...)
        if(plotknots){
            abline(v=knots,col=col.knots,lty=lty.knots,lwd=lwd.knots,...)
            lines(dens,type="l",lty=lty,col=col,lwd=lwd)
        }
    }
    if(which=="coef" | which=="all"){
        burnin<-x$control$burnin
        if(is.null(burnin)) burnin<-x$call$control$burnin
        if(is.null(burnin)) burnin<-dim(x$history$coefficients)[2]
        if(is.null(xlab)) xlab1<-"x" else xlab1<-xlab
        if(is.null(ylab)) ylab1<-"Posterior density" else ylab1<-ylab
        coefs<-x$history$coefficients
        coefnames<-colnames(coefs)
        if(length(coefnames)>1) par(ask=TRUE)
        for(i in 1:dim(coefs)[2]){
            if(is.null(main)) main1<-paste("Coefficient of",coefnames[i]) else main1<-main
            betai<-coefs[burnin:(dim(coefs)[1]),i]
            plot(density(betai),lty=lty,col=col,lwd=lwd,main=main1,xlab=xlab1,ylab=ylab1,xlim=xlim,ylim=ylim,...)
        }     
    }
    par(ask=oldask)
}

splinesurvtkplot<-function(x,...)
{
    library(tkrplot)
	which <- "h"
	tt <- tktoplevel()
	tktitle(tt)<-"SplineSurv"
	iter <- 0
	tcliter<-tclVar(iter)
	tclwhich<-tclVar(which)
	maxiter=x$control$maxiter
	res<-round(maxiter/100)
	img <- tkrplot(tt, function() plot(x=x,which=which,iter=iter,...))	
	setiter<-function(...){
		thisiter<-as.numeric(tclvalue(tcliter))
		if(iter != thisiter){
			assign("iter",thisiter,inherits=TRUE)
			tkrreplot(img)
		}	
	}
	setwhich_h<-function(...) {
		assign("which","h",inherits=TRUE)
		tkrreplot(img)
	}
	setwhich_s<-function(...) {
		assign("which","s",inherits=TRUE)
		tkrreplot(img)
	}
	setwhich_f<-function(...){
		assign("which","f",inherits=TRUE)		
		tkrreplot(img)
	}
	setpost<-function(...){
		assign("iter",0,inherits=TRUE)
		tclvalue(tcliter)<-0
		tkrreplot(img)
	}

	iter_scale <- tkscale(tt, command=setiter, from=0, to=maxiter, resolution=res, showvalue=T,orient="horiz",variable=tcliter, length=400)
	ff<-tkframe(tt,relief="ridge",borderwidth=2,width=150,height=100)
	which_b_h <- tkradiobutton(ff, command=setwhich_h,text="hazard",variable=tclwhich, value="h" )
	which_b_s <- tkradiobutton(ff, command=setwhich_s,text="survival",variable=tclwhich, value="s" )
	which_b_f <- tkradiobutton(ff, command=setwhich_f,text="frailty",variable=tclwhich, value="f" )
	post_b <- tkbutton(tt, command=setpost, text="Posterior" )
	tkgrid(img,columnspan=2)
	tkgrid(which_b_h,which_b_s,which_b_f)
	tkgrid(ff,columnspan=2)
	tkgrid(post_b,iter_scale)
	tkgrid.configure(iter_scale,sticky="e")
	tkgrid.configure(post_b,sticky="sw")

}



################# predict method

predict.splinesurv<-function(object,type=c("hazard","lp","risk","frailty"),x=NULL,newdata=NULL,iter=NULL,...)
{
    type<-match.arg(type)
    fit<-object; haz<-1
    ntimes<-100
    if(type=="hazard" | (type=="risk" & !is.null(x))){
        if(type=="hazard" & !is.null(newdata)) stop("newdata not allowed for type \"hazard\"")
        if(is.null(x) | is.character(x))   times<-seq(from=min(fit$data$time),to=max(fit$data$time),length=ntimes) else times<-x
        hazard<-fit$hazard
        hazard$x<-times; hazard$haspar<-haspar(hazard); hazard$hasspline<-hasspline(hazard); hazard$spline.norm<-FALSE
        if(is.null(iter)){
            hazard$spline.par<-fit$posterior.mean$hazard.spline.par
            hazard$param.par<-fit$posterior.mean$hazard.param.par
            hazard$weight<-fit$posterior.mean$hazard.weight
        }else{
            hazard$spline.par<-fit$history$hazard.spline.par[iter,]
            hazard$param.par<-fit$history$hazard.param.par[iter,]
            hazard$weight<-fit$history$hazard.weight[iter,]
        }
        hazard<-makesplinebasis(hazard,quick=TRUE)
        hazard<-evalparametric(hazard)
        hazard<-evalspline(hazard,quick=TRUE)
        haz<-hazard$y
        if(type=="hazard") return(data.frame(time=times,hazard=haz))
    } 
    if(type=="lp" | type=="risk")
    {
        if(is.null(newdata)) {
            vars<-model.matrix(fit$terms,model.frame(fit))[,-1,drop=FALSE]
        }else{
            temp<-cbind(data.frame(i=0,j=0,time=0,delta=0,newdata))
            vars<-model.matrix(fit$terms,model.frame(fit,data=temp))[,-1,drop=FALSE]
        }
        if(is.null(iter))
            coef<-fit$posterior.mean$coefficients
        else
            coef<-fit$history$coefficients[iter,]
        lp<-as.numeric(vars%*%coef)
        if(type=="lp") return(data.frame(lp))
    }
    if(type=="risk")
    {
        pred<-exp(lp)
        if(is.null(newdata)) {
            Ji<-table(fit$data$i)
            if(is.null(iter))
                frailty<-rep(fit$posterior.mean$frailty,Ji)
            else
                frailty<-rep(fit$history$frailty[iter,],Ji)
            pred<-frailty*pred
        }
        risk<-outer(haz,pred,"*")
        colnames(risk)<-if(is.null(newdata)) rownames(fit$data) else rownames(newdata)
        if(dim(risk)[1]==1){
             risk<-as.data.frame(t(risk))
             colnames(risk)<-"risk"
        }else{
            risk<-data.frame(time=times,risk)
        }
        return(risk)
    }
    if(type=="frailty")
    {
        frailty<-fit$frailty
        if(is.null(x)) x<-seq(from=min(frailty$spline.knots),to=max(frailty$spline.knots),length=ntimes)
        frailty$x<-x; frailty$haspar<-haspar(frailty); frailty$hasspline<-hasspline(frailty); frailty$spline.norm<-TRUE
        if(is.null(iter)){
            frailty$spline.par<-fit$posterior.mean$frailty.spline.par
            frailty$param.par<-fit$posterior.mean$frailty.param.par
            frailty$weight<-fit$posterior.mean$frailty.weight
        }else{
            frailty$spline.par<-fit$history$frailty.spline.par[iter,]
            frailty$param.par<-fit$history$frailty.param.par[iter,]
            frailty$weight<-fit$history$frailty.weight[iter,]
        }
        frailty<-makesplinebasis(frailty,quick=TRUE)
        frailty<-evalparametric(frailty)
        frailty<-evalspline(frailty,quick=TRUE)
        density<-frailty$y        
        return(data.frame(x=x,density=density))
    }
}

}}}

{{{ #Main
##############################################################
# \section{Main} MAIN FUNCTION
##############################################################


splinesurv.agdata<-function(x,hazard=NULL,frailty=NULL,regression=NULL,control=NULL,coda=FALSE,initial=NULL,verbose=3,usec=TRUE,...)
{
        
    if(verbose>=1) cat("Initializing...\n")
    
    agdata<-x
    rm(x)
    call<-match.call()
    m<-length(unique(agdata$i))
    Ji<-table(agdata$i)
    
    if(verbose>=2) cat("\tSetting initial parameters...\n")
    
    # Parse input (control)
    control.in<-control
    control.default<-list(
        burnin=500, # Length of the burn-in period
        maxiter=1000, # Max number of iterations
        tun.auto=TRUE, # Auto-calibrate tuning parameters
        tun.int=100 # Interval for calibration of the acceptance rate
    )
    control<-control.default
    controlnames<-names(control)
    innames<-names(control.in)
    if(!is.null(control.in)){
        for(n in innames) eval(parse(text=paste("control$",match.arg(n,controlnames),"<-control.in$",n,sep="")))
    }
    if(control$burnin>control$maxiter) {
        stop("Burnin cannot be greater than maxiter")
    }
     
     # Parse input (frailty)   
    frailty.in<-frailty
    frailty.default<-list(
        type="spline",
        spline.ord=4,
        spline.knotspacing="equal",
        spline.nknots=NULL,
        spline.knots=NULL,
        spline.par=NULL,
        spline.min=-5,
        spline.penalty="2diff",
        spline.penaltyfactor=1,
        spline.meanpenalty=100,
        spline.priorvar=0.1,
        spline.hyper=c(0.01,0.01),
        spline.tun=1,
        spline.accept=0,       
        param.dist="none",
        param.par=NULL,
        param.priorvar=0.1,
        param.hyper=c(0.01,0.01),
        param.tun=1,
        param.accept=0,
        weight=0.5,
        weight.priorvar=0.1,
        weight.hyper=c(1,2),
        weight.tun=0.01,
        weight.accept=0,
        accept=0
    )
    frailty<-frailty.default
    frailtynames<-names(frailty)
    if(!is.null(frailty.in)){
        for(n in names(frailty.in)) eval(parse(text=paste("frailty$",match.arg(n,frailtynames),"<-frailty.in$",n,sep="")))
    }
    frailty$type<-match.arg(frailty$type,c("spline","parametric","both"))
    frailty$hasspline<-frailty$type=="spline" | frailty$type=="both"
    frailty$haspar<-frailty$type=="parametric" | frailty$type=="both"
    frailty$spline.knotspacing<-match.arg(frailty$spline.knotspacing,c("equal"))
    frailty$spline.penalty<-match.arg(frailty$spline.penalty,c("2diff","2deriv","log2deriv","none"))
    frailty$param.dist<-match.arg(frailty$param.dist,c("none","gamma","lognormal"))
    if(m==1 & frailty$haspar) stop("parametric component not allowed for single cluster")
    if(frailty$haspar & frailty$param.dist=="none") {
        warning("no distribution specified for frailty parametric component -- setting to gamma")
        frailty$param.dist<-"gamma"
    }
    frailty$spline.norm<-TRUE
    frailty$name<-"frailty"
    
     # Parse input (frailty)   
    hazard.in<-hazard
    hazard.default<-list(
        type="spline",
        spline.ord=4,
        spline.knotspacing="quantile",
        spline.nknots=NULL,
        spline.knots=NULL,
        spline.par=NULL,
        spline.min=-5,
        spline.penalty="2deriv",
        spline.penaltyfactor=1,
        spline.priorvar=0.1,
        spline.hyper=c(0.01,0.01),
        spline.tun=1,      
        spline.accept=0, 
        param.dist="none",
        param.par=NULL,
        param.priorvar=0.1,
        param.hyper=c(0.01,0.01),
        param.tun=1,
        param.accept=0,
        weight=0.5,
        weight.priorvar=0.1,
        weight.hyper=c(1,2),
        weight.tun=0.01,
        weight.accept=0
    )
    hazard<-hazard.default
    haznames<-names(hazard)
    if(!is.null(hazard.in)){
        for(n in names(hazard.in)) eval(parse(text=paste("hazard$",match.arg(n,haznames),"<-hazard.in$",n,sep="")))
    }
    hazard$type<-match.arg(hazard$type,c("spline","parametric","both"))
    hazard$hasspline<-hazard$type=="spline" | hazard$type=="both"
    hazard$haspar<-hazard$type=="parametric" | hazard$type=="both"
    hazard$spline.knotspacing<-match.arg(hazard$spline.knotspacing,c("quantile","equal","mindiff"))
    hazard$spline.penalty<-match.arg(hazard$spline.penalty,c("2diff","2deriv","log2deriv","none"))
    hazard$param.dist<-match.arg(hazard$param.dist,c("none","exponential","weibull","lognormal")) 
    if(hazard$haspar & hazard$param.dist=="none") {
        warning("no distribution specified for hazard parametric component -- setting to weibull")
        hazard$param.dist<-"weibull"
    }
    if(!hazard$haspar) hazard$weight<-1
    if(!hazard$hasspline) hazard$weight<-0
    if(!frailty$haspar) frailty$weight<-1
    if(!frailty$hasspline) frailty$weight<-0
    
    hazard$spline.norm<-FALSE
    hazard$name<-"hazard"
    hazard$x<-agdata$time
     
    # Parse input (regression)
    reg.in<-regression
    reg.default<-list(
        priorvar=0.1,
        hyper=c(0.01,0.01),
        tun=1,
        accept=0
    )
    regression<-reg.default
    regnames<-names(regression)
    if(!is.null(reg.in)) for(n in names(reg.in)) eval(parse(text=paste("regression$",match.arg(n,regnames),"<-reg.in$",n,sep="")))

     # Automatic number of knots
    if(is.null(hazard$spline.nknots)) hazard$spline.nknots<-min(round(sum(Ji)/4),35)
    if(is.null(frailty$spline.nknots)) frailty$spline.nknots<-min(round(m/4),35)
    
    
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
    frailty$x<-Ui   
    beta<-coxfit$coef
    regression$m<-m
    regression$Ji<-Ji
    regression$covariates<-as.matrix(agdata[,-(1:4)],sum(Ji),length(beta))
    regression$time<-agdata$time
    regression$status<-agdata$delta
    regression$cluster<-as.integer(agdata$i)
    regression<-updateregression(regression,beta)
        
    rm(coxfit)
    
    # Parametric fits
    if(verbose>=2 & (frailty$haspar | hazard$haspar)) cat("\tFitting parametric components...\n")
    
    hazard<-fitparametric(hazard,agdata)
    frailty<-fitparametric(frailty,Ui)

    # Spline knots
    if(verbose>=2 & (frailty$hasspline | hazard$hasspline)) cat("\tComputing spline knots...\n")

    hazard<-makeknots(hazard,agdata$time[agdata$delta==1],bounds=c(min(agdata$time),max(agdata$time)))
    frailty<-makeknots(frailty,Ui,bounds=c(0,2*max(Ui)))
        
    # Evaluate the splines and integrals
    if(verbose>=2 & (frailty$hasspline | hazard$hasspline)) cat("\tConstructing spline basis functions...\n")
    
    hazard<-makesplinebasis(hazard, usec=usec)
    frailty<-makesplinebasis(frailty, usec=usec)
       
    if(verbose>=2 & (frailty$hasspline | hazard$hasspline)) cat("\tInitializing penalty matrices...\n")

    # Penalty matrices
    hazard<-makepenalty(hazard, usec=usec)
    frailty<-makepenalty(frailty, usec=usec)    
        

    if(verbose>=2 & (frailty$hasspline | hazard$hasspline))  cat("\tObtaining initial values for spline parameters...\n")

    {{{# Initial values for the theta vectors
    
    if(hazard$haspar & hazard$hasspline){
        oldhazweight<-hazard$weight
        hazard$weight<-1;
        hazard<-weightcurve(hazard)
    }
    if(hazard$haspar & hazard$hasspline){
        oldfrailweight<-frailty$weight
        frailty$weight<-1;
        frailty<-weightcurve(frailty)
    }

    # Initial values for hazard parameters 
    if(hazard$hasspline){
	    theta.haz<-rep(0,hazard$spline.nknots+hazard$spline.ord)
        if(usec){
            par<-as.double(theta.haz); status<-as.double(regression$status);
            lp<-as.double(regression$lp); frailrep<-as.double(rep(frailty$x,Ji));
            hazParY<-as.double(if(hazard$haspar) hazard$param.y else rep(0,length(lp))); 
            hazParYcum<-as.double(if(hazard$haspar) hazard$param.ycum else rep(0,length(lp))); 
            weight<-hazard$weight; B<-as.double(hazard$spline.basis);
            C<-as.double(hazard$spline.basiscum);
            P<-as.double(hazard$spline.penaltyfactor * hazard$spline.penaltymatrix);
            penaltyType<-as.integer(pmatch(hazard$spline.penalty,c("none","2diff","2deriv","log2deriv"))-1);
            splinemin<-as.double(hazard$spline.min)
            sigma2<-100; sigma2target<-hazard$spline.priorvar
            # compute initial values by slowly decreasing the prior variance
            while(sigma2>sigma2target){
                sigma2<-as.double(sigma2/10)
                opt.theta.haz<-optim(par,fn=cmklik.spline.haz,gr=cmkgr.spline.haz,
                status=status,lp=lp,frailrep=frailrep,hazParY=hazParY,hazParYcum=hazParYcum,weight=weight,B=B,C=C,P=P,penaltyType=penaltyType,sigma2=sigma2,min=splinemin,
                method="BFGS",control=list(fnscale=-1))
                par<-as.double(opt.theta.haz$par)
            }
            opt.theta.haz<-optim(par,fn=cmklik.spline.haz,gr=cmkgr.spline.haz,
            status=status,lp=lp,frailrep=frailrep,hazParY=hazParY,hazParYcum=hazParYcum,weight=weight,B=B,C=C,P=P,penaltyType=penaltyType,sigma2=sigma2,min=splinemin,
            method="BFGS",control=list(fnscale=-1, maxit=1),hessian=TRUE)
            rm(par,status,lp,hazParY,hazParYcum,weight,B,C,P,penaltyType,splinemin,sigma2)
            gcout<-gc()
        }else{
            hazard<-updatespline(hazard,theta.haz)
            gcout<-gc()
            sigma2<-100; sigma2target<-hazard$spline.priorvar
            # compute initial values by slowly decreasing the prior variance
            while(sigma2>sigma2target){
                sigma2<-sigma2/10;
                hazard$spline.priorvar<-sigma2
                opt.theta.haz<-optim(hazard$spline.par,
                        fn=mklik.spline.haz,
                        gr=mkgr.spline.haz,
                        method="BFGS",
                        control=list(fnscale=-1),
                        hazard=hazard,
                        frailty=frailty,
                        regression=regression,
                        hessian=FALSE)
                hazard<-updatespline(hazard,opt.theta.haz$par)
            }
            opt.theta.haz<-optim(hazard$spline.par,
                    fn=mklik.spline.haz,
                    gr=mkgr.spline.haz,
                    method="BFGS",
                    control=list(fnscale=-1),
                    hazard=hazard,
                    frailty=frailty,
                    regression=regression,
                    hessian=TRUE)
        }
        gcout<-gc()
        hazard<-updatespline(hazard,opt.theta.haz$par)

    }

        # Initial values for frailty parameters
    if(frailty$hasspline){
	    frailty$spline.fixedind<-max(1,which.min((frailty$spline.knots-mean(frailty$x))^2)-frailty$spline.ord)
	    theta.frail<-rep(0,frailty$spline.nknots+frailty$spline.ord-1)
        #if(usec){
        #    frailParY<-as.double(if(frailty$haspar) frailty$param.y else rep(0,length(frailty$x))); 
        #    fixedind<-as.integer(frailty$spline.fixedind)
        #    weight<-as.double(frailty$weight)
        #    B<-as.double(frailty$spline.basis);
        #    E<-as.double(frailty$spline.basisexp);
        #    M<-as.double(frailty$spline.meanpenalty);
        #    P<-as.double(frailty$spline.penaltyfactor * frailty$spline.penaltymatrix);
        #    penaltyType<-as.integer(pmatch(frailty$spline.penalty,c("none","2diff","2deriv","log2deriv"))-1);
        #    splinemin<-as.double(frailty$spline.min)
        #    sigma2<-as.double(frailty$spline.priorvar)
        #    opt.theta.frail<-optim(theta.frail,
        #        fn=cmklik.spline.frail,
        #        fixedind=fixedind, frailParY=frailParY,weight=weight,B=B,E=E,M=M,P=P,penaltyType=penaltyType,sigma2=sigma2, min=splinemin,
        #        method="BFGS", control=list(fnscale=-1), hessian=TRUE)
        #}
        #else
        {
            opt.theta.frail<-optim(theta.frail,
                    fn=mklik.spline.frail,
                    method="BFGS",
                    control=list(fnscale=-1),
                    hazard=hazard,
                    frailty=frailty,
                    regression=regression,
                    hessian=TRUE)    
        }
        gcout<-gc()
        frailty<-updatespline(frailty,opt.theta.frail$par)
    }

    gcout<-gc()
    if(hazard$hasspline & hazard$haspar){
        hazard$weight<-oldhazweight
        hazard<-weightcurve(hazard);
    }
    if(frailty$hasspline & frailty$haspar){
        frailty$weight<-oldfrailweight
        frailty<-weightcurve(frailty)
    }

    gcout<-gc() 
    # Evaluate variances and hessians for candidate generation
    frailty$tun<-diff(range(Ui))^2/6
    hess.coef<-mkhess.coef(regression$coefficients,hazard,frailty,regression)
    Sigma.coef<-solve(-hess.coef)
    regression$candcov<-Sigma.coef
    regression$cholcandcov<-chol(Sigma.coef,pivot=TRUE)
    if(hazard$hasspline){
        hess.haz<-opt.theta.haz$hess
        Sigma.haz<-inverthessian(hess.haz)
        hazard$spline.candcov<-Sigma.haz
        hazard$spline.cholcandcov<-chol(Sigma.haz,pivot=TRUE)
        rm(hess.haz,Sigma.haz)
    }
    if(frailty$hasspline){
        hess.frail<-opt.theta.frail$hess
        Sigma.frail<-inverthessian(hess.frail)
        frailty$spline.candcov<-Sigma.frail
        frailty$spline.cholcandcov<-chol(Sigma.frail,pivot=TRUE)
        rm(hess.frail,Sigma.frail)
    }
    if(hazard$haspar){ 
        temphaz<-list(haspar=hazard$haspar,hasspline=hazard$hasspline,weight=hazard$weight,spline.y=hazard$spline.y,spline.ycum=hazard$spline.ycum,name=hazard$name,param.dist=hazard$param.dist,x=hazard$x,y=hazard$y,ycum=hazard$ycum,param.y=hazard$param.y,param.ycum=hazard$param.ycum,param.par=hazard$param.par,param.priorvar=hazard$param.priorvar)
        eps<-1e-5
        par<-temphaz$param.par
        hess<-matrix(0,length(par),length(par))
        for(i in 1:length(par)){
            for(j in 1:length(par)){
            par1<-par;par2<-par;par3<-par;
            if(i==j) {par1[i]<-par[i]+eps;par2<-par1;par3[i]<-par[i]+2*eps}
            else {par1[i]<-par[i]+eps; par2[j]<-par[j]+eps; par3<-par1; par3[j]<-par[j]+eps}
            g1<-(mklik.param.haz(par1,temphaz,frailty,regression)-mklik.param.haz(par,temphaz,frailty,regression))/eps
            g2<-(mklik.param.haz(par3,temphaz,frailty,regression)-mklik.param.haz(par2,temphaz,frailty,regression))/eps
            hess[i,j]<-(g2-g1)/eps
            }
        }
        Sigma.par.haz<-inverthessian(hess)
        hazard$param.candcov<-Sigma.par.haz
        hazard$param.cholcandcov<-chol(Sigma.par.haz,pivot=TRUE)
        rm(Sigma.par.haz,hess,temphaz)
    }
    if(frailty$haspar){
        Sigma.par.frail<-inverthessian(numHess.par(frailty$param.par,mklik.param.frail,hazard=hazard,frailty=frailty,regression=regression))
        frailty$param.candcov<-Sigma.par.frail
        frailty$param.cholcandcov<-chol(Sigma.par.frail,pivot=TRUE)
        rm(Sigma.par.frail)
    }
    
    
    }}}

    #browser()
    gcout<-gc()
    # Store initial values in parameter history
    history<-inithistory(hazard,frailty,regression,control)
    avg.tunhist<-NULL
    avg.accepthist<-NULL
    
    ######################
     main<-function() {}
        
    #browser()
    if(verbose>=1) cat("Starting MCMC...\n")
    
    iter<-1
    
    if(verbose>=3) cat(iter," ")

    
    while(iter<control$maxiter)
    {
        if(usec){
            # C version of the main loop
            gcout<-gc()
            nexttunint<-iter-iter%%control$tun.int+control$tun.int
            enditer <- min(nexttunint, control$maxiter)
            out<-.Call("SplineSurvMainLoop",hazard,frailty,regression,history,iter,enditer,verbose) 
            iter<-enditer
        }else{

            # R version of the main loop

            iter<-iter+1

            if(verbose>=3) cat(iter," ")
                
            # MH update of frailties
            frailty<-mh.frail(hazard,frailty,regression)
            
            # MH update of regression parameters
            regression<-mh.coef(hazard,frailty,regression)

            # MH update of baseline parameters
            hazard<-mh.hazard.spline(hazard,frailty,regression)
            
            # MH update of frailty density parameters
            frailty<-mh.frailty.spline(hazard,frailty,regression)
                    
            # MH update of parametric baseline parameters
            hazard<-mh.hazard.param(hazard,frailty,regression)
            
            # MH update of parametric frailty parameters
            frailty<-mh.frailty.param(hazard,frailty,regression)

            # MH update of weights
            hazard<-mh.weight("hazard",hazard,frailty,regression)
            frailty<-mh.weight("frailty",hazard,frailty,regression)

            # Update of the sigmas / taus
            hazard<-updatepostvar.curve(hazard)
            frailty<-updatepostvar.curve(frailty)
            regression<-updatepostvar.coef(regression)
                    
            {{{# Update parameter history
            history$frailty[iter,]<-frailty$x
            history$coefficients[iter,]<-regression$coefficients
            if(hazard$hasspline) history$hazard.spline.par[iter,]<-hazard$spline.par
            if(frailty$hasspline) history$frailty.spline.par[iter,]<-frailty$spline.par
            if(hazard$haspar) history$hazard.param.par[iter,]<-hazard$param.par
            if(frailty$haspar) history$frailty.param.par[iter,]<-frailty$param.par
            if(hazard$hasspline & hazard$haspar) history$hazard.weight[iter]<-hazard$weight
            if(frailty$hasspline & frailty$haspar) history$frailty.weight[iter]<-frailty$weight
            history$priorvar[iter,]<-c(regression$priorvar,hazard$spline.priorvar,frailty$spline.priorvar,
                hazard$param.priorvar,frailty$param.priorvar,hazard$weight.priorvar,frailty$weight.priorvar)
            history$accept[iter,]<-c(regression$accept,hazard$spline.accept,frailty$spline.accept,
                hazard$param.accept,frailty$param.accept,hazard$weight.accept,frailty$weight.accept,frailty$accept)
            }}}
        }

        {{{  # Periodic calibration check
        if(iter%%control$tun.int==0 & iter<control$maxiter){
  
            if(verbose==1 | verbose==2) cat(iter," ")
            
            calinds<-(iter-control$tun.int+1):iter   

            # Calibration of the tuning parameters for acceptance rate
            if(control$tun.auto & iter<=control$burnin){
               if(verbose>=3) cat("\n Calibration ...\n")
                      
                #browser()
                          
                tunnames<-c("regression$tun",
                    "hazard$spline.tun",
                    "frailty$spline.tun",
                    "hazard$param.tun",
                    "frailty$param.tun",
                    "hazard$weight.tun",
                    "frailty$weight.tun",
                    "frailty$tun")             
                alltun<-rep(0,length(tunnames))
                for(i in 1:length(alltun)) eval(parse(text=paste("alltun[",i,"]<-",tunnames[i])))
                avg.tunhist<-rbind(avg.tunhist,alltun)
                avg.accepthist<-rbind(avg.accepthist,
                                      apply(history$accept[calinds,],2,mean))
                for(g in 1:length(alltun)){
                    if(all(avg.accepthist[,g]>.25)) alltun[g]<-alltun[g]*2
                    if(all(avg.accepthist[,g]<.25)) alltun[g]<-alltun[g]/2
                    if(any(avg.accepthist[,g]>.25) & any(avg.accepthist[,g]<.25)){
                        fit<-lm(y~x,data.frame(x=avg.tunhist[,g],y=avg.accepthist[,g]))
                        out<-try(max(1e-3,nlm(accrate.predict.lm,alltun[g],m=fit)$est))
                        if(inherits(out,"try-error"))
                            alltun[g]<-avg.tunhist[which.min((avg.accepthist-.25)^2)]
                        else
                            alltun[g]<-out
                    }
                }
                for(i in 1:length(alltun)) eval(parse(text=paste(tunnames[i],"<-alltun[",i,"]")))

                if(verbose>=4){
                    outmat<-cbind(avg.accepthist,alltun); rownames(outmat)<-tunnames; colnames(outmat)<-c("acceptance","tuning")
                    print(outmat)
                }
            }
            # Print full iteration info
            if(verbose>=5){
                cat("Frailties: ",submean(history$frailty,calinds),"\n")
                cat("Coefficients: ",submean(history$coefficients,calinds),"\n")
                cat("Hazard.spline: ",submean(history$hazard.spline.par,calinds),"\n")
                cat("Frailty.spline: ",submean(history$frailty.spline.par,calinds),"\n")
                cat("Hazard.param: ",submean(history$hazard.param.par,calinds),"\n")
                cat("Frailty.param: ",submean(history$frailty.param.par,calinds),"\n")
                cat("Prior variances: ",submean(history$priorvar,calinds),"\n")
            }
            
            #browser()
            
        }
        }}}
        
    }
    gcout<-gc()
    if(verbose>0) cat("Done!\n")
    hazard<-makeoutputcurve(hazard)
    frailty<-makeoutputcurve(frailty)
   
    gcout<-gc()
   {{{ # Construct output
    if(control$burnin<iter){
        sub<-(1:(iter))>(control$burnin)
        posterior.mean<-list(coefficients=submean(history$coefficients,sub),
                            frailty=submean(history$frailty,sub),
                            hazard.spline.par=submean(history$hazard.spline.par,sub),
                            hazard.param.par=submean(history$hazard.param.par,sub),
                            hazard.weight=submean(history$hazard.weight,sub),
                            frailty.spline.par=submean(history$frailty.spline.par,sub),
                            frailty.param.par=submean(history$frailty.param.par,sub),
                            frailty.weight=submean(history$frailty.weight,sub)
                        )
        names(posterior.mean$coefficients)<-varnames
    }else{
        postmean=NULL
    }
    if(coda){
        library(coda)
        history$frailty<-mcmc(history$frailty); history$coefficients<-mcmc(history$coefficients)
        history$hazard.spline.par<-mcmc(history$hazard.spline.par); history$frailty.spline.par<-as.mcmc(history$frailty.spline.par)
        history$hazard.param.par<-mcmc(history$hazard.param.par); history$frailty.param.par<-mcmc(history$frailty.param.par)
        history$priorvar<-mcmc(history$priorvar); history$accept<-mcmc(history$accept)
    }
    rownames(history$frailty)<-rownames(history$coefficients)<-rownames(history$splinepar.haz)<-rownames(history$splinepar.frail)<-NULL
    control$iter<-iter
    }}}

    gcout<-gc()
    out<-list(call=call,history=history,posterior.mean=posterior.mean,hazard=hazard,frailty=frailty,control=control,data=agdata)
    class(out)<-"splinesurv"
    return(out)
}
}}}
