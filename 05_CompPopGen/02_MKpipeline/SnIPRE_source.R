# G(n) from equation 1.4 in dissertation

Gnint <- function(s,n=30){
int <- 0
  f <- function(x,g){(1-x)^(n-1) * ( 1- exp(-2*g*x))/(2*g*x)}
for(i in 1:length(s)){
  int[i] <- integrate(f,g= s[i], 0,1)$value}
return(int)
}

# F(n) from equation 1.3 in dissertation

Fnint <- function(s, n=30){
int <- 0
  f <- function(x,g){(1-x^(n) - (1-x)^(n))/(1-x)* ( 1- exp(-2*g*x))/(2*g*x)}
for(i in 1:length(s)){
  int[i] <- integrate(f,g= s[i], 0,1)$value}
return(int)
}

# root of function is estimate f (g = estimate of gamma)

funcf <- function(frac, g, r.eff, npop=30, nout=1){
  denom <- Ln(npop)+Ln(nout)
  f <- r.eff - log(frac*(2*g)/(1-exp(-2*g))*(Fnint(g,npop)+Fnint(g, nout))/denom)
  return(f)
}
            
#functions that find roots of funcf       

frac.est <- function(r.eff, g, npop = rep(30, length(r.eff)), 
                     nout = rep(1, length(r.eff))){
  frac <- 0
  for(i in 1:length(r.eff)){
    frac[i] <-  uniroot(funcf, lower = 0, upper = 100, r.eff = r.eff[i],g = g[i],
                                   npop=npop[i], nout=nout[i])$root
  }
  return(frac)
}

# L(n) from equation 1.2 in disseration  

Ln <- function(n){
  v <- 0
  t <- 0
  for(i in 1:length(n)){
  v <- c(1:(n[i]-1))
  t[i] <- sum(1/v)
  if(n[i]==1) t[i] <- 0}
  return(t)
}

# estimate mutation rate

theta <- function(beta, betaG, m, n){
  th <- exp(beta+betaG)/(Ln(m)+Ln(n))
  return(th)
}

# solution to "func" is the selection effect 

func <- function(g, ef, tau, n=30,m =1){
    L.n <- sum(1/(c(1:(n-1))))
  f <- ef - log((tau +Gnint(g,m)+Gnint(g,n))/(Fnint(g, n))*((L.n)/(tau+1/m+1/n)))
return(f)
}

# uses the unit root function to find root of "func"

LLest2gamm <- function(effect,tau.est = rep(10, length(effect)), 
                       n = rep(30, length(effect)), m = rep(1, length(effect))){
gest <- 0
neg <- which(tau.est < 0)
tau.est[neg] <- 0
for(i in 1:length(effect)){
    worked <- try(gest[i] <- uniroot(func, lower=-350, upper=300,
                                     ef=effect[i],tau=tau.est[i], 
                                     n=n[i],m = m[i])$root)
    if (class(worked)=="try-error") gest[i] <- NA
                }
return(gest)
}

# function calls estimation of mutation rate, and estimates tau 

tau.theta <- function(beta, betaG, betaF, betaFG, m=1, n=30){
  th <- theta(beta, betaG, m, n)
  t <- exp(beta+betaG +betaF+betaFG)/th  - 1/m - 1/n
  return(list(tau.est = t, theta.est=th))
}


## Main function call

SnIPRE <-function(mydata){
 data <- mydata
 PS <- data$PS
 PR <- data$PR
 FR <- data$FR
 FS <- data$FS

 n <- length(FS)
 TS <- data$Tsil
 TR <- data$Trepl
 nout = data$nout
 npop = data$npop

 Ivec <- matrix(1, nrow = n)  # makes one vector of appropriate size subset
 d.mu <-as.numeric( matrix(c(PS,PR,FS,FR), ncol =1))
 d.replacement <- as.numeric(c(0,1,0,1)%x%Ivec)
 d.fixed <- as.numeric(c(0,0,1,1)%x%Ivec)
 d.gene <- as.vector(rep(1,4)%x%c(1:n))
 d.TS <- as.vector(rep(1,4)%x%c(TS))
 d.TR <- as.vector(rep(1,4)%x%c(TR))

 count <- as.vector(d.mu)
 R <- as.vector(d.replacement)
 F <- as.vector(d.fixed)
 G <- as.vector(d.gene)
 RF <- as.vector(R*F)
 TR <- as.vector(d.TR)
 TS <- as.vector(d.TS)

 modGEN <- glmer(count~ 1+ R + F + RF +(1+R+F+RF|G),offset =  log(TS*(1-R)+TR*R), family = poisson)
                                        # sel effect
 se.RFG = se.ranef(modGEN)$G[,4]
 re.RFG = ranef(modGEN)$G$RF
 lbound <- fixef(modGEN)[4]+re.RFG - 1.96*se.RFG
 ubound <- fixef(modGEN)[4]+re.RFG + 1.96*se.RFG

 beta <- fixef(modGEN)[1]
 betaG <- ranef(modGEN)$G[,1]
 betaF <- fixef(modGEN)[3]
 betaFG <- ranef(modGEN)$G[,3]

 mydata$SnIPRE.lbound <- lbound
 mydata$SnIPRE.ubound <- ubound

 negC <- which(ubound<=0)
 posC <- which(lbound>=0)

 mydata$SnIPRE.class <- "neut"
 mydata$SnIPRE.class[negC] <- "neg"
 mydata$SnIPRE.class[posC] <- "pos"
 mydata$SnIPRE.est <- fixef(modGEN)[4]+re.RFG
                                        # replacement effect

re.RG = ranef(modGEN)$G$R
se.RG = se.ranef(modGEN)$G[,2]
Rlbound <- fixef(modGEN)[2]+re.RG - 1.96*se.RG
Rubound <- fixef(modGEN)[2]+re.RG + 1.96*se.RG

mydata$SnIPRE.Rlbound <- Rlbound
mydata$SnIPRE.Rubound <- Rubound

negR <- which(Rubound<=0)
posR <- which(Rlbound>=0)
params <- tau.theta(beta, betaG,betaF,betaFG, npop, nout)
mydata$SnIPRE.tau = params$tau.est
mydata$SnIPRE.theta = params$theta.est
mydata$SnIPRE.Rclass <- "neut"
mydata$SnIPRE.Rclass[negR] <- "neg"
mydata$SnIPRE.Rclass[posR] <- "pos"
mydata$SnIPRE.Rest <- fixef(modGEN)[2]+ranef(modGEN)$G[,2]
mydata$SnIPRE.gamma[mydata$SnIPRE.est<=-4.4] <- -Inf 
mydata$SnIPRE.gamma[mydata$SnIPRE.est>-4.4] <-  LLest2gamm(mydata$SnIPRE.est[mydata$SnIPRE.est>-4.4],
                                                        mydata$SnIPRE.tau[mydata$SnIPRE.est>-4.4], n = npop, m= nout)
na.set <- which(is.na(mydata$SnIPRE.gamma) == TRUE)
good.set <- which(mydata$SnIPRE.est > - 4.4)
use <- setdiff(good.set, na.set)
mydata$SnIPRE.f <- NA
mydata$SnIPRE.f.lb <- NA
mydata$SnIPRE.f.ub <- NA
g.ests <- mydata$SnIPRE.gamma
g.zeros <- which(mydata$SnIPRE.class == "neut")
g.ests[g.zeros] <- .000001
mydata$SnIPRE.f[use] <- frac.est(mydata$SnIPRE.Rest[use], g.ests[use],
                                          npop[use], nout[use])
mydata$SnIPRE.f.lb[use] <-frac.est(mydata$SnIPRE.Rlbound[use], g.ests[use], npop[use],nout[use])
mydata$SnIPRE.f.ub[use] <-frac.est(mydata$SnIPRE.Rubound[use], g.ests[use], npop[use],nout[use])
mydata$SnIPRE.f.class <- "neut"
mydata$SnIPRE.f.class[which(mydata$SnIPRE.f.ub <1)] <- "neg"
mydata$SnIPRE.f.class[which(mydata$SnIPRE.f.lb >1)] <- "pos"
return(list(new.dataset = mydata, model = modGEN))
}