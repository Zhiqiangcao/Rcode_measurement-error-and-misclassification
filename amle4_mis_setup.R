###multivariate amle for logistic regression (keep taylor expansion at 4th order)
###along with misclassification based on model (4)


###Score function of amle4 for logistic regression with misclassification#######################
###Parameter description
##p: number of variables measured with error;
##betav: start points for parameter estimation
##m: conditional mean of x|w; n*p matrix
##tau2: conditional variance of x|w
##z: error free covariates matrix
##zo: observed binary misclassification variable
##c11, c12, c21, c22: estimation of \Delta^{-1}\Delta_e
##y: response variable

score_amle4_mis = function(betav){
  b00 = betav[1]; b01 = betav[p+2] 
  beta0 = matrix(betav[2:(p+1)], ncol=1)
  beta1 = matrix(betav[(p+3):(2+2*p)], ncol=1)
  lam0 = exp(b00 + m%*%beta0) 
  u0 = lam0 / (1 + lam0)
  lam1 = exp(b01 + m%*%beta1) 
  u1 = lam1 / (1 + lam1)
  zo0 = 1*(zo==0); zo1 = 1*(zo==1)
  A0 = exp(y*(b00 + m%*%beta0) - log(1 + lam0))
  A1 = exp(y*(b01 + m%*%beta1) - log(1 + lam1))
  #first derivative
  g1_u0 = y ^ 2 - 2 * u0 * y + 2 * u0 ^ 2 - u0
  g1_u1 = y ^ 2 - 2 * u1 * y + 2 * u1 ^ 2 - u1
  g2_u0 = y ^ 4 - 4 * u0 * y ^ 3 + (12 * u0 ^ 2 - 6 * u0) * y ^ 2 + (-24 * u0 ^ 3 + 24 * u0 ^ 2 - 4 * u0) * y + (24 * u0 ^ 4 - 36 * u0 ^ 3 + 14 * u0 ^ 2 - u0)
  g2_u1 = y ^ 4 - 4 * u1 * y ^ 3 + (12 * u1 ^ 2 - 6 * u1) * y ^ 2 + (-24 * u1 ^ 3 + 24 * u1 ^ 2 - 4 * u1) * y + (24 * u1 ^ 4 - 36 * u1 ^ 3 + 14 * u1 ^ 2 - u1)
  gamma1_u0 = 0.5*as.numeric(t(beta0) %*% tau2 %*% beta0)
  gamma1_u1 = 0.5*as.numeric(t(beta1) %*% tau2 %*% beta1)
  ##the following terms are some constant terms about elements in tau2
  #for zo=0
  s1112 = 4*tau2[1,1]*tau2[1,2]; 
  s2221 = 4*tau2[2,2]*tau2[2,1]; 
  s1212 = 2*(tau2[1,1]*tau2[2,2]+2*tau2[1,2]^2);
  b40t0 = tau2[1,1]^2*beta0[1]^4 + tau2[2,2]^2*beta0[2]^4
  b31t0 = s1112*beta0[2]*beta0[1]^3 + s2221*beta0[1]*beta0[2]^3
  b22t0 = s1212*beta0[1]^2*beta0[2]^2 
  gamma2_u0 = 1/8*(b40t0 + b31t0 + b22t0)
  #for zo=1
  b40t1 = tau2[1,1]^2*beta1[1]^4 + tau2[2,2]^2*beta1[2]^4
  b31t1 = (s1112*beta1[2])*beta1[1]^3 + (s2221*beta1[1])*beta1[2]^3 
  b22t1 = s1212*beta1[1]^2*beta1[2]^2
  gamma2_u1 = 1/8*(b40t1 + b31t1 + b22t1)
  B0 = 1 + g1_u0*gamma1_u0 + g2_u0*gamma2_u0
  B1 = 1 + g1_u1*gamma1_u1 + g2_u1*gamma2_u1
  A01 = c11*A0*B0; A02 = c12*A1*B1; 
  A03 = c21*A0*B0; A04 = c22*A1*B1;
  pg1_u0 = (2*u0^2 - 2*u0)*y + (-4*u0^3 + 5*u0^2 - u0)
  pg1_u1 = (2*u1^2 - 2*u1)*y + (-4*u1^3 + 5*u1^2 - u1)
  pg2_u0 = (4*u0^2 - 4*u0)*y^3 + (-24*u0^3 + 30*u0^2 - 6*u0)*y^2 + (72*u0^4 - 120*u0^3 + 52*u0^2 - 4*u0)*y + 
           (-96*u0^5 + 204*u0^4 - 136*u0^3 + 29*u0^2 - u0)
  pg2_u1 = (4*u1^2 - 4*u1)*y^3 + (-24*u1^3 + 30*u1^2 - 6*u1)*y^2 + (72*u1^4 - 120*u1^3 + 52*u1^2 - 4*u1)*y + 
           (-96*u1^5 + 204*u1^4 - 136*u1^3 + 29*u1^2 - u1)
  cc = diag(p); pgamma1_u0 = pgamma1_u1 = rep(0,p)
  for(i in 1:p){
    ci = matrix(cc[i,],nrow=1) #c1,c2,...,cp
    pgamma1_u0[i] = as.numeric(ci%*%tau2%*%beta0) #P\gamma_1(\beta0), 1*p vector
    pgamma1_u1[i] = as.numeric(ci%*%tau2%*%beta1) #p\gamma_1(\beta1), 1*p vector
  }
  #partial diervative of gamm2_u0 and gamm2_u1
  g31_0 = 4*tau2[1,1]^2*beta0[1]^3 + 3*(s1112*beta0[2])*beta0[1]^2 + 
          (s2221*beta0[2]^3) + 2*beta0[1]*(s1212*beta0[2]^2)
   
  g32_0 = 4*tau2[2,2]^2*beta0[2]^3 + 3*(s2221*beta0[1])*beta0[2]^2 + 
          (s1112*beta0[1]^3) + 2*beta0[2]*(s1212*beta0[1]^2)
  pgamma2_u0 = 1/8*c(g31_0,g32_0) #p \gamma_2(\beta0), 1*p vector
  #for zo=1
  g31_1 = 4*tau2[1,1]^2*beta1[1]^3 + 3*(s1112*beta1[2])*beta1[1]^2 + 
          (s2221*beta1[2]^3) + 2*beta1[1]*(s1212*beta1[2]^2)
  g32_1 = 4*tau2[2,2]^2*beta1[2]^3 + 3*(s2221*beta1[1])*beta1[2]^2 + 
          (s1112*beta1[1]^3) + 2*beta1[2]*(s1212*beta1[1]^2) 
  pgamma2_u1 = 1/8*c(g31_1,g32_1) #p \gamma_2(\beta_1), 1*p vector
  #gamma1_u0 = 0.5*as.numeric(t(beta0) %*% tau2 %*% beta0)
  #gamma1_u1 = 0.5*as.numeric(t(beta1) %*% tau2 %*% beta1)
  #for zo=0
  PA01_1 = A01*(y-u0) + c11*A0*(pg1_u0*gamma1_u0 + pg2_u0*gamma2_u0)
  PA03_1 = A03*(y-u0) + c21*A0*(pg1_u0*gamma1_u0 + pg2_u0*gamma2_u0)
  #for zo=1
  PA02_1 = A02*(y-u1) + c12*A1*(pg1_u1*gamma1_u1 + pg2_u1*gamma2_u1)
  PA04_1 = A04*(y-u1) + c22*A1*(pg1_u1*gamma1_u1 + pg2_u1*gamma2_u1)
  #score for intercept
  s00 = (PA01_1/(A01+A02))*zo0 + (PA03_1/(A03+A04))*zo1
  s01 = (PA02_1/(A01+A02))*zo0 + (PA04_1/(A03+A04))*zo1
  #score for each variable measured with error
  s00p = s01p = matrix(0,n,p)
  for(i in 1:p){
    PA01i = PA01_1*m[,i] + c11*A0*(pgamma1_u0[i]*g1_u0 + pgamma2_u0[i]*g2_u0)
    PA03i = PA03_1*m[,i] + c21*A0*(pgamma1_u0[i]*g1_u0 + pgamma2_u0[i]*g2_u0)
    PA02i = PA02_1*m[,i] + c12*A1*(pgamma1_u1[i]*g1_u1 + pgamma2_u1[i]*g2_u1)
    PA04i = PA04_1*m[,i] + c22*A1*(pgamma1_u1[i]*g1_u1 + pgamma2_u1[i]*g2_u1)
    s00p[,i] = (PA01i/(A01+A02))*zo0 + (PA03i/(A03+A04))*zo1
    s01p[,i] = (PA02i/(A01+A02))*zo0 + (PA04i/(A03+A04))*zo1
  }
  sfc = cbind(s00,s00p,s01,s01p)
  sf = apply(sfc,2,mean)
  return(sf)
}

###Fisher information of amle4 for logistic regression with misclassification#######################
###Parameter description
##p: number of variables measured with error;
##betav: start points for parameter estimation
##m: conditional mean of x|w; n*p matrix
##tau2: conditional variance of x|w
##z: error free covariates matrix
##zo: observed binary misclassification variable
##c11, c12, c21, c22: estimation of \Delta^{-1}\Delta_e
##y: response variable

var_amle4_mis = function(betav, y, m, zo, tau2, c11, c12, c21, c22, p, n){
  b00 = betav[1]; b01 = betav[p+2] 
  beta0 = matrix(betav[2:(p+1)], ncol=1)
  beta1 = matrix(betav[(p+3):(2+2*p)], ncol=1)
  lam0 = exp(b00 + m%*%beta0) 
  u0 = lam0 / (1 + lam0)
  lam1 = exp(b01 + m%*%beta1) 
  u1 = lam1 / (1 + lam1)
  zo0 = 1*(zo==0); zo1 = 1*(zo==1)
  A0 = exp(y*(b00 + m%*%beta0) - log(1 + lam0))
  A1 = exp(y*(b01 + m%*%beta1) - log(1 + lam1))
  #first derivative
  g1_u0 = y ^ 2 - 2 * u0 * y + 2 * u0 ^ 2 - u0
  g1_u1 = y ^ 2 - 2 * u1 * y + 2 * u1 ^ 2 - u1
  g2_u0 = y ^ 4 - 4 * u0 * y ^ 3 + (12 * u0 ^ 2 - 6 * u0) * y ^ 2 + (-24 * u0 ^ 3 + 24 * u0 ^ 2 - 4 * u0) * y + (24 * u0 ^ 4 - 36 * u0 ^ 3 + 14 * u0 ^ 2 - u0)
  g2_u1 = y ^ 4 - 4 * u1 * y ^ 3 + (12 * u1 ^ 2 - 6 * u1) * y ^ 2 + (-24 * u1 ^ 3 + 24 * u1 ^ 2 - 4 * u1) * y + (24 * u1 ^ 4 - 36 * u1 ^ 3 + 14 * u1 ^ 2 - u1)
  gamma1_u0 = 0.5*as.numeric(t(beta0) %*% tau2 %*% beta0)
  gamma1_u1 = 0.5*as.numeric(t(beta1) %*% tau2 %*% beta1)
  ##the following terms are some constant terms about elements in tau2
  #for zo=0
  s1112 = 4*tau2[1,1]*tau2[1,2]; 
  s2221 = 4*tau2[2,2]*tau2[2,1]; 
  s1212 = 2*(tau2[1,1]*tau2[2,2] + 2*tau2[1,2]^2);
  b40t0 = tau2[1,1]^2*beta0[1]^4 + tau2[2,2]^2*beta0[2]^4
  b31t0 = s1112*beta0[2]*beta0[1]^3 + s2221*beta0[1]*beta0[2]^3
  b22t0 = s1212*beta0[1]^2*beta0[2]^2 
  gamma2_u0 = 1/8*(b40t0 + b31t0 + b22t0)
  #for zo=1
  b40t1 = tau2[1,1]^2*beta1[1]^4 + tau2[2,2]^2*beta1[2]^4
  b31t1 = (s1112*beta1[2])*beta1[1]^3 + (s2221*beta1[1])*beta1[2]^3 
  b22t1 = s1212*beta1[1]^2*beta1[2]^2
  gamma2_u1 = 1/8*(b40t1 + b31t1 + b22t1)
  B0 = 1 + g1_u0*gamma1_u0 + g2_u0*gamma2_u0
  B1 = 1 + g1_u1*gamma1_u1 + g2_u1*gamma2_u1
  A01 = c11*A0*B0; A02 = c12*A1*B1; 
  A03 = c21*A0*B0; A04 = c22*A1*B1;
  pg1_u0 = (2 * u0 ^ 2 - 2 * u0) * y + (-4 * u0 ^ 3 + 5 * u0 ^ 2 - u0)
  pg1_u1 = (2 * u1 ^ 2 - 2 * u1) * y + (-4 * u1 ^ 3 + 5 * u1 ^ 2 - u1)
  pg2_u0 = (4 * u0 ^ 2 - 4 * u0) * y ^ 3 + (-24 * u0 ^ 3 + 30 * u0 ^ 2 - 6 * u0) * y ^ 2 + (72 * u0 ^ 4 - 120 * u0 ^ 3 + 52 * u0 ^ 2 - 4 * u0) * y + 
    (-96 * u0 ^ 5 + 204 * u0 ^ 4 - 136 * u0 ^ 3 + 29 * u0 ^ 2 - u0)
  pg2_u1 = (4 * u1 ^ 2 - 4 * u1) * y ^ 3 + (-24 * u1 ^ 3 + 30 * u1 ^ 2 - 6 * u1) * y ^ 2 + (72 * u1 ^ 4 - 120 * u1 ^ 3 + 52 * u1 ^ 2 - 4 * u1) * y + 
    (-96 * u1 ^ 5 + 204 * u1 ^ 4 - 136 * u1 ^ 3 + 29 * u1 ^ 2 - u1)
  cc = diag(p); pgamma1_u0 = pgamma1_u1 = rep(0,p)
  for(i in 1:p){
    ci = matrix(cc[i,],nrow=1) #c1,c2,...,cp
    pgamma1_u0[i] = as.numeric(ci%*%tau2%*%beta0) #p \gamma_1(\beta_0)
    pgamma1_u1[i] = as.numeric(ci%*%tau2%*%beta1) #p \gamma_1(\beta_1)
  }
  #partial diervative of gamm2_u0 and gamm2_u1
  g31_0 = 4*tau2[1,1]^2*beta0[1]^3 + 3*(s1112*beta0[2])*beta0[1]^2 + 
    (s2221*beta0[2]^3) + 2*beta0[1]*(s1212*beta0[2]^2)
  
  g32_0 = 4*tau2[2,2]^2*beta0[2]^3 + 3*(s2221*beta0[1])*beta0[2]^2 + 
    (s1112*beta0[1]^3) + 2*beta0[2]*(s1212*beta0[1]^2)
  pgamma2_u0 = 1/8*c(g31_0,g32_0) #p \gamma_2(\beta0), 1*p vector
  #for zo=1
  g31_1 = 4*tau2[1,1]^2*beta1[1]^3 + 3*(s1112*beta1[2])*beta1[1]^2 + 
    (s2221*beta1[2]^3) + 2*beta1[1]*(s1212*beta1[2]^2)
  g32_1 = 4*tau2[2,2]^2*beta1[2]^3 + 3*(s2221*beta1[1])*beta1[2]^2 + 
    (s1112*beta1[1]^3) + 2*beta1[2]*(s1212*beta1[1]^2) 
  pgamma2_u1 = 1/8*c(g31_1,g32_1) #p \gamma_2(\beta_1), 1*p vector
  #gamma1_u0 = 0.5*as.numeric(t(beta0) %*% tau2 %*% beta0)
  #gamma1_u1 = 0.5*as.numeric(t(beta1) %*% tau2 %*% beta1)
  #for zo=0
  PA01_1 = A01*(y-u0) + c11*A0*(pg1_u0*gamma1_u0 + pg2_u0*gamma2_u0)
  PA03_1 = A03*(y-u0) + c21*A0*(pg1_u0*gamma1_u0 + pg2_u0*gamma2_u0)
  #for zo=1
  PA02_1 = A02*(y-u1) + c12*A1*(pg1_u1*gamma1_u1 + pg2_u1*gamma2_u1)
  PA04_1 = A04*(y-u1) + c22*A1*(pg1_u1*gamma1_u1 + pg2_u1*gamma2_u1)
  #score for each variable measured with error
  PA01_p = PA03_p = matrix(0,n,p)
  PA02_p = PA04_p = matrix(0,n,p)
  for(i in 1:p){
    PA01_p[,i] = PA01_1*m[,i] + c11*A0*(pgamma1_u0[i]*g1_u0 + pgamma2_u0[i]*g2_u0)
    PA03_p[,i] = PA03_1*m[,i] + c21*A0*(pgamma1_u0[i]*g1_u0 + pgamma2_u0[i]*g2_u0)
    PA02_p[,i] = PA02_1*m[,i] + c12*A1*(pgamma1_u1[i]*g1_u1 + pgamma2_u1[i]*g2_u1)
    PA04_p[,i] = PA04_1*m[,i] + c22*A1*(pgamma1_u1[i]*g1_u1 + pgamma2_u1[i]*g2_u1)
  }
  
  ###for the second derivative
  ppg1_u0 = (-4 * u0 ^ 3 + 6 * u0 ^ 2 - 2 * u0) * y + (12 * u0 ^ 4 - 22 * u0 ^ 3 + 11 * u0 ^ 2 - u0)
  ppg2_u0 = (-8 * u0 ^ 3 + 12 * u0 ^ 2 - 4 * u0) * y ^ 3 + (72 * u0 ^ 4 - 132 * u0 ^ 3 + 66 * u0 ^ 2 - 6 * u0) * y ^ 2 + (-288 * u0 ^ 5 + 648 * u0 ^ 4 - 464 * u0 ^ 3 + 108 * u0 ^ 2 - 4 * u0) * y + (480 * u0 ^ 6 - 1296 * u0 ^ 5 + 1224 * u0 ^ 4 - 466 * u0 ^ 3 + 59 * u0 ^ 2 - u0)
  ppg1_u1 = (-4 * u1 ^ 3 + 6 * u1 ^ 2 - 2 * u1) * y + (12 * u1 ^ 4 - 22 * u1 ^ 3 + 11 * u1 ^ 2 - u1)
  ppg2_u1 = (-8 * u1 ^ 3 + 12 * u1 ^ 2 - 4 * u1) * y ^ 3 + (72 * u1 ^ 4 - 132 * u1 ^ 3 + 66 * u1 ^ 2 - 6 * u1) * y ^ 2 + (-288 * u1 ^ 5 + 648 * u1 ^ 4 - 464 * u1 ^ 3 + 108 * u1 ^ 2 - 4 * u1) * y + (480 * u1 ^ 6 - 1296 * u1 ^ 5 + 1224 * u1 ^ 4 - 466 * u1 ^ 3 + 59 * u1 ^ 2 - u1)
  um1um_0 = u0*(u0-1); um1um_1 = u1*(u1-1)
  P1PA01_1 = A01*(y-u0)^2 + A01*um1um_0 + 2*c11*A0*(y-u0)*(pg1_u0*gamma1_u0 + 
                                                             pg2_u0*gamma2_u0) + c11*A0*(ppg1_u0*gamma1_u0 + ppg2_u0*gamma2_u0)
  P1PA03_1 = A03*(y-u0)^2 + A03*um1um_0 + 2*c21*A0*(y-u0)*(pg1_u0*gamma1_u0 + 
                                                             pg2_u0*gamma2_u0) + c21*A0*(ppg1_u0*gamma1_u0 + ppg2_u0*gamma2_u0)
  P1PA02_1 = A02*(y-u1)^2 + A02*um1um_1 + 2*c12*A1*(y-u1)*(pg1_u1*gamma1_u1 + 
                                                             pg2_u1*gamma2_u1) + c12*A1*(ppg1_u1*gamma1_u1 + ppg2_u1*gamma2_u1)
  P1PA04_1 = A04*(y-u1)^2 + A04*um1um_1 + 2*c22*A1*(y-u1)*(pg1_u1*gamma1_u1 + 
                                                             pg2_u1*gamma2_u1) + c22*A1*(ppg1_u1*gamma1_u1 + ppg2_u1*gamma2_u1)
  pt = 2*(1+p)
  fall = matrix(0,pt,pt)  #second derivative matrix
  #second derivative for the intercept term corresponding to zo=0
  Bijp1 = (P1PA01_1/(A01+A02) - (PA01_1/(A01+A02))^2)*zo0
  Bijp3 = (P1PA03_1/(A03+A04) - (PA03_1/(A03+A04))^2)*zo1 
  fall[1,1] = mean(Bijp1+Bijp3)
  #the first row for p covariates measured with error and pz error free covariates
  temp1_0 = (y-u0)*g1_u0 + pg1_u0; temp2_0 = (y-u0)*g2_u0 + pg2_u0
  for(i in 2:(p+1)){
    PPPA01ij = P1PA01_1*m[,i-1] + c11*A0*(temp1_0*pgamma1_u0[i-1] + temp2_0*pgamma2_u0[i-1])
    PPPA03ij = P1PA03_1*m[,i-1] + c21*A0*(temp1_0*pgamma1_u0[i-1] + temp2_0*pgamma2_u0[i-1])
    Bijp1 = (PPPA01ij/(A01+A02) - PA01_1*PA01_p[,i-1]/(A01+A02)^2)*zo0
    Bijp3 = (PPPA03ij/(A03+A04) - PA03_1*PA03_p[,i-1]/(A03+A04)^2)*zo1
    fall[1,i] = mean(Bijp1 + Bijp3)
  }
  fall[1,2+p] = mean(-(PA01_1*PA02_1/(A01+A02)^2)*zo0 - (PA03_1*PA04_1/(A03+A04)^2)*zo1)
  for(i in (3+p):(2+2*p)){
    Bijp1 = -(PA01_1*PA02_p[,i-2-p]/(A01+A02)^2)*zo0
    Bijp3 = -(PA03_1*PA04_p[,i-2-p]/(A03+A04)^2)*zo1
    fall[1,i] =  mean(Bijp1 + Bijp3)
  }
  ###construct derivative of pgamma2_u0 and pgamma2_u1
  g311_0 = 12*tau2[1,1]^2*beta0[1]^2+6*(s1112*beta0[2])*beta0[1]+2*(s1212*beta0[2]^2)
  g312_0 = 3*s1112*beta0[1]^2+3*s2221*beta0[2]^2+4*beta0[1]*s1212*beta0[2]
 
  g321_0 = 3*s2221*beta0[2]^2+3*s1112*beta0[1]^2+4*beta0[2]*s1212*beta0[1]
  g322_0 = 12*tau2[2,2]^2*beta0[2]^2+6*(s2221*beta0[1])*beta0[2]+2*(s1212*beta0[1]^2)
  ppgamma2_0 = 1/8*matrix(c(g311_0,g312_0,g321_0,g322_0), byrow = T,p,p) #pp \gamma_2(\beta_0)
  #for zo=1
  g311_1 = 12*tau2[1,1]^2*beta1[1]^2+6*(s1112*beta1[2])*beta1[1]+2*(s1212*beta1[2]^2)
  g312_1 = 3*s1112*beta1[1]^2+3*s2221*beta1[2]^2+4*beta1[1]*s1212*beta1[2]
  
  g321_1 = 3*s2221*beta1[2]^2+3*s1112*beta1[1]^2+4*beta1[2]*s1212*beta1[1]
  g322_1 = 12*tau2[2,2]^2*beta1[2]^2+6*(s2221*beta1[1])*beta1[2]+2*(s1212*beta1[1]^2)
  ppgamma2_1 = 1/8*matrix(c(g311_1,g312_1,g321_1,g322_1), byrow = T,p,p) #pp \gamma_2(\beta_0)
  ppgamma1_0 = ppgamma1_1 = matrix(0,p,p)
  cc = diag(p)
  for(i in 1:p){
    ci = matrix(cc[i,],nrow=1)
    for(j in 1:p){
      cj = matrix(cc[j,],nrow=1)
      ppgamma1_0[i,j] = as.numeric(ci%*%tau2%*%t(cj))
      ppgamma1_1[i,j] = ppgamma1_0[i,j]
    }
  }
  #pgamma1_u0[i] = as.numeric(ci%*%tau2%*%beta0) #p \gamma_1(\beta_0)
  #pgamma1_u1[i] = as.numeric(ci%*%tau2%*%beta1) #p \gamma_1(\beta_1)
  #temp1_0 = (y-u0)*g1_u0 + pg1_u0; temp2_0 = (y-u0)*g2_u0 + pg2_u0
  #i stands for row, j stands for column
  for(i in 2:(p+1)){
    tem1 = A0*(temp1_0*pgamma1_u0[i-1]+temp2_0*pgamma2_u0[i-1])
    for(j in 1:pt){
      if(j < i) fall[i,j] = fall[j,i]
      else if(j <= (p+1)){
        tem2 = A0*(temp1_0*pgamma1_u0[j-1]+temp2_0*pgamma2_u0[j-1])
        tem3 = A0*(ppgamma1_0[i-1,j-1]*g1_u0+ppgamma2_0[i-1,j-1]*g2_u0)
        PPPA01ij = P1PA01_1*m[,i-1]*m[,j-1]+c11*tem1*m[,j-1]+c11*tem2*m[,i-1]+c11*tem3
        PPPA03ij = P1PA03_1*m[,i-1]*m[,j-1]+c21*tem1*m[,j-1]+c21*tem2*m[,i-1]+c21*tem3
        Bijp1 = (PPPA01ij/(A01+A02) - PA01_p[,i-1]*PA01_p[,j-1]/(A01+A02)^2)*zo0 
        Bijp3 = (PPPA03ij/(A03+A04) - PA03_p[,i-1]*PA03_p[,j-1]/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp1+Bijp3)
      }
      else if(j == (2+p)){
        Bijp1 = - (PA01_p[,i-1]*PA02_1/(A01+A02)^2)*zo0 
        Bijp3 = - (PA03_p[,i-1]*PA04_1/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp1+Bijp3)
      }
      else{
        Bijp1 = - (PA01_p[,i-1]*PA02_p[,j-2-p]/(A01+A02)^2)*zo0 
        Bijp3 = - (PA03_p[,i-1]*PA04_p[,j-2-p]/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp1+Bijp3)
      }
    }
  }
  #for row 2+p
  temp1_1 = (y-u1)*g1_u1 + pg1_u1; temp2_1 = (y-u1)*g2_u1 + pg2_u1
  i = 2+p
  for(j in 1:pt){
    if(j < i) fall[i,j] = fall[j,i]
    else if(j == (2+p)){
      Bijp2 = (P1PA02_1/(A01+A02) - (PA02_1/(A01+A02))^2)*zo0
      Bijp4 = (P1PA04_1/(A03+A04) - (PA04_1/(A03+A04))^2)*zo1 
      fall[i,j] = mean(Bijp2+Bijp4)
    }
    else if(j<=2+2*p){   #c11*A0*(temp1_0*pgamma1_u0[i-1] + (temp2_0)*pgamma2_u0[i-1]
      tem1 = A1*(temp1_1*pgamma1_u1[j-2-p]+temp2_1*pgamma2_u1[j-2-p])
      PPPA02ij = P1PA02_1*m[,j-2-p]+c12*tem1
      PPPA04ij = P1PA04_1*m[,j-2-p]+c22*tem1
      Bijp2 = (PPPA02ij/(A01+A02) - PA02_1*PA02_p[,j-2-p]/(A01+A02)^2)*zo0 
      Bijp4 = (PPPA04ij/(A03+A04) - PA04_1*PA04_p[,j-2-p]/(A03+A04)^2)*zo1
      fall[i,j] = mean(Bijp2+Bijp4)
    }
  }
  #for row (3+p):(2+2*p)
  for(i in (3+p):(2+2*p)){
    tem1 = A1*(temp1_1*pgamma1_u1[i-2-p]+temp2_1*pgamma2_u1[i-2-p])
    for(j in 1:pt){
      if(j < i) fall[i,j] = fall[j,i]
      else{
        tem2 = A1*(temp1_1*pgamma1_u1[j-2-p]+temp2_1*pgamma2_u1[j-2-p])
        tem3 = A1*(ppgamma1_1[i-2-p,j-2-p]*g1_u1+ppgamma2_1[i-2-p,j-2-p]*g2_u1)
        PPPA02ij = P1PA02_1*m[,i-2-p]*m[,j-2-p]+c12*tem1*m[,j-2-p]+c12*tem2*m[,i-2-p]+c12*tem3
        PPPA04ij = P1PA04_1*m[,i-2-p]*m[,j-2-p]+c22*tem1*m[,j-2-p]+c22*tem2*m[,i-2-p]+c22*tem3
        Bijp2 = (PPPA02ij/(A01+A02) - PA02_p[,i-2-p]*PA02_p[,j-2-p]/(A01+A02)^2)*zo0 
        Bijp4 = (PPPA04ij/(A03+A04) - PA04_p[,i-2-p]*PA04_p[,j-2-p]/(A03+A04)^2)*zo1 
        fall[i,j] = mean(Bijp2+Bijp4)
      }
    } 
  }
  f1 = -fall
  var_call = (1/n)*ginv(f1); ##generalized inverse
  return (var_call)
}