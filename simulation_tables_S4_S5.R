###reproduce simulation results in Tables S4-S5, all methods are converged N times
###author: Zhiqiang CAO
rm(list = ls(all = TRUE))

library(nleqslv)
library(MASS)
library(mvtnorm)

#input program settings
#Set directory to the "R Code" folder
# setwd(".../R Code") # change working directory to the current directory
setwd("C:/Users/user/Dropbox/research/measurment_error/paper6_mem_mis_pa/simulation/Rcode")
source("amle4_mis_setup.R")

###generate simulated data (true x follows standardized gammma)
dataset_gamma = function (n, pt, pi0, pi1) {
  x1 = rgamma(n,shape=1,rate=1); x2=rgamma(n,shape = 2,rate=2)
  x1 = scale(x1); x2 = scale(x2);
  datae = rmvnorm(n = n, mean = c(0, 0), sigma = cove) 
  e1 = datae[, 1]
  e2 = datae[, 2]
  #surrgorate variable for true variable x
  w1 = aa1 + bb1 * x1 + e1
  w2 = aa2 + bb2 * x2 + e2  
  z = rbinom(n, 1, prob = pt)
  b00 = trubet[1]; b10 = trubet[2]
  b20 = trubet[3]; 
  b01 = trubet[4]; b11 = trubet[5]
  b21 = trubet[6];
  y = numeric(n)
  x10 = x1[z==0]; x11 = x1[z==1] 
  x20 = x2[z==0]; x21 = x2[z==1]
  n0 = sum(z==0); n1 = sum(z==1)
  zo = numeric(n)
  zo[z==1] = rbinom(n1, 1, 1-pi1)
  zo[z==0] = rbinom(n0, 1, pi0)
  oddrat0 = exp(b00 + b10 * x10 + b20* x20)
  pr0 = oddrat0 / (1 + oddrat0)
  y[z==0] = rbinom(n0, 1, pr0)
  oddrat1 = exp(b01 + b11 * x11 + b21* x21)
  pr1 = oddrat1 / (1 + oddrat1)
  y[z==1] = rbinom(n1, 1, pr1)
  dim(y) = c(n, 1)
  dim(w1) = c(n, 1)
  dim(w2) = c(n, 1)
  dim(x1) = c(n, 1)
  dim(x2) = c(n, 1)
  dim(z) = c(n, 1)
  dim(zo) = c(n, 1)
  dat = cbind(y, w1, w2, x1, x2, z, zo)
  return(dat)
}

##generate simulated data (true error model is multiplication error model)
dataset_mult = function (n, pt, pi0, pi1) {
  #pt: binary variable
  datax1 = rmvnorm(n = n, mean = c(0, 0), sigma = covt)
  x1 = datax1[, 1]
  x2 = datax1[, 2]
  datae = rmvnorm(n = n, mean = c(0, 0), sigma = cove) 
  e1 = datae[, 1]
  e2 = datae[, 2]
  #multiplication error model
  w1 = x1*exp(e1)
  w2 = x2*exp(e2)  
  z = rbinom(n, 1, prob = pt)
  b00 = trubet[1]; b10 = trubet[2]
  b20 = trubet[3]; 
  b01 = trubet[4]; b11 = trubet[5]
  b21 = trubet[6];
  y = numeric(n)
  x10 = x1[z==0]; x11 = x1[z==1] 
  x20 = x2[z==0]; x21 = x2[z==1]
  n0 = sum(z==0); n1 = sum(z==1)
  zo = numeric(n)
  zo[z==1] = rbinom(n1, 1, 1-pi1)
  zo[z==0] = rbinom(n0, 1, pi0)
  oddrat0 = exp(b00 + b10 * x10 + b20* x20)
  pr0 = oddrat0 / (1 + oddrat0)
  y[z==0] = rbinom(n0, 1, pr0)
  oddrat1 = exp(b01 + b11 * x11 + b21* x21)
  pr1 = oddrat1 / (1 + oddrat1)
  y[z==1] = rbinom(n1, 1, pr1)
  dim(y) = c(n, 1)
  dim(w1) = c(n, 1)
  dim(w2) = c(n, 1)
  dim(x1) = c(n, 1)
  dim(x2) = c(n, 1)
  dim(z) = c(n, 1)
  dim(zo) = c(n, 1)
  dat = cbind(y, w1, w2, x1, x2, z, zo)
  return(dat)
}


#sample size
n = 3000 
#simulation times
N = 1000

#set regression coefficients (choice=1 for Table 1 and choice=2 for Table 2 of manuscript)
choice = 1
if(choice == 1){
  trubet = c(0, 0.3, 0.3, -0.1, 0.2, 0.2) 
}else if(choice == 2){
  trubet = c(0, 0.5, 0.5, -0.1, 0.4, 0.4) 
}

#intercept terms in additive error model
aa1 = aa2 = 0  
#scale-bias terms in additive error model
#In Table S5, bb1 = bb2 = 1
bb1 = bb2 = 0.4  
#variacne-covariance of true variables
covt = matrix(c(1, 0, 0, 1), 2, 2)  
#correlation coefficient between error variables 
rho_uu = 0.75
##variacne-covariance of error variables
cove = matrix(c(0.8, rho_uu*sqrt(0.8*0.85), rho_uu*sqrt(0.8*0.85), 0.85), 2, 2)
#set misclassification rates (three cases) 
case = 3
if(case == 1){
  po = 0.22
  pi0 = 0.1; pi1 = 0.1    #Pr(V=0|V_o=0)=0.981; Pr(V=1|V_o=1)=0.614
}else if(case == 3){
  po = 0.44
  pi0 = 0.2; pi1 = 0.2    #Pr(V=0|V_o=0)=0.857; Pr(V=1|V_o=1)=0.727
}

pt = (po-pi0)/(1-pi0-pi1);
delta = matrix(c((1-pi0)*(1-pt)+pi1*pt,0,0,pi0*(1-pt)+(1-pi1)*pt),2,2)
delta_e = matrix(c((1-pi0)*(1-pt),pi1*pt,pi0*(1-pt),(1-pi1)*pt),byrow = TRUE,2,2)
ccm = solve(delta)%*%delta_e
c11 = ccm[1,1]; c12 = ccm[1,2]
c21 = ccm[2,1]; c22 = ccm[2,2]

#measurement error
sigwx = matrix(c(bb1 * covt[1, 1], bb1 * covt[1, 2], bb2 * covt[2, 1], bb2 *
                   covt[2, 2]), byrow = T, ncol = 2)
sigxw = t(sigwx)
msig_11 = bb1 ^ 2 * covt[1, 1] + cove[1, 1]
msig_12 = bb1 * bb2 * covt[1, 2] + cove[1, 2]
msig_22 = bb2 ^ 2 * covt[2, 2] + cove[2, 2]
sigww = matrix(c(msig_11, msig_12, msig_12, msig_22), 2, 2)
rho2 = solve(sigww)%*%sigwx
tau2 = covt - sigxw %*% solve(sigww) %*% t(sigxw)
rho2f = sigxw %*% solve(sigww)
mt = c(0, 0)
mmu = matrix(rep(mt, n), nrow = 2)
mux = matrix(c(0,0),2,1);
muw = matrix(c(aa1,aa2),2,1)
mu1 = c(1,t(mux)-t(muw)%*%rho2)
mu2 = cbind(matrix(rep(0,2),2,1),rho2)
R0 = rbind(mu1,mu2)
R = kronecker(R0,ccm)
p = 2; taub2 = max(abs(trubet))^2*max(tau2);
qna = qnorm(0.975)


#main function
betae1 = betae2  = betae3 = matrix(0, N, 6)
se1 = se2 = se3 = matrix(0, N, 6)
cova1 = cova2 = cova3 = matrix(0, N, 6)
flag1 = flag2 = flag3 = numeric(N)
#some conditons to exclude extreme estimation results
extrfac = 10 * max(trubet)
seed0 = seq(1, N) + 110007
varfac = 100
i = k = 1;
while (i<=N) {
  cat("iter=",i,"\n")
  set.seed(seed0[i]+k)
  #In Table S4: when choice=1 and case=1; using set.seed(100000+k);
  #             when choice=1 and case=3; using set.seed(seed0[i]+k)
  #             when choice=2;            using set.seed(seed0+k)
  #In Table S5: using set.seed(seed0[i]+k)
  dtaset = dataset_gamma(n, pt, pi0, pi1) #to generate simulated data in Table S4
  #dtaset = dataset_mult(n, pt, pi0, pi1) #to generate simulated data in Table S5
  y = dtaset[, 1]
  w1 = dtaset[, 2]; w2 = dtaset[, 3]
  x1 = dtaset[, 4]; x2 = dtaset[, 5]
  z = dtaset[,6]; zo = dtaset[,7]
  zo0 = 1*(zo==0); zo1 = 1*(zo==1)
  w10 = w1*zo0; w20 = w2*zo0
  w11 = w1*zo1; w21 = w2*zo1
  #naive estimation
  naive = glm(y ~ zo0 + w10 + w20 + zo1 + w11 + w21 - 1, family = binomial(link = 'logit')) #without intercept
  betav = naive$coef
  est_temp1 = betav
  cond1 = max(abs(est_temp1)) <= extrfac
  #RCE
  wm = cbind(zo0,zo1,w1*zo0,w1*zo1,w2*zo0,w2*zo1)
  mwm = wm%*%R
  zm0 = mwm[,1]; zm1 = mwm[,2]
  m10 = mwm[,3]; m20 = mwm[,5]; #when Zo = 0
  m11 = mwm[,4]; m21 = mwm[,6]; #when Zo = 1
  rc = glm(y ~ zm0 + m10 + m20 + zm1 + m11 + m21 - 1, family = binomial(link = 'logit')) #without intercept
  est_temp2 = rc$coef
  cond2 = max(abs(est_temp2)) <= extrfac
  #ALE
  dataw = cbind(w1, w2)
  mw = apply(dataw, 2, mean)
  mmu1 = matrix(rep(mw, n), nrow = 2)
  m = rho2f %*% (t(dataw) - mmu1) + mmu
  m = t(m); m1 = m[, 1]; m2 = m[, 2]
  amle4 = nleqslv(betav, score_amle4_mis)
  flag = amle4$termcd == 1 #converged or not
  est_temp3 = amle4$x
  cond3 = flag & max(abs(est_temp3)) <= extrfac
  cc1 = cond1 & cond2 & cond3
  if(cc1==TRUE){
    #variance of naive 
    sum_naive = summary(naive)
    var_temp1 = (sum_naive$coefficients[,2])^2
    cond1 = max(var_temp1)>varfac
    #variance of RCE
    sum_rc = summary(rc)
    var_temp2 = (sum_rc$coefficients[,2])^2
    cond2 = max(var_temp2)>varfac
    #variance of ALE
    var_temp3 = var_amle4_mis(est_temp3, y, m, zo, tau2, c11, c12, c21, c22, p, n)
    diagm = diag(var_temp3)
    cond3 = is.nan(diagm[1]) | is.nan(diagm[2]) | is.nan(diagm[3]) | is.nan(diagm[4]) |is.nan(diagm[5]) | is.nan(diagm[6]) | any(diagm<0)| max(diagm)>varfac
    cc2 = !cond1 & !cond2 & !cond3 
    if(cc2 == TRUE){
      #naive estimation
      betae1[i, ] = est_temp1
      se1[i,] = sqrt(var_temp1)
      low1 = betae1[i,] - qna * se1[i,]
      high1 = betae1[i,] + qna * se1[i,]
      if (low1[1] <= trubet[1] & trubet[1] <= high1[1]) cova1[i, 1] = 1
      if (low1[2] <= trubet[2] & trubet[2] <= high1[2]) cova1[i, 2] = 1
      if (low1[3] <= trubet[3] & trubet[3] <= high1[3]) cova1[i, 3] = 1
      if (low1[4] <= trubet[4] & trubet[4] <= high1[4]) cova1[i, 4] = 1
      if (low1[5] <= trubet[5] & trubet[5] <= high1[5]) cova1[i, 5] = 1
      if (low1[6] <= trubet[6] & trubet[6] <= high1[6]) cova1[i, 6] = 1
      #RCE
      betae2[i, ] = est_temp2
      se2[i,] = sqrt(var_temp2)
      low2 = betae2[i, ] - qna * se2[i,]
      high2 = betae2[i, ] + qna * se2[i,]
      if (low2[1] <= trubet[1] & trubet[1] <= high2[1]) cova2[i, 1] = 1
      if (low2[2] <= trubet[2] & trubet[2] <= high2[2]) cova2[i, 2] = 1
      if (low2[3] <= trubet[3] & trubet[3] <= high2[3]) cova2[i, 3] = 1
      if (low2[4] <= trubet[4] & trubet[4] <= high2[4]) cova2[i, 4] = 1
      if (low2[5] <= trubet[5] & trubet[5] <= high2[5]) cova2[i, 5] = 1
      if (low2[6] <= trubet[6] & trubet[6] <= high2[6]) cova2[i, 6] = 1
      #ALE
      betae3[i, ] = est_temp3
      se3[i, ] = sqrt(diag(var_temp3))
      low3 = betae3[i, ] - qna * se3[i, ]
      high3 = betae3[i, ] + qna * se3[i, ]
      if (low3[1] <= trubet[1] & trubet[1] <= high3[1]) cova3[i, 1] = 1
      if (low3[2] <= trubet[2] & trubet[2] <= high3[2]) cova3[i, 2] = 1
      if (low3[3] <= trubet[3] & trubet[3] <= high3[3]) cova3[i, 3] = 1
      if (low3[4] <= trubet[4] & trubet[4] <= high3[4]) cova3[i, 4] = 1
      if (low3[5] <= trubet[5] & trubet[5] <= high3[5]) cova3[i, 5] = 1
      if (low3[6] <= trubet[6] & trubet[6] <= high3[6]) cova3[i, 6] = 1
      i = i+1; k = k+1
    }
    else {i = i; k = k+1}
  }
  else {i = i; k = k+1}
}


#naive estimate
ebetae1 = apply(betae1, 2, mean)
bias1 = ebetae1 - trubet
rb1 = (bias1/trubet)*100  #relative bias (%)
ese1 = apply(se1, 2, mean)
ecova1 = apply(cova1, 2, mean)

resu1 = c(bias1[1], rb1[1], ese1[1], ecova1[1], 
          bias1[2], rb1[2], ese1[2], ecova1[2], 
          bias1[3], rb1[3], ese1[3], ecova1[3],
          bias1[4], rb1[4], ese1[4], ecova1[4], 
          bias1[5], rb1[5], ese1[5], ecova1[5], 
          bias1[6], rb1[6], ese1[6], ecova1[6])
#RCE
ebetae2 = apply(betae2, 2, mean)
bias2 = ebetae2 - trubet
rb2 = (bias2/trubet)*100  #relative bias (%)
ese2 = apply(se2, 2, mean)
ecova2 = apply(cova2, 2, mean)
resu2 = c(bias2[1], rb2[1], ese2[1], ecova2[1],
          bias2[2], rb2[2], ese2[2], ecova2[2],
          bias2[3], rb2[3], ese2[3], ecova2[3],
          bias2[4], rb2[4], ese2[4], ecova2[4],
          bias2[5], rb2[5], ese2[5], ecova2[5],
          bias2[6], rb2[6], ese2[6], ecova2[6])
#ALE
ebetae3 = apply(betae3, 2, mean)
bias3 = ebetae3 - trubet
rb3 = (bias3/trubet)*100  #relative bias (%)
ese3 = apply(se3, 2, mean)
ecova3 = apply(cova3, 2, mean)
resu3 = c(bias3[1], rb3[1], ese3[1], ecova3[1],
          bias3[2], rb3[2], ese3[2], ecova3[2],
          bias3[3], rb3[3], ese3[3], ecova3[3],
          bias3[4], rb3[4], ese3[4], ecova3[4],
          bias3[5], rb3[5], ese3[5], ecova3[5],
          bias3[6], rb3[6], ese3[6], ecova3[6])

resut = rbind(resu1, resu2, resu3)
colnames(resut) = c("Bias00", "RB00", "SE00", "CP00", 
                    "Bias10", "RB10", "SE10", "CP10",
                    "Bias20", "RB20", "SE20", "CP20",
                    "Bias01", "RB01", "SE01", "CP01", 
                    "Bias11", "RB11", "SE11", "CP11",
                    "Bias21", "RB21", "SE21", "CP21")
rownames(resut) = c("naive","rce","ale")
res = round(resut,4)
ress = res[,-c(1:4,13:16)]
#keep results
write.csv(ress,paste0("choice_",choice,"_case_",case,"_n_",n,".csv"))
