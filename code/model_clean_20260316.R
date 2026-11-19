library("rstan")
library(dplyr)
library(tidyr)
#library(tidybayes)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(ggplot2)
library(gridExtra)
#setwd("C:/Users/Ying/OneDrive - UBC (1)/Documents/Research/mortality paper")
setwd("C:/Users/jxlei/OneDrive - UBC/Documents/Research/mortality paper")
load("Rdata0206.RData")

####Extract predicted (underreported) true cases####
#p4 = readRDS("final0721/final0721/p4_829_ge_full_new.rds")
p4 = readRDS("p4_gen_829_0306.rds")
risk4 = summary(p4)$summary
risk4 = risk4[grep("lambda", rownames(risk4))[-1], "50%"]
risk4 = matrix(risk4,nrow=49,ncol=7,byrow = T)
p3 = readRDS("p3_829_ge_full_new.rds")
risk3 = summary(p3)$summary
risk3 = risk3[grep("lambda", rownames(risk3))[-1], "50%"]
risk3 = matrix(risk3,nrow=49,ncol=5,byrow = T)
p2 = readRDS("p2_829_ge_full_new.rds")
risk2 = summary(p2)$summary
risk2 = risk2[grep("lambda", rownames(risk2))[-1], "50%"]
risk2 = matrix(risk2,nrow=49,ncol=5,byrow = T)
p1 = readRDS("p1_66_ge_full_new.rds")
risk1 = summary(p1)$summary
risk1 = risk1[grep("lambda", rownames(risk1))[c(-1,-2)], "50%"]
risk1 = matrix(risk1,nrow=49,ncol=6,byrow = T)
risk=cbind(risk1,risk2,risk3,risk4)
rm(risk1,risk2,risk3,risk4,p1,p2,p3,p4)
########Data preparation#########
data_full = list(
  N = 49,
  T = 21,
  num=num,
  N_edges=107,
  nodes=nodes,
  node1 = nodes$node1,
  node2 = nodes$node2,
  eigen_d2_w=eigen_d2_w,
  cases = as.matrix(cases[,7:27]),
  risk = risk[,3:23],
  lag1_risk = risk[,2:22],
  lag2_risk = risk[,1:21],
  
  pop = pop$pop_20_age65,
  log_pop=log(pop$pop_20_age65),
  death_natural = as.matrix(death_natural[,7:27]),
  dose2_65 = vacc_65[,4:24],#-colMeans(vacc_65[,18:24]),
  lag1_dose2_65 = vacc_65[,3:23],#-colMeans(vacc_65[,17:23]),
  # lag2_dose2_65 = vacc_65[,16:22],#-colMeans(vacc_65[,17:23]),
  # dose2_65_ave=colMeans(vacc_65[,1:24]),
  #dose2_65 = vacc_65[,18:24],
  lpopd = log(pop$density)-mean(log(pop$density)),
  mintry = pop$minrty-mean(pop$minrty),
  #  election = presidental_election$Trump,
  
  age_dist1=age_distri$percent65.74-mean(age_distri$percent65.74),
  age_dist2=age_distri$percent85-mean(age_distri$percent85),
  #dose1_65=dose1_65[,18:24]-colMeans(dose1_65[,18:24]),
  mcc=mcc$X65_or_older-mean(mcc$X65_or_older),
  cardio=cardio$CD_hosp_rate-mean(cardio$CD_hosp_rate),
  
  svi1=svi_index$svi1-mean(svi_index$svi1),#ranking for Socioeconomic Status theme summary
  svi2=svi_index$svi2-mean(svi_index$svi2),#ranking for Household Characteristics theme summary
  svi3=svi_index$svi3-mean(svi_index$svi3),#ranking for Racial and Ethnic Minority Status theme
  svi4=svi_index$svi4-mean(svi_index$svi4),#ranking for Housing Type/ Transportation theme
  svi_tot=svi_index$svi_tot-mean(svi_index$svi_tot)
)

#######Modeling ####
result_main <- stan(file = 'model specifications/mortality_simplified_v11.stan', 
                    data = data_full,
                    chains = 4, iter =10000,seed = 1234,warmup = 5000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
result_sum_main = summary(result_main)$summary
saveRDS(result_main,"result_main_21m.rds")

###########Test covariates##########
#SVI4(2 significant), svi3, in 
#lpopd, cardio, svi1, age_dist, out
result_main <- stan(file = 'model specifications/mortality_simplified_v11.stan', 
                    data = data_full,
                    chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),#thin=10,
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, #b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, #b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, #b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, #b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
#MV with varying coeff
result_main1 <- stan(file = 'model specifications/mortality_simplified_v11.stan', 
                    data = data_full,
                    chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
#16438sec
saveRDS(result_main,"result_main_21m_svi34_timefix1.rds")
result_main=readRDS("result_main_21m_svi34_timefix.rds")

result_main1 <- stan(file = 'model specifications/mortality_simplified_v11.stan', 
                     data = data_full,
                     chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                     control = list(adapt_delta=0.9,max_treedepth=15),
                     init = list(
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=rep(0,21),b_svi4=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=rep(0,21),b_svi4=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=rep(0,21),b_svi4=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=rep(0,21),b_svi4=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

########Sen 1: Full####
result_sen1 <- stan(file = 'model specifications/mortality_simplified_v11_sen1_full_sd_HN05.stan', 
                     data = data_full,
                     chains = 4, iter =10000,seed = 1234,warmup = 5000,save_warmup=F,
                     control = list(adapt_delta=0.9,max_treedepth=15),
                     init = list(
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
#time-invariant sigma
result_sen1 <- stan(file = 'model specifications/mortality_simplified_v11_single_sigma.stan', 
                    data = data_full,
                    chains = 4, iter =15000,seed = 1234,warmup = 10000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),thin = 5,
                    init = list(
                      list(sigma_s=0.1, alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=0.1, alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=0.1, alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=0.1, alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

#time-varying intercept
result_sen1 <- stan(file = 'model specifications/mortality_simplified_v11_sen1_full_sd_HN05_tvb0.stan', 
                    data = data_full,
                    chains = 4, iter =10000,seed = 1234,warmup = 7000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=12),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=rep(-5,21),b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=rep(-5,21),b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=rep(-5,21),b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=rep(-5,21),b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

########Sen 2: Empty####
#The largest R-hat is 1.1; need a re-run
result_sen2 <- stan(file = 'model specifications/mortality_simplified_v11_sen2_empty.stan', 
                     data = data_full,
                     chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                     control = list(adapt_delta=0.9,max_treedepth=15),
                     init = list(
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5,b_infect=rep(0,21),f_infect=0.5,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

########Sen 3: addtive RE####
W_t <- diag(0, 21)
W_t[abs(row(W_t) - col(W_t)) == 1] <- 1
node1_t = node2_t = c()
for (i in 1:21){
  node1_t = cbind(node1_t,t(rep(i,length(which(W_t[i,] == 1)))))
  node2_t = cbind(node2_t,t(which(W_t[i,] == 1)))
}
nodes_t = data.frame(node1_t = t(node1_t), node2_t = t(node2_t))
nodes_t = subset(nodes_t,node1_t<node2_t)
num_t=c(1,rep(2,19),1)

inv_sqrt_diag_t <- diag(1/sqrt(num_t))
eigen_d2_w_t <- eigen(
  inv_sqrt_diag_t %*% W_t %*% inv_sqrt_diag_t,
  symmetric = TRUE,
  only.values = TRUE
)$values

data_full_st = list(
  N = 49,
  T = 21,
  num=num,
  N_edges=107,
  nodes=nodes,
  eigen_d2_w=eigen_d2_w,
  
  num_t=num_t,
  N_edges_t=20,
  nodes_t=nodes_t,
  eigen_d2_w_t=eigen_d2_w_t,
  risk = risk[,3:23],
  lag1_risk = risk[,2:22],

  pop = pop$pop_20_age65,
  log_pop=log(pop$pop_20_age65),
  death_natural = as.matrix(death_natural[,7:27]),
  dose2_65 = vacc_65[,4:24],#-colMeans(vacc_65[,18:24]),
  lag1_dose2_65 = vacc_65[,3:23],#-colMeans(vacc_65[,17:23]),
  # lag2_dose2_65 = vacc_65[,16:22],#-colMeans(vacc_65[,17:23]),
  # dose2_65_ave=colMeans(vacc_65[,1:24]),
  #dose2_65 = vacc_65[,18:24],
  lpopd = log(pop$density)-mean(log(pop$density)),
  mintry = pop$minrty-mean(pop$minrty),
  #  election = presidental_election$Trump,
  
  age_dist1=age_distri$percent65.74-mean(age_distri$percent65.74),
  age_dist2=age_distri$percent85-mean(age_distri$percent85),
  #dose1_65=dose1_65[,18:24]-colMeans(dose1_65[,18:24]),
  mcc=mcc$X65_or_older-mean(mcc$X65_or_older),
  cardio=cardio$CD_hosp_rate-mean(cardio$CD_hosp_rate),
  
  svi1=svi_index$svi1-mean(svi_index$svi1),#ranking for Socioeconomic Status theme summary
  svi2=svi_index$svi2-mean(svi_index$svi2),#ranking for Household Characteristics theme summary
  svi3=svi_index$svi3-mean(svi_index$svi3),#ranking for Racial and Ethnic Minority Status theme
  svi4=svi_index$svi4-mean(svi_index$svi4),#ranking for Housing Type/ Transportation theme
  svi_tot=svi_index$svi_tot-mean(svi_index$svi_tot)
)
result_sen3 <- stan(file = 'model specifications/mortality_simplified_v11_sen3_additive.stan', 
                     data = data_full_st,
                     chains = 4, iter =10000,seed = 1234,warmup = 3000,save_warmup=F,
                     control = list(adapt_delta=0.9,max_treedepth=15),
                     init = list(
                       list(sigma_s=0.01, sigma_t=0.01,alpha_s=0.8,alpha_t=0.8,
                            b0=-5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1)),
                       list(sigma_s=0.01, sigma_t=0.01, alpha_s=0.8,alpha_t=0.8,
                            b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1)),
                       list(sigma_s=0.01, sigma_t=0.01, alpha_s=0.8,alpha_t=0.8,
                            b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1)),
                       list(sigma_s=0.01, sigma_t=0.01, alpha_s=0.8,alpha_t=0.8,
                            b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1))))

result_tv <- stan(file = 'model specifications/mortality_simplified_v11_additive_tvsigma_s.stan', 
                    data = data_full_st,
                    chains = 4, iter =12000,seed = 1234,warmup = 5000,save_warmup=F,
                    control = list(adapt_delta=0.99,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), sigma_t=0.01,alpha_s=0.8,alpha_t=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49), 
                           ar_wn_t = rnorm(21, 0, 0.1)),
                      list(sigma_s=rep(0.01,21), sigma_t=0.01, alpha_s=0.8,alpha_t=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49), 
                           ar_wn_t = rnorm(21, 0, 0.1)),
                      list(sigma_s=rep(0.01,21), sigma_t=0.01, alpha_s=0.8,alpha_t=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49), 
                           ar_wn_t = rnorm(21, 0, 0.1)),
                      list(sigma_s=rep(0.01,21), sigma_t=0.01, alpha_s=0.8,alpha_t=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49), 
                           ar_wn_t = rnorm(21, 0, 0.1))))

########Sen 4: iid RE####
result_sen4 <- stan(file = 'model specifications/mortality_simplified_v11_sen4_iidre.stan', 
                    data = data_full,
                    chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma=0.1, 
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=0.1, 
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=0.1, 
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=0.1, 
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

result_sen4 <- stan(file = 'model specifications/mortality_simplified_v11_sen4_idre_tvsigma.stan', 
                    data = data_full,
                    chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma=rep(0.1,21), 
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=rep(0.1,21), 
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=rep(0.1,21), 
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=rep(0.1,21), 
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

########Sen 5: same month infection####
#The largest R-hat is 1.08; need a re-run
result_sen5 <- stan(file = 'model specifications/mortality_simplified_v11_sen5_same.stan', 
                     data = data_full,
                     chains = 4, iter =10000,seed = 1234,warmup = 5000,save_warmup=F,
                     control = list(adapt_delta=0.9,max_treedepth=15),
                     init = list(
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5,b_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.5,b_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.2,b_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.8,b_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

########Sen 6: lag 2 months infection####
#The largest R-hat is 1.1; need a re-run
result_sen6 <- stan(file = 'model specifications/mortality_simplified_v11_sen6_lag2.stan', 
                    data = data_full,
                    chains = 4, iter =10000,seed = 1234,warmup = 5000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

########Sen 7: observed infection####
#The largest R-hat is 1.09; need a re-run
risk_obs = as.matrix(cases[5:27])/pop$pop_20_age65

data_obs = list(
  N = 49,
  T = 21,
  num=num,
  N_edges=107,
  nodes=nodes,
  eigen_d2_w=eigen_d2_w,
  cases = as.matrix(cases[,7:27]),
  risk = risk_obs[,3:23],
  lag1_risk = risk_obs[,2:22],

  pop = pop$pop_20_age65,
  log_pop=log(pop$pop_20_age65),
  death_natural = as.matrix(death_natural[,7:27]),

  mcc=mcc$X65_or_older-mean(mcc$X65_or_older),

  svi1=svi_index$svi1-mean(svi_index$svi1),#ranking for Socioeconomic Status theme summary
  svi2=svi_index$svi2-mean(svi_index$svi2),#ranking for Household Characteristics theme summary
  svi3=svi_index$svi3-mean(svi_index$svi3),#ranking for Racial and Ethnic Minority Status theme
  svi4=svi_index$svi4-mean(svi_index$svi4),#ranking for Housing Type/ Transportation theme
  svi_tot=svi_index$svi_tot-mean(svi_index$svi_tot)
)

result_sen7 <- stan(file = 'model specifications/mortality_simplified_v11_sen1_full_sd_HN05.stan', 
                    data = data_obs,
                    chains = 4, iter =10000,seed = 1234,warmup = 5000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0, 
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0, 
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
#result_sen7=readRDS("result_main_21m_sen7.rds")

#Sen 8: marginal NB (Poisson-Gamma mixture)####
result_sen8 <- stan(file = 'model specifications/mortality_simplified_v11_sen8_NB.stan', 
                    data = data_full,
                    chains = 4, iter =10000,seed = 1234,warmup = 9000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,xi=500,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,xi=500,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,xi=500,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,xi=500,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

#Interaction termo: infection vs mcc####
result_mcc_infe <- stan(file = 'model specifications/mortality_simplified_v11_mccxinf_tv.stan', 
                     data = data_full,
                     chains = 4, iter =10000,seed = 1234,warmup = 9000,save_warmup=F,
                     control = list(adapt_delta=0.9,max_treedepth=15),
                     init = list(
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                       list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                            b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                            b_mcc_infect=rep(0,21),
                            b_mcc=0, b_svi3=0,b_svi4=0,
                            b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                            ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
#also try only same month or only lag one month

####Looking at age and interactions####
result_age_main <- stan(file = 'model specifications/mortality_simplified_v11_age.stan', 
                        data = data_full,
                        chains = 4, iter =10000,seed = 1234,warmup = 9000,save_warmup=F,
                        control = list(adapt_delta=0.9,max_treedepth=15),
                        init = list(
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5,b_infect=rep(0,21),f_infect=0.5,
                               b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
                               b_agex1=rep(0,21),b_agex2=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                               b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
                               b_agex1=rep(0,21),b_agex2=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                               b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
                               b_agex1=rep(0,21),b_agex2=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                               b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
                               b_agex1=rep(0,21),b_agex2=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
result_main_21m_age0 = readRDS("result_main_21m_age.rds")
result_main_21m_mcc = readRDS("result_main_21m_interaction_mcctv.rds")
# result_age_inte <- stan(file = 'model specifications/mortality_simplified_v11_age.stan', 
#                         data = data_full,
#                         chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
#                         control = list(adapt_delta=0.9,max_treedepth=15),
#                         init = list(
#                           list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
#                                b0=-5,b_infect=rep(0,21),#f_infect=0.5,
#                                b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
#                                ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
#                           list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
#                                b0=-5.5,b_infect=rep(0,21),#f_infect=0.5,
#                                b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
#                                ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
#                           list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
#                                b0=-5.2,b_infect=rep(0,21),#f_infect=0.5,
#                                b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
#                                ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
#                           list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
#                                b0=-5.8,b_infect=rep(0,21),#f_infect=0.5,
#                                b_mcc=0, b_age1=rep(0,21),b_age2=rep(0,21),
#                                ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

#################Model selection process####
result_1 <- stan(file = 'model specifications/mortality_simplified_v11_selection1.stan', 
                        data = data_full,
                        chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                        control = list(adapt_delta=0.9,max_treedepth=15),
                        init = list(
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5,b_infect=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5.5,b_infect=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5.2,b_infect=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                          list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                               b0=-5.8,b_infect=rep(0,21),
                               ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
result_2 <- stan(file = 'model specifications/mortality_simplified_v11_selection2.stan', 
                 data = data_full_st,
                 chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                 control = list(adapt_delta=0.9,max_treedepth=15),
                 init = list(
                   list(sigma_s=0.01, sigma_t=0.01,alpha_s=0.8,alpha_t=0.8,
                        b0=-5,b_infect=rep(0,21),
                        ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1)),
                   list(sigma_s=0.01, sigma_t=0.01,alpha_s=0.8,alpha_t=0.8,
                        b0=-5.5,b_infect=rep(0,21),
                        ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1)),
                   list(sigma_s=0.01, sigma_t=0.01,alpha_s=0.8,alpha_t=0.8,
                        b0=-5.2,b_infect=rep(0,21),
                        ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1)),
                   list(sigma_s=0.01, sigma_t=0.01,alpha_s=0.8,alpha_t=0.8,
                        b0=-5.8,b_infect=rep(0,21),
                        ar_wn_s = rnorm(49, 0, 0.1), ar_wn_t = rnorm(21, 0, 0.1))))

result_3 <- stan(file = 'model specifications/mortality_simplified_v11_selection3.stan', 
                    data = data_full,
                    chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma=0.1, 
                           b0=-5,b_infect=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=0.1, 
                           b0=-5.5,b_infect=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=0.1, 
                           b0=-5.2,b_infect=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma=0.1, 
                           b0=-5.8,b_infect=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

result_4 <- stan(file = 'model specifications/mortality_simplified_v11_selection4.stan', 
                 data = data_full,
                 chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                 control = list(adapt_delta=0.9,max_treedepth=15),
                 init = list(
                   list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                        b0=-5,b_infect=rep(0,21),
                        ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                   list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                        b0=-5.5,b_infect=rep(0,21),
                        ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                   list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                        b0=-5.2,b_infect=rep(0,21),
                        ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                   list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                        b0=-5.8,b_infect=rep(0,21),
                        ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

result_6 <- stan(file = 'model specifications/mortality_simplified_v11_selection6.stan', 
                    data = data_full,
                    chains = 4, iter =7000,seed = 1234,warmup = 6000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect1=0.3,f_infect2=0.3,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

####Check vaccine####
result_vacc <- stan(file = 'model specifications/mortality_simplified_v11_vacc.stan', 
                    data = data_full,
                    chains = 4, iter =40000,seed = 1234,warmup = 30000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),thin=10,
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_vac=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_vac=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_vac=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.8,rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_vac=rep(0,21),
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

####Try Gini index####
data_full = list(
  N = 49,
  T = 21,
  num=num,
  N_edges=107,
  nodes=nodes,
  eigen_d2_w=eigen_d2_w,
  cases = as.matrix(cases[,7:27]),
  risk = risk[,3:23],
  lag1_risk = risk[,2:22],
  lag2_risk = risk[,1:21],
  
  pop = pop$pop_20_age65,
  log_pop=log(pop$pop_20_age65),
  death_natural = as.matrix(death_natural[,7:27]),

  mcc=mcc$X65_or_older-mean(mcc$X65_or_older),

  gini=gini$Gini-mean(gini$Gini),
  earning=(earning$Mean_earnings - mean(earning$Mean_earnings))/10000 #annual earning in 10k
)

result_gini <- stan(file = 'model specifications/mortality_simplified_v11_gini.stan', 
                    data = data_full,
                    chains = 4, iter =10000,seed = 1234,warmup = 9000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
#association was negative, relatively stable at -2
#gini index is positively associated with income, corr=0.52
#adjust both
result_gini_inc <- stan(file = 'model specifications/mortality_simplified_v11_gini.stan', 
                    data = data_full,
                    chains = 4, iter =20000,seed = 1234,warmup = 19000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,b_earn=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,b_earn=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,b_earn=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, b_gini=0,b_earn=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))

####Try Period 4 with degenerative model####
p4 = readRDS("final0721/final0721/p4_829_deg_full_new.rds")
risk4 = summary(p4)$summary
risk4 = risk4[grep("lambda", rownames(risk4))[c(-1)], "50%"]
risk4 = matrix(risk4,nrow=49,ncol=7,byrow = T)
p3 = readRDS("p3_829_ge_full_new.rds")
risk3 = summary(p3)$summary
risk3 = risk3[grep("lambda", rownames(risk3))[-1], "50%"]
risk3 = matrix(risk3,nrow=49,ncol=5,byrow = T)
p2 = readRDS("p2_829_ge_full_new.rds")
risk2 = summary(p2)$summary
risk2 = risk2[grep("lambda", rownames(risk2))[-1], "50%"]
risk2 = matrix(risk2,nrow=49,ncol=5,byrow = T)
p1 = readRDS("p1_66_ge_full_new.rds")
risk1 = summary(p1)$summary
risk1 = risk1[grep("lambda", rownames(risk1))[c(-1,-2)], "50%"]
risk1 = matrix(risk1,nrow=49,ncol=6,byrow = T)
risk=cbind(risk1,risk2,risk3,risk4)
rm(risk1,risk2,risk3,risk4,p1,p2,p3,p4)
data_full = list(
  N = 49,
  T = 21,
  num=num,
  N_edges=107,
  nodes=nodes,
  eigen_d2_w=eigen_d2_w,
  cases = as.matrix(cases[,7:27]),
  risk = risk[,3:23],
  lag1_risk = risk[,2:22],
  lag2_risk = risk[,1:21],
  
  pop = pop$pop_20_age65,
  log_pop=log(pop$pop_20_age65),
  death_natural = as.matrix(death_natural[,7:27]),
  
  mcc=mcc$X65_or_older-mean(mcc$X65_or_older)
)
result_maindeg <- stan(file = 'model specifications/mortality_simplified_v11.stan', 
                    data = data_full,
                    chains = 4, iter =20000,seed = 1234,warmup = 19000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0, 
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), alpha_s=0.9,rho=0.9,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.6,
                           b_mcc=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
saveRDS(result_maindeg,"result_main_21m_degen.rds")

####iCAR####
result_icar <- stan(file = 'model specifications/mortality_simplified_v11_iCAR.stan', 
                    data = data_full,
                    chains = 4, iter =30000,seed = 1234,warmup = 29000,save_warmup=F,
                    control = list(adapt_delta=0.9,max_treedepth=15),
                    init = list(
                      list(sigma_s=rep(0.01,21), rho=0.8,
                           b0=-5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), rho=0.8,
                           b0=-5.5,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), rho=0.8,
                           b0=-5.2,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49)),
                      list(sigma_s=rep(0.01,21), rho=0.8,
                           b0=-5.8,b_infect=rep(0,21),f_infect=0.5,
                           b_mcc=0, b_svi3=0,b_svi4=0,
                           b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,
                           ar_wn_s = matrix(rnorm(49*21, 0, 0.1), nrow=49))))
result_icar = summary(result_icar)$summary
saveRDS(result_icar,"result_main_21m_iCAR.rds")
