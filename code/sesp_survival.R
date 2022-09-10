library(tidyverse)
library(ggplot2)
library(MASS)
library(glmnet)
library(mgcv)
library(gee)
library(gt)
library(gtsummary)
library(geepack)
library(nleqslv)
library(webshot)

library(ggsci)
library(reshape2)
  

## シュミレーション設定----
# 初期化
rm(list = ls(all.names = TRUE))　

# set.seed(12345)
# サンプルサイズ
n = 1000

# シミュレーションの繰り返し数
kk_T <- 50
# kk_T <- 500

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

# 共変量行列の次元
# 100以上でglmはエラーになる
nval <- 1

results.corrected <- data.frame(matrix(NaN,nrow = kk_T, ncol=2))
colnames(results.corrected) <- c("beta1", "beta2")

results.naive <- data.frame(matrix(NaN,nrow = kk_T, ncol=2))
colnames(results.naive) <- c("beta1", "beta2")

results.true <- data.frame(matrix(NaN,nrow = kk_T, ncol=2))
colnames(results.true) <- c("beta1", "beta2")

sigmoid <- function(x){
  return( 1 / (1 + exp(-x)) )
}


for (i in 1:kk_T) {
  setTxtProgressBar(pb,i)
  
  ## 共変量の分布の設定、乱数の発生----
  mu<- rep(0,nval)
  
  # 相関パラメータ
  rho <- 0
  
  Sigma <- diag(nval)
  
  XX <- mvrnorm(n, mu, Sigma)
  
  ## アウトカム生成の設定----
  
  ## 回帰パラメータ
  beta_0 <- 1.2
  
  # ベースラインハザード
  j_num <- 5
  lambda_0 <- seq(0, 4, by=1) *0.05  + 0.1
  
  Y <- numeric()
  for(k in 1:length(lambda_0)){
    # ハザード
    lambda <- exp(XX %*% beta_0) %*% lambda_0[1:k]

    lambda.sum <- rowSums(lambda)

    # 生存関数
    p_survival <- exp(-lambda.sum)
    
    Y_j <- rep(0, n)
    for(j in 1:n){
      Y_j[j] <- 1 - rbinom(1, size=1, prob=p_survival[j])
    }
    
    Y <- cbind(Y, Y_j)
  }
  
  colnames(Y)<-c("j0","j1","j2","j3","j4")
  
  
  ##　誤判別を含むアウトカムデータの生成----
  
  # 感度・特異度に関する設定
  theta0 <- -1
  theta1 <- 2
  theta2 <- -2
  
  SEx <- sigmoid(theta0 + theta1 * 1 + theta2 * XX ) # mean 0.65
  SPx <- 1 - sigmoid(theta0 + theta1 * 0 + theta2 * XX ) # mean 0.648
  
  # Y_starの生成

  Y_star <- numeric()
  for(k in 1:length(lambda_0)){
    # ハザード
    lambda <- exp(XX %*% beta_0) %*% lambda_0[1:k]
    
    lambda.sum <- rowSums(lambda)
    
    # 生存関数
    p_survival <- exp(-lambda.sum) # median 0.904
    
    p_Y_star <- SPx*p_survival + (1-SEx)*(1-p_survival) # median 0.688 
    Y_j <- rep(0, n)
    for(j in 1:n){
      Y_j[j] <- rbinom(1, size=1, prob=1-p_Y_star[j])
    }
    
    Y_star <- cbind(Y_star, Y_j)
  }
  
  colnames(Y_star)<-c("j0","j1","j2","j3","j4")
  
  
  ## パラメータ推定----
  
  # アウトカムの補正
  corrected_outcome_fn <- function(Y){
    return((Y - (1-SPx)) / (SEx + SPx -1))
  }
  
 CO <- numeric()
  for (j in 1:dim(Y_star)[2]) {
    CO_j <- Y_star[,j] %>% corrected_outcome_fn()
    CO <- cbind(CO, CO_j)
  }
 
 # カプランマイヤー推定量の構成
 
 # 引数：Y (n x d)
 # 返り値：surv (d x 1)
 km_fn <- function(Y){
   d <- ncol(Y)
   if(is.null(d)){d <- 1}
   
   lambda_sum <- numeric()
   for (j in 1:d) {
     if(j==1){lambda_sum <- 0}
     else{
     lambda <- (sum(Y[,j] - Y[,j-1])) / (sum(1-Y[,j-1]))
     lambda_sum <- append(lambda_sum, lambda_sum[length(lambda_sum)] + lambda)
     }
   }
   surv <- exp(-lambda_sum)
   
   return(surv)
 }
 
 
 # 擬似値の生成
 
 # 引数：Y (n x d)
 # 返り値：pseudo_obs(n x d)
 pseudo_obs_fn <- function(Y){
   n <- nrow(Y)
   if(is.null(n)){n<-length(Y)}
   d <- ncol(Y)
   if(is.null(d)){d<-1}
   
   pseudo_obs <- numeric()
   for (j in 1:n) {
     pso_j <- n*km_fn(Y) - (n-1)*km_fn(Y[-j,])
     pseudo_obs <- rbind(pseudo_obs, pso_j)
   }
  
  rownames(pseudo_obs) <- 1:n
  return(pseudo_obs)
 }
 
 ## 擬似値のプロット
 
 pseudo_obs <- pseudo_obs_fn(CO)
 x <- pseudo_obs %>% head(4)

 colnames(x) <- 1:5
 rownames(x) <- c("id1", "id2", "id3", "id4")
 y <- melt(x)
 colnames(y) <- c("id","time", "PO")

 g <- ggplot(y, aes(x = time, y = PO, color = id))
 g <- g + geom_line()
 g <- g + scale_color_nejm()
 g <- g + labs(title="Pseudo-obs(corrected Y)")
 plot(g)
 # ggsave(file = "./results/plot/0910/po_plot_corrected.png", plot = g)
 
 
  
 # パラメータ推定----
 
 outcome.naive <- pseudo_obs_fn(Y_star)[,5]
 outcome.corrected <- pseudo_obs_fn(CO)[,5]
 outcome.true <- pseudo_obs_fn(Y)[,5]
  
  # 推定方程式
  # ロジットリンク
  fn <- function(X, Y, beta){
    t(X) %*% (sigmoid(X %*% beta)*(1-sigmoid(X %*% beta))*(Y - sigmoid(X %*% beta)))
  }
  
  XX_1 <- cbind(rep(1,n), XX)
  
  fn1.naive <- function(beta) fn(XX_1, outcome.naive, beta)
  ans.naive <- nleqslv(c(0.1, 0.1), fn1.naive)
  
  fn1.corrected <- function(beta) fn(XX_1, outcome.corrected, beta)
  ans.corrected <- nleqslv(c(0.1, 0.1), fn1.corrected)
  
  fn1.true <- function(beta) fn(XX_1, outcome.true, beta)
  ans.true <- nleqslv(c(0.1, 0.1), fn1.true)
  
  
  ### 結果の保存
  results.corrected[i,] <- ans.corrected$x
  results.naive[i,] <- ans.naive$x
  results.true[i,] <- ans.true$x
}


##　ボックスプロット----
# results.corr %>% boxplot()
# results.bias1 %>% boxplot()
# results.bias10 %>% boxplot()
# results.bias100 %>% boxplot()
# results.bias1000 %>% boxplot()

## result_dfの作成----
# result_df <- data.frame(
cor_df <- data.frame(results.corrected, method='corrected')
naive_df <- data.frame(results.naive, method='naive')
true_df <- data.frame(results.true, method='true')
temp_df<- rbind(cor_df,naive_df) 
temp_df <- rbind(temp_df, true_df)

# min_df <- temp_df %>% filter(miss_type=='min') %>% select(-miss_type)

theme_gtsummary_journal(journal = "jama")
# theme_gtsummary_mean_sd()
result_df <- temp_df %>% tbl_summary(by=method, digits = everything()~2) %>% 
  modify_header(label="Parameter") %>% 
  modify_caption("Logit Link")

## テーブル出力----
# result_df %>% as_gt() %>% gtsave("./results/table/0910/res_table_med.png")


## アウトプット
# export_data <- results.corr

# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/ite_sesp/0709.csv")

# 
# df <- read_csv2("~/Projects/gamma_DR_ver0.1/results/0423/g0_m01.csv")
# df <- read_csv2("~/gamma_DR_ver0.1/results/tian/misbin_nval50.csv")
# df %>% summary()
# df %>% boxplot()

