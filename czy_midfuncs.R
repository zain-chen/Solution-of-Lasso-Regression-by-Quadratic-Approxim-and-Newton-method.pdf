GCV.lambda <- function(x, y, lambdas){
  # 此方法用GCV方法选出最优的lambdas（只能选出一个）
  # 传入的lambda为一个范围向量
  p <- dim(x)[2]
  GCV <- vector(length = length(lambdas))
  q = 1
  for(lambda in lambdas){
    beta.gcv <- LASSO.MM(beta0, x, y, lambda, d, thd1, thd2)               # 对输入的lambda拟合lasso
    P <- x %*% solve(t(x)%*%x + n*diag(as.vector(lambda/abs(beta.gcv)))) %*% t(x)   # 计算帽子矩阵
    d <- sum(diag(P))                                                               # 计算帽子矩阵的迹
    e2 <- sum((y-x%*%beta.gcv)^2)                                                   # 计算lasso估计的RSS
    gcv <- e2 / (1-d/n)^2                                                           # 计算GCV
    GCV[q] <- gcv
    q = q + 1
  }
  bestlambda <- lambdas[which.min(GCV)]                                             # 找出使GCV最小的lambda
  return(bestlambda)
}

LASSO.LQA <- function(beta0, x, y, lambda, thd1, M, thd2){
  # 用LQA算法计算lasso的数值解
  # beta0为初始系数向量（p维），n为样本量，p为变量维数，lambda为调节参数的待选范围，thd1为将β分量压缩为0的阈值，thd2为停止迭代的阈值
  # 在迭代过程中将β小于阈值的分量压缩为零
  n <- dim(x)[1]
  p <- dim(x)[2]
  sigma.func <- function(beta){                   # Σ矩阵
    beta.copy = beta
    for(i in 1:p){
      if(beta.copy[i] == 0) beta.copy[i] <- M     # 为使矩阵逆存在，将1/0用一个很大的数M代替
      else beta.copy[i] <- 1/abs(beta.copy[i])
    }
    sigma <- diag(as.vector(lambda*beta.copy))    # 对角线为lambda_i / |beta|
    return(sigma)
  }
  
  beta <- beta0                                   # 初值为最小二乘估计
  k = 0
  repeat{
    sigma <- sigma.func(beta)                        
    newbeta <- solve( t(x)%*%x + n*sigma ) %*% t(x) %*% y  # 计算beta_k+1
    for(beta.i in newbeta){                                # 当其分量足够小时将其压缩为0
      if(abs(beta.i) < thd1) beta.i <- 0
    }
    if(t(newbeta-beta)%*%(newbeta-beta) < thd2) break      # 当两次迭代的向量差距(差的二范数)足够小时停止迭代
    beta <- newbeta
    k = k + 1
    #cat("\r", k)                                           # 输出迭代次数
  }
  return(beta) 
}

LASSO.LQA2 <- function(beta0, x, y, lambda, thd1, M, thd2){
  # 此函数在迭代过程中不将β小于阈值的分量压缩为零，而在收敛后再将小于阈值的分量压缩为零
  # 其余部分与LASSO.LQA相同
  n <- dim(x)[1]
  p <- dim(x)[2]
  sigma.func <- function(beta){
    sigma <- diag(as.vector(lambda/abs(beta)))
    return(sigma)
  }
  beta <- beta0
  k = 0
  repeat{
    sigma <- sigma.func(beta)                        
    newbeta <- solve( t(x)%*%x + n*sigma ) %*% t(x) %*% y
    if(t(newbeta-beta)%*%(newbeta-beta) < thd2) break          
    beta <- newbeta
    k = k + 1
    #cat("\r", k)
  }
  for(beta.i in newbeta){                                  # 收敛后判断：当其分量足够小时将其压缩为0，再返回lasso估计
    if(abs(beta.i) < thd1) beta.i <- 0
  }
  return(beta) 
}

LASSO.MM <- function(beta0, x, y, lambda, d, thd1, thd2){
  # 此函数在迭代过程中将β小于阈值的分量压缩为零，但不将1/|β_k|替换为M，而是在分母上加上一个较小的数d，使得矩阵的对角线元素不会太大
  # 稳定性更好
  n <- dim(x)[1]
  p <- dim(x)[2]
  sigma.func <- function(beta){                   
    beta.copy = beta
    for(i in 1:p){    
      beta.copy[i] <- 1/(abs(beta.copy[i]) + d)           # 在分母上加上一个较小数d以保证矩阵可逆
    }
    sigma <- diag(as.vector(lambda*beta.copy))
    return(sigma)
  }
  beta <- beta0
  k = 0
  repeat{
    sigma <- sigma.func(beta)                        
    newbeta <- solve( t(x)%*%x + n*sigma ) %*% t(x) %*% y  
    for(beta.i in newbeta){                                
      if(abs(beta.i) < thd1) beta.i <- 0                           # 当其分量足够小时将其压缩为0
    }
    if(t(newbeta-beta)%*%(newbeta-beta) < thd2) break      
    beta <- newbeta
    k = k + 1
    #cat("\r", k)                                           
  }
  return(beta) 
}

Adaptive.Lasso <- function(beta0, x, y, lambda0, d, thd1, thd2){
  # 用MM算法计算adaptive lasso的数值解，因为相对LQA更稳定
  # 这里用每次迭代的1/|β_k|作为计算β_k+1时的权重而非一直用LS，是考虑到LS可能估计不好，而每次迭代的估计会愈来愈精确
  # 指数姑且取1
  n <- dim(x)[1]
  p <- dim(x)[2]
  sigma.func <- function(beta){                   
    lambda <- lambda0/abs(beta)                        #用每次迭代的1/|β_k|代替固定的一个lambda值
    beta.copy = beta
    for(i in 1:p){    
      beta.copy[i] <- 1/(abs(beta.copy[i]) + d)           
    }
    sigma <- diag(as.vector(lambda*beta.copy))
    return(sigma)
  }
  beta <- beta0
  k = 0
  for(l in 1:70){                                       # 因为adaptive lasso在实际中更难收敛，故设定迭代上限70次
    sigma <- sigma.func(beta)                        
    newbeta <- solve( t(x)%*%x + n*sigma ) %*% t(x) %*% y  
    for(beta.i in newbeta){                                
      if(abs(beta.i) < thd1) beta.i <- 0
    }
    if(t(newbeta-beta)%*%(newbeta-beta) < thd2) break      
    beta <- newbeta
    k = k + 1
    #cat("\r", k)                                           
  }
  return(beta)
}

# 最优子集选择：使用AIC与BIC
Best.Subset <- function(x, y, criterion){
  # 用最优子集选择进行变量选择
  # criterion是可选判断准则，此函数提供AIC和BIC
  # x[,i]代表列值，即Xi
  n <- dim(x)[1]
  p <- dim(x)[2]
  models <- list()
  for(i in 1:p){
    # 拟合所有包含i个变量的模型，用RSS最小准则选出最优的model_i
    choice <- combn(p, i)                         # 所有组合
    RSS1 <- vector(length = choose(p,i))          # 存放所有组合拟合结果
    for(j in choose(p,i)){
      x1 <- x[,choice[,j]]
      e1 <- y - x1 %*% solve(t(x1)%*%x1) %*% t(x1) %*% y
      rss1 <- t(e1) %*% e1
      RSS1[j] <- rss1                             # 计算所有组合的rss并存放（同时得到了全模型RSS，后面会用到）
    }
    bestchoice <- choice[,which.min(RSS1)]        # 按rss最小选出最优model_i
    models[[i]] <- bestchoice[1:i]
  }
  # 计算p个最优model_i的RSS和β方差估计，为计算AIC、BIC做准备
  RSS2 <- vector(length = p)
  D2 <- 1:p
  q = 1
  for(model in models){
    x2 <- x[,model]
    d = length(model)
    e2 <- y - x2 %*% solve(t(x2)%*%x2) %*% t(x2) %*% y    # 残差
    rss2 <- t(e2) %*% e2                                  # 部分变量模型的RSS
    RSS2[q] <- rss2
    q = q + 1
  }
  SIG2 <- rep(RSS2[p]/(n-p), p)
  if(criterion=='AIC'){
    AIC <- (RSS2 + 2*D2*SIG2)/(n*SIG2)
    bestmodel <- models[[which.min(AIC)]]
  }
  
  if(criterion=='BIC'){
    BIC <- (RSS2 + log(n)*D2*SIG2)/n
    bestmodel <- models[[which.min(BIC)]]
  }
  
  return(bestmodel)
}

# 模拟函数
simu <- function(e.sd, x.sigma, lambdas, d)
{
  # 只设计了epsilon的标准差，X的协方差矩阵，选择最优lambda的范围和MM估计中的参数d作为可改变变量
  # 因为接下来的不同条件的实验只会改变这些参数
  # 其他参数指定为全局变量
  correct.lasso.lqa <- c()
  incorrect.lasso.lqa <- c()
  correct.lasso.mm <- c()
  incorrect.lasso.mm <- c()
  correct.bic <- rep(0, rep.num)
  incorrect.bic <- rep(0, rep.num)
  correct.aic <- rep(0, rep.num)
  incorrect.aic <- rep(0, rep.num)
  correct.adaptivelasso <- c()
  incorrect.adaptivelasso <- c()
  for(i in 1:rep.num){
    # 更新
    e <- rnorm(n, 0, e.sd)
    x <- rmvnorm(n, rep(0,p), x.sigma)
    colnames(x) <- c("x1","x2","x3","x4","x5","x6")
    y <- x %*% realbeta + e
    beta0 <- solve(t(x)%*%x)%*%t(x)%*%y
    opt.lam <- GCV.lambda(x,y,lambdas)
    # mm
    correct.lasso.mm[i] <- sum(abs(LASSO.MM(beta0, x, y, opt.lam, d, thd1, thd2))[4:6]<0.005)
    incorrect.lasso.mm[i] <- sum(abs(LASSO.MM(beta0, x, y, opt.lam, d, thd1, thd2))[1:3]<0.005)
    # adaptive lasso
    correct.adaptivelasso[i] <- sum(abs(Adaptive.Lasso(beta0, x, y, opt.lam, d, thd1, 0.1))[4:6]<0.005)
    incorrect.adaptivelasso[i] <- sum(abs(Adaptive.Lasso(beta0, x, y, opt.lam, d, thd1, 0.1))[1:3]<0.005)
    # aic bic
    for(j in 1:3){
      if(!(j %in% Best.Subset(x,y,criterion = 'AIC'))) incorrect.aic[i] = incorrect.aic[i] + 1
    }
    for(j in 4:6){
      if(!(j %in% Best.Subset(x,y,criterion = 'AIC'))) correct.aic[i] = correct.aic[i] + 1
    }
    for(j in 1:3){
      if(!(j %in% Best.Subset(x,y,criterion = 'BIC'))) incorrect.bic[i] = incorrect.bic[i] + 1
    }
    for(j in 4:6){
      if(!(j %in% Best.Subset(x,y,criterion = 'BIC'))) correct.bic[i] = correct.bic[i] + 1
    }
  }
  # 把结果存放在dataframe里做成表格展示
  result.co <- colMeans(data.frame(correct.lasso.mm, correct.aic, correct.bic, correct.adaptivelasso))
  result.in <- colMeans(data.frame(incorrect.lasso.mm, incorrect.aic, incorrect.bic, incorrect.adaptivelasso))
  result <- rbind.data.frame(result.co,result.in)
  colnames(result) <- c("lasso","aic","bic","adaptive")
  rownames(result) <- c("correct", "incorrect")
  return(result)
}
