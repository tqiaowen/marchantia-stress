TFelnet <- function(y, x, goi_n, resdir, intercept=FALSE, K=5, report=FALSE, cmode) {
  set.seed(2019)
  #do elastic net regression for gene of interest profile
  #find optimal lambda and s by cross validation and return coefficients for that conditions
  #if data is not centralized intercept should be set to TRUE
  source('deps/dircreater.r')
  library(elasticnet)
  library(ggplot2)

  if(report){
    if (!(substr(goi_n, 1,2)=='Mp') & dir.exists(file.path(resdir, goi_n))) {
      print('Warning! filename does already exist, maybe check for doubled gene names')
    }
    dircreater(file.path(resdir, goi_n))
  }
  #center x values rowwise (per condition)
  #x <- t(scale(t(x), center=TRUE, scale=FALSE))
  # set lambda [lambda2 the quadratic punishment]
  lambdas <- c(0,0.001,.01,.05,.1,.5,1,1.5,2,10,100)
  
  if (cmode == "lasso"){
    lambdas <- c(0)
  }
  #print(paste(x, y))
  #print(paste("lambdas defined:", lambdas))
  #initalise matrix to store values for fit with lowest prediction error for each lambda
  cv_res <- matrix(nrow=4, ncol=length(lambdas))
  rownames(cv_res) <- c('lambda', 's', 'cv', 'cv_er')
  cv_res[1,] <- lambdas
  #do the crossvalidation
  for (i in seq(1, length(lambdas))) {
    cv_lm <- cv.enet(x=x, K=K, y=y, s=seq(0.1,1, 0.1), mode='fraction', lambda=lambdas[i], intercept=intercept)
    idx_min <- which.min(cv_lm$cv)
    cv_res[2:4,i] <- c(cv_lm$s[idx_min], cv_lm$cv[idx_min], cv_lm$cv.error[idx_min])
  }
  #print("cv done")
  #print(cv_res)
  #only keep res with <0.91
  #select lambda2 with lowest cv and build model
  if (cmode == "elnet"){
    cv_res <- cv_res[, cv_res['s',]<1 & cv_res['s',] >0.09]
  }
  #(cv_res)
  #plot cv vs lambda
  if(report) {
    p<- ggplot(data.frame(t(cv_res)), aes(x=lambda, y=cv)) + 
      geom_point() +  scale_x_log10() + scale_y_log10() +
      geom_errorbar(aes(ymin=cv-cv_er, ymax=cv+cv_er), width=.2)
    ggsave(file.path(resdir, goi_n, 'cvmin_lamb_s.pdf'))
  }
  #select lambda2 with lowest cv and build model
  if (cmode == "lasso"){
    lambda2 <- 0
  } else if (cmode == "elnet"){
    lambda2 <- cv_res[1, which.min(cv_res['cv',])]
  }
  #print("lambda2 with lowest cv selected")  
  
  lm <- enet(x=x, y=y, lambda=lambda2, intercept=intercept)
  #print("model built")
  if(report) {
  save(lm, file=file.path(resdir, goi_n, 'lm.obj'))
  }
  if(report) {
    pdf(file.path(resdir, goi_n, 'betavs_s.pdf'))
    plot(lm)
    abline(v=cv_res['s', which.min(cv_res['cv',])], col='red', lty='dashed')
    dev.off()
  }
  
  #select s, lambda and cv with lowest cv and build model
  if (cmode == "lasso"){
    mins <- cv_res["s",]
    lambda_used <- cv_res["lambda",]
    cv_used <- cv_res["cv",]
    print(paste("mins:", mins, "lambda_used:", lambda_used, "cv_used:", cv_used))
  } else if (cmode == "elnet"){
    mins <- cv_res['s', which.min(cv_res['cv',])]
    lambda_used <- cv_res['lambda', which.min(cv_res['cv',])]
    cv_used <- cv_res['cv', which.min(cv_res['cv',])]
  }
  
  
  #select s showing lowest cv and generate coefficients (select lambda 1)
  beta_lm <- predict(lm, s=mins, type='coefficients', mode ='fraction')
  #print("beta predicted")
  x_lm <- predict(lm, newx= x, s=mins, type='fit', mode ='fraction')
  #print("lm predicted")
  regulr2 <- function(obs, fit) {
    return(1-(sum((obs-fit)^2))/var(obs))
  }
  relvar <- regulr2(y,x_lm$fit)
  beta_lm_df <- data.frame(predicted = rep(goi_n, length(beta_lm$coefficients)), Gene.ID = names(beta_lm$coefficients), coefficients=beta_lm$coefficients, 
                           relvar=rep(relvar, length(beta_lm$coefficients)), 
                           s=rep(mins, length(beta_lm$coefficients)), lambda=rep(lambda_used, length(beta_lm$coefficients)),
                           cv=rep(cv_used, length(beta_lm$coefficients)))
  beta_lm_df <- beta_lm_df[order(abs(beta_lm_df$coefficients), decreasing=TRUE),]
  #write.table(beta_lm_df, file = file.path(resdir, goi_n, 'coeff.txt'))
  #p2 <- ggplot(beta_lm_df[abs(beta_lm_df$coefficients) >0.05,], aes(x=Gene.Name, y=coefficients)) + geom_bar(stat='identity')+
  #  theme(axis.text.x = element_text(angle = 45))
  #ggsave(file.path(resdir, goi_n, 'coeff.pdf'))
  return(beta_lm_df)
} 