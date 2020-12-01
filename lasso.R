library(pROC)
library(glmnet)

set.seed(0)

group_info  # group of the sample
BRCA_HC_matrix  # count matrix of the significant genes

# Discretization
roc = apply(as.matrix(t(BRCA_HC_matrix)), 2, function(x){roc(group_info,x, quiet = T)$auc})
roc_cutoff = apply(as.matrix(t(BRCA_HC_matrix)), 2, function(x){y = roc(group_info,x, quiet = T) ; 
return(y$thresholds[which.max(y$sensitivities + y$specificities)])})

x1 = as.matrix(t(BRCA_HC_matrix))
z1 = t(as.matrix(t(BRCA_HC_matrix)))
z1_cutoff = z1 >= roc_cutoff
z1[z1_cutoff] = 1
z1[!z1_cutoff] = 0
z1 = t(z1)

# LASSO
set.seed(0)
cv.out1 = cv.glmnet(z1, group_info, alpha =1, nfolds = 5,  family = "binomial", type.measure = "auc")
plot(cv.out1)
lam1 =cv.out1$lambda.1se

# 5-fold cross validation
bs_round = 0
lasso_auc_list = c()

set.seed(0)
while(bs_round < 100){
  bs_round = bs_round + 1
  if(bs_round %% 10 == 0){
    print(bs_round)
  }
  
  sample_random = sample(170, 35)  
  
  x_trn = x1[-sample_random,]
  y_trn = group_info[-sample_random]
  
  roc_trn = apply(x_trn, 2, function(x){roc(y_trn,x, quiet = T)$auc})
  roc_cutoff_trn = apply(x_trn, 2, function(x){y = roc(y_trn,x, quiet = T) ; 
  return(y$thresholds[which.max(y$sensitivities + y$specificities)])})
  
  z_trn = t(x_trn)
  z_trn_cutoff = z_trn >= roc_cutoff_trn
  z_trn[z_trn_cutoff] = 1
  z_trn[!z_trn_cutoff] = 0
  z_trn = t(z_trn)
  
  z_val = t(x1[sample_random,])
  z_val_cutoff = z_val >= roc_cutoff_trn
  z_val[z_val_cutoff] = 1
  z_val[!z_val_cutoff] = 0
  z_val = t(z_val)
  
  bs.out1 = glmnet(z_trn, y_trn, alpha =1, lambda = lam1, family = "binomial")  
  
  bs.pred1 = predict(bs.out1, s=lam1, newx =z_val, type="link")
  
  bsroc_1 <- roc(group_id_bi[sample_random],bs.pred1[,1])
  
  lasso_auc_list = c(lasso_auc_list1, bsroc_1$auc)
  if(bs_round == 1){
    plot.roc(bsroc_1, print.auc=F, add=F, print.thres=F, col = rgb(0, 0, 0, 20, maxColorValue=255))
    }else{
      plot.roc(bsroc_1, print.auc=F, add=T, print.thres=F, col = rgb(0, 0, 0, 20, maxColorValue=255))
    }
}

print(median(lasso_auc_list))


# Validation Set
BRCA_HC_matrix_val  # validation count matrix 
group_info_val  # group of the validation sample

x2 = as.matrix(t(BRCA_HC_matrix_val))
z2 = t(as.matrix(t(BRCA_HC_matrix_val)))
z2_cutoff = z2 >= roc_cutoff
z2[z2_cutoff] = 1
z2[!z2_cutoff] = 0
z2 = t(z2)

bs.out = glmnet(z1, group_info, alpha =1, lambda = lam1,  family = "binomial")   
bs.pred = predict(bs.out, s=lam1, newx =z2, type="link")
bsroc_val <- roc(group_info_val,bs.pred[,1])

plot.roc(bsroc_val, print.auc=T, add=F, print.thres=F, print.auc.x = .35, print.auc.y = 0.1)
