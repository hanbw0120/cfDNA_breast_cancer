cfDNA_up = c()
cfDNA_down = c()
RNA_up = c()
RNA_down = c()

c1 = table(cfDNA_up %in% RNA_down)["TRUE"]/length(unique(cfDNA_up, RNA_down))
c2 = table(cfDNA_up %in% RNA_up)["TRUE"]/length(unique(cfDNA_up, RNA_up))
c3 = table(cfDNA_down %in% RNA_up)["TRUE"]/length(unique(cfDNA_down, RNA_up))
c4 = table(cfDNA_down %in% RNA_down)["TRUE"]/length(unique(cfDNA_down, RNA_down))

# 1000 times
set.seed(1)
permutaion_times = 1000
matrix_out = matrix(0,permutaion_times,4)
colnames(matrix_out) = c("c1", "c2", "c3", "c4")
rownames(matrix_out) = 1:permutaion_times


for(i in 1:permutaion_times){
  cf_all = c(cfDNA_up, cfDNA_down)  

  random_index = sample(length(cf_all), length(cfDNA_up))
  cf_up_re = cf_all[random_index]
  cf_down_re = cf_all[-random_index]
  
  c1re = table(cf_up_re %in% RNA_down)["TRUE"]/length(unique(cf_up_re, RNA_down))
  c2re = table(cf_up_re %in% RNA_up)["TRUE"]/length(unique(cf_up_re, RNA_up))
  
  random_index = sample(length(cf_all), length(cfDNA_up))
  cf_up_re = cf_all[random_index]
  cf_down_re = cf_all[-random_index]
  
  c3re = table(cf_down_re %in% RNA_up)["TRUE"]/length(unique(cf_up_re, RNA_up))
  c4re = table(cf_down_re %in% RNA_down)["TRUE"]/length(unique(cf_up_re, RNA_down))
  
  matrix_out[i,] = c(c1re, c2re, c3re, c4re)
}

cs_list = list()
cj_list = list()
cjp_list = list()
for(i in 1:4){
  ci = c(c1, c2, c3, c4)[i]
  cs = scale(matrix_out[,i])
  
  ci_scale = attr(cs,"scaled:scale")
  ci_center = attr(cs,"scaled:center")
  cj = (ci - ci_center) / ci_scale
  
  cs_list[[i]] = cs
  cj_list[[i]] = cj
  cjp_list[[i]] = (1-pnorm( abs(cj)) ) * 2
  
}

# Z-score
print(cs_list)

# Distribution
print(cj_list)

# p-value
print(cjp_list)
