

# Compute marginal probability
prob = lapply(Y,function(y){apply(y,1,function(x){mean(x>0)})})


group.colors <- c('1'= "#FF3000",'2'= "#FF9900",'3'= "#FFDB6D",'4'= "#009999",'5'= "#333BFF",'6'= "#9633FF")

df = data.frame(p = unlist(prob))
df$Region = factor(rep(rownames(EC8_new[[1]]), M), levels = rownames(EC8_new[[1]]))
df$Mouse = c(rep('1', R),rep('2', R),rep('3', R),rep('4', R),rep('5', R),rep('6', R))

ggplot(df, mapping = aes(x = Region, y= p, color = Mouse, group = Mouse)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_color_manual(values=group.colors) +
  ylab('P(y>0)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Compute mean projection strength
es = lapply(Y,function(y){apply(y,2,function(x){x/sum(x)})})
es = lapply(es,function(y){apply(y,1,function(x){mean(x)})})


group.colors <- c('1'= "#FF3000",'2'= "#FF9900",'3'= "#FFDB6D",'4'= "#009999",'5'= "#333BFF",'6'= "#9633FF")

df = data.frame(p = unlist(es))
df$Region = factor(rep(rownames(EC8_new[[1]]), M), levels = rownames(EC8_new[[1]]))
df$Mouse = c(rep('1', R),rep('2', R),rep('3', R),rep('4', R),rep('5', R),rep('6', R))

ggplot(df, mapping = aes(x = Region, y= p, color = Mouse, group = Mouse)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_color_manual(values=group.colors) +
  ylab('E(y/n)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Compute conditional prob

jprob = lapply(Y,function(y){matrix(apply(apply(y,2,function(x){(x>0)%*%t(x>0)}),1,mean),R,R)})
cprob = list()
for (m in 1:M){
  cprob[[m]] = jprob[[m]]/matrix(prob[[m]],R,R,byrow = TRUE)
}

#Mouse 1
pp1 = data.frame(cprob[[1]], row.names = rownames(EC8_new[[1]]))
names(pp1) = rownames(EC8_new[[1]])
print(pp1)

pp1 = data.frame(pp1, "Region B" = factor(rownames(EC8_new[[1]]), levels = rownames(EC8_new[[1]])))
names(pp1)[1:R] = rownames(EC8_new[[1]])
pp1 <-  pivot_longer(pp1,
                     cols = !Region.B,
                     names_to = "Region.A", 
                     values_to = "CP"
)
pp1$Region.A = factor(pp1$Region.A, levels = rownames(EC8_new[[1]]))

ggplot(pp1, aes(x = Region.A, y = Region.B, fill = CP)) +
  geom_tile() +
  labs(x = "Region A",
       y = "Region B",
       fill= "P(B|A)") +
  theme_bw() +
  scale_fill_gradient2(high = "red", mid = "yellow", low="white", midpoint = 0.5) +
  geom_text(aes(label = round(CP,3)), color = "black", size = 4) +
  coord_fixed()

