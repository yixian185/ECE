#SOM

dat <- read.csv('f00.csv', row.names = 1)

dat<-t(dat)
library(kohonen)
som_grid <- somgrid(xdim = 10, ydim = 10, topo = 'hexagonal')

set.seed(100)
som_model <- supersom(dat, grid = som_grid, rlen = 10000)

plot(som_model, type = 'changes')
plot(som_model, type = 'count')
plot(som_model, type = 'dist.neighbours')
plot(som_model, type = 'codes')

coolBlueHotRed <- function(n, alpha = 0.7) rainbow(n, end=4/6, alpha=alpha)[n:1]

color_by <- apply(som_model$data[[1]], 1, mean)
unit_colors <- aggregate(color_by, by = list(som_model$unit.classif), FUN = mean, simplify = TRUE)
unit_dat <- data.frame(value = rep(0, 100))
unit_dat[unit_colors$Group.1,'value'] <- unit_colors$x

plot(som_model, type = 'property', property = unit_dat[[1]], palette.name = coolBlueHotRed, 
     shape = 'round', keepMargins = TRUE, border = NA)

plot(som_model, type = 'property', property = unit_dat[[1]], palette.name = coolBlueHotRed, 
     shape = 'straight', keepMargins = TRUE, border = NA)

som_model$grid

som_model_class <- data.frame(gene_name = rownames(som_model$data[[1]]), code_class = som_model$unit.classif)
som_model_class <- cbind(som_model_class, dat)
head(som_model_class) 
som_model$grid

som.hc <- cutree(hclust(object.distances(
  som_model, "codes")), 2)
#plot(som_model)
add.cluster.boundaries(som_model, som.hc)

dev.off()
rm(list=ls())



#RDA

library(vegan)
library(ggplot2)

cera_table <- read.csv('cge.csv', header = T,row.names=1)
cera_matrix <- cera_table

envdat_raw <- read.csv('c4.csv', header = T, row.names = 1)
h.envdat_raw<-decostand(envdat_raw,method='hellinger')
envdat <- h.envdat_raw[match(row.names(cera_matrix), row.names(envdat_raw)),]


decorana(cera_matrix)

res <- rda(cera_matrix ~ ., envdat)

res
plot(res)

xxxx <- summary(res)
aa <- xxxx$concont$importance
aa <- round(aa, 4)

aa

library(dplyr)
clu<-read.csv('ccluster.csv')


pdat <- res$CCA
samples<-data.frame(sample = row.names(pdat$u),RDA1 = pdat$u[,1],RDA2 = pdat$u[,2],cluster=clu$cluster)
species<-data.frame(spece = row.names(pdat$v),RDA1 = pdat$v[,1],RDA2 = pdat$v[,2])
envi<-data.frame(en = row.names(pdat$biplot),RDA1 = pdat$biplot[,1],RDA2 = pdat$biplot[,2])


library(ggrepel)
p <- ggplot() + 
  geom_hline(aes(yintercept = 0), colour="gray88", linetype="dashed") + 
  geom_vline(aes(xintercept = 0), colour="gray88", linetype="dashed")  +

  geom_segment(data = species,aes(x=0, xend= RDA1, y=0, yend= RDA2 ), arrow = arrow(length = unit(0.3, "cm")), colour = 'navy',size=1) +
  geom_text(data = species,aes(x = RDA1*1.1, y = RDA2*1.1, label = spece), size = 3, colour = 'navy') +

  geom_segment(data = envi,aes(x=0, xend= RDA1, y=0, yend= RDA2 ), arrow = arrow(length = unit(0.3, "cm")), colour = 'brown',size=1) +
  geom_text(data = envi,aes(x = RDA1*1.1, y = RDA2*1.1, label = en), size = 3, colour = 'brown', check_overlap = F) +
  geom_point(data = samples, aes(x=RDA1, y=RDA2,color=factor(cluster)),size = 6)+

  geom_text_repel(data = samples,aes(x = RDA1*1.1, y = RDA2*1.1, label = sample),size = 3, colour = 'black')+
  theme_bw() +
  theme(panel.grid.major=element_line(color=NA), panel.grid.minor = element_blank(),
        panel.border = element_rect(color='black',size=1.5))+
        

  xlab(paste('RDA1 (', aa[2,1]*100, '%)', sep ='')) + ylab(paste('RDA2 (', aa[2,2]*100, '%)', sep ='')) +
  theme(text=element_text(family='Arial',face='bold',size=18))+
  theme(axis.text.x=element_text(vjust=1,size=15,face = "bold",color='black'))+
  theme(axis.text.y=element_text(vjust=1,size=15,face = "bold",color='black'))
  #stat_ellipse(aes(clu2$RDA1,clu2$RDA2,color=clu2$cluster,group=clu2$cluster ),size = 1,level = 0.6,linetype='dashed',show.legend = F)

print(p)
dev.off()

#UPSET

library(vegan)
library(rdacca.hp)
library(UpSetVP)
gen<-read.csv('cepi1.csv',row.names=1)
env<-read.csv('capcoa.csv',row.names=1)

mod <- rdacca.hp(gen,env,method = 'RDA', var.part = TRUE, type = 'adjR2', scale = FALSE)

upset_vp(mod, plot.hp = TRUE, order.part = 'effect', nVar = 40)
upset_vp(mod, plot.hp = TRUE, order.part = 'degree', nVar = 40)

barplot_hp(mod, col.fill = 'var', 
           col.color = c('#F39B7FFF', '#8491B4FF', '#91D1C2FF', '#DC0000FF'))

dev.off()

#RF

struc<-read.csv('v19h.csv',row.names=1)
library(tidyverse)
library(randomForest)
struc2 <- struc %>%
  mutate(Structure = as.factor(Structure))
set.seed(100)
cera_rf= randomForest(Structure ~ ., data = struc2, importance=TRUE, proximity=TRUE)
cera_rf

imp_cera <- as_tibble(round(importance(cera_rf), 2),rownames = "Env") %>%
  arrange(desc(MeanDecreaseAccuracy))
imp_cera

ncol(struc2)
mycera= struc2[,-1]
set.seed(100)

result<-rfcv(mycera, struc2$Structure, cv.fold=10, scale = "log", step = 0.9)
result1<-result
result1$n.var

with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

error.cv<-data.frame(num=result$n.var,error.1=result$error.cv)
for (i in 100:104){
  print(i)
  set.seed(i)
  result= rfcv(mycera, struc2$Structure, cv.fold=10, scale = "log", step = 0.9)
  error.cv = cbind(error.cv, result$error.cv)
}

n.var <- error.cv$num
error.cv <- error.cv[,-1]
colnames(error.cv)<- paste('err',1:5,sep='.')
err.mean <-apply(error.cv,1,mean)
allerr<-data.frame(num=n.var,err.mean=err.mean,error.cv)
head(allerr[,1:6])

optimal = 2
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=1, colour="black"),
                   axis.line.y=element_line(size=1, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=12),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=12),
                   text=element_text(family="sans", size=12))
max(allerr$num)

ggplot() +
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey',size=0.5) +
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey',size=0.5) +
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey',size=0.5) +
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey',size=0.5) +
  geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey',size=0.5) +
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black',size=1) +
  geom_vline(xintercept = optimal, colour='black', lwd=0.5, linetype=2) +
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 5, 10, 20, 30, 50)) + 
  labs( x='Number of Envs ', y='Cross-validation error rate') +
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("Optimal = ", optimal, sep=""))+
  main_theme

result$n.var

write.csv(imp_cera[1:9,], "biomarks.csv")
library(ggsci)
library(RColorBrewer)
col<-brewer.pal(9,"YlGnBu")
col
imp_cera[1:8,] %>%
  select(Env,MeanDecreaseGini) %>%
  arrange(MeanDecreaseGini) %>%
  mutate(Env = forcats::fct_inorder(Env)) %>%
  ggplot(aes(x = Env, y = MeanDecreaseGini))+

  geom_bar(aes(fill=Env),stat = "identity")+
  #scale_fill_lancet()+
  scale_fill_manual(values = col)+

  labs(x = "", y = "Mean decrease accuracy")+
  coord_flip()+

  main_theme ->p2
p2

