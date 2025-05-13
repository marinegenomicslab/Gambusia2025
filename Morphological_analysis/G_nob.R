#welcome to my code 
#load in packages
library(readxl)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library (ggfortify)
library(geomorph)
library(writexl)
#install.packages("geomorph")

#just start here for geomorph stuff
FallLandmarks <- readmulti.tps(c("BalF.TPS","BitF.TPS","DiaF.TPS","type.TPS"),specID="imageID")
MallLandmarks<-readmulti.tps(c("BalM.TPS","BitM.TPS","DiaM.TPS"),specID="imageID")
#divide based on sex
#run procrustes
Fprocrustes <- gpagen(FallLandmarks) 
Mprocrustes <- gpagen(MallLandmarks)

#run PCA on each sex
#female first
Fpca<-gm.prcomp(Fprocrustes$coords)
Fpcadata <- as.data.frame(Fpca$x) #turn into data frame
Fpcadata<-Fpcadata[-c(6),]#remove the extra 13fish
summ<-summary(Fpca)

#find eigenvectors
Feigenvectors <- Fpca$rotation
#males
Mpca <- gm.prcomp(Mprocrustes$coords)
summ<-summary(Mpca)
Mpcadata <- as.data.frame(Mpca$x) #turn into data frame
Meigenvectors <- Mpca$rotation


#load in site and sex data so I can make a PCA

feminfo<-read_xlsx("g_nob.xlsx", sheet = "female")
maleinfo<-read_xlsx("g_nob.xlsx", sheet = "male")

#add the info into female and male
Fpcadata$site<-feminfo$site
Fpcadata$sex<-feminfo$sex

Mpcadata$site<-maleinfo$site
Mpcadata$sex<-maleinfo$sex


#female PCA
Fpc<-ggplot()+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=0.8,colour="lightgray")+
  geom_vline(aes(xintercept=0),linetype="dashed",size=0.8,colour="lightgray")+
  geom_point(data=Fpcadata,
             mapping = aes(x=Comp1,y=Comp2, colour = Fpcadata$site),
             alpha=0.8,size=3.5)+
  geom_segment(data=Fpcadata,
               mapping=aes(x=0,y=0,xend=Comp1,yend=Comp2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="lightgray",
               linewidth=0.8)+
  scale_color_manual(values = c("Balmorhae" = "lightpink1","Bitter_lake" = "steelblue2", "Diamond_y" = "lightgreen", "Type" = "green4"))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.15, 0.85),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.line = element_line(colour="black"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
xlab("PC1 32.5%")+
  ylab("PC2 17.8%")+
  labs(color="Site")
Fpc


 
#male PCA
Mpc<-ggplot()+
  geom_abline(intercept=0,slope=0,linetype="dashed",size=0.8,colour="lightgray")+
  geom_vline(aes(xintercept=0),linetype="dashed",size=0.8,colour="lightgray")+
  geom_point(data=Mpcadata,
             mapping = aes(x=Comp1,y=Comp2, colour = Mpcadata$site),
             alpha=0.8,size=3.5)+
  geom_segment(data=Mpcadata,
               mapping=aes(x=0,y=0,xend=Comp1,yend=Comp2),
               arrow=arrow(length=unit(0.015,"npc"),
                           type="closed"),
               colour="lightgray",
               linewidth=0.8)+
  scale_color_manual(values = c("Balmorhae" = "lightpink1","Bitter_lake" = "steelblue2", "Diamond_y" = "lightgreen"))+
theme(panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_blank(),
      legend.position = c(0.85, 0.85),
      legend.title = element_text(size=12),
      legend.text = element_text(size=12),
      axis.line = element_line(colour="black"),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12))+
  xlab("PC1 39.8%")+
  ylab("PC2 22.7%")+
  labs(color="Site")
Mpc

 
#start with the females
#make mean shape
femalemean<-mshape(Fprocrustes$coords)
malemean<-mshape(Mprocrustes$coords)

#Diamond_y female 9-24
#subset
gp1.mn<-Fprocrustes$coords[1:10,,31]
plotRefToTarget(femalemean,gp1.mn,mag = 2)
Links<-define.links(femalemean, ptsize = 1, links = NULL)
plotRefToTarget(femalemean,gp1.mn,mag = 2, link = Links)
#BitterLake Female BL 26-24
gp2.mn<-Fprocrustes$coords[1:10,,22]
plotRefToTarget(femalemean,gp2.mn,mag = 2)
Links<-define.links(femalemean, ptsize = 1, links = NULL)
plotRefToTarget(femalemean,gp2.mn,mag = 2, link = Links)
#Balhmorae Female B_14-24
gp3.mn<-Fprocrustes$coords[1:10,,15]
plotRefToTarget(femalemean,gp3.mn,mag = 2)
Links<-define.links(femalemean, ptsize = 1, links = NULL)
plotRefToTarget(femalemean,gp3.mn,mag = 2, link = Links)

#next for the males
#DiamondY- 08-24 male
gp4.mn<-Mprocrustes$coords[1:10,,18]
plotRefToTarget(malemean,gp4.mn,mag = 2)
Links<-define.links(malemean, ptsize = 1, links = NULL)
plotRefToTarget(malemean,gp4.mn,mag = 2, link = Links)
#BitterLake- BL-20-24
gp5.mn<-Mprocrustes$coords[1:10,,11]
plotRefToTarget(malemean,gp5.mn,mag = 2)
Links<-define.links(malemean, ptsize = 1, links = NULL)
plotRefToTarget(malemean,gp5.mn,mag = 2, link = Links)
#Balhmorae- B-17-24
gp6.mn<-Mprocrustes$coords[1:10,,5]
plotRefToTarget(malemean,gp6.mn,mag = 2)
Links<-define.links(malemean, ptsize = 1, links = NULL)
plotRefToTarget(malemean,gp6.mn,mag = 2, link = Links)







