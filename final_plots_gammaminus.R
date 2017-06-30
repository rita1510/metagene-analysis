

library(psych)
library(stats)
library(fields)
library("chemometrics")
library(RcppRoll)
library(ggplot2)
#==================================Non_normalized Damage===========================================================



WTCHG_350045_236_F2_sense_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/sense_WTCHG_350045_236_AsiSI_hg38_gammaminus94"
WTCHG_350045_236_F2_anti_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/anti_WTCHG_350045_236_AsiSI_hg38_gammaminus94"


WTCHG_350045_236_F2_sense <- as.matrix(read.table(WTCHG_350045_236_F2_sense_file))
WTCHG_350045_236_F2_anti <- as.matrix(read.table(WTCHG_350045_236_F2_anti_file))


dim(WTCHG_350045_236_F2_sense)


RPM_WTCHG_350045_236_F2_sense1<-as.vector(apply(WTCHG_350045_236_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_350045_236_F2_anti1<-as.vector(apply(WTCHG_350045_236_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_350045_236_F2_sense2<- matrix(RPM_WTCHG_350045_236_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_350045_236_F2_anti2<- matrix(RPM_WTCHG_350045_236_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_350045_236_F2_sense3<-rowMeans(RPM_WTCHG_350045_236_F2_sense2)
RPM_WTCHG_350045_236_F2_anti3<-rowMeans(RPM_WTCHG_350045_236_F2_anti2)


RPM_WTCHG_350045_236_F2_sense4<-roll_mean(RPM_WTCHG_350045_236_F2_sense3, n = 5)
RPM_WTCHG_350045_236_F2_anti4<-roll_mean(RPM_WTCHG_350045_236_F2_anti3, n = 5)

length(RPM_WTCHG_350045_236_F2_sense4)

sems_WTCHG_350045_236_F2_sense3 <-apply(RPM_WTCHG_350045_236_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_350045_236_F2_sense2)[2])
sems_WTCHG_350045_236_F2_anti3 <-apply(RPM_WTCHG_350045_236_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_350045_236_F2_anti2)[2])


WTCHG_350045_240_F2_sense_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/sense_WTCHG_350045_240_AsiSI_hg38_gammaminus94"
WTCHG_350045_240_F2_anti_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/anti_WTCHG_350045_240_AsiSI_hg38_gammaminus94"



WTCHG_350045_240_F2_sense <- as.matrix(read.table(WTCHG_350045_240_F2_sense_file))
WTCHG_350045_240_F2_anti <- as.matrix(read.table(WTCHG_350045_240_F2_anti_file))


RPM_WTCHG_350045_240_F2_sense1<-as.vector(apply(WTCHG_350045_240_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_350045_240_F2_anti1<-as.vector(apply(WTCHG_350045_240_F2_anti,2,mean,trim=0.0))

RPM_WTCHG_350045_240_F2_sense2<- matrix(RPM_WTCHG_350045_240_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_350045_240_F2_anti2<- matrix(RPM_WTCHG_350045_240_F2_anti1[1:4001], ncol=1, byrow=TRUE)



RPM_WTCHG_350045_240_F2_sense3<-rowMeans(RPM_WTCHG_350045_240_F2_sense2)
RPM_WTCHG_350045_240_F2_anti3<-rowMeans(RPM_WTCHG_350045_240_F2_anti2)


RPM_WTCHG_350045_240_F2_sense4<-roll_mean(RPM_WTCHG_350045_240_F2_sense3, n = 5)
RPM_WTCHG_350045_240_F2_anti4<-roll_mean(RPM_WTCHG_350045_240_F2_anti3, n = 5)

length(RPM_WTCHG_350045_240_F2_sense4)
sems_WTCHG_350045_240_F2_sense3 <-apply(RPM_WTCHG_350045_240_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_350045_240_F2_sense2)[2])
sems_WTCHG_350045_240_F2_anti3 <-apply(RPM_WTCHG_350045_240_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_350045_240_F2_anti2)[2])

(a<-(-2000+length(RPM_WTCHG_350045_236_F2_sense4)-1))
distance_from_end=c(-2000:a)

df= data.frame(x=distance_from_end,RPM_WTCHG_350045_236_F2_sense4,RPM_WTCHG_350045_236_F2_anti4)

(h<-ggplot(df,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_350045_236_F2_sense4, colour="sense"))
  +geom_line(aes(y = RPM_WTCHG_350045_236_F2_anti4, colour="antisense"))
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  #+xlim(-500,1200)
  + ggtitle("Y1P OH+ yH2AX-, uncut 94 genic AsiSI sites"))
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Tyr_nonnormalized_gammaminus86sites_damage.png", width=10, height=7)



(a<-(-2000+length(RPM_WTCHG_350045_236_F2_sense4)-1))
distance_from_end=c(-2000:a)

df= data.frame(x=distance_from_end,RPM_WTCHG_350045_236_F2_sense4,RPM_WTCHG_350045_236_F2_anti4,RPM_WTCHG_350045_240_F2_sense4,RPM_WTCHG_350045_240_F2_anti4)

(h<-ggplot(df,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_350045_240_F2_sense4, colour="sense"))
  +geom_line(aes(y = RPM_WTCHG_350045_240_F2_anti4, colour="antisense"))
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  #+xlim(-500,1200)
  + ggtitle("Total OH+ yH2AX-, uncut 94 genic sites"))
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Total_nonnormalized_gammaminus86sites_damage.png", width=10, height=7)



#==================================Non_normalized Damage===========================================================



WTCHG_350045_236_F2_sense_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/sense_WTCHG_350045_236_AsiSI_hg38_intersect_94"
WTCHG_350045_236_F2_anti_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/anti_WTCHG_350045_236_AsiSI_hg38_intersect_94"


WTCHG_350045_236_F2_sense <- as.matrix(read.table(WTCHG_350045_236_F2_sense_file))
WTCHG_350045_236_F2_anti <- as.matrix(read.table(WTCHG_350045_236_F2_anti_file))


dim(WTCHG_350045_236_F2_sense)


RPM_WTCHG_350045_236_F2_sense1<-as.vector(apply(WTCHG_350045_236_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_350045_236_F2_anti1<-as.vector(apply(WTCHG_350045_236_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_350045_236_F2_sense2<- matrix(RPM_WTCHG_350045_236_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_350045_236_F2_anti2<- matrix(RPM_WTCHG_350045_236_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_350045_236_F2_sense3<-rowMeans(RPM_WTCHG_350045_236_F2_sense2)
RPM_WTCHG_350045_236_F2_anti3<-rowMeans(RPM_WTCHG_350045_236_F2_anti2)


RPM_WTCHG_350045_236_F2_sense4<-roll_mean(RPM_WTCHG_350045_236_F2_sense3, n = 5)
RPM_WTCHG_350045_236_F2_anti4<-roll_mean(RPM_WTCHG_350045_236_F2_anti3, n = 5)

length(RPM_WTCHG_350045_236_F2_sense4)

sems_WTCHG_350045_236_F2_sense3 <-apply(RPM_WTCHG_350045_236_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_350045_236_F2_sense2)[2])
sems_WTCHG_350045_236_F2_anti3 <-apply(RPM_WTCHG_350045_236_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_350045_236_F2_anti2)[2])


WTCHG_350045_240_F2_sense_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/sense_WTCHG_350045_240_AsiSI_hg38_intersect_94"
WTCHG_350045_240_F2_anti_file<-"/Volumes/Maxtor/Rita/15052017/Y1P_Damage_new/metagenes_final/merge/anti_WTCHG_350045_240_AsiSI_hg38_intersect_94"



WTCHG_350045_240_F2_sense <- as.matrix(read.table(WTCHG_350045_240_F2_sense_file))
WTCHG_350045_240_F2_anti <- as.matrix(read.table(WTCHG_350045_240_F2_anti_file))


RPM_WTCHG_350045_240_F2_sense1<-as.vector(apply(WTCHG_350045_240_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_350045_240_F2_anti1<-as.vector(apply(WTCHG_350045_240_F2_anti,2,mean,trim=0.0))

RPM_WTCHG_350045_240_F2_sense2<- matrix(RPM_WTCHG_350045_240_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_350045_240_F2_anti2<- matrix(RPM_WTCHG_350045_240_F2_anti1[1:4001], ncol=1, byrow=TRUE)



RPM_WTCHG_350045_240_F2_sense3<-rowMeans(RPM_WTCHG_350045_240_F2_sense2)
RPM_WTCHG_350045_240_F2_anti3<-rowMeans(RPM_WTCHG_350045_240_F2_anti2)


RPM_WTCHG_350045_240_F2_sense4<-roll_mean(RPM_WTCHG_350045_240_F2_sense3, n = 5)
RPM_WTCHG_350045_240_F2_anti4<-roll_mean(RPM_WTCHG_350045_240_F2_anti3, n = 5)

length(RPM_WTCHG_350045_240_F2_sense4)
sems_WTCHG_350045_240_F2_sense3 <-apply(RPM_WTCHG_350045_240_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_350045_240_F2_sense2)[2])
sems_WTCHG_350045_240_F2_anti3 <-apply(RPM_WTCHG_350045_240_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_350045_240_F2_anti2)[2])

(a<-(-2000+length(RPM_WTCHG_350045_236_F2_sense4)-1))
distance_from_end=c(-2000:a)

df= data.frame(x=distance_from_end,RPM_WTCHG_350045_236_F2_sense4,RPM_WTCHG_350045_236_F2_anti4)

(h<-ggplot(df,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_350045_236_F2_sense4, colour="sense"))
  +geom_line(aes(y = RPM_WTCHG_350045_236_F2_anti4, colour="antisense"))
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  #+xlim(-500,1200)
  + ggtitle("Y1P OH+ yH2AX+, cut 94 genic AsiSI sites"))
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Tyr_nonnormalized_gammaminus86sites_damage.png", width=10, height=7)



(a<-(-2000+length(RPM_WTCHG_350045_236_F2_sense4)-1))
distance_from_end=c(-2000:a)

df= data.frame(x=distance_from_end,RPM_WTCHG_350045_236_F2_sense4,RPM_WTCHG_350045_236_F2_anti4,RPM_WTCHG_350045_240_F2_sense4,RPM_WTCHG_350045_240_F2_anti4)

(h<-ggplot(df,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_350045_240_F2_sense4, colour="sense"))
  +geom_line(aes(y = RPM_WTCHG_350045_240_F2_anti4, colour="antisense"))
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  #+xlim(-500,1200)
  + ggtitle("Total OH+ yH2AX+, cut 94 genic sites"))
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Total_nonnormalized_gammaminus86sites_damage.png", width=10, height=7)







