

library(psych)
library(stats)
library(fields)
library("chemometrics")
library(RcppRoll)
library(ggplot2)

#=========================Y1P and Total +OH====================================

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
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Tyr_nonnormalized_gammaminus94sites_damage.png", width=10, height=7)



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
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Total_nonnormalized_gammaminus94sites_damage.png", width=10, height=7)



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
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Tyr_nonnormalized_gammaplus94sites_damage.png", width=10, height=7)



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
#ggsave(h, file="~/Desktop/Y1P_Damage_new/metagenes_final/plots_final_gammaminus/Total_nonnormalized_gammaplus94sites_damage.png", width=10, height=7)












#===========================gammaminus damage and non damage conditions=========================================================

gamma_minus_dist_file<-"~/Desktop/ELIFE/tss_from_asisi_distance_gammaminus.txt"
(gamma_minus_dist<-as.vector(read.table(gamma_minus_dist_file)[1]))


WTCHG_417720_289_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_289_AsiSI_hg38_gammaminus94"
WTCHG_417720_289_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_289_AsiSI_hg38_gammaminus94"


WTCHG_417720_289_F2_sense <- as.matrix(read.table(WTCHG_417720_289_F2_sense_file))
WTCHG_417720_289_F2_anti <- as.matrix(read.table(WTCHG_417720_289_F2_anti_file))





RPM_WTCHG_417720_289_F2_sense1<-as.vector(apply(WTCHG_417720_289_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_289_F2_anti1<-as.vector(apply(WTCHG_417720_289_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_417720_289_F2_sense2<- matrix(RPM_WTCHG_417720_289_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_289_F2_anti2<- matrix(RPM_WTCHG_417720_289_F2_anti1[1:4001], ncol=1, byrow=TRUE)



RPM_WTCHG_417720_289_F2_sense3<-rowMeans(RPM_WTCHG_417720_289_F2_sense2)
RPM_WTCHG_417720_289_F2_anti3<-rowMeans(RPM_WTCHG_417720_289_F2_anti2)


RPM_WTCHG_417720_289_F2_sense4<-roll_mean(RPM_WTCHG_417720_289_F2_sense3, n = 5)
RPM_WTCHG_417720_289_F2_anti4<-roll_mean(RPM_WTCHG_417720_289_F2_anti3, n = 5)

length(RPM_WTCHG_417720_289_F2_sense4)

sems_WTCHG_417720_289_F2_sense3 <-apply(RPM_WTCHG_417720_289_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_289_F2_sense2)[2])
sems_WTCHG_417720_289_F2_anti3 <-apply(RPM_WTCHG_417720_289_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_289_F2_anti2)[2])

WTCHG_417720_290_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_290_AsiSI_hg38_gammaminus94"
WTCHG_417720_290_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_290_AsiSI_hg38_gammaminus94"


WTCHG_417720_290_F2_sense <- as.matrix(read.table(WTCHG_417720_290_F2_sense_file))
WTCHG_417720_290_F2_anti <- as.matrix(read.table(WTCHG_417720_290_F2_anti_file))


dim(WTCHG_417720_290_F2_sense)


RPM_WTCHG_417720_290_F2_sense1<-as.vector(apply(WTCHG_417720_290_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_290_F2_anti1<-as.vector(apply(WTCHG_417720_290_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_417720_290_F2_sense2<- matrix(RPM_WTCHG_417720_290_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_290_F2_anti2<- matrix(RPM_WTCHG_417720_290_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_417720_290_F2_sense3<-rowMeans(RPM_WTCHG_417720_290_F2_sense2)
RPM_WTCHG_417720_290_F2_anti3<-rowMeans(RPM_WTCHG_417720_290_F2_anti2)


RPM_WTCHG_417720_290_F2_sense4<-roll_mean(RPM_WTCHG_417720_290_F2_sense3, n = 5)
RPM_WTCHG_417720_290_F2_anti4<-roll_mean(RPM_WTCHG_417720_290_F2_anti3, n = 5)

length(RPM_WTCHG_417720_290_F2_sense4)

sems_WTCHG_417720_290_F2_sense3 <-apply(RPM_WTCHG_417720_290_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_290_F2_sense2)[2])
sems_WTCHG_417720_290_F2_anti3 <-apply(RPM_WTCHG_417720_290_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_290_F2_anti2)[2])

WTCHG_417720_293_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_293_AsiSI_hg38_gammaminus94"
WTCHG_417720_293_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_293_AsiSI_hg38_gammaminus94"


WTCHG_417720_293_F2_sense <- as.matrix(read.table(WTCHG_417720_293_F2_sense_file))
WTCHG_417720_293_F2_anti <- as.matrix(read.table(WTCHG_417720_293_F2_anti_file))


dim(WTCHG_417720_293_F2_sense)


RPM_WTCHG_417720_293_F2_sense1<-as.vector(apply(WTCHG_417720_293_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_293_F2_anti1<-as.vector(apply(WTCHG_417720_293_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_417720_293_F2_sense2<- matrix(RPM_WTCHG_417720_293_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_293_F2_anti2<- matrix(RPM_WTCHG_417720_293_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_417720_293_F2_sense3<-rowMeans(RPM_WTCHG_417720_293_F2_sense2)
RPM_WTCHG_417720_293_F2_anti3<-rowMeans(RPM_WTCHG_417720_293_F2_anti2)


RPM_WTCHG_417720_293_F2_sense4<-roll_mean(RPM_WTCHG_417720_293_F2_sense3, n = 5)
RPM_WTCHG_417720_293_F2_anti4<-roll_mean(RPM_WTCHG_417720_293_F2_anti3, n = 5)

length(RPM_WTCHG_417720_293_F2_sense4)

sems_WTCHG_417720_293_F2_sense3 <-apply(RPM_WTCHG_417720_293_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_293_F2_sense2)[2])
sems_WTCHG_417720_293_F2_anti3 <-apply(RPM_WTCHG_417720_293_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_293_F2_anti2)[2])

WTCHG_417720_294_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_294_AsiSI_hg38_gammaminus94"
WTCHG_417720_294_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_294_AsiSI_hg38_gammaminus94"


WTCHG_417720_294_F2_sense <- as.matrix(read.table(WTCHG_417720_294_F2_sense_file))
WTCHG_417720_294_F2_anti <- as.matrix(read.table(WTCHG_417720_294_F2_anti_file))


dim(WTCHG_417720_294_F2_sense)


RPM_WTCHG_417720_294_F2_sense1<-as.vector(apply(WTCHG_417720_294_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_294_F2_anti1<-as.vector(apply(WTCHG_417720_294_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_417720_294_F2_sense2<- matrix(RPM_WTCHG_417720_294_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_294_F2_anti2<- matrix(RPM_WTCHG_417720_294_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_417720_294_F2_sense3<-rowMeans(RPM_WTCHG_417720_294_F2_sense2)
RPM_WTCHG_417720_294_F2_anti3<-rowMeans(RPM_WTCHG_417720_294_F2_anti2)


RPM_WTCHG_417720_294_F2_sense4<-roll_mean(RPM_WTCHG_417720_294_F2_sense3, n = 5)
RPM_WTCHG_417720_294_F2_anti4<-roll_mean(RPM_WTCHG_417720_294_F2_anti3, n = 5)

length(RPM_WTCHG_417720_294_F2_sense4)

sems_WTCHG_417720_294_F2_sense3 <-apply(RPM_WTCHG_417720_294_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_294_F2_sense2)[2])
sems_WTCHG_417720_294_F2_anti3 <-apply(RPM_WTCHG_417720_294_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_294_F2_anti2)[2])
tss_minus<-rep(1, dim(gamma_minus_dist)[1])
df_dist_minus=data.frame(x=gamma_minus_dist[1], y=tss_minus)

(a<-(-2000+length(RPM_WTCHG_417720_289_F2_sense4)-1))
distance_from_end=c(-2000:a)


dftot= data.frame(x=distance_from_end,RPM_WTCHG_417720_289_F2_sense4,RPM_WTCHG_417720_289_F2_anti4,RPM_WTCHG_417720_293_F2_sense4,RPM_WTCHG_417720_293_F2_anti4)

(h<-ggplot(dftot,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist_minus,aes(x=gamma_minus_dist[1], y=tss_minus),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-2000,2000)
  + ggtitle("Total yH2AX-, uncut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Total_nonnormalized_gammaminus94sites_damage.png", width=10, height=7)

dfy1p= data.frame(x=distance_from_end,RPM_WTCHG_417720_290_F2_sense4,RPM_WTCHG_417720_290_F2_anti4,RPM_WTCHG_417720_294_F2_sense4,RPM_WTCHG_417720_294_F2_anti4)


(h<-ggplot(dfy1p,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist_minus,aes(x=gamma_minus_dist[1], y=tss_minus),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-2000,2000)
  + ggtitle("Y1P yH2AX-, uncut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Y1P_nonnormalized_gammaminus94sites_damage.png", width=10, height=7)



dftot= data.frame(x=distance_from_end,RPM_WTCHG_417720_289_F2_sense4,RPM_WTCHG_417720_289_F2_anti4,RPM_WTCHG_417720_293_F2_sense4,RPM_WTCHG_417720_293_F2_anti4)

(h<-ggplot(dftot,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist_minus,aes(x=gamma_minus_dist[1], y=tss_minus),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-500,500)
  + ggtitle("Total yH2AX-, uncut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Total_nonnormalized_gammaminus94sites_damage_pm500.png", width=10, height=7)


dfy1p= data.frame(x=distance_from_end,RPM_WTCHG_417720_290_F2_sense4,RPM_WTCHG_417720_290_F2_anti4,RPM_WTCHG_417720_294_F2_sense4,RPM_WTCHG_417720_294_F2_anti4)


(h<-ggplot(dfy1p,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist_minus,aes(x=gamma_minus_dist[1], y=tss_minus),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-500,500)
  + ggtitle("Y1P yH2AX-, uncut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Y1P_nonnormalized_gammaminus94sites_damage_pm500.png", width=10, height=7)





#========================================gammaplus===========================================
gamma_plus_dist_file<-"~/Desktop/ELIFE/tss_from_asisi_distance_gammaplus.txt"
(gamma_plus_dist<-as.vector(read.table(gamma_plus_dist_file)[1]))



WTCHG_417720_289_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_289_AsiSI_hg38_intersect_94"
WTCHG_417720_289_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_289_AsiSI_hg38_intersect_94"

WTCHG_417720_289_F2_sense <- as.matrix(read.table(WTCHG_417720_289_F2_sense_file))
WTCHG_417720_289_F2_anti <- as.matrix(read.table(WTCHG_417720_289_F2_anti_file))


dim(WTCHG_417720_289_F2_sense)

WTCHG_417720_289_F2_sense[,2001:2501]

RPM_WTCHG_417720_289_F2_sense1<-as.vector(apply(WTCHG_417720_289_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_289_F2_anti1<-as.vector(apply(WTCHG_417720_289_F2_anti,2,mean,trim=0.0))



RPM_WTCHG_417720_289_F2_sense2<- matrix(RPM_WTCHG_417720_289_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_289_F2_anti2<- matrix(RPM_WTCHG_417720_289_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_417720_289_F2_sense3<-rowMeans(RPM_WTCHG_417720_289_F2_sense2)
RPM_WTCHG_417720_289_F2_anti3<-rowMeans(RPM_WTCHG_417720_289_F2_anti2)


RPM_WTCHG_417720_289_F2_sense4<-roll_mean(RPM_WTCHG_417720_289_F2_sense3, n = 5)
RPM_WTCHG_417720_289_F2_anti4<-roll_mean(RPM_WTCHG_417720_289_F2_anti3, n = 5)

length(RPM_WTCHG_417720_289_F2_sense4)

sems_WTCHG_417720_289_F2_sense3 <-apply(RPM_WTCHG_417720_289_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_289_F2_sense2)[2])
sems_WTCHG_417720_289_F2_anti3 <-apply(RPM_WTCHG_417720_289_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_289_F2_anti2)[2])

WTCHG_417720_290_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_290_AsiSI_hg38_intersect_94"
WTCHG_417720_290_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_290_AsiSI_hg38_intersect_94"


WTCHG_417720_290_F2_sense <- as.matrix(read.table(WTCHG_417720_290_F2_sense_file))
WTCHG_417720_290_F2_anti <- as.matrix(read.table(WTCHG_417720_290_F2_anti_file))


dim(WTCHG_417720_290_F2_sense)


RPM_WTCHG_417720_290_F2_sense1<-as.vector(apply(WTCHG_417720_290_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_290_F2_anti1<-as.vector(apply(WTCHG_417720_290_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_417720_290_F2_sense2<- matrix(RPM_WTCHG_417720_290_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_290_F2_anti2<- matrix(RPM_WTCHG_417720_290_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_417720_290_F2_sense3<-rowMeans(RPM_WTCHG_417720_290_F2_sense2)
RPM_WTCHG_417720_290_F2_anti3<-rowMeans(RPM_WTCHG_417720_290_F2_anti2)


RPM_WTCHG_417720_290_F2_sense4<-roll_mean(RPM_WTCHG_417720_290_F2_sense3, n = 5)
RPM_WTCHG_417720_290_F2_anti4<-roll_mean(RPM_WTCHG_417720_290_F2_anti3, n = 5)

length(RPM_WTCHG_417720_290_F2_sense4)

sems_WTCHG_417720_290_F2_sense3 <-apply(RPM_WTCHG_417720_290_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_290_F2_sense2)[2])
sems_WTCHG_417720_290_F2_anti3 <-apply(RPM_WTCHG_417720_290_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_290_F2_anti2)[2])


WTCHG_417720_293_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_293_AsiSI_hg38_intersect_94"
WTCHG_417720_293_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_293_AsiSI_hg38_intersect_94"


WTCHG_417720_293_F2_sense <- as.matrix(read.table(WTCHG_417720_293_F2_sense_file))
WTCHG_417720_293_F2_anti <- as.matrix(read.table(WTCHG_417720_293_F2_anti_file))


dim(WTCHG_417720_293_F2_sense)


RPM_WTCHG_417720_293_F2_sense1<-as.vector(apply(WTCHG_417720_293_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_293_F2_anti1<-as.vector(apply(WTCHG_417720_293_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_417720_293_F2_sense2<- matrix(RPM_WTCHG_417720_293_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_293_F2_anti2<- matrix(RPM_WTCHG_417720_293_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_417720_293_F2_sense3<-rowMeans(RPM_WTCHG_417720_293_F2_sense2)
RPM_WTCHG_417720_293_F2_anti3<-rowMeans(RPM_WTCHG_417720_293_F2_anti2)


RPM_WTCHG_417720_293_F2_sense4<-roll_mean(RPM_WTCHG_417720_293_F2_sense3, n = 5)
RPM_WTCHG_417720_293_F2_anti4<-roll_mean(RPM_WTCHG_417720_293_F2_anti3, n = 5)

length(RPM_WTCHG_417720_293_F2_sense4)

sems_WTCHG_417720_293_F2_sense3 <-apply(RPM_WTCHG_417720_293_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_293_F2_sense2)[2])
sems_WTCHG_417720_293_F2_anti3 <-apply(RPM_WTCHG_417720_293_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_293_F2_anti2)[2])

WTCHG_417720_294_F2_sense_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/sense_WTCHG_417720_294_AsiSI_hg38_intersect_94"
WTCHG_417720_294_F2_anti_file<-"~/Desktop/ELIFE/metagenes_ELIFE/merge/anti_WTCHG_417720_294_AsiSI_hg38_intersect_94"


WTCHG_417720_294_F2_sense <- as.matrix(read.table(WTCHG_417720_294_F2_sense_file))
WTCHG_417720_294_F2_anti <- as.matrix(read.table(WTCHG_417720_294_F2_anti_file))


dim(WTCHG_417720_294_F2_sense)


RPM_WTCHG_417720_294_F2_sense1<-as.vector(apply(WTCHG_417720_294_F2_sense,2,mean,trim=0.0))
RPM_WTCHG_417720_294_F2_anti1<-as.vector(apply(WTCHG_417720_294_F2_anti,2,mean,trim=0.0))


RPM_WTCHG_417720_294_F2_sense2<- matrix(RPM_WTCHG_417720_294_F2_sense1[1:4001], ncol=1, byrow=TRUE)
RPM_WTCHG_417720_294_F2_anti2<- matrix(RPM_WTCHG_417720_294_F2_anti1[1:4001], ncol=1, byrow=TRUE)




RPM_WTCHG_417720_294_F2_sense3<-rowMeans(RPM_WTCHG_417720_294_F2_sense2)
RPM_WTCHG_417720_294_F2_anti3<-rowMeans(RPM_WTCHG_417720_294_F2_anti2)


RPM_WTCHG_417720_294_F2_sense4<-roll_mean(RPM_WTCHG_417720_294_F2_sense3, n = 5)
RPM_WTCHG_417720_294_F2_anti4<-roll_mean(RPM_WTCHG_417720_294_F2_anti3, n = 5)

length(RPM_WTCHG_417720_294_F2_sense4)

sems_WTCHG_417720_294_F2_sense3 <-apply(RPM_WTCHG_417720_294_F2_sense2,1,SD)/sqrt(dim(RPM_WTCHG_417720_294_F2_sense2)[2])
sems_WTCHG_417720_294_F2_anti3 <-apply(RPM_WTCHG_417720_294_F2_anti2,1,SD)/sqrt(dim(RPM_WTCHG_417720_294_F2_anti2)[2])

(a<-(-2000+length(RPM_WTCHG_417720_289_F2_sense4)-1))
distance_from_end=c(-2000:a)

tss<-rep(1, dim(gamma_plus_dist)[1])
df_dist=data.frame(x=gamma_plus_dist[1], y=tss)



dftot= data.frame(x=distance_from_end,RPM_WTCHG_417720_289_F2_sense4,RPM_WTCHG_417720_289_F2_anti4,RPM_WTCHG_417720_293_F2_sense4,RPM_WTCHG_417720_293_F2_anti4)


(h<-ggplot(dftot,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist,aes(x=gamma_plus_dist[1], y=tss),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-2000,2000)
  + ggtitle("Total yH2AX+, cut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Total_nonnormalized_intersect_94sites_damage.png", width=10, height=7)

dfy1p= data.frame(x=distance_from_end,RPM_WTCHG_417720_290_F2_sense4,RPM_WTCHG_417720_290_F2_anti4,RPM_WTCHG_417720_294_F2_sense4,RPM_WTCHG_417720_294_F2_anti4)


(h<-ggplot(dfy1p,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist,aes(x=gamma_plus_dist[1], y=tss),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-2000,2000)
  + ggtitle("Y1P yH2AX+, cut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Y1P_nonnormalized_intersect_94sites_damage.png", width=10, height=7)


dftot= data.frame(x=distance_from_end,RPM_WTCHG_417720_289_F2_sense4,RPM_WTCHG_417720_289_F2_anti4,RPM_WTCHG_417720_293_F2_sense4,RPM_WTCHG_417720_293_F2_anti4)

(h<-ggplot(dftot,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_289_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_293_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist,aes(x=gamma_plus_dist[1], y=tss),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-500,500)
  + ggtitle("Total yH2AX+, cut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Total_nonnormalized_intersect_94sites_damage_pm500.png", width=10, height=7)


dfy1p= data.frame(x=distance_from_end,RPM_WTCHG_417720_290_F2_sense4,RPM_WTCHG_417720_290_F2_anti4,RPM_WTCHG_417720_294_F2_sense4,RPM_WTCHG_417720_294_F2_anti4)


(h<-ggplot(dfy1p,aes(x = distance_from_end))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_sense4, colour="sense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_290_F2_anti4, colour="antisense WT U20S"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_sense4, colour="sense AsiSI-ER U2OS"))
  +geom_line(aes(y = RPM_WTCHG_417720_294_F2_anti4, colour="antisense AsiSI-ER U2OS"))
  #+geom_point(data=df_dist,aes(x=gamma_plus_dist[1], y=tss),shape="I")
  +theme(legend.title=element_blank())
  + theme_bw()
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  +ylab(label="signal average")
  +xlab("distance from AsiSI site")
  +ylim(-0.5,1.5)
  +xlim(-500,500)
  + ggtitle("Y1P yH2AX+, cut 94 genic AsiSI sites"))
ggsave(h, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/Y1P_nonnormalized_intersect_94sites_damage_pm500.png", width=10, height=7)
ggplot_build(h)$data

#=====================================boxplots unpaired [0 500]/[1500 2000]====================


RPM_WTCHG_417720_289_F2_sense_sum_reg1<-as.vector(apply(WTCHG_417720_289_F2_sense[,2001:2501],1,sum,trim=0.0))
RPM_WTCHG_417720_289_F2_sense_sum_reg2<-as.vector(apply(WTCHG_417720_289_F2_sense[,3501:4001],1,sum,trim=0.0))
RPM_WTCHG_417720_289_F2_sense_sum_ratio<-log2(RPM_WTCHG_417720_289_F2_sense_sum_reg1/RPM_WTCHG_417720_289_F2_sense_sum_reg2)

RPM_WTCHG_417720_290_F2_sense_sum_reg1<-as.vector(apply(WTCHG_417720_290_F2_sense[,2001:2501],1,sum,trim=0.0))
RPM_WTCHG_417720_290_F2_sense_sum_reg2<-as.vector(apply(WTCHG_417720_290_F2_sense[,3501:4001],1,sum,trim=0.0))
RPM_WTCHG_417720_290_F2_sense_sum_ratio<-log2(RPM_WTCHG_417720_290_F2_sense_sum_reg1/RPM_WTCHG_417720_290_F2_sense_sum_reg2)


RPM_WTCHG_417720_293_F2_sense_sum_reg1<-as.vector(apply(WTCHG_417720_293_F2_sense[,2001:2501],1,sum,trim=0.0))
RPM_WTCHG_417720_293_F2_sense_sum_reg2<-as.vector(apply(WTCHG_417720_293_F2_sense[,3501:4001],1,sum,trim=0.0))
RPM_WTCHG_417720_293_F2_sense_sum_ratio<-log2(RPM_WTCHG_417720_293_F2_sense_sum_reg1/RPM_WTCHG_417720_293_F2_sense_sum_reg2)

RPM_WTCHG_417720_294_F2_sense_sum_reg1<-as.vector(apply(WTCHG_417720_294_F2_sense[,2001:2501],1,sum,trim=0.0))
RPM_WTCHG_417720_294_F2_sense_sum_reg2<-as.vector(apply(WTCHG_417720_294_F2_sense[,3501:4001],1,sum,trim=0.0))
RPM_WTCHG_417720_294_F2_sense_sum_ratio<-log2(RPM_WTCHG_417720_294_F2_sense_sum_reg1/RPM_WTCHG_417720_294_F2_sense_sum_reg2)

(a2_289<-RPM_WTCHG_417720_289_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_289_F2_sense_sum_ratio)])
(a2_293<-RPM_WTCHG_417720_293_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_293_F2_sense_sum_ratio)])
(a2_290<-RPM_WTCHG_417720_290_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_290_F2_sense_sum_ratio)])
(a2_294<-RPM_WTCHG_417720_294_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_294_F2_sense_sum_ratio)])
W_289_293_u<-wilcox.test(a2_289, a2_293, paired=F)
p_289_293_u<-format(W_289_293_u$p.value, scientific = TRUE, digits=2)
W_290_294_u<-wilcox.test(a2_290, a2_294, paired=F)
p_290_294_u<-format(W_290_294_u$p.value, scientific = TRUE, digits=2)

df1 <- data.frame(type=c(rep("Total WT (74)",length(a2_289)),rep("Total AsiSI-ER (58)",length(a2_293)),rep("Y1P WT (76)",length(a2_290)),rep("Y1P AsiSI-ER (62)",length(a2_294))), val=c(a2_289,a2_293,a2_290,a2_294))
col=c("#7e51a0","#008689","#c77cff","#00bfc4")
(p <- ggplot(df1, aes(factor(type,levels = c("Total WT (74)","Total AsiSI-ER (58)","Y1P WT (76)","Y1P AsiSI-ER (62)")), val)) 
+theme(legend.title=element_blank())
+ theme_bw()
+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

(p<-p + geom_boxplot(notch = TRUE, fill=col)+xlab(label="Data")+ylab(label="log2(Signal Sum Ratio)")+ ggtitle(paste0("Unpaired data p_total=",p_289_293_u," p_Y1P=",p_290_294_u)))

ggsave(p, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/signal_boxplots_1500-2000_unpaired_col.png", width=10, height=7)
#=====================================boxplots paired [0 500]/[1500 2000]====================

a1_289<-RPM_WTCHG_417720_289_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_289_F2_sense_sum_ratio)]
(a1_293<-RPM_WTCHG_417720_293_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_289_F2_sense_sum_ratio)])
(b1_289<-a1_289[is.finite(a1_293)])
(b1_293<-a1_293[is.finite(a1_293)])

W_289_293<-wilcox.test(b1_289, b1_293, paired=T)
p_289_293<-format(W_289_293$p.value, scientific = TRUE, digits=2)


(a1_290<-RPM_WTCHG_417720_290_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_290_F2_sense_sum_ratio)])
(a1_294<-RPM_WTCHG_417720_294_F2_sense_sum_ratio[is.finite(RPM_WTCHG_417720_290_F2_sense_sum_ratio)])
(b1_290<-a1_290[is.finite(a1_294)])
(b1_294<-a1_294[is.finite(a1_294)])

W_290_294<-wilcox.test(b1_290, b1_294, paired=T)
p_290_294<-format(W_290_294$p.value, scientific = TRUE, digits=2)

df_paired <- data.frame(type=c(rep("Total WT (53)",length(b1_289)),rep("Total AsiSI-ER (53)",length(b1_293)),rep("Y1P WT (57)",length(b1_290)),rep("Y1P AsiSI-ER (57)",length(b1_294))), val=c(b1_289,b1_293,b1_290,b1_294))
col=c("#7e51a0","#008689","#c77cff","#00bfc4")
p <- ggplot(df_paired, aes(factor(type,levels = c("Total WT (53)","Total AsiSI-ER (53)","Y1P WT (57)","Y1P AsiSI-ER (57)")), val))+theme(legend.title=element_blank())

(p<-p + geom_boxplot(notch = TRUE, fill=col)+xlab(label="Data")+ylab(label="log2(Signal Sum Ratio) - Paired")+ ggtitle(paste0("Paired data p_total=",p_289_293," p_Y1P=",p_290_294))
  +theme(legend.title=element_blank())
+ theme_bw()
+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))



ggsave(p, file="~/Desktop/ELIFE/metagenes_ELIFE/plots/signal_boxplots_1500-2000_paired_col.png", width=10, height=7)




