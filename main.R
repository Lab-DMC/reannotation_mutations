source("/code/functions_four_methods.R")

install.packages("openxlsx")
library(openxlsx)


gene_tier12<-read.table("/data/Census_all_15_03_2024.tsv",header=T,sep="\t",stringsAsFactors = F)

skcm_maf<-read.table("/data/skcm_maf.txt",sep="\t",row.names=1,header=T,stringsAsFactors=F)
nsmpl<-length(unique(skcm_maf[,16]))



qresult_seq_length<-read.table("/data/length_transcripts_sequence_mc3_all.txt",sep="\t",row.names=NULL,header=T,stringsAsFactors=F)
tr_pr_length<-qresult_seq_length[,3]
names(tr_pr_length)<-qresult_seq_length[,1]



#Entropy analysis 
skcmentr_nSMs<-entropy_cluster(skcm_maf[which(skcm_maf[,51]!="synonymous_variant"),],38,1,5,6,56,nsmpl,116,37,16,51,"missense_variant",qresult_seq_length,0.20,rm_mut=5,1000,gene_tier12[,1],signif=T)
skcmentr_SMs<-entropy_cluster(skcm_maf[which(skcm_maf[,51]=="synonymous_variant"),],38,1,5,6,56,nsmpl,116,37,16,51,"synonymous_variant",qresult_seq_length,0.20,rm_mut=5,1000,gene_tier12[,1],signif=T)
skcmentr_all<-entropy_cluster(skcm_maf,38,1,5,6,56,nsmpl,116,37,16,51,"synonymous_variant",qresult_seq_length,0.20,rm_mut=5,1000,gene_tier12[,1],signif=T)
list_entr_skcm<-list(skcmentr_nSMs,skcmentr_SMs,skcmentr_all)
names(list_entr_skcm)<-c("nSMs_entr","SMs_entr","All_entr")
  
 #Concentration analysis
skcmconc_nSMs<-concentration_cluster(skcm_maf[which(skcm_maf[,51]!="synonymous_variant"),],38,1,6,nsmpl,116,37,56,16,51,"missense_variant",tr_pr_length,40,thr_mut=70,rm_mut=5,1000,gene_tier12[,1],signif=T,qresult_seq_length,5)
skcmconc_SMs<-concentration_cluster(skcm_maf[which(skcm_maf[,51]=="synonymous_variant"),],38,1,6,nsmpl,116,37,56,16,51,"synonymous_variant",tr_pr_length,40,thr_mut=70,rm_mut=5,1000,gene_tier12[,1],signif=T,qresult_seq_length,5)
skcmconc_all<-concentration_cluster(skcm_maf,38,1,6,nsmpl,116,37,56,16,51,"synonymous_variant",tr_pr_length,40,thr_mut=70,rm_mut=5,1000,gene_tier12[,1],signif=T,qresult_seq_length,5)
list_conc_skcm<-list(skcmconc_nSMs,skcmconc_SMs,skcmconc_all)
names(list_conc_skcm)<-c("nSMs_conc","SMs_conc","All_conc")

#Hotspot 12 analysis
skcmnSMs_hot<-window_clustering(skcm_maf[which(skcm_maf[,51]!="synonymous_variant"),],6,5,38,1,56,37,116,11,16,51,"missense_variant","Transcript_ID",5, nsmpl,40,qresult_seq_length,1000,gene_tier12[,1])
skcmSMs_hot<-window_clustering(skcm_maf[which(skcm_maf[,51]=="synonymous_variant"),],6,5,38,1,56,37,116,11,16,51,"synonymous_variant","Transcript_ID",5, nsmpl,40,qresult_seq_length,1000,gene_tier12[,1])
skcmall_hot<-window_clustering(skcm_maf,6,5,38,1,56,37,116,11,16,51,"synonymous_variant","Transcript_ID",5,nsmpl,40,qresult_seq_length,1000,gene_tier12[,1])
skcm_aggrhot<-list(skcmnSMs_hot[[3]],skcmSMs_hot[[3]],skcmall_hot[[3]])
names(skcm_aggrhot)<-c("nSMs_hot12","SMs_hot12","All_hot12")
library(openxlsx)
  
results_skcm<-c(list_entr_skcm,list_conc_skcm,skcm_aggrhot)
names(results_skcm)<-c(names(list_entr_skcm),names(list_conc_skcm),names(skcm_aggrhot))
write.xlsx(results_skcm,"/results/Three_methods_skcm.xlsx",rowNames=T)
