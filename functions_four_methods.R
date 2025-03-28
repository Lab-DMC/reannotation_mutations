entropy_cluster<-function(maf_file,pos_tr,pos_gn,pos_chr,pos_start,pos_codon,nsmpl_df,cds_pos,prot_pos,pos_barcode,pos_variant,type_mut,pr_length,entr_thr,rm_mut=15,num_perm,signif=T,canc_genes)
{
  tr_gn_table<-unique(maf_file[,c(pos_tr,pos_gn)])
  sign_dataset<-vector()
  tr_unique<-unique(maf_file[,pos_tr])
  range<-vector()
  range_norm<-vector()
  n_type<-vector()
  tot_m<-vector()
  p_value<-vector()
  entropy_gene<-vector()
  nsaml_vec<-vector()
  nsaml_perc<-vector()
  perc_mut<-vector()
  mut<-vector()
  chrm_min<-vector()
  chrm_max<-vector()
  chr<-vector()
  mut_freq<-vector()
  mut_cds<-vector()
  mut_aa<-vector()
  vtype_mut<-vector()
  cod_mut<-vector()
  chr_pos_mut<-vector()
  entr_v_save<-NULL
  c<-1
  for(i in 1:length(tr_unique))
  {
  
    tr_maf<-maf_file[which(maf_file[,pos_tr]==tr_unique[i]),]
    nsample<-length(unique(tr_maf[,pos_barcode]))
    if(nsample>=rm_mut)
    {
      pos_null<-which(tr_maf[,cds_pos]=="")
      if(length(pos_null)>0)
        tr_maf<-tr_maf[-pos_null,]
      
      tr_maf1<-tr_maf
      
      prot_l<-pr_length[which(pr_length[,1]==tr_unique[i]),3]
      delins<-grep("delins",tr_maf1[,cds_pos])
      if(length(delins)>0)
        tr_maf1<-tr_maf1[-delins,]
      mut_pos_pr<-gsub("([A-Z.-]+)", "", tr_maf1[,cds_pos])
      mut_pos_pr<-as.numeric(gsub(">", "", mut_pos_pr))
      pos_na<-which(is.na(mut_pos_pr)==T)
      if(length(pos_na)>0)
      {
        mut_pos_pr<- mut_pos_pr[-pos_na]
        tr_maf1<-tr_maf1[-pos_na,]
      }
      mut_pr_freq<-table(mut_pos_pr)/nrow(tr_maf1)
      entr_v<-(-1)*sum(mut_pr_freq*log2(mut_pr_freq))/log2(prot_l)
      entr_v_save<-c(entr_v_save,entr_v)
      if(entr_v<=entr_thr)
      {
        entropy_gene[c]<-entr_v
        p_value[c]<-entr_randistr1(length(mut_pos_pr),prot_l,entropy_gene[c],num_perm)
        names(p_value)[c]<-tr_unique[i]
        names(entropy_gene)[c]<-tr_unique[i]
        range[c]<-max(mut_pos_pr)-min(mut_pos_pr)
        leng_gene<-pr_length[which(pr_length[,1]==tr_unique[i]),3]
        range_norm[c]<-range[c]/leng_gene
        names(range)[c]<-tr_unique[i]
        names(range_norm)[c]<-tr_unique[i]
        tb_type<-table(tr_maf1[,pos_variant])
        freq_mut<-table(tr_maf1[,cds_pos])
        pos_max<-which(freq_mut==max(freq_mut))[1]
        mut_freq[c]<-freq_mut[pos_max]
        mut_cds[c]<-unique(tr_maf1[which(tr_maf1[,cds_pos]==names(freq_mut)[pos_max]),cds_pos])
        mut_aa[c]<-unique(tr_maf1[which(tr_maf1[,cds_pos]==names(freq_mut)[pos_max]),prot_pos])
        chr_pos_mut[c]<-unique(tr_maf1[which(tr_maf1[,cds_pos]==names(freq_mut)[pos_max]),pos_start])
        vtype_mut[c]<-unique(tr_maf1[which(tr_maf1[,cds_pos]==names(freq_mut)[pos_max]),pos_variant])
        cod_mut[c]<-unique(tr_maf1[which(tr_maf1[,cds_pos]==names(freq_mut)[pos_max]),pos_codon])
        if( length(which(names(tb_type)==type_mut))==0)
        {
          perc_type<-0
          n_type[c]<-0
        }
        else
        {
          
          perc_type<-tb_type[which(names(tb_type)==type_mut)]/sum(tb_type)
          if(length(perc_type)>1)
          n_type[c]<-tb_type[which(names(tb_type)==type_mut)]
        }
        max_chrm<-max(tr_maf1[,pos_start])
        min_chrm<-min(tr_maf1[,pos_start])
        info_tr2<-tr_maf1[tr_maf1[,pos_start]<=max_chrm,]
        info_tr2<-tr_maf1[tr_maf1[,pos_start]>=min_chrm,]
        
        nsaml_vec[c]<-length(unique(info_tr2[,pos_barcode]))
        nsaml_perc[c]<- nsaml_vec[c]/nsmpl_df*100
        names(nsaml_perc)[c]<-tr_unique[i]
        names(nsaml_vec)[c]<-tr_unique[i]
        
        tot_m[c]<-sum(tb_type)
        perc_mut[c]<-perc_type
        names(perc_mut)<-tr_unique[i]
        chr[c]<-unique(tr_maf1[,pos_chr])
        names(chr)<-tr_unique[i]
        chrm_min[c]<-min(tr_maf1[,pos_start],na.rm=T)
        names(chrm_min)<-tr_unique[i]
        chrm_max[c]<-max(tr_maf1[,pos_start],na.rm=T)
        names(chrm_max)<-tr_unique[i]
        c<-c+1
      }
    }
  }
  if(length(p_value)==0)
  {
    return("no gene passed the threshold")
  }
  else
  {
    
    adj_pvalue<-p.adjust(p_value, method = "BH", n = length(p_value))
    adj_pvalue<-round(adj_pvalue,digits=5)
    fin_res<-cbind(names(entropy_gene),entropy_gene,adj_pvalue,chr,nsaml_perc,chr_pos_mut, mut_freq,mut_cds,mut_aa, vtype_mut,cod_mut)
    gene_symb<-vector()
    for(t in 1:nrow(fin_res))
    {
      gene_symb[t]<-tr_gn_table[which(tr_gn_table[,1]==rownames(fin_res)[t]),2]
    }
    canc_gene<-rep("unknown",length(gene_symb))
    pos_canc_gene<-which(gene_symb%in%canc_genes==T)
    if(length(pos_canc_gene)>0)
      canc_gene[pos_canc_gene]<-"Cancer gene"
    fin_res<-cbind(gene_symb,canc_gene,fin_res)
    if(signif==T)
    {
      pos_sig<-which(adj_pvalue<=0.05)
      
      if(length(pos_sig)==0)
        return("No significant genes")
      else
      {
        if(length(pos_sig)==1)
        {
          nam<-rownames(fin_res)[pos_sig]
          fin_res<-fin_res[pos_sig,]
          fin_res<-as.data.frame(t(fin_res),stringsAsFactors=F)
          rownames(fin_res)<-nam
        }
        if(length(pos_sig)>1)
          fin_res<-fin_res[pos_sig,]
      }
    }
   
    colnames(fin_res)<-c("Gene", "Canc_gene","Transcript","Entropy score","Adj_pvalue","Chr","Nsample%","Chr_pos","#Mut_freq","Mut_freq_cDNA","Mut_freq_AA","Type_mut","Orig/Mut_cod")
    fin_res[,4]<-round(as.numeric(fin_res[,4]),digits=4)
    fin_res[,7]<-round(as.numeric(fin_res[,7]),digits=4)
    fin_res<-fin_res[order(as.numeric(fin_res[,9]),decreasing=T),]
    return(fin_res)
    
  }
}
entr_randistr1<-function(n_mut,pr_length,real_entr,nperm)
{
  c<-1
  entr_mut<-vector()
  i<-1
  for(i in 1:nperm)
  {
    pos_mut<-sample(1:pr_length,n_mut,replace=T)
    if(length(unique(pos_mut))>1)
    {
      mut_pr_freq<-table(pos_mut)/length(pos_mut)
      entr_mut[c]<-(-1)*sum(mut_pr_freq*log2(mut_pr_freq))/log2(pr_length)
      c<-c+1
    }
  }
  summary(entr_mut)
  num_extr1<-(length(which(entr_mut<real_entr))+1)/(nperm+1)
  return(num_extr1)
}

concentration_cluster<-function(maf_file,pos_tr,pos_gn,pos_chrm,nsmpl_df,cds_pos,aa_pos,pos_codon,pos_barcode,pos_variant,name_mut,length_pr,thr_conc,thr_mut,rm_mut=5,perm_numb,cancer_gene,signif=T,leng_seq,chr_col)
{
  tr_gn_table<-unique(maf_file[,c(pos_tr,pos_gn)])
  range_vec<-vector()
  num_mut<-vector()
  mut_in_range<-vector()
  min_chrm<-vector()
  max_chrm<-vector()
  chr<-vector()
  mut_clust_perc<-vector()
  sample_clust_perc<-vector()
  n_nucl<-vector()
  sil_perc<-vector()
  mut_freq<-vector()
  mut_cds<-vector()
  mut_aa<-vector()
  vtype_mut<-vector()
  cod_mut<-vector()
  chr_pos_mut<-vector()
  p_value<-vector()
  final_res<-vector()
  nsample_save<-vector()
  nsample_save_perc<-vector()
  unique_tr<-as.character(unique(maf_file[,pos_tr]))
  c<-1
  for(i in 1:length(unique_tr))
  {
    
    info_tr<-maf_file[which(maf_file[,pos_tr]==unique_tr[i]),]
    nsample<-length(unique(info_tr[,pos_barcode]))
    if(nsample>=rm_mut)
    {
      pos_null<-which(info_tr[,cds_pos]=="")
      if(length(pos_null)>0)
        info_tr<-info_tr[-pos_null,]
      delins<-grep("delins",info_tr[,cds_pos])
      if(length(delins)>0)
        info_tr<-info_tr[-delins,]
      pos_rm<-mut_outliers(info_tr[,cds_pos])
      if(length(pos_rm)>0)
      {
        info_tr1<-info_tr[-pos_rm,]
      } else {
        info_tr1<-info_tr
      }
      mut_pos<-gsub("([A-Z.-]+)", "", info_tr1[,cds_pos])
      mut_pos<-as.numeric(gsub(">", "", mut_pos))
      
      if(length(unique(mut_pos))>1)
      {
        range_mut<-(max(mut_pos)-min(mut_pos))/length_pr[which(names(length_pr)==unique_tr[i])]
        range_mut<-range_mut*100
        n_mut<-length(mut_pos)
        mut_inarea<-nrow(info_tr1)/nrow(info_tr)*100
        if(range_mut<=thr_conc && mut_inarea>=thr_mut)
        {
          p_value[c]<-conc_randistr(n_mut,length_pr[which(names(length_pr)==unique_tr[i])],range_mut,perm_numb)
          range_vec[c]<-range_mut
          num_mut[c]<-paste(nrow(info_tr1)," out of ",nrow(info_tr))
          mut_in_range[c]<-nrow(info_tr1)/nrow(info_tr)*100
          
          n_nucl[c]<-max(mut_pos)-min(mut_pos)
          min_chrm[c]<-min(info_tr1[,pos_chrm])
          max_chrm[c]<-max(info_tr1[,pos_chrm])
          chr[c]<-unique(info_tr1[,chr_col])
          info_tr2<-info_tr1[info_tr1[,pos_chrm]<=max_chrm[c],]
          info_tr2<-info_tr2[info_tr2[,pos_chrm]>=min_chrm[c],]
          nsample_save[c]<-length(unique(info_tr2[,pos_barcode]))
          nsample_save_perc[c]<-nsample_save[c]/nsmpl_df*100
          freq_mut<-table(info_tr2[,cds_pos])
          pos_max<-which(freq_mut==max(freq_mut))[1]
          mut_freq[c]<-freq_mut[pos_max]
          mut_cds[c]<-unique(info_tr2[which(info_tr2[,cds_pos]==names(freq_mut)[pos_max]),cds_pos])
          mut_aa[c]<-unique(info_tr2[which(info_tr2[,cds_pos]==names(freq_mut)[pos_max]),aa_pos])
          chr_pos_mut[c]<-unique(info_tr2[which(info_tr2[,cds_pos]==names(freq_mut)[pos_max]),pos_chrm])
          vtype_mut[c]<-unique(info_tr2[which(info_tr2[,cds_pos]==names(freq_mut)[pos_max]),pos_variant])
          cod_mut[c]<-unique(info_tr2[which(info_tr2[,cds_pos]==names(freq_mut)[pos_max]),pos_codon])
          sil_freq_t<-table(info_tr2[,pos_variant])
          pos_var_sel<-which(names(sil_freq_t)==name_mut)
          if(length(pos_var_sel)==0)
            sil_perc[c]=0 else
            sil_perc[c]<-sil_freq_t[pos_var_sel]/sum(sil_freq_t)
          names(mut_in_range)[c]<-unique_tr[i]
          names(num_mut)[c]<-unique_tr[i]
          names(range_vec)[c]<-unique_tr[i]
          names(p_value)[c]<-unique_tr[i]
          names(nsample_save)[c]<-unique_tr[i]
          names(min_chrm)[c]<-unique_tr[i]
          names(max_chrm)[c]<-unique_tr[i]
          names(chr)[c]<-unique_tr[i]
          names(sil_perc)[c]<-unique_tr[i]
          names(n_nucl)[c]<-unique_tr[i]
          names(nsample_save_perc)[c]<-unique_tr[i]
          names(mut_freq)[c]<-unique_tr[i]
          names(mut_cds)[c]<-unique_tr[i]
          names(mut_aa)[c]<-unique_tr[i]
          names(vtype_mut)[c]<-unique_tr[i]
          names(cod_mut)[c]<-unique_tr[i]
          names(chr_pos_mut)[c]<-unique_tr[i]
          c<-c+1
        }
      }
    }
  }
  if(length(p_value)==0)
  {
    return("No transcripts passed the thresholds")
  }
  else
  {
    adj_pvalue<-p.adjust(p_value, method = "BH", n = length(p_value))    
    fin_res<-cbind(names(p_value),adj_pvalue,chr,min_chrm,max_chrm,n_nucl,range_vec,num_mut,mut_in_range,nsample_save,nsample_save_perc,chr_pos_mut, mut_freq,mut_cds,mut_aa, vtype_mut,cod_mut)
    gene_symb<-vector()
    for(t in 1:nrow(fin_res))
    {
      gene_symb[t]<-tr_gn_table[which(tr_gn_table[,1]==rownames(fin_res)[t]),2]
    }
    canc_gene<-rep("unknown",length(gene_symb))
    pos_canc_gene<-which(gene_symb%in%cancer_gene==T)
    if(length(pos_canc_gene)>0)
      canc_gene[pos_canc_gene]<-"Cancer gene"
    fin_res<-cbind(gene_symb,canc_gene,fin_res)
    if(signif==F)
      fin_res<-fin_res
    else
    {
      pos_sig<-which(adj_pvalue<=0.05)
      if(length(pos_sig)==0)
        return("No significant genes") else {
        if(length(pos_sig)==1)
        {
          nam<-rownames(fin_res)[pos_sig]
          fin_res<-fin_res[pos_sig,]
          fin_res<-as.data.frame(t(fin_res),stringsAsFactors=F)
          rownames(fin_res)<-nam
        }
        if(length(pos_sig)>1)
          fin_res<-fin_res[pos_sig,]
      }
    }
    colnames(fin_res)<-c("Gene","Canc_gene","Transcript","Adj_pvalue","Chr","Chr_min","Chr_max","Nnucl","Range_mut%","Nmut_cluster","Nmut_cluster%","Nsample","Nsample%","Chr_pos","#Mut_freq","Mut_freq_cDNA","Mut_freq_AA","Type_mut","Orig/Mut_cod")
    fin_res[,4]<-round( as.numeric(fin_res[,4]),digits=5)
    fin_res[,9]<-round( as.numeric(fin_res[,9]),digits=4)
    fin_res[,11]<-round( as.numeric(fin_res[,11]),digits=4)
    fin_res[,13]<-round( as.numeric(fin_res[,13]),digits=4)
    fin_res<-fin_res[order(as.numeric(fin_res[,15]),decreasing=T),]
    return(fin_res)
  }
}
conc_randistr<-function(n_mut,pr_length,real_val,nperm)
{
  c<-1
  range_mut<-vector()
  for(i in 1:nperm)
  {
    pos_mut<-sample(1:pr_length,n_mut,replace=T)
    pos_rm<-mut_outliers(pos_mut)
    if(length(pos_rm)!=0)
      pos_mut<-pos_mut[-pos_rm]
    if(length(unique(pos_mut))>1)
    {
      range_mut[c]<-((max(pos_mut)-min(pos_mut))/pr_length)*100
      c<-c+1
    }
  }
  length(range_mut)
  num_extr<-(length(which(range_mut<real_val))+1)/(nperm+1)
  return(num_extr)
}
window_clustering<-function(maf_info,pos_chrm,pos_chr,pos_transcript,pos_gn,pos_codon,pos_aa,cds_pos,wind_size,pos_barcode,pos_variant,name_mut,col_tr_name,thr_samples,tot_samples,thr_mut,tr_leng,nperm,cancer_gene)
{ 
  tr_gn_table<-unique(maf_info[,c(pos_transcript,pos_gn)])
  genes_name<-unique(maf_info[,pos_transcript])
  gene_list_mut<-list()
  n_tr<-NULL
  s<-1
  for(i in 1:length(genes_name))
  {
    gene_info<-maf_info[which(maf_info[,pos_transcript]==genes_name[i]),]
    wind_to_save<-list()
    n_mut<-nrow(gene_info)
    c<-1
    if(thr_samples<=length(unique(gene_info[,pos_barcode])))
    {
      n_tr<-c(n_tr,i)
      chr_pos<-unique(gene_info[,pos_chrm])
      j<-1
      for(j in 1:length(chr_pos))
      {
        
        pos_to_cons<-which(gene_info[,pos_chrm]==chr_pos[j])[1]
        wind<-subset(gene_info, gene_info[,pos_chrm]>=gene_info[pos_to_cons,pos_chrm] & gene_info[,pos_chrm]<=(as.numeric(gene_info[pos_to_cons,pos_chrm])+wind_size))
        wind_lim<-max(as.numeric(wind[,pos_chrm]))-min(as.numeric(wind[,pos_chrm]))
        freq_m<-nrow(wind)/n_mut*100
        freq_m_abs<-nrow(wind)
        freq_s_tr_abs<-length(unique(gene_info[,pos_barcode]))
        freq_s_tr_perc<-length(unique(gene_info[,pos_barcode]))/tot_samples*100
        freq_s_abs<-length(unique(wind[,pos_barcode]))
        freq_s<-freq_s_abs/tot_samples*100
        perc_syn<-table(wind[,pos_variant])
        pos_sil<-which(names(perc_syn)==name_mut)
        if(length(pos_sil)==0)
          sil_freq<-0
        else
          sil_freq<-perc_syn[pos_sil]/sum(perc_syn)
        if(freq_m>thr_mut)
        {
          if(c==1)
          {
            wind<-wind[,c(pos_chr,pos_chrm,pos_transcript,pos_variant)]
            wind<-wind[!duplicated(wind),]
            wind<-cbind(wind,freq_s,freq_m,freq_m_abs,freq_s_abs,freq_s_tr_abs,freq_s_tr_perc,n_mut,wind_lim,sil_freq)
            wind_to_save[[c]]<-wind
            c<-c+1
          }
          else
          {
            trial<-merge(wind_to_save[[c-1]],wind)
            if(nrow(trial)!=nrow(wind))
            {
              wind<-wind[,c(pos_chr,pos_chrm,pos_transcript,pos_variant)]
              wind<-wind[!duplicated(wind),]
              wind<-cbind(wind,freq_s,freq_m,freq_m_abs,freq_s_abs,freq_s_tr_abs,freq_s_tr_perc,n_mut,wind_lim,sil_freq)
              wind_to_save[[c]]<-wind
              c<-c+1
            }
          }
        }
      }
      if(length(wind_to_save)!=0)
      {
        gene_list_mut[[s]]<-wind_to_save
        s<-s+1
      }
    }
  }
  data_fr_res<-data.frame()
  if(length(gene_list_mut)==0)
    return(data_fr_res) else {
    list1 <- unlist(gene_list_mut, recursive = FALSE)
    data_fr_res<-do.call(rbind,list1)
    freq_m_col<-which(colnames(data_fr_res)=="freq_m")
    tr_col<-which(colnames(data_fr_res)==col_tr_name)
    n_mut_col<-which(colnames(data_fr_res)=="n_mut")
    wind_lim<-which(colnames(data_fr_res)=="wind_lim")
    sil_freq_n<-which(colnames(data_fr_res)=="sil_freq")
    pos_sm_tr_abs<-which(colnames(data_fr_res)=="freq_s_abs")
    pos_sm_tr_perc<-which(colnames(data_fr_res)=="freq_s")
    datafr_1<- data_fr_res[,c(freq_m_col,tr_col,n_mut_col,wind_lim,sil_freq_n,pos_sm_tr_abs,pos_sm_tr_perc)]
    data_fr2<-unique(datafr_1)
   
    p_value_hot<-vector()
    for(i in 1:nrow(data_fr2))
    {
      print(paste("Permutation hotspot",i," out of",nrow(data_fr2)))
      tr_l<-tr_leng[which(tr_leng[,1]==data_fr2[i,2]),3]
      p_value_hot[i]<-window_randistr(data_fr2[i,3],tr_l,data_fr2[i,1],data_fr2[i,4],nperm)
    }
    adj_pvalue<-p.adjust(p_value_hot, method = "BH", n = length(p_value_hot))
    datafr_2<-cbind(data_fr2,adj_pvalue)
    hot_procdata<-hotspot_processing(list(data_fr_res,datafr_2),cancer_gene,3,6,tr_gn_table)
    hot_agg<-hotspot_unit_results(hot_procdata)
    colnames(hot_agg)<- c("Ensembl_tr","chrm","min_chrm","max_chrm","%samples_hotspot","%Mutations_hotspot", "#mutations_hotspot",
                          "#samples_hotspot", "#samples_transcript","%samples_transcript","#Mutations_transcript","size_window","%Reference mutation","Cancer_gene","Symbol","pvalue")
    range_mutPerc<-vector()
    nmut_cluster<-vector()
    mut_freq<-vector()
    mut_cds<-vector()
    mut_aa<-vector()
    vtype_mut<-vector()
    cod_mut<-vector()
    chr_pos_mut<-vector()
    nsmpl_window<-vector()
    nsmpl_window_perc<-vector()
    for(i in 1:nrow(hot_agg))
    {
      maf_tr<-maf_info[which(maf_info[,pos_transcript]==hot_agg[i,1]),]
      maf_tr_wind<-maf_tr[which(maf_tr[,pos_chrm]>=hot_agg[i,3]),]
      maf_tr_wind<-maf_tr_wind[which(maf_tr_wind[,pos_chrm]<=hot_agg[i,4]),]
      range_mutPerc[i]<-nrow(maf_tr_wind)/nrow(maf_tr)*100
      nmut_cluster[i]<-nrow(maf_tr_wind)
      freq_mut<-table(maf_tr_wind[,cds_pos])
      pos_max<-which(freq_mut==max(freq_mut))[1]
      mut_freq[i]<-freq_mut[pos_max]
      mut_cds[i]<-unique(maf_tr_wind[which(maf_tr_wind[,cds_pos]==names(freq_mut)[pos_max]),cds_pos])
      mut_aa[i]<-unique(maf_tr_wind[which(maf_tr_wind[,cds_pos]==names(freq_mut)[pos_max]),pos_aa])
      chr_pos_mut[i]<-unique(maf_tr_wind[which(maf_tr_wind[,cds_pos]==names(freq_mut)[pos_max]),pos_chrm])
      vtype_mut[i]<-unique(maf_tr_wind[which(maf_tr_wind[,cds_pos]==names(freq_mut)[pos_max]),pos_variant])
      cod_mut[i]<-unique(maf_tr_wind[which(maf_tr_wind[,cds_pos]==names(freq_mut)[pos_max]),pos_codon])
    } 
     
     hot_agg<-cbind(hot_agg,range_mutPerc,mut_freq,mut_cds,mut_aa,chr_pos_mut,vtype_mut,cod_mut,as.numeric(hot_agg[,4])-as.numeric(hot_agg[,3]))
     hot_agg<-hot_agg[,c(15,14,1,16,2,3,4,23,17,7,6,8,5,21,18:20,22,23)]
     colnames(hot_agg)<-c("Gene","Canc_gene","Transcript","Adj_pvalue","Chr","Chr_min","Chr_max",	"Nnucl","Range_mut%","Nmut_cluster","Nmut_cluster%","Nsample","Nsample%","Chr_mut",	"#Mut_freq","Mut_freq_cDNA","Mut_freq_AA","Type_mut", "Orig/Mut_cod")
     hot_agg[,4]<-round( as.numeric(hot_agg[,4]),digits=5)
     hot_agg[,9]<-round( as.numeric(hot_agg[,9]),digits=4)
     hot_agg[,11]<-round( as.numeric(hot_agg[,11]),digits=4)
     hot_agg[,13]<-round( as.numeric(hot_agg[,13]),digits=4)
     hot_agg<-hot_agg[order(hot_agg[,15],decreasing=T),]
     return(list(data_fr_res,datafr_2,hot_agg))
  }
}

hotspot_processing<-function(hotspot_out,cancer_gene,tr_pos=3,pos_fm=6,tr_gn_table)
{
  sign_hotspots<-hotspot_out[[2]][which(hotspot_out[[2]][,8]<=0.05),]
  sign_hotmaf<-hotspot_out[[1]][hotspot_out[[1]][,tr_pos]%in%sign_hotspots[,2],]
  p_value_adapted<-vector()
  list_res<-list()
  for(i in 1:nrow(hotspot_out[[2]]))
  {
    list_res[[i]]<-subset(hotspot_out[[1]],hotspot_out[[1]][,tr_pos]==hotspot_out[[2]][i,2] & hotspot_out[[1]][,pos_fm]==hotspot_out[[2]][i,1])
    list_res[[i]]<-cbind(list_res[[i]],rep(hotspot_out[[2]][i,8],nrow(list_res[[i]])))
    colnames(list_res[[i]])[ncol(list_res[[i]])]<-"adj_pvalue"
  }
  data_fr_res<-do.call(rbind,list_res)
  pos_sign<-which(data_fr_res[,ncol(list_res[[i]])]<=0.05)
  data_fr_res_sign<-data_fr_res[pos_sign,]
  colnames(data_fr_res_sign)
  sign_genes<-vector()
  for(i in 1:nrow(data_fr_res_sign))
  {
    sign_genes[i]<-unique(tr_gn_table[which(tr_gn_table[,1]==data_fr_res_sign[i,tr_pos]),2])
  }
  canc_gene<-rep("unknown",length(sign_genes))
  pos_canc_gene<-which(sign_genes%in%cancer_gene==T)
  if(length(pos_canc_gene)>0)
    canc_gene[pos_canc_gene]<-"Cancer gene"
  data_fr_res_sign<-cbind( data_fr_res_sign,canc_gene,sign_genes)
  colnames(data_fr_res_sign)[ncol(data_fr_res_sign)-1]<-"Cancer_gene"
  colnames(data_fr_res_sign)[ncol(data_fr_res_sign)]<-"Symbol_name"
  return(data_fr_res_sign)
}
hotspot_unit_results<-function(res_hotspot)
{
  library(dplyr)
  c<-1
  hotspot<-data.frame()
  numb_mut<-nrow(res_hotspot)
  colnames(res_hotspot)
  while(numb_mut>0)
  {
    tr<-res_hotspot[1,3]
    res_hotspot1<-res_hotspot[which(res_hotspot[,3]==tr),]
    res_hotspot1<-res_hotspot1[order(res_hotspot1[,2],decreasing=F),]
    pos_mut<-as.numeric(res_hotspot1[1,2])
    wind_sz<-as.numeric(res_hotspot1[1,12])+1
    sub_set<-subset(res_hotspot1,res_hotspot1[1,2]>=pos_mut,res_hotspot1[1,2]<=pos_mut+wind_sz)
    min_chrm<-pos_mut
    max_chrm<-pos_mut+wind_sz
    chrm<-res_hotspot1[1,1]
    genesymb<-as.character(res_hotspot1[1,16])
    cancer<-as.character(res_hotspot1[1,15])
    freq_s<-res_hotspot1[1,5:13]
    adj_pvalue<-sub_set[which(sub_set[,12]==max(sub_set[,12]))[1],14]
    hot_res<-cbind(tr,chrm,min_chrm,max_chrm,freq_s,cancer,genesymb,adj_pvalue,stringsAsFactors=F)
    hotspot<-rbind(hotspot,hot_res,stringsAsFactors=F)
    res_hotspot<-anti_join(res_hotspot,sub_set)
    numb_mut<-nrow(res_hotspot)
    c<-c+1
  }
  return(hotspot)
}

window_randistr<-function(n_mut,gene_length,thrs_hot,wind_size,nperm)
{
  c<-1
  hot_pos<-vector()
  hot_perm<-vector()
  i<-1
  for(i in 1:nperm)
  {
    pos_mut<-sample(1:gene_length,n_mut,replace=T)
    pos_mut<-pos_mut[order(pos_mut,decreasing = F)]
    tab_mut<-table(pos_mut)
    tab_mut<-cbind(tab_mut,as.numeric(names(tab_mut)))
    j<-1
    for(j in 1:nrow(tab_mut))
    {
      wind<-subset(tab_mut[,2], tab_mut[,2]>=tab_mut[j,2] & tab_mut[,2]<=(tab_mut[j,2]+wind_size))
      freq_m<-sum(tab_mut[tab_mut[,2]%in%wind,1])/n_mut*100
      if(freq_m>=thrs_hot)
      {
        hot_pos[j]<-1
      }
      else
      {
        hot_pos[j]<-0
      }
    }
    if(sum(hot_pos)>0)
      hot_perm[i]=1
    else
      hot_perm[i]=0
  }
  num_extr1<-(length(which(hot_perm==1))+1)/(nperm+1)
  return(num_extr1)
}
mut_outliers<-function(amino_pos)
{
  mut_pos<-gsub("([A-Z.-]+)", "", amino_pos)
  mut_pos1<-as.numeric(gsub(">", "", mut_pos))
  pos_na<-which(is.na(mut_pos1)==T)
  lowerq = quantile(mut_pos1,na.rm=T)[2]
  upperq = quantile(mut_pos1,na.rm=T)[4]
  iqr = upperq - lowerq
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  if(length(pos_na)>0)
  {
    pos_rm<-which(mut_pos1>extreme.threshold.upper)
    pos_rm<-c(pos_rm,which(mut_pos1<extreme.threshold.lower))
    pos_rm<-c(pos_rm,pos_na)
    return(unique(pos_rm))
  }
  else
  {
    pos_rm<-which(mut_pos1>extreme.threshold.upper)
    pos_rm<-c(pos_rm,which(mut_pos1<extreme.threshold.lower))
    return(pos_rm)
  }
  
}
