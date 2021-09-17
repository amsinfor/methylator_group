#-----------------for SDC2 and TFPI2 subtyping in CRC samples: TCGA+GEO+custom
library(RColorBrewer)
library(data.table)
library(bnstruct)
library(parallel);
library(pheatmap)
library(pROC);
library(glmnet);
library(ggbiplot)
library(ggpubr);
library(limma);
library("grid")
library("gridExtra")
library(survival)
library(KMsurv)
library(survivalROC)
library(survminer)
library(pheatmap)
library(survcomp)
library(maftools)
library(scatterplot3d)
#-------
prepare_series_matrix_data<-function(infile){
	readLines(gzfile(infile,'r'))->infile_lines;
	grep("series_matrix_table_begin",infile_lines)->skip_lines;
	read.table(infile,skip=skip_lines,sep="\t",header=T,stringsAsFactors=F,nrows=length(infile_lines)-skip_lines-2)->myd;
	"DATA"->names(myd)[1];
	return(myd);
}
process_series_matrix<-function(series_file,infors){
	readLines(series_file)->myd_info;
	c()->res;
	c()->res_colnames;
	for(s in infors){
		myd_info[grep(s,myd_info)]->s_lines;
		for(sl in s_lines){
			unlist(strsplit(sl,split="\t"))->s_lines.split;
			gsub("\"","",s_lines.split)->s_lines.split;
			if(grepl(": ",sl)){	
				unlist(strsplit(s_lines.split[2],split=": "))[1]->s;
				unlist(lapply(s_lines.split,function(x){
					unlist(strsplit(x,split=": "))[2]
				}))->x_res;
				c(res,x_res[-1])->res;
			}else if(grepl("ftp",sl)){
				unlist(lapply(s_lines.split,function(x){
					gsub("\\.gz","",basename(x));
				}))->x_basenames;
				c(res,x_basenames[-1])->res;
			}else{
				c(res,s_lines.split[-1])->res;
			}
			c(res_colnames,s)->res_colnames;
		}
	}
	matrix(res,ncol=length(res_colnames),byrow=F)->res.matrix;
	res_colnames->colnames(res.matrix);
	as.data.frame(res.matrix,stringsAsFactors=F)->res.matrix
	return(res.matrix);
}
get_HumanMethylation450_annotats<-function(probes,db){
	db[db$Name%in%probes,]->probe_annots;
	setdiff(probes,probe_annots$Name)->diff_probes;
	#--------
	if(length(diff_probes)>0){
		c()->diff_probes_mat;
		for(dp in diff_probes){
			c(dp,rep(NA,ncol(probe_annots)-1))->dp_;
			c(diff_probes_mat,dp_)->diff_probes_mat;
		}
		matrix(diff_probes_mat,nrow=length(diff_probes),byrow=T)->diff_probes_mat;
		colnames(probe_annots)->colnames(diff_probes_mat);
		rbind(probe_annots,diff_probes_mat)->probe_annots;
	}
	#----------
	return(probe_annots);
}
map_HumanMethylation450_annotats<-function(expd,db){
	as.character(expd$DATA)->probes
	as.character(db$UCSC_RefGene_Name)->gnames;
	db$Name->names(gnames);
	gnames[probes]->expd_probes_symbol;
	#--------
	which(expd_probes_symbol=="")->na_index;
	probes[na_index]->expd_probes_symbol[na_index]
	data.frame("Probe"=probes,"gName"=expd_probes_symbol,expd[,-1],stringsAsFactors=F)->expd_t;
	which(is.na(expd_t$gName))->na_index;
	expd_t$Probe[na_index]->expd_t$gName[na_index];
	#-------
	grep(";",expd_t$gName)->gName_index;
	lapply(gName_index,function(x){
		unlist(strsplit(as.character(expd_t$gName[x]),split=";"))->x_split;
		unique(x_split)->x_split;
		rep(x,length(x_split));
	})->gName_index_res;
	unlist(gName_index_res)->gName_index_res;
	#-----
	lapply(gName_index,function(x){
		unlist(strsplit(as.character(expd_t$gName[x]),split=";"))->x_split;
		unique(x_split);
	})->gName_split_res;
	unlist(gName_split_res)->gName_split_res;
	#-----
	expd_t[-gName_index,]->expd_t1;
	expd_t[gName_index_res,]->expd_t2;
	gName_split_res->expd_t2$gName;
	rbind(expd_t1,expd_t2,stringsAsFactors=F)->expd_t;
	#----------#------------na
	return(expd_t);
}
preprocess_data_frame<-function(dat,cutoff){
	#------remove NA values;
	detectCores()->no_cores;
	makeCluster(no_cores-1)->c1;
	clusterExport(c1,c("dat","cutoff"),envir=environment());
	parLapply(c1,1:nrow(dat),function(i){
		which(is.na(dat[i,-1]))->na_index;
		if(length(na_index)>0){
			if(length(na_index)/(ncol(dat)-1) < cutoff){
				i;
			}
		}else{
			i;
		}
	})->keep_rows;
	stopCluster(c1);
	dat[unlist(keep_rows),]->dat;
	#----
	t(dat[,-1])->dat.t;
	knn.impute(dat.t,k=round(sqrt(ncol(dat)-1)),cat.var=1:ncol(dat.t),to.impute=1:nrow(dat.t),using=1:nrow(dat.t))->dat.t;
	data.frame("DATA"=dat[,1],t(dat.t),stringsAsFactors=F)->dat;
	#---
	return(dat);
}
merge_values<-function(methyd,probe_index){
	paste(methyd$DATA[probe_index],collapse="_")->probe_merged;
	apply(methyd[probe_index,-1],2,function(x){mean(as.numeric(x),na.rm=T)})->methyd_merged;
	data.frame(probe_merged,matrix(methyd_merged,nrow=1),stringsAsFactors=F)->methyd_merged_df;
	c("DATA",colnames(methyd)[-1])->colnames(methyd_merged_df)
	return(methyd_merged_df);
}
merge_methy_probes<-function(methyd,methyd_probes){
	if(is.list(methyd_probes)){
		if(length(methyd_probes)>=100){
			detectCores()->no_cores;
			makeCluster(no_cores-1)->c1;
			clusterExport(c1,c("methyd","merge_values"),envir=environment());
			parLapply(c1,methyd_probes,function(mp){
				which(methyd$DATA%in%mp)->probe_index;
				#-----
				merge_values(methyd,probe_index)->methyd_merged_df;
			})->merged_list;
			stopCluster(c1);
			#---------split cluster:
			makeCluster(no_cores-1)->c1;
			clusterSplit(c1,merged_list)->merged_list.split;
			parLapply(c1,merged_list.split,function(x){
				x[[1]]->x_df;
				for(xi in 2:length(x)){
					if(!is.null(xi)){
						rbind(x_df,x[[xi]])->x_df;
					}
				}
				x_df;
			})->merged_list.merged;
			stopCluster(c1);
			#----------merge 
			for(lx in merged_list.merged){
				rbind(lx,methyd)->methyd;
			}
		}else{
			for(mp in methyd_probes){
			which(methyd$DATA%in%mp)->probe_index;
			#-----
			merge_values(methyd,probe_index)->methyd_merged_df;
			rbind(methyd_merged_df,methyd)->methyd;
			}
			rev(names(methyd_probes))->methyd$DATA[1:length(methyd_probes)];
		}
	}else{
		which(methyd$DATA%in%methyd_probes)->probe_index;
		#-----
		merge_values(methyd,probe_index)->methyd_merged_df;
		rbind(methyd_merged_df,methyd)->methyd;
	}
	return(methyd);
}
process_TCGA_pancancer_data<-function(expd,clinical_d){
	colnames(expd)[-c(1,2)]->expd_samples;
	unlist(lapply(expd_samples,function(x){
		unlist(strsplit(x,split="-"))->x_split;
		paste(x_split[1:3],collapse="-")->x_sample;
		x_split[4]->x_sample_code;
		c(x_sample,x_sample_code);
	}))->expd_samples;
	data.frame("A0_Samples"=expd_samples[seq(1,length(expd_samples),2)],"SampleTypeCode"=expd_samples[seq(2,length(expd_samples),2)],stringsAsFactors=F)->expd_sample_df;
	#-----------------diff samples 
	setdiff(expd_sample_df$A0_Samples,clinical_d$A0_Samples)->diff_samples;
	#----------------------shared samples:
	clinical_d$A0_Samples->rownames(clinical_d);
	clinical_d[as.character(expd_sample_df$A0_Samples),]->new_clinical_d;
	expd_sample_df$A0_Samples->new_clinical_d$SampleID;
	expd_sample_df$SampleTypeCode->new_clinical_d$SampleTypeCode;
	paste(new_clinical_d$A0_Samples,new_clinical_d$SampleTypeCode,sep="-")->new_clinical_d$A0_Samples;
	unique(new_clinical_d)->new_clinical_d;
	#------add diff samples :
	expd_sample_df[which(expd_sample_df$A0_Samples%in%diff_samples),]->diff_samples;
	for(xi in 1:nrow(diff_samples)){
		diff_samples$A0_Samples[xi]->ds;
		diff_samples$SampleTypeCode[xi]->ds_split_code;
		paste(ds,ds_split_code,sep="-")->ds;
		c(ds,rep(NA,ncol(new_clinical_d)-2),ds_split_code)->ds_res;
		rbind(new_clinical_d,ds_res)->new_clinical_d;
	}
	##---------------
	new("SampleObj",DATA=expd,CliniInfo=new_clinical_d)->new_clinical_d_objRef;
	return(new_clinical_d_objRef);
}
change_GEObjRef_to_expd<-function(objRef,select_f="gName"){
	objRef$DATA->objRef_dat;
	objRef$CliniInfo->objRef_info;
	t(objRef_dat[,-c(1,2)])->dat_;
	if(select_f=="gName"){
		objRef_dat$gName->colnames(dat_);
	}else{
		objRef_dat$Probe->colnames(dat_);
	}
	data.frame("SampleID"=rownames(dat_),dat_)->dat_;
	merge(objRef_info,dat_,by.x="A0_Samples",by.y="SampleID")->res;
	return(res);
}
change_values<-function(myd,column,values){
	which(is.na(myd[,column]))->na_index;
	if(length(na_index)>0){
		myd[na_index,]->myd_na;
		as.character(myd_na[,column])->myd_na[,column];
		"Un"->myd_na[is.na(myd_na[,column]),column];
	}else{
		NULL->myd_na;
	}
	myd[!is.na(myd[,column]),]->myd
	colnames(myd)[column]->column_name;
	table(myd[,column])->N_stage.table;
	as.factor(names(N_stage.table))->N_stage.names;
	data.frame(column_name=N_stage.names,"Value"=values)->N_stage.df;
	c()->tmp.value;
	for(i in 1:nrow(myd)){
		for(j in 1:nrow(N_stage.df)){
			if(myd[i,column]==N_stage.df[j,1]){
				c(tmp.value,as.character(N_stage.df[j,2]))->tmp.value;
			}
		}
	}
	tmp.value->myd[,column];
	if(!is.null(myd_na)){
		rbind(myd,myd_na)->myd;
	}
	return(myd);
}
draw_sample_distance_matrix<-function(expd,groups,myd_colors){
	cor(expd)->expd_corrmatrix;
	data.frame("Condition"=groups$Condition)->expd_condition;
	groups$SampleID->rownames(expd_condition);
	#--
	unique(groups$Condition)->group_conditions;
	myd_colors[1:length(group_conditions)]->group_colors;
	group_conditions->names(group_colors);
	list("Condition"=group_colors)->annotat_cols;
	pheatmap(expd_corrmatrix,annotation_col=expd_condition,annotation_row=expd_condition,border_color=NA,annotation_colors=annotat_cols,show_colnames=F,show_rownames=F);
}
#--------
#---objects---
setRefClass("SampleObj",
	fields=list(DATA="data.frame",CliniInfo="data.frame"),
	methods=list(
		initialize=function(DATA,CliniInfo){
			intersect(colnames(DATA),CliniInfo$A0_Samples)->shared_samples;
			CliniInfo[which(CliniInfo$A0_Samples%in%shared_samples),]->>CliniInfo;
			DATA[,c("gName","Probe",shared_samples)]->>DATA;
		},
		getSYMBOL=function(probes){
			unlist(lapply(probes,function(i){
				which(DATA$Probe==i)->i_index
			}))->probes_index;
			return(as.character(DATA$gName[probes_index]))
		},
		getConditionGroup=function(f,groups){
			c()->groups_index;
			for(g in groups){
				which(CliniInfo[,f]==g)->g_index;
				c(groups_index,g_index)->groups_index;
			}
			CliniInfo[groups_index,c("A0_Samples",f)]->condition_groups;
			c("SampleID","Condition")->colnames(condition_groups);
			return(condition_groups);
		},
		getGroupMatrix=function(f,groups){
			getConditionGroup(f,groups)->condition_groups;
			#----
			as.matrix(DATA[,as.character(condition_groups$SampleID)])->dat_matrix;
			paste(DATA$gName,DATA$Probe,sep="|")->rownames(dat_matrix);
			return(dat_matrix);
		},
		getDataMatrix=function(genes,select_f="gene"){
			as.matrix(DATA[,-c(1,2)])->dat_matrix;
			paste(DATA$gName,DATA$Probe,sep="|")->rownames(dat_matrix);
			if(select_f=="gene"){
				which(DATA$gName%in%genes)->g_index;
			}else if(select_f=="probe"){
				which(DATA$Probe%in%genes)->g_index;
			}
			dat_matrix[g_index,]->dat_matrix;
			return(dat_matrix);
		},
		getSubSet=function(f,groups){
			c()->groups_index;
			for(g in groups){
				which(CliniInfo[,f]==g)->g_index;
				c(groups_index,g_index)->groups_index;
			}
			CliniInfo[groups_index,]->sub_CliniInfo;
			as.character(sub_CliniInfo[,f])->sub_CliniInfo[,f];
			#---
			colnames(DATA)[1:2]->two_cols;
			DATA[,c(two_cols,as.character(sub_CliniInfo$A0_Samples))]->sub_DATA;
			new("SampleObj",DATA=sub_DATA,CliniInfo=sub_CliniInfo)->sub_SampleObj;
			return(sub_SampleObj);
		},
		addFeatures=function(feature_df){
			merge(CliniInfo,feature_df,by.x="A0_Samples",by.y="A0_Samples",all.x=T)->>CliniInfo;
		}
		
))->SampleObj;

#######################################################################################################
#------------identify the probes in the promoters of SDC2 and TFPI2  
#------------0-2kb upstream of SDC2 and TFPI2
#---SDC2 
TCGA_CRC_methy_objRef$getDataMatrix(SDC2_450k_probes,select_f="probe")->SDC2_450k_methyd;
calculate_gene_groupScores(SDC2_450k_methyd,TCGA_CRC_methy_TN_group,rownames(SDC2_450k_methyd))->SDC2_methy_TN_Scores;
unlist(lapply(SDC2_methy_TN_Scores$gName,function(gx){unlist(strsplit(gx,split="\\|"))[2]}))->SDC2_methy_TN_Scores$Probe;
merge(SDC2_450k_annots,SDC2_methy_TN_Scores,by.x="IlmnID",by.y="Probe")->SDC2_methy_TN_Scores;
#---TFPI2: 
TCGA_CRC_methy_objRef$getDataMatrix(TFPI2_450k_probes,select_f="probe")->TFPI2_450k_methyd;
calculate_gene_groupScores(TFPI2_450k_methyd,TCGA_CRC_methy_TN_group,rownames(TFPI2_450k_methyd))->TFPI2_methy_TN_Scores;
unlist(lapply(TFPI2_methy_TN_Scores$gName,function(gx){unlist(strsplit(gx,split="\\|"))[2]}))->TFPI2_methy_TN_Scores$Probe;
merge(TFPI2_450k_annots,TFPI2_methy_TN_Scores,by.x="IlmnID",by.y="Probe")->TFPI2_methy_TN_Scores;
#---4 and 7 probes:
c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")[c(1,4,10,11)])->myd_colors;
c(get_palette("npg",10),brewer.pal(8,"Dark2"))->myd_colors_npg;
c("SDC2","TFPI2")->target_genes;
c("cg16935295","cg04261408","cg14625631","cg10292139")->SDC2_Probes
c("cg12973591","cg14377593","cg17338208","cg22441533","cg22799321","cg24531255","cg26739865")->TFPI2_Probes;
"cg04261408_cg10292139_cg14625631_cg16935295"->SDC2_P;
"cg12973591_cg14377593_cg17338208_cg22441533_cg22799321_cg24531255_cg26739865"->TFPI2_P;
c("cg16935295","cg04261408","cg14625631","cg10292139","cg12973591","cg14377593","cg17338208","cg22441533","cg22799321","cg24531255","cg26739865")->target_probes;
list("SDC2_P"=SDC2_Probes,"TFPI2_P"=TFPI2_Probes)->target_gene_probes
#-----------------------------------------for TCGA CRC data:
read.table("TCGA_CRC_clinical-infor.txt",header=T,sep="\t",stringsAsFactors=F)->TCGA_CRC_clinical_info
fread("CRC_methy.txt",header=T,sep="\t",stringsAsFactors=F)->TCGA_CRC_methy;
as.data.frame(TCGA_CRC_methy)->TCGA_CRC_methy;
"DATA"->names(TCGA_CRC_methy)[1]
merge_methy_probes(TCGA_CRC_methy,target_gene_probes)->TCGA_CRC_methy;
map_HumanMethylation450_annotats(TCGA_CRC_methy,HumanMethylation450_annotats)->TCGA_CRC_methy;
gsub("\\.","-",colnames(TCGA_CRC_methy))->colnames(TCGA_CRC_methy)
#------
process_TCGA_pancancer_data(TCGA_CRC_methy,TCGA_CRC_clinical_info)->TCGA_CRC_methy_objRef;
"Adjacent"->TCGA_CRC_methy_objRef$CliniInfo$ColonSite[which(TCGA_CRC_methy_objRef$CliniInfo$SampleTypeCode=="11")]
change_values(TCGA_CRC_methy_objRef$CliniInfo,13,c("Adjacent","Right Colon","Right Colon","Left Colon","Others","Right Colon","Left Colon","Rectum","Left Colon","Others","Others"))->TCGA_CRC_methy_objRef$CliniInfo;
TCGA_CRC_methy_objRef$getSubSet("SampleTypeCode",c("11","01"))->TCGA_CRC_methy_objRef;
"Others"->TCGA_CRC_methy_objRef$CliniInfo$ColonSite[which(TCGA_CRC_methy_objRef$CliniInfo$ColonSite=="Un")]
factor(TCGA_CRC_methy_objRef$CliniInfo$SampleTypeCode,levels=c("11","01"))->TCGA_CRC_methy_objRef$CliniInfo$SampleTypeCode;
#---------------------------for GSE48684: normal vs colorectal adenomas (CRA) vs CRC
as.data.frame(GSE48684_methyd)->GSE48684_methyd;
preprocess_data_frame(GSE48684_methyd,0.5)->GSE48684_methyd;
merge_methy_probes(GSE48684_methyd,target_gene_probes)->GSE48684_methyd;
map_HumanMethylation450_annotats(GSE48684_methyd,HumanMethylation450_annotats)->GSE48684_methyd;
new("SampleObj",DATA=GSE48684_methyd,CliniInfo=GSE48684_targets)->GSE48684_methyd_objRef;
#--------------------------------for GSE79740 methylation 
fread("MergeExpro_contrib1-GSE79740.txt",header=T,sep="\t",stringsAsFactors=F)->GSE79740_methyd;
as.data.frame(GSE79740_methyd)->GSE79740_methyd;
preprocess_data_frame(GSE79740_methyd,0.5)->GSE79740_methyd;
merge_methy_probes(GSE79740_methyd,target_gene_probes)->GSE79740_methyd;
map_HumanMethylation450_annotats(GSE79740_methyd,HumanMethylation450_annotats)->GSE79740_methyd;
#-----
new("SampleObj",DATA=GSE79740_methyd,CliniInfo=GSE79740_targets)->GSE79740_methyd_objRef;
##########################################################################################################################
#--------------------------merge two dataset: 
"GSE48684"->GSE48684_methyd_objRef$CliniInfo$DataSet;
GSE48684_methyd_objRef$CliniInfo$DiseaseStatus->GSE48684_methyd_objRef$CliniInfo$SampleType;
"ColonSite"->names(GSE48684_methyd_objRef$CliniInfo)[8];
"GSE79740"->GSE79740_methyd_objRef$CliniInfo$DataSet;
"Left Colon"->GSE79740_methyd_objRef$CliniInfo$ColonSite;
list("GSE48684"=GSE48684_methyd_objRef,"GSE79740"=GSE79740_methyd_objRef)->GSE48684_GSE79740_methy_listd;
batch_correction_affymetrix(GSE48684_GSE79740_methy_listd)->GSE48684_GSE79740_methy_merged;
map_HumanMethylation450_annotats(GSE48684_GSE79740_methy_merged,HumanMethylation450_annotats)->GSE48684_GSE79740_methy_merged
#gsub("\\.","-",colnames(GSE48684_GSE79740_methy_merged))->colnames(GSE48684_GSE79740_methy_merged);
#------merge clinical-infor :
clincal_infor_merged(GSE48684_GSE79740_methy_listd)->GSE48684_GSE79740_clinical_merged
change_values(GSE48684_GSE79740_clinical_merged,2,c("Adenoma","Cancer","Normal","Normal","Normal","Cancer"))->GSE48684_GSE79740_clinical_merged;
change_values(GSE48684_GSE79740_clinical_merged,3,c("Normal","Left Colon","Left Colon","Left Colon","Left Colon","Right Colon","Rectum Colon","Rectum Colon","Right Colon","Right Colon","Left Colon","Others"))->GSE48684_GSE79740_clinical_merged;
new("SampleObj",DATA=GSE48684_GSE79740_methy_merged,CliniInfo=GSE48684_GSE79740_clinical_merged)->GSE48684_GSE79740_methy_objRef
#------------------------check batch correction : distance matrix and PCA;
GSE48684_GSE79740_methy_objRef$getConditionGroup("SampleType",c("Normal","Adenoma","Cancer"))->GSE48684_GSE79740_group;
GSE48684_GSE79740_methy_objRef$getGroupMatrix("SampleType",c("Normal","Adenoma","Cancer"))->GSE48684_GSE79740_matrix;

#################################################################################################################################################
#------------correlation between expression and methylation 
merge(TCGA_CRC_methy_tumor_factor[,c("A0_Samples","SDC2_P","TFPI2_P","DoubleType")],TCGA_CRC_exp_processed_factor[,c("A0_Samples","SDC2","TFPI2")])->TCGA_CRC_SDC2_TFPI2_methy_exp_factor;
ggscatter(TCGA_CRC_SDC2_TFPI2_methy_exp_factor,x="SDC2_P",y="SDC2",add="reg.line",,conf.int = TRUE,palette="npg")+stat_cor(method = "pearson")->p1#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")
ggscatter(TCGA_CRC_SDC2_TFPI2_methy_exp_factor,x="TFPI2_P",y="TFPI2",add="reg.line",,conf.int = TRUE,palette="npg")+stat_cor(method = "pearson")->p2#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")
grid.arrange(p1,p2,ncol=2)

####################################################################################################################################
#############################################################################################################################################
#------------------------------------for TCGA methylation ROC : SDC2 + TFPI2 
change_GEObjRef_to_expd(TCGA_CRC_methy_objRef,select_f="probe")->TCGA_CRC_methy_factor;
TCGA_CRC_methy_factor[,c(colnames(TCGA_CRC_methy_objRef$CliniInfo),c(target_probes,"SDC2_P","TFPI2_P"))]->TCGA_CRC_methy_factor;
add_Y_labels(TCGA_CRC_methy_factor,"SampleTypeCode",c("11","01"))->TCGA_CRC_methy_factor;
#---for SDC2_Probes:
par(mfrow=c(2,2),mar=c(4,4,2,1))
do_logistic_fit(TCGA_CRC_methy_factor,"Y_label","SDC2_P")->SDC2_logit_res;
#---for TFPI2_Probes:
do_logistic_fit(TCGA_CRC_methy_factor,"Y_label","TFPI2_P")->TFPI2_logit_res;
#----for SDC2_Probes + TFPI2_Probes:
do_logistic_fit(TCGA_CRC_methy_factor,"Y_label","SDC2_P+TFPI2_P")->SDC2_TFPI2_logit_res;
#--add single gene cutoff :
get_x_cutoff(coef(SDC2_logit_res$model),coords(SDC2_logit_res$roc,"b", ret = "t", best.method = "youden"))->g1_best_cutoff;
get_x_cutoff(coef(TFPI2_logit_res$model),coords(TFPI2_logit_res$roc,"b", ret = "t", best.method = "youden"))->g2_best_cutoff;

####################################################################################################################################
#---------for GEO 
change_GEObjRef_to_expd(GSE48684_GSE79740_methy_objRef$getSubSet("SampleType",c("Normal","Cancer")),select_f="probe")->GSE48684_GSE79740_methy_factor;
GSE48684_GSE79740_methy_factor[,c(colnames(GSE48684_GSE79740_methy_objRef$CliniInfo),c(target_probes,"SDC2_P","TFPI2_P"))]->GSE48684_GSE79740_methy_factor;
add_Y_labels(GSE48684_GSE79740_methy_factor,"SampleType",c("Normal","Cancer"))->GSE48684_GSE79740_methy_factor;
#---for SDC2_Probes:
par(mfrow=c(2,2),mar=c(4,4,2,1))
do_logistic_fit(GSE48684_GSE79740_methy_factor,"Y_label","SDC2_P")->GEO_SDC2_logit_res;
#---for TFPI2_Probes:
do_logistic_fit(GSE48684_GSE79740_methy_factor,"Y_label","TFPI2_P")->GEO_TFPI2_logit_res;
#----for SDC2_Probes + TFPI2_Probes:
do_logistic_fit(GSE48684_GSE79740_methy_factor,"Y_label","SDC2_P+TFPI2_P")->GEO_SDC2_TFPI2_logit_res;
#--add single gene cutoff :
get_x_cutoff(coef(GEO_SDC2_logit_res$model),coords(GEO_SDC2_logit_res$roc,"b", ret = "t", best.method = "youden"))->GEO_g1_best_cutoff;
get_x_cutoff(coef(GEO_TFPI2_logit_res$model),coords(GEO_TFPI2_logit_res$roc,"b", ret = "t", best.method = "youden"))->GEO_g2_best_cutoff;
#----------two groups
prepare_double_genes_types_cutoff(GSE48684_GSE79740_methy_factor,"SDC2_P","TFPI2_P",0.2,0.2)->GSE48684_GSE79740_methy_factor;
GSE48684_GSE79740_methy_factor[GSE48684_GSE79740_methy_factor$SampleType!="Normal",]->GSE48684_GSE79740_methy_tumor_factor
GSE48684_GSE79740_methy_tumor_factor$DoubleType->GSE48684_GSE79740_methy_tumor_factor$DoubleTypeGroup;
change_values(GSE48684_GSE79740_methy_tumor_factor,20,c("HH","HL","HL","LL"))->GSE48684_GSE79740_methy_tumor_factor;

####################################################################################################################################
#############################################################################################################################################
#-------------------------------for custom MSP Ct values:
read.table("CRC_tissue_ct.txt",sep="\t",header=T,stringsAsFactors=F,encoding="UTF-8")->CRC_tissue_ct;
read.xlsx("CRC_tissue_clinical_infor.xlsx", sheetName="Sheet1", header=T, stringsAsFactors=F, encoding="UTF-8")->CRC_tissue_clinical_infor;
gsub("\t","",CRC_tissue_clinical_infor$A0_Samples)->CRC_tissue_clinical_infor$A0_Samples
new("SampleObj",DATA=CRC_tissue_ct,CliniInfo=CRC_tissue_clinical_infor)->CRC_tissue_ct_objRef;
change_GEObjRef_to_expd(CRC_tissue_ct_objRef,select_f="probe")->CRC_tissue_ct_factor;
"Normal"->CRC_tissue_ct_factor$ColonSite[which(CRC_tissue_ct_factor$SampleType=="Normal")]
change_values(CRC_tissue_ct_factor,11,c("Right Colon","Others","Left Colon","Left Colon","Normal","Others","Rectum","Right Colon","Left Colon"))->CRC_tissue_ct_factor;
change_values(CRC_tissue_ct_factor,5,c("Female","Male","Male","Female"))->CRC_tissue_ct_factor;
"Others"->CRC_tissue_ct_factor$ColonSite[which(CRC_tissue_ct_factor$ColonSite=="Un")]
#--double group:
CRC_tissue_ct_factor$DoubleType->CRC_tissue_ct_factor$DoubleTypeGroup;
change_values(CRC_tissue_ct_factor,31,c("LL","HL","HL","HH"))->CRC_tissue_ct_factor;
############################################################################################
#_------------------plot sankey: 
table(CRC_tissue_ct_factor$SampleType,CRC_tissue_ct_factor$DoubleTypeGroup)->x
x/apply(x,1,sum)
#-for CRC vs colon site:
table(CRC_tissue_ct_tumor_factor$DoubleTypeGroup,CRC_tissue_ct_tumor_factor$ColonSite)->x;
x/apply(x,1,sum)
table(CRC_tissue_ct_tumor_factor$DoubleTypeGroup,CRC_tissue_ct_tumor_factor$MSI_Status)->x;
x/apply(x,1,sum)
table(CRC_tissue_ct_tumor_factor$DoubleTypeGroup,CRC_tissue_ct_tumor_factor$Age>60)->x;
x/apply(x,2,sum)



#------
add_Y_labels(CRC_tissue_ct_factor,"SampleType",c("Normal","Tumor"))->CRC_tissue_ct_factor;
#---for SDC2_Probes:
par(mfrow=c(2,2),mar=c(4,4,2,1))
do_logistic_fit(CRC_tissue_ct_factor,"Y_label","SDC2_ct")->SDC2_ct_logit_res;
#---for TFPI2_Probes:
do_logistic_fit(CRC_tissue_ct_factor,"Y_label","TFPI2_ct")->TFPI2_ct_logit_res;
#----for SDC2_Probes + TFPI2_Probes:
do_logistic_fit(CRC_tissue_ct_factor,"Y_label","SDC2_ct+TFPI2_ct")->SDC2_TFPI2_ct_logit_res;
#--add single gene cutoff :
get_x_cutoff(coef(SDC2_ct_logit_res$model),coords(SDC2_ct_logit_res$roc,"b", ret = "t", best.method = "youden"))->SDC2_ct_best_cutoff;
get_x_cutoff(coef(TFPI2_ct_logit_res$model),coords(TFPI2_ct_logit_res$roc,"b", ret = "t", best.method = "youden"))->TFPI2_ct_best_cutoff;
prepare_double_genes_types_cutoff(CRC_tissue_ct_factor,"SDC2_ct","TFPI2_ct",38,38)->CRC_tissue_ct_factor;
CRC_tissue_ct_factor[CRC_tissue_ct_factor$SampleType!="Normal",]->CRC_tissue_ct_tumor_factor
CRC_tissue_ct_tumor_factor$DoubleType->CRC_tissue_ct_tumor_factor$DoubleTypeGroup;
change_values(CRC_tissue_ct_tumor_factor,31,c("LL","HL","HL","HH"))->CRC_tissue_ct_tumor_factor;
table(CRC_tissue_ct_tumor_factor$DoubleTypeGroup,CRC_tissue_ct_tumor_factor$ColonSite)
#----------------------------plot points:
ggscatter(CRC_tissue_ct_tumor_factor, x = "SDC2_ct", y = "TFPI2_ct",color="DoubleType",conf.int = TRUE,palette="npg",shape="ColonSite")+stat_cor(method = "pearson")#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")
ggboxplot(CRC_tissue_ct_tumor_factor,x="DoubleType",y="Age")

##############################################################################################################################
##-------------------------------------------boxplot for SDC2 and TFPI2: + ROC
prepare_ggboxplot_data<-function(objRef_dat,genes,group,group_f=NULL){
	if(is.null(group_f)){
		names(table(objRef_dat$CliniInfo[,group]))->group_f;
	}
	#---
	objRef_dat$getSubSet(group,group_f)->objRef_dat;
	objRef_dat$getConditionGroup(group,group_f)->dat_group;
	objRef_dat$getDataMatrix(genes)->dat_matrix;
	#--------------
	rownames(dat_matrix)->gene_names;
	unlist(lapply(gene_names,function(gx){
		dat_matrix[gx,];
	}))->gene_values;
	data.frame("gName"=rep(gene_names,each=nrow(dat_group)),"SampleID"=rep(colnames(dat_matrix),length(gene_names)),"Expression"=gene_values,stringsAsFactors=F)->res;
	#data.frame("SampleID"=rep(dat_group$SampleID,length(gene_names)),"Condition"=rep(dat_group$Condition,length(gene_names)),stringsAsFactors=F)->res2;
	#---
	merge(res,dat_group,by.x="SampleID",by.y="SampleID")->res;
	factor(res$Condition,levels=group_f)->res$Condition;
	return(res);
}
prepare_ggboxplot_data(TCGA_CRC_methy_objRef,c("SDC2_P","TFPI2_P"),"SampleTypeCode",c("11","01"))->TCGA_CRC_methy_boxplot_data;
ggboxplot(TCGA_CRC_methy_boxplot_data, x = "gName", y = "Expression",color = "Condition", palette = c("gray","black"),add = "jitter",size=0.05)+stat_compare_means(aes(group=Condition),label = "p")->p1;
ggbarplot(TCGA_CRC_methy_boxplot_data, x = "gName", y = "Expression",fill = "Condition",color="black",add = "mean_se",palette = c("gray","black"),position= position_dodge())+stat_compare_means(aes(group=Condition),label="p")#jco
#-----for GEO;
prepare_ggboxplot_data(GSE48684_GSE79740_methy_objRef,c("SDC2_P","TFPI2_P"),"SampleType",c("Normal","Cancer"))->GSE48684_GSE79740_methy_boxplot_data;
ggboxplot(GSE48684_GSE79740_methy_boxplot_data, x = "gName", y = "Expression",color = "Condition", palette = c("gray","black"),add = "jitter",size=0.05)+stat_compare_means(aes(group=Condition),label = "p.signif")->p2;
ggbarplot(GSE48684_GSE79740_methy_boxplot_data, x = "gName", y = "Expression",fill = "Condition",color="black",add = "mean_se",palette = c("gray","black"),position= position_dodge())+stat_compare_means(aes(group=Condition),label="p.signif")#jco
#-----for CRC tissue;
prepare_ggboxplot_data(CRC_tissue_ct_objRef,c("SDC2_ct","TFPI2_ct"),"SampleType",c("Normal","Tumor"))->CRC_tissue_ct_boxplot_data;
ggboxplot(CRC_tissue_ct_boxplot_data, x = "gName", y = "Expression",color = "Condition", palette = c("gray","black"),add = "jitter",size=0.05)+stat_compare_means(aes(group=Condition),label = "p.signif")
ggbarplot(CRC_tissue_ct_boxplot_data, x = "gName", y = "Expression",fill = "Condition",color="black",add = "mean_se",palette = c("gray","black"),position= position_dodge())+stat_compare_means(aes(group=Condition),label="p.signif")#jco
###########################################################################################################################################
draw_multiClass_roc_v3<-function(dat,Y_label,variables,myd_colors,gtitle,ci_value=0.95){
	#---coords:
	c("sensitivity","specificity","ppv","npv","accuracy")->coords_ret;
	#-------------------
	c()->v_best
	c()->v_ci;
	#------------------
	variables[1]->v1;
	as.formula(paste(Y_label,v1,sep="~"))->v1_formula;
	glm(v1_formula, data = dat, family = "binomial")->model_glm
	predict(model_glm, newdata = dat, type = "response")->v1.test_prob
	pROC::roc(dat[,Y_label]~v1.test_prob,plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->v1.test_roc#
	plot(v1.test_roc,main=gtitle, col = myd_colors[1], print.thres = "best", print.auc = T,ci.type="shape");
	#----summary: youden index
	which.max(v1.test_roc$sensitivities+v1.test_roc$specificities)[1]->v_max;
	#----coords:
	coords(v1.test_roc,ret=coords_ret)->v1_coords;
	c(v_best,round(as.numeric(v1_coords[v_max,]),2),round(v1.test_roc$auc,2),round(v1.test_roc$ci[c(1,3)],2),paste(round(v1.test_roc$ci[c(1,3)],2),collapse="~"))->v_best;
	#---for sensi/speci ci 95%:
	v1_coords$sensitivity[v_max]/100->v_sensi_max;
	v1_coords$specificity[v_max]/100->v_speci_max;
	length(which(dat$Y_label==1))->v_sensi_size;
	length(which(dat$Y_label==0))->v_speci_size;
	c(v_sensi_max-qnorm((1+ci_value)/2)*sqrt(v_sensi_max*(1-v_sensi_max)/v_sensi_size),v_sensi_max+qnorm((1+ci_value)/2)*sqrt(v_sensi_max*(1-v_sensi_max)/v_sensi_size))->v1_sensi_ci;
	c(v_speci_max-qnorm((1+ci_value)/2)*sqrt(v_speci_max*(1-v_speci_max)/v_speci_size),v_speci_max+qnorm((1+ci_value)/2)*sqrt(v_speci_max*(1-v_speci_max)/v_speci_size))->v1_speci_ci;
	#---ppv, npv
	v1_coords$ppv[v_max]/100->v_ppv_max;
	v1_coords$npv[v_max]/100->v_npv_max;
	length(which(dat[,v1]==1))->v_ppv_size;
	length(which(dat[,v1]==0))->v_npv_size;
	c(v_ppv_max-qnorm((1+ci_value)/2)*sqrt(v_ppv_max*(1-v_ppv_max)/v_ppv_size),v_ppv_max+qnorm((1+ci_value)/2)*sqrt(v_ppv_max*(1-v_ppv_max)/v_ppv_size))->v1_ppv_ci;
	c(v_npv_max-qnorm((1+ci_value)/2)*sqrt(v_npv_max*(1-v_npv_max)/v_npv_size),v_npv_max+qnorm((1+ci_value)/2)*sqrt(v_npv_max*(1-v_npv_max)/v_npv_size))->v1_npv_ci;
	#--accuracy
	v1_coords$accuracy[v_max]/100->v_accu;
	v_ppv_size+v_npv_size->v_size;
	c(v_accu-qnorm((1+ci_value)/2)*sqrt(v_accu*(1-v_accu)/v_size),v_accu+qnorm((1+ci_value)/2)*sqrt(v_accu*(1-v_accu)/v_size))->v1_accu_ci;
	c(v_ci,v1_sensi_ci,v1_speci_ci,v1_ppv_ci,v1_npv_ci,v1_accu_ci)->v_ci;
	#print(sqrt(v_sensi_max*(1-v_sensi_max)/v_sensi_size));flush.console();
	#--------
	2->i;
	for(vx in variables[-1]){
		as.formula(paste(Y_label,vx,sep="~"))->vx_formula;
		glm(vx_formula, data = dat, family = "binomial")->model_glm
		predict(model_glm, newdata = dat, type = "response")->vx.test_prob
		pROC::roc(dat[,Y_label]~vx.test_prob, plot = F,ci=T,print.auc = F,ci.method="boot",boot.n=100,percent=T,algorithm = 6)->vx.test_roc
		plot(vx.test_roc,main=gtitle, col = myd_colors[i],print.thres = "best", print.auc = T,ci.type="shape",add=T);
		i+1->i;
		#-----summary 
		which.max(vx.test_roc$sensitivities+vx.test_roc$specificities)[1]->vx_max;
		#----coords:
		coords(vx.test_roc,ret=coords_ret)->vx_coords;
		c(v_best,round(as.numeric(vx_coords[vx_max,]),2),round(vx.test_roc$auc,2),round(vx.test_roc$ci[c(1,3)],2),paste(round(vx.test_roc$ci[c(1,3)],2),collapse="~"))->v_best;
		#---sens, speci: ci95%
		vx_coords$sensitivity[vx_max]/100->vx_sensi_max;
		vx_coords$specificity[vx_max]/100->vx_speci_max;
		c(vx_sensi_max-qnorm((1+ci_value)/2)*sqrt(vx_sensi_max*(1-vx_sensi_max)/v_sensi_size),vx_sensi_max+qnorm((1+ci_value)/2)*sqrt(vx_sensi_max*(1-vx_sensi_max)/v_sensi_size))->vx_sensi_ci;
		c(vx_speci_max-qnorm((1+ci_value)/2)*sqrt(vx_speci_max*(1-vx_speci_max)/v_speci_size),vx_speci_max+qnorm((1+ci_value)/2)*sqrt(vx_speci_max*(1-vx_speci_max)/v_speci_size))->vx_speci_ci;
		#---ppv, npv
		vx_coords$ppv[v_max]/100->vx_ppv_max;
		vx_coords$npv[v_max]/100->vx_npv_max;
		length(which(dat[,vx]==1))->vx_ppv_size;
		length(which(dat[,vx]==0))->vx_npv_size;
		c(vx_ppv_max-qnorm((1+ci_value)/2)*sqrt(vx_ppv_max*(1-vx_ppv_max)/vx_ppv_size),vx_ppv_max+qnorm((1+ci_value)/2)*sqrt(vx_ppv_max*(1-vx_ppv_max)/vx_ppv_size))->vx_ppv_ci;
		c(vx_npv_max-qnorm((1+ci_value)/2)*sqrt(vx_npv_max*(1-vx_npv_max)/vx_npv_size),vx_npv_max+qnorm((1+ci_value)/2)*sqrt(vx_npv_max*(1-vx_npv_max)/vx_npv_size))->vx_npv_ci;
		#--accuracy
		vx_coords$accuracy[v_max]/100->vx_accu;
		vx_ppv_size+vx_npv_size->vx_size;
		c(vx_accu-qnorm((1+ci_value)/2)*sqrt(vx_accu*(1-vx_accu)/vx_size),vx_accu+qnorm((1+ci_value)/2)*sqrt(vx_accu*(1-vx_accu)/vx_size))->vx_accu_ci;
		c(v_ci,vx_sensi_ci,vx_speci_ci,vx_ppv_ci,vx_npv_ci,vx_accu_ci)->v_ci;
	}
	legend("bottomright",legend=variables,lty=1,lwd=1,col=myd_colors[1:length(variables)]);
	#---
	matrix(v_best,nrow=length(variables),byrow=T)->v_best;
	c("best_sensi","best_spec","ppv","npv","accuracy","auc","auc95cil","auc95ciu","auc_ci")->colnames(v_best);
	matrix(v_ci,nrow=length(variables),byrow=T)->v_ci;
	c("se.low","se.high","sp.low","sp.high","ppv.low","ppv.high","npv.low","npv.high","accuracy.low","accuracy.high")->colnames(v_ci);
	data.frame("Feature"=variables,v_best,v_ci,stringsAsFactors=F)->v_best;
	return(v_best);
}
#---------------------------------------------------------------------------------------------------------------------------------------------
#----------------the optimal beta values 
#sampling 100 times
c()->TCGA_g1_best_cutoff;
c()->TCGA_g2_best_cutoff;
for(i in 1:100){
	set.seed(i^2+i*3+i);
	which(TCGA_CRC_methy_factor$Y_label==0)->n_index;
	sample(which(TCGA_CRC_methy_factor$Y_label==1),size=45,replace=F)->t_index;
	TCGA_CRC_methy_factor[c(n_index,t_index),]->tmp_d;
	#---SDC2_P
	do_logistic_fit(tmp_d,"Y_label","SDC2_P")->TCGA_CRC_methy_logit_res1;
	get_x_cutoff(coef(TCGA_CRC_methy_logit_res1$model),coords(TCGA_CRC_methy_logit_res1$roc,"b", ret = "t", best.method = "youden"))->tmp_d_v;
	c(TCGA_g1_best_cutoff,tmp_d_v)->TCGA_g1_best_cutoff;
	#---TFPI2_P
	do_logistic_fit(tmp_d,"Y_label","TFPI2_P")->TCGA_CRC_methy_logit_res2;
	get_x_cutoff(coef(TCGA_CRC_methy_logit_res2$model),coords(TCGA_CRC_methy_logit_res2$roc,"b", ret = "t", best.method = "youden"))->tmp_d_v;
	c(TCGA_g2_best_cutoff,tmp_d_v)->TCGA_g2_best_cutoff;
}
writeLines(paste(TCGA_g1_best_cutoff+0.02,"",sep=""),"test.txt")
writeLines(paste(TCGA_g2_best_cutoff-0.07,"",sep=""),"test.txt")
draw_multiClass_roc_v3(tmp_d,"Y_label",c("SDC2_P","TFPI2_P"),get_palette("aaas",5),"Normal vs CRC")->test_;
#--
writeLines(paste(GEO_g1_best_cutoff-0.02,"",sep=""),"test.txt")
writeLines(paste(GEO_g2_best_cutoff-0.02,"",sep=""),"test.txt")

#######################################################################################################################################
#----------------------------------------------------------------------
#-------------ROC 
par(mfrow=c(1,3))
change_methy_values(TCGA_CRC_methy_factor,"SDC2_P",0.2)->test_;
do_logistic_fit(test_,"Y_label","SDC2_P")->SDC2_TFPI2_logit_res;
change_methy_values(test_,"TFPI2_P",0.2)->test_;
do_logistic_fit(test_,"Y_label","TFPI2_P")->SDC2_TFPI2_logit_res;
do_logistic_fit(test_,"Y_label","SDC2_P+TFPI2_P")->SDC2_TFPI2_logit_res;
unlist(lapply(1:nrow(test_),function(i){if(test_$SDC2_P[i]==0 && test_$TFPI2_P[i]==0){0}else{1}}))->x
table(test_$Y_label,x)

#---GEO
change_methy_values(GSE48684_GSE79740_methy_factor,"SDC2_P",0.2)->test_;
do_logistic_fit(test_,"Y_label","SDC2_P")->GEO_SDC2_TFPI2_logit_res;
change_methy_values(test_,"TFPI2_P",0.2)->test_;
do_logistic_fit(test_,"Y_label","TFPI2_P")->GEO_SDC2_TFPI2_logit_res;
do_logistic_fit(test_,"Y_label","SDC2_P+TFPI2_P")->SDC2_TFPI2_logit_res;
unlist(lapply(1:nrow(test_),function(i){if(test_$SDC2_P[i]==0 && test_$TFPI2_P[i]==0){0}else{1}}))->x
table(test_$Y_label,x)

#----for CRC tissue ct values:
change_methy_values(CRC_tissue_ct_factor,"SDC2_ct",38,reverse=T)->test_;
do_logistic_fit(test_,"Y_label","SDC2_ct")->SDC2_TFPI2_ct_logit_res;
change_methy_values(test_,"TFPI2_ct",38,reverse=T)->test_;
do_logistic_fit(test_,"Y_label","TFPI2_ct")->SDC2_TFPI2_ct_logit_res;
do_logistic_fit(test_,"Y_label","SDC2_ct+TFPI2_ct")->SDC2_TFPI2_logit_res;
unlist(lapply(1:nrow(test_),function(i){if(test_$SDC2_ct[i]==0 && test_$TFPI2_ct[i]==0){0}else{1}}))->x
table(test_$Y_label,x)


##########################################################################################################
#--------------compare methylation to MSI score 
read.table("TCGA-MSI-mantis_score.txt",header=T,sep="\t",stringsAsFactors=F)->TCGA_pancancer_MANTIS_Score;
substr(TCGA_pancancer_MANTIS_Score$Tumor.Filename,1,15)->TCGA_pancancer_MANTIS_Score$A0_Samples;
TCGA_pancancer_MANTIS_Score[,c("A0_Samples","MANTIS.Score")]->TCGA_pancancer_MANTIS_Score;
unlist(lapply(TCGA_pancancer_MANTIS_Score$MANTIS.Score,function(x){
	if(is.na(x)){
		"Un";
	}else if(x>=0.4){"MSI-H"}else{"MSS"};
}))->TCGA_pancancer_MANTIS_Score$MSI_Status;
merge(TCGA_CRC_methy_tumor_factor,TCGA_pancancer_MANTIS_Score,by.x="A0_Samples",by.y="A0_Samples",all.x=T)->TCGA_CRC_methy_tumor_factor
#------------------ for TMB analysis : not good 
read.table("immune-landscape.txt",header=T,sep="\t",stringsAsFactors=F)->TCGA_pancancer_immuneladscape;
"A0_Samples"->names(TCGA_pancancer_immuneladscape)[1];
paste(TCGA_pancancer_immuneladscape$A0_Samples,"01",sep="-")->TCGA_pancancer_immuneladscape$A0_Samples
merge(TCGA_CRC_methy_tumor_factor,TCGA_pancancer_immuneladscape,by.x="A0_Samples",by.y="A0_Samples",all.x=T)->TCGA_CRC_immune_factor
TCGA_CRC_immune_factor$DoubleType->TCGA_CRC_immune_factor$DoubleTypeGroup;
change_values(TCGA_CRC_immune_factor,107,c("HH","HL","HL","LL"))->TCGA_CRC_immune_factor
#-----for SDC2/TFPI2 methylation :
c("SDC2_P","TFPI2_P","Nonsilent.Mutation.Rate","MANTIS.Score")->selected_geneSet;
#subset(TCGA_CRC_immune_factor,MSI_Status=="MSI-H")->test_;
ggscatter(TCGA_CRC_immune_factor,x="SDC2_P",y="MANTIS.Score",color="MSI_Status",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(5,9)],add="reg.line")+stat_cor(method = "pearson")#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")
ggscatter(TCGA_CRC_immune_factor,x="TFPI2_P",y="MANTIS.Score",color="MSI_Status",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(5,9)],add="reg.line")+stat_cor(method = "pearson")#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")
#ggboxplot(test_,x="DoubleTypeGroup",y="MANTIS.Score",color="DoubleTypeGroup",conf.int = TRUE,palette=rev(brewer.pal(9,"Greys")),yscale="log2",add="jitter")+stat_compare_means()
#------------------------------------------summary TMB:
ggboxplot(TCGA_CRC_immune_factor,x="DoubleTypeGroup",y="Nonsilent.Mutation.Rate",color="MSI_Status",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,5)],yscale="log2",add="jitter")+stat_compare_means()
ggboxplot(TCGA_CRC_immune_factor,x="DoubleTypeGroup",y="Silent.Mutation.Rate",color="MSI_Status",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,5)],yscale="log2",add="jitter")+stat_compare_means()

ggboxplot(TCGA_CRC_immune_factor,x="DoubleTypeGroup",y="A1_OS",color="DoubleTypeGroup",conf.int = TRUE,palette="npg",add="jitter")+stat_compare_means()

#----------------------for 3D scatter
prepare_colors<-function(dat,f,myd_col){
	names(table(dat[,f]))->f_names;
	myd_col[1:length(f_names)]->myd_col;
	f_names->names(myd_col);
	myd_col[as.character(dat[,f])]->f_cols;
	return(f_cols);
}
draw_scatterplot3d<-function(dat,x,y,z,f,myd_col){
	#----
	names(table(dat[,f]))->f_names;
	myd_col[1:length(f_names)]->myd_col;
	f_names->names(myd_col);
	myd_col[as.character(dat[,f])]->f_cols;
	#----
	scatterplot3d(dat[,x],dat[,y],dat[,z],pch=20,color=f_cols,xlab=x,ylab=y,zlab=z);
	legend("topleft",pch=rep(20,length(f_names)),legend=f_names,col=myd_col,cex=1.5);

}
draw_scatterplot3d(TCGA_CRC_immune_factor,"SDC2_P","TFPI2_P","MANTIS.Score","MSI_Status",myd_colors_npg)
#--------------for OS/PFI survival:
add_status<-function(expd,alive_symbol){
	expd[which(!is.na(expd$A2_Event)),]->expd;
	expd[which(!is.na(expd$A1_OS)),]->expd;
	as.numeric(expd$A1_OS)->expd$A1_OS
	lapply(expd$A2_Event,function(xi){
		if(xi==alive_symbol){
			0
		}else{
			1
		}
	})->x_status;
	unlist(x_status)->expd$Status;
	return(expd);
}
draw_survial_curve_custom_v2<-function(myd,column,bk,myd_colors,show_table=NULL){
	myd[myd[,column]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,column]),]->myd.rm;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~myd.rm[,column])->myd.fit;
	survdiff(formula=myd.surv~myd.rm[,column],rho=0)->myd.diff;
	table(myd.rm[,column])->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	#--------------------------
	plot(myd.fit,col=myd_colors,xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,lty=3,axes=F,main=paste("KM(",colnames(myd.rm)[column],")",sep=""),xlim=c(-100,max_xlim*1.1));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk),pos=0);
	rug(x=seq(0,max_xlim,bk)+bk/2,ticksize=-0.01,side=1,pos=0);
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20),pos=0);
	rug(x=seq(0,0.9,0.2)+0.1,ticksize=-0.01,side=2,pos=0);
	#----add point shapes:
	names(myd.fit$strata)->strata_names;
	myd.fit$time->strata_times;
	myd.fit$surv->strata_surv;
	1->strata_start;
	for(i in 1:length(strata_names)){
		myd.fit$strata[i]+strata_start-1->strata_end;
		strata_times[strata_start:strata_end]->i_times;
		strata_surv[strata_start:strata_end]->i_surv;
		points(x=i_times,y=i_surv,pch=3,col=myd_colors[i],cex=0.8);
		strata_end+1->strata_start
	}
	#--------
	#abline(h=seq(0.2,1,0.2),col=brewer.pal(9,"Greys")[3],lty=3)
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(N=",myd.table,")",sep="")),pch=3,pt.cex=1.2,col=myd_colors,text.col=myd_colors,bty="n",cex=1.2);
	if(pvalue<1e-5){
		legend(x=bk/2,y=0.15,legend="p < 1e-5",bty="n",cex=1.2)
	}else{
		legend(x=bk/2,y=0.15,legend=paste("p=",round(pvalue,5),sep=""),bty="n",cex=1.2)
	}
	return(c(pvalue,myd.table));
}
par(mfrow=c(1,2))
add_status(TCGA_CRC_immune_factor,"Alive")->test_;
draw_survial_curve_custom_v2(test_,107,365,myd_colors)#p=0.3829106
#---for PFS:
TCGA_CRC_immune_factor[,c("A0_Samples","NewEvent","NewEventTime","DoubleTypeGroup")]->test_;
test_$NewEvent->test_$A2_Event;
test_$NewEventTime->test_$A1_OS;
add_status(test_,"0")->test_;
draw_survial_curve_custom_v2(test_,4,365,myd_colors)#p=0.822071
#------------------------------------------------------------------------------------------------------------
prepare_graphprim_column<-function(expd_fator,genes,f,f_groups=NULL){
	table(expd_fator[,f])->f_table;
	if(is.null(f_groups)){
		names(f_table)->f_groups;
	}
	f_table[f_groups]->f_table;
	max(f_table)->f_group_max;
	#------------
	lapply(genes,function(gx){
		c()->gx_values;
		for(fx in f_groups){
			which(expd_fator[,f]==fx)->fx_index;
			c(gx_values,expd_fator[fx_index,gx],rep(NA,f_group_max-f_table[fx]))->gx_values;
		}
		matrix(gx_values,ncol=length(f_groups),byrow=F)->gx_matrix;
		f_groups->colnames(gx_matrix);
		gx_matrix;
	})->genes_res;
	genes->names(genes_res)
	return(genes_res);
}
prepare_graphprim_column(TCGA_CRC_immune_factor,c("Age"),f="DoubleTypeGroup")->test_;
write.csv(test_[[1]],"results/TCGA_CRC_age_doublegroups.csv",quote=F,row.names=F)
prepare_graphprim_column(CRC_tissue_ct_tumor_factor,c("Age"),f="DoubleTypeGroup")->test_;
write.csv(test_[[1]],"results/CRC_tissue_ct_age_doublegroups.csv",quote=F,row.names=F)
#----------------------------for age :
ggboxplot(TCGA_CRC_immune_factor,x="DoubleTypeGroup",y="Age",color="DoubleTypeGroup",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,7,4)],add="jitter")+stat_compare_means()#p<0.05
ggboxplot(CRC_tissue_ct_tumor_factor,x="DoubleType",y="Age",color="DoubleType",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,7,4)],add="jitter")+stat_compare_means()#
#------------------------------------------------------------------------------------------------------------
#---------for colon site with MSI-H
table(TCGA_CRC_immune_factor$ColonSite,TCGA_CRC_immune_factor$MSI_Status)->test_;
table(CRC_tissue_ct_tumor_factor$ColonSite,CRC_tissue_ct_tumor_factor$MSI_Status)->test_;
#draw_barplot(test_,myd_colors_npg)
balloonplot(test_,main="TCGA CRC tissue source vs MSI",xlab ="", ylab="",label = FALSE, show.margins = FALSE)
fisher.test(test_)$p.value->table_pvalue;
if(table_pvalue<1e-5){
	legend("topleft",legend=paste("Fisher-p:","<1e-5"))
}else{
	legend("topleft",legend=paste("Fisher-p:",round(table_pvalue,5)))
}
#------------------------------relationship of delta Ct with age 
as.numeric(CRC_tissue_ct_tumor_factor$Age)->CRC_tissue_ct_tumor_factor$Age;
2^(-CRC_tissue_ct_tumor_factor$SDC2_ct+CRC_tissue_ct_tumor_factor$ACTB_ct)->CRC_tissue_ct_tumor_factor$SDC2_delta_ct;
2^(-CRC_tissue_ct_tumor_factor$TFPI2_ct+CRC_tissue_ct_tumor_factor$ACTB_ct)->CRC_tissue_ct_tumor_factor$TFPI2_delta_ct;
ggscatter(CRC_tissue_ct_tumor_factor,x="Age",y="SDC2_delta_ct",yscale="log2",add="reg.line",cor.coef=TRUE,cor.method="spearman",conf.int=TRUE)->p1;
ggscatter(CRC_tissue_ct_tumor_factor,x="Age",y="TFPI2_delta_ct",yscale="log2",add="reg.line",cor.coef=TRUE,cor.method="spearman",conf.int=TRUE)->p2;#p1+rotate()
grid.arrange(p1,p2,ncol=2)

#################################################################################################################################
#---------------------------------------------------compare with age:
ggscatter(TCGA_CRC_immune_factor,x="SDC2_P",y="Age",conf.int = TRUE,add="reg.line")+stat_cor(method = "pearson")#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")
ggscatter(TCGA_CRC_immune_factor,x="TFPI2_P",y="Age",conf.int = TRUE,add="reg.line")+stat_cor(method = "pearson")#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")

#ggscatter(CRC_tissue_ct_tumor_factor[CRC_tissue_ct_tumor_factor$SDC2_ct!=45,],x="SDC2_ct",y="Age",conf.int = TRUE,add="reg.line")+stat_cor(method = "pearson")#add = "reg.line",add.params = list(color = "blue",fill = "lightgray")
draw_scatterplot3d(TCGA_CRC_immune_factor,"SDC2_P","TFPI2_P","Age","SampleTypeCode",c("white","black"))
draw_scatterplot3d(CRC_tissue_ct_tumor_factor,"SDC2_ct","TFPI2_ct","Age","SampleType",c("black"))

####################################################################################################
#---------------------------MMR_genes: methylation;"MLH1"  "MSH2"  "MSH6"  "PMS2"  "EXO1"  "POLE"  "POLD1"
HumanMethylation450_annotats[grep("MLH1",HumanMethylation450_annotats$UCSC_RefGene_Name),]->MLH1_annotats;
MLH1_annotats[grep("TSS",MLH1_annotats$UCSC_RefGene_Group),]->MLH1_annotats;
as.character(MLH1_annotats$IlmnID)->MLH1_probes;

#----
change_GEObjRef_to_expd(TCGA_CRC_methy_objRef,select_f="probe")->TCGA_CRC_methy_MLH1_factor;
TCGA_CRC_methy_MLH1_factor[,c(colnames(TCGA_CRC_methy_objRef$CliniInfo),MLH1_probes)]->TCGA_CRC_methy_MLH1_factor;
merge(TCGA_CRC_methy_MLH1_factor,TCGA_CRC_immune_factor[,c("A0_Samples","DoubleType","DoubleTypeGroup")],by.x="A0_Samples",by.y="A0_Samples",all.x=T)->TCGA_CRC_methy_MLH1_factor;
"Normal"->TCGA_CRC_methy_MLH1_factor$DoubleType[which(is.na(TCGA_CRC_methy_MLH1_factor$DoubleType))];
"Normal"->TCGA_CRC_methy_MLH1_factor$DoubleTypeGroup[which(is.na(TCGA_CRC_methy_MLH1_factor$DoubleTypeGroup))];
#---merge by sample types:
TCGA_CRC_methy_MLH1_factor[,c("A0_Samples","DoubleTypeGroup")]->MLH1_TN_group;
c("SampleID","Condition")->colnames(MLH1_TN_group)
calculate_gene_groupScores(MLH1_methy_matrix,MLH1_TN_group,rownames(MLH1_methy_matrix))->MLH1_methy_groupScores;
write.csv(MLH1_methy_groupScores,"results/MLH1_methy_groupScores.csv",quote=F,row.names=F)
##_--------------------------merge by doubletype:
HumanMethylation450_annotats[grep("MLH1",HumanMethylation450_annotats$UCSC_RefGene_Name),]->MLH1_annotats;
MLH1_annotats[grep("TSS",MLH1_annotats$UCSC_RefGene_Group),]->MLH1_annotats;
as.character(MLH1_annotats$IlmnID)->MLH1_probes;

TCGA_CRC_methy_objRef$getDataMatrix(MLH1_probes,select_f="probe")->MLH1_methy_matrix;
MLH1_methy_matrix[grep("MLH1",rownames(MLH1_methy_matrix)),]->MLH1_methy_matrix;
calculate_gene_probeScores(MLH1_methy_matrix,MLH1_TN_group,rownames(MLH1_methy_matrix))->test_;
ggboxplot(test_,x="Condition",y="Expression")

#---heatmap:
preorder_samples_by_factor(TCGA_CRC_methy_MLH1_factor,MLH1_probes,"DoubleTypeGroup",c("Normal","LL","HL","HH"))->TCGA_CRC_methy_MLH1_sort;
#,"Normal","Colon polypus","Coloproctitis","Diverticulum","Hemorrhoids","Esophagitis","Gastric polyps","Gastritis","Other","Non-advanced adenoma","Advanced adenoma","CRC","Anal cancer","Appendix cancer","Liver cancer","Stomach cancer"
TCGA_CRC_methy_MLH1_sort[,c(13,72)]->annot_df;
TCGA_CRC_methy_MLH1_sort$A0_Samples->rownames(annot_df);
TCGA_CRC_methy_MLH1_factor[,MLH1_probes]->heatmap_dat;
TCGA_CRC_methy_MLH1_sort$A0_Samples->rownames(heatmap_dat);
list("DoubleTypeGroup"=generate_color(annot_df$DoubleTypeGroup,myd_colors),"ColonSite"=generate_color(annot_df$ColonSite,myd_colors_npg))->col_colors;
pheatmap(t(heatmap_dat),cluster_cols=F,show_colnames=F,annotation_col=annot_df,annotation_colors=col_colors,cellheight=20)

#########################################################################################################################
#-----------------------------------------SNV: 
read.maf(maf="TCGA.CRC.mutect.somatic-maftools.maf",clinicalData="TCGA_CRC_immune_factor.txt")->TCGA_CRC_snv_maf;
#--------
preorder_samples_GeneC_OS(TCGA_CRC_immune_factor,c("HH","HL","LL"))->myd_exp_cluster.sort;
as.character(myd_exp_cluster.sort$Tumor_Sample_Barcode)->SNV.sample_order;
list("DoubleTypeGroup"=generate_color(myd_exp_cluster.sort$DoubleTypeGroup,myd_colors_npg),"ColonSite"=generate_color(myd_exp_cluster.sort$ColonSite,myd_colors_npg))->SNV.annot_IRGCluster;
plot.new()
oncoplot(maf=TCGA_CRC_snv_maf,genes=MMR_genes,bgCol="white",borderCol=NA,removeNonMutated=F,clinicalFeatures=c("DoubleTypeGroup","ColonSite"), sortByAnnotation =T,draw_titv = T,annotationColor=SNV.annot_IRGCluster)#
fisher.test(matrix(c(87,2,1,353,34,7),ncol=2,byrow=F))#p=0.075
#-----------------------------------------------for DoubleTypeGroup:
subsetMaf(maf = TCGA_CRC_snv_maf,clinQuery = "DoubleTypeGroup %in% 'HH'")->TCGA_CRC_snv_maf_HH;
plotmafSummary(TCGA_CRC_snv_maf_HH,top = 20);
subsetMaf(maf = TCGA_CRC_snv_maf, clinQuery = "DoubleTypeGroup %in% 'HL'")->TCGA_CRC_snv_maf_HL;
plotmafSummary(TCGA_CRC_snv_maf_HL,top = 20);
subsetMaf(maf = TCGA_CRC_snv_maf, clinQuery = "DoubleTypeGroup %in% 'LL'")->TCGA_CRC_snv_maf_LL;
plotmafSummary(TCGA_CRC_snv_maf_LL,top = 20);

#-----------------for CIMP-H : 
c("BRAF","PIK3CA","KRAS","TP53","APC")->CIMP_genes;
subsetMaf(maf = TCGA_CRC_snv_maf,genes = CIMP_genes,clinQuery = "DoubleTypeGroup %in% 'HH'")->TCGA_CRC_snv_maf_HH;
subsetMaf(maf = TCGA_CRC_snv_maf,genes = CIMP_genes,clinQuery = "DoubleTypeGroup %in% 'HL'")->TCGA_CRC_snv_maf_HL;
subsetMaf(maf = TCGA_CRC_snv_maf,genes = CIMP_genes,clinQuery = "DoubleTypeGroup %in% 'LL'")->TCGA_CRC_snv_maf_LL;
#-
plot.new()
oncoplot(maf=TCGA_CRC_snv_maf,genes=CIMP_genes,bgCol="white",borderCol=NA,removeNonMutated=F,clinicalFeatures=c("DoubleTypeGroup","ColonSite"), sortByAnnotation =T,draw_titv = T,annotationColor=SNV.annot_IRGCluster)#
#---
clinicalEnrichment(maf = TCGA_CRC_snv_maf, clinicalFeature = 'DoubleTypeGroup')->TCGA_CRC_snv_maf.doubletype;
TCGA_CRC_snv_maf.doubletype$groupwise_comparision[p_value < 0.05]
plotEnrichmentResults(enrich_res = TCGA_CRC_snv_maf.doubletype, pVal = 0.05)
#------------------------MMR+CIMP
plot.new()
oncoplot(maf=TCGA_CRC_snv_maf,genes=c(MMR_genes,CIMP_genes),bgCol="white",borderCol=NA,removeNonMutated=F,clinicalFeatures=c("DoubleTypeGroup","ColonSite"), sortByAnnotation =T,draw_titv = T,annotationColor=SNV.annot_IRGCluster)#



################################################################################################################################################################
##############################################################################################################################################
#--------------------------------------------------------------------------differentially expressed genes in HH,HL and LL
do_anova_mutliClass<-function(expd,groups){
	as.character(groups$SampleID)->group_samples;
	rownames(expd)->gene_names;
	lapply(1:nrow(expd),function(i){
		data.frame("y"=expd[i,group_samples],groups)->i_df;
		aov(y~Condition,data=i_df)->i_df_aov;
		summary.lm(i_df_aov)->i_df_aov_summary;
		summary(i_df_aov)[[1]][["Pr(>F)"]][1]->i_pvalue;
		#----
		i_df_aov_summary$coefficients->i_coefs;
		gsub("[(|)]","",rownames(i_coefs))->rownames(i_coefs);
		i_coefs[,"Estimate"]->i_coefs_esti;
		i_coefs[,"Pr(>|t|)"]->i_coefs_pvalue;
		i_df_aov_summary$fstatistic[1]->i_fstatistic;
		i_df_aov_summary$sigma->i_rse;
		i_df_aov_summary$r.squared->i_rsquared;
		#--------
		c(i_coefs_esti,i_coefs_pvalue,i_fstatistic,i_rse,i_rsquared,i_pvalue)->i_res;
		c(paste(rownames(i_coefs),"estimate",sep="_"),paste(rownames(i_coefs),"pvalue",sep="_"),"fstatistic","RSE","RSquared","P.Value")->names(i_res);
		i_res;
	})->res;
	#--
	names(res[[1]])->tmp_colnames;
	unlist(res)->res;
	matrix(res,nrow=length(gene_names),byrow=T)->res_mat;
	tmp_colnames->colnames(res_mat);
	#--------
	data.frame("gName"=gene_names,res_mat,stringsAsFactors=F)->res_mat;
	p.adjust(res_mat$P.Value)->res_mat$FDR;
	#----
	return(res_mat);
}
#-----------------------------------------compare 3 groups : HH, HL and LL specific genes
#--do anova 
TCGA_CRC_object_processed$getConditionGroup("DoubleTypeGroup",c("HH","HL","LL"))->TCGA_CRC_expd_group#,
TCGA_CRC_object_processed$getGroupMatrix("DoubleTypeGroup",c("HH","HL","LL"))->TCGA_CRC_expd_matrix
do_anova_mutliClass(TCGA_CRC_expd_matrix,TCGA_CRC_expd_group)->TCGA_CRC_expd_diff_anova;
TCGA_CRC_expd_diff_anova[which(TCGA_CRC_expd_diff_anova$FDR<0.05),]->TCGA_CRC_expd_diff_anova_filter;
#------------------------------------- select group specific genes:
select_group_specific_genes(TCGA_CRC_expd_matrix,TCGA_CRC_expd_group,method="ranktest")->TCGA_CRC_expd_anova_res;
#-----show the results:
table(TCGA_CRC_expd_anova_res$Group,TCGA_CRC_expd_anova_res$gName)->TCGA_CRC_expd_anova_res_table;
barplot(table(TCGA_CRC_expd_anova_res$Group,TCGA_CRC_expd_anova_res$gName),beside=T,col=myd_colors_npg[1:3],border=NA,las=2)
legend("topright",fill=myd_colors_npg[1:3],legend=rownames(table(TCGA_CRC_expd_anova_res$Group,TCGA_CRC_expd_anova_res$gName)),border=NA)
#-
#list("gName"=generate_color(TCGA_CRC_methy_heatmap$ProbeAnnot$gName,brewer.pal(8,"Set1")))->col_colors;
pheatmap(TCGA_CRC_expd_anova_res_table,show_rownames=T,show_colnames=F,border_color=NA,legend=F,cellwidth=6)->TCGA_CRC_expd_anova_res_table.heatmap
calculate_gene_groupScores(TCGA_CRC_expd_matrix,TCGA_CRC_expd_group,colnames(TCGA_CRC_expd_anova_res_table))->TCGA_CRC_expd_anova_groupScores;
pheatmap(t(TCGA_CRC_expd_anova_groupScores[TCGA_CRC_expd_anova_res_table.heatmap$tree_col$order,-1]),cluster_cols=F,scale="column",show_rownames=T,show_colnames=F,border_color=NA,legend=F,cellwidth=6)->TCGA_CRC_expd_anova_groupScores.pheatmap;
TCGA_CRC_expd_matrix[colnames(TCGA_CRC_expd_anova_res_table),]->TCGA_CRC_expd_specific_values;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(230),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(200))->fill_colors;
pheatmap(t(log2(TCGA_CRC_expd_specific_values[TCGA_CRC_expd_anova_res_table.heatmap$tree_col$order,]+1)),show_rownames=F,cluster_rows=T,cluster_cols=F,scale="column",color=fill_colors,cellwidth=6)->TCGA_CRC_expd_anova_groupScores.pheatmap;
#----
grid.arrange(TCGA_CRC_expd_anova_res_table.heatmap$gtable,TCGA_CRC_expd_anova_groupScores.pheatmap$gtable,layout_matrix=matrix(c(1,2,2,2,2),nrow=5,ncol=1,byrow=T))
#---------------------------use dotplot for KEGG and GO:
draw_GO_KEGG_dotplot<-function(enrich_df,show_c=10){
	ifelse(nrow(enrich_df)<show_c,nrow(enrich_df),show_c)->show_c;
	enrich_df[1:show_c,]->enrich_df;
	enrich_df[order(enrich_df$term_genes_found),]->enrich_df;
	factor(enrich_df$term,levels=enrich_df$term)->enrich_df$term;
	#-------
	enrich_df$term_genes_found/enrich_df$input_size->enrich_df$GeneRatio;
	-log10(enrich_df$hyp_pval_adj)->enrich_df$Log10_pvalue;
	ggplot(enrich_df,aes(x=term,y=GeneRatio))->p1
	p1+geom_point(aes(color=Log10_pvalue,size=term_genes_found))+coord_flip()+scale_color_gradient(low ="blue",high ="red")+scale_size_continuous(range = c(1,10))->p1;
	return(p1);
}
read.table("GO-KEGG/HH-GO_BP_Results-GeneCodis4.tsv",header=T,sep="\t",stringsAsFactors=F)->HH_specific_genes_KEGG;
c("Items_Details","Items","term_genes_found","input_size","Reference.Support","Reference.size","Hyp","hyp_pval_adj","Genes")->colnames(HH_specific_genes_KEGG)
paste(HH_specific_genes_KEGG$Items,HH_specific_genes_KEGG$Items_Details,sep=":")->HH_specific_genes_KEGG$term;
draw_GO_KEGG_dotplot(HH_specific_genes_KEGG,10)->HH_specific_genes_KEGG_p1;
#-
read.table("GO-KEGG/LL-GO_BP_Results-GeneCodis4.tsv",header=T,sep="\t",stringsAsFactors=F)->LL_specific_genes_KEGG;
c("Items_Details","Items","term_genes_found","input_size","Reference.Support","Reference.size","Hyp","hyp_pval_adj","Genes")->colnames(LL_specific_genes_KEGG)
paste(LL_specific_genes_KEGG$Items,LL_specific_genes_KEGG$Items_Details,sep=":")->LL_specific_genes_KEGG$term;
draw_GO_KEGG_dotplot(LL_specific_genes_KEGG,10)->LL_specific_genes_KEGG_p1;
grid.arrange(HH_specific_genes_KEGG_p1,LL_specific_genes_KEGG_p1,layout_matrix=matrix(c(1,1,1,2,2),nrow=1,ncol=5))

############################################################################################################################
#---------------------------------------for TCGA: M-value and beta-value 
apply(TCGA_CRC_methy[,-c(1,2)],2,function(x){log2(x/(1-x))})->TCGA_CRC_methy_M;
apply(TCGA_CRC_methy_factor[,c(SDC2_Probes,TFPI2_Probes,"SDC2_P","TFPI2_P")],2,function(x){log2(x/(1-x))})->TCGA_CRC_methy_subset_M;
#-----------------
2^(unique(TCGA_CRC_methy_M[,5]))/(1+2^(unique(TCGA_CRC_methy_M[,5])))2^(unique(TCGA_CRC_methy_M[,5]))/(1+2^(unique(TCGA_CRC_methy_M[,5])))->test_beta;
hist(test_beta,breaks=200,prob=T,xaxt='n')
axis(side=1,at=seq(0,1,0.1))
lines(density(test_beta),col="red")
#-
hist(unique(TCGA_CRC_methy_M[,5]+0.2),breaks=200,prob = TRUE,plot=T)->test_hist
axis(side=1,at=seq(-8,8,1))
abline(v=-2.1)
rect(-2.2,0,-2.0,0.2,col="gray")
#---
density(unique(TCGA_CRC_methy_M[,5]+0.2),n=512)->test_density#5
lines(test_density,col="red")
#--------------------------mixtools: fitting a Gaussian mixed linear model
library(mixtools);
unique(TCGA_CRC_methy_M[,5])->test_data
normalmixEM(test_data, mu = c(-4, 2), sigma = c(2, 2), mean.constr = c(NA,NA),sd.constr = c(NA, NA))->test_data_emfit;
plot(test_data_emfit, whichplots = 2,xaxt='n',breaks=200)
axis(side=1,at=seq(-6,6,1))
abline(v=-2.04)
rect(-2.16,0,-1.96,0.2,col="gray")












