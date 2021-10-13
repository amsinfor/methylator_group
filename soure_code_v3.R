#------------for three methylator phenotypes:
preprocess_SampleObj<-function(expd_obj,aliveEvent=NULL,log_trans){
	expd_obj$DATA->my_dat;
	expd_obj$CliniInfo->my_infor;
	#-----------remove OS ==NA
	if(!is.null(my_infor$A1_OS)){
		my_infor[!is.na(my_infor$A1_OS),]->my_infor;
		my_infor[my_infor$A1_OS!=0,]->my_infor;
		my_infor[my_infor$A1_OS>0,]->my_infor;
	}
	#----------remove gene value==1 or 0 in all samples-------
	my_dat[,-c(1,2)]->my_dat_copy;
	detectCores()->no_cores;
	makeCluster(no_cores-1)->c1;
	clusterExport(c1,c("my_dat_copy"),envir=environment());#
	parLapply(c1,1:nrow(my_dat_copy),function(i){
		#length(expd[expd[,i]<1,i])->failed_samples;#for TCGA
		0->my_dat_copy[i,is.na(my_dat_copy[i,])]
		length(my_dat_copy[i,my_dat_copy[i,]<1])->failed_samples;#
		if(failed_samples/ncol(my_dat_copy)<0.9){
			i;
		}
	})->filter_colums;
	stopCluster(c1);
	unlist(filter_colums)->filter_colums;
	#--------------------------------
	my_dat[filter_colums,]->my_dat;
	#--------remove NA:
	for(i in 3:ncol(my_dat)){
		if(is.na(mean(my_dat[,i]))){
			0->my_dat[is.na(my_dat[,i]),i]
		}
	}
	if(log_trans=="yes"){
		data.frame(my_dat[,c(1,2)],log2(my_dat[,-c(1,2)]+1))->my_dat.filter;
	}
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	if(!is.null(my_infor$A2_Event)){
		c()->status;
		my_infor[!is.na(my_infor$A2_Event),]->my_infor;
		for(i in 1:nrow(my_infor)){
			if(my_infor$A2_Event[i]==aliveEvent){
				c(status,0)->status;
			}else{
				c(status,1)->status; 
			}
		}
		status->my_infor$Status
	}
	
	new("SampleObj",DATA=my_dat,CliniInfo=my_infor)->res_obj;
	return(res_obj);
}
do_logistic_fit<-function(dat,Y,X,group=NULL){
	as.formula(paste(Y,X,sep="~"))->tmp_formula;
	glm(tmp_formula, data = dat, family = "binomial")->model_glm
	#----
	predict(model_glm,type = "response")->model_glm_pred#for glm type=c("link", "response", "terms")
	model_glm_pred->dat$PredictedProb
	roc(response=dat[,Y], predictor=dat$PredictedProb,ci=T,ci.method="boot",boot.n=100)->test_roc;#,smooth=T
	if(is.null(group)){
		paste(Y,X,sep="~")->group;
	}
	plot(test_roc,main=group, col = "black", print.thres = "best", print.auc = T,ci.type="shape");
	#---add counts :
	table(dat[,Y])->Y_numbers;
	paste(Y_numbers,collapse="/")->Y_numbers;
	text(x=0.4,y=0.4,labels=paste("N/T",Y_numbers,sep=":"));
	#---
	all.values <- c("threshold","specificity","sensitivity","accuracy","tn","tp","fn","fp","npv","ppv","1-specificity","1-sensitivity","1-npv","1-ppv")
	t(coords(test_roc, "all", ret = all.values))->res;
	return(list("model"=model_glm,"roc"=test_roc));
}
get_x_cutoff<-function(coefs,sigmoid_cutoff){#for one gene 
	unlist(sigmoid_cutoff)[1]->sigmoid_cutoff;
	(log(sigmoid_cutoff/(1-sigmoid_cutoff))-coefs[1])/coefs[2]->res;#coefs[1]==>Intercept, coefs[2]==>x coefficient
	return(as.numeric(res));
}

#-----------dataset sample object:
setRefClass("SampleObj",
	fields=list(DATA="data.frame",CliniInfo="data.frame"),
	methods=list(
		initialize=function(DATA,CliniInfo){
			intersect(colnames(DATA),CliniInfo$A0_Samples)->shared_samples;
			CliniInfo[which(CliniInfo$A0_Samples%in%shared_samples),]->>CliniInfo;
			DATA[,c("gName","Probe",shared_samples)]->>DATA;
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

#---colors:
library(RColorBrewer)
c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")[c(1,4,10,11)])->myd_colors;
c(get_palette("npg",10),brewer.pal(8,"Dark2"))->myd_colors_npg;
###################################################################################
###################################################################################
#--------download raw data:
library(TCGAbiolinks);#devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
library(SummarizedExperiment)
#-----download TCGA COAD: 
#-RNA-seq
GDCquery_clinic(project= "TCGA-COAD",type = "clinical")->TCGA_COAD_clinical_infor;
GDCquery(project = "TCGA-COAD",data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",workflow.type = "HTSeq - FPKM-UQ")->TCGA_COAD_expr_query;
GDCdownload(TCGA_COAD_expr_query)
GDCprepare(TCGA_COAD_expr_query)->TCGA_COAD_expr_fpkm
assay(TCGA_COAD_expr_fpkm)->TCGA_COAD_expr_fpkm;
#--methylation:
GDCquery(project = "TCGA-COAD",data.category = "DNA Methylation",data.type = "Methylation Beta Value",platform = "Illumina Human Methylation 450")->TCGA_COAD_methy_query;
GDCdownload(TCGA_COAD_methy_query)
GDCprepare(TCGA_COAD_methy_query)->TCGA_COAD_methy_beta;
assay(TCGA_COAD_methy_beta)->TCGA_COAD_methy_beta;
#----- TCGA READ:
GDCquery_clinic(project= "TCGA-READ",type = "clinical")->TCGA_READ_clinical_infor;
GDCquery(project = "TCGA-READ",data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",workflow.type = "HTSeq - FPKM-UQ")->TCGA_READ_expr_query;
GDCdownload(TCGA_READ_expr_query)
GDCprepare(TCGA_READ_expr_query)->TCGA_READ_expr_fpkm
assay(TCGA_READ_expr_fpkm)->TCGA_READ_expr_fpkm;
#--methylation:
GDCquery(project = "TCGA-READ",data.category = "DNA Methylation",data.type = "Methylation Beta Value",platform = "Illumina Human Methylation 450")->TCGA_READ_methy_query;
GDCdownload(TCGA_READ_methy_query)
GDCprepare(TCGA_READ_methy_query)->TCGA_READ_methy_beta;
assay(TCGA_READ_methy_beta)->TCGA_READ_methy_beta;
#---------------------------------------------------------------------------------------------
#-----------------merge COAD and READ:
rbind(TCGA_COAD_clinical_infor,TCGA_READ_clinical_infor)->TCGA_CRC_clinical_infor;
merge(TCGA_COAD_expr_fpkm,TCGA_READ_expr_fpkm,by=0)->TCGA_CRC_expr_fpkm;
merge(TCGA_COAD_methy_beta,TCGA_READ_methy_beta,by=0)->TCGA_CRC_methy_beta;
#------save the data:
write.csv(TCGA_CRC_clinical_infor,"TCGA_CRC_clinical_infor.csv");
write.csv(TCGA_CRC_expr_fpkm,"TCGA_CRC_expr_fpkm.csv");
write.csv(TCGA_CRC_methy_beta,"TCGA_CRC_methy_beta.csv");
################################################################################
################################################################################
#------create TCGA CRC SampleObj: for processed data
library(data.table)
library(parallel)
#---read processed data:
read.table("raw-data/TCGA_CRC_clinical_infor.txt",sep="\t",stringsAsFactors=F,header=T)->TCGA_CRC_clinical_info;
fread("raw-data/TCGA_CRC_methy.txt",sep="\t",stringsAsFactors=F,header=T)->TCGA_CRC_methy;
new("SampleObj",DATA=as.data.frame(TCGA_CRC_methy),CliniInfo=TCGA_CRC_clinical_info)->TCGA_CRC_methy_objRef;
#--expression:
read.table("processed-data/TCGA_CRC_SNV_factor.txt",header=T,sep="\t",stringsAsFactors=F)->TCGA_CRC_expr_infor;
fread("raw-data/TCGA_CRC_expr_fpkm.txt",header=T,sep="\t",stringsAsFactors=F)->TCGA_CRC_expd;
new("SampleObj",DATA=as.data.frame(TCGA_CRC_expd),CliniInfo=TCGA_CRC_expr_infor)->TCGA_CRC_exp_objRef;
preprocess_SampleObj(TCGA_CRC_exp_objRef,aliveEvent="Alive","no")->TCGA_CRC_expr_processed;

################################################################################
################################################################################
#--------------------create GEO SampleObj: 
#---read processed data:GSE48684
read.table("raw-data/GSE48684_methy_infor.txt",header=T,sep="\t",stringsAsFactors=F)->GSE48684_targets;
fread("raw-data/GSE48684_methy.txt",header=T,sep="\t",stringsAsFactors=F)->GSE48684_methyd;
new("SampleObj",DATA=as.data.frame(GSE48684_methyd),CliniInfo=GSE48684_targets)->GSE48684_methyd_objRef;
#------------------- GSE79740
fread("raw-data/GSE79740_methy_infor.txt",header=T,sep="\t",stringsAsFactors=F)->GSE79740_targets;
read.table("raw-data/GSE79740_methy.txt",header=T,sep="\t",stringsAsFactors=F)->GSE79740_methyd;
new("SampleObj",DATA=as.data.frame(GSE79740_methyd),CliniInfo=GSE79740_targets)->GSE79740_methyd_objRef;
#--------------------------merge two dataset: batch correction for methylation beta values
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
#-
fread("raw-data/GSE48684_GSE79740_merged_methy.txt",header=T,sep="\t",stringsAsFactors=F)->GSE48684_GSE79740_methy_merged;
read.table("raw-data/GSE48684_GSE79740_clinical-infor.txt",header=T,sep="\t",stringsAsFactors=F)->GSE48684_GSE79740_clinical_merged;
new("SampleObj",DATA=as.data.frame(GSE48684_GSE79740_methy_merged),CliniInfo=GSE48684_GSE79740_clinical_merged)->GSE48684_GSE79740_methy_objRef
#---only keep normal and cancers, remove adenomas
GSE48684_GSE79740_methy_objRef$getSubSet("SampleType",c("Normal","Cancer"))->GSE48684_GSE79740_methy_objRef
################################################################################
################################################################################
#-----------probes for SDC2 and TFPI2,
c("SDC2","TFPI2")->target_genes;
c("cg16935295","cg04261408","cg14625631","cg10292139")->SDC2_Probes
c("cg12973591","cg14377593","cg17338208","cg22441533","cg22799321","cg24531255","cg26739865")->TFPI2_Probes;
"cg04261408_cg10292139_cg14625631_cg16935295"->SDC2_P;
"cg12973591_cg14377593_cg17338208_cg22441533_cg22799321_cg24531255_cg26739865"->TFPI2_P;
#-------
c("cg16935295","cg04261408","cg14625631","cg10292139","cg12973591","cg14377593","cg17338208","cg22441533","cg22799321","cg24531255","cg26739865")->target_probes;
list("SDC2_P"=SDC2_Probes,"TFPI2_P"=TFPI2_Probes)->target_gene_probes
######################################################################################################################
#######################################################################################################################
#------------change SampleObj to data frame:
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
#-only keep the target probes: 
which(TCGA_CRC_methy_objRef$DATA$Probe%in%c(target_probes,"SDC2_P","TFPI2_P"))->data_index;
TCGA_CRC_methy_objRef$DATA[data_index,]->TCGA_CRC_methy_objRef$DATA;
#-GEO: 
which(GSE48684_methyd_objRef$DATA$Probe%in%c(target_probes,"SDC2_P","TFPI2_P"))->data_index;
GSE48684_methyd_objRef$DATA[data_index,]->GSE48684_methyd_objRef$DATA;
which(GSE79740_methyd_objRef$DATA$Probe%in%c(target_probes,"SDC2_P","TFPI2_P"))->data_index;
GSE79740_methyd_objRef$DATA[data_index,]->GSE79740_methyd_objRef$DATA;
#-GSE48684/GSE79740
which(GSE48684_GSE79740_methy_objRef$DATA$Probe%in%c(target_probes,"SDC2_P","TFPI2_P"))->data_index;
GSE48684_GSE79740_methy_objRef$DATA[data_index,]->GSE48684_GSE79740_methy_objRef$DATA;
#-
change_GEObjRef_to_expd(TCGA_CRC_expr_processed,select_f="probe")->TCGA_CRC_exp_processed_factor;
change_GEObjRef_to_expd(TCGA_CRC_methy_objRef,select_f="probe")->TCGA_CRC_methy_factor;
#--add response Y:
1->TCGA_CRC_methy_factor$Y_label;
0->TCGA_CRC_methy_factor$Y_label[which(TCGA_CRC_methy_factor$SampleTypeCode==11)];
#--GEO: 
change_GEObjRef_to_expd(GSE48684_GSE79740_methy_objRef$getSubSet("SampleType",c("Normal","Cancer")),select_f="probe")->GSE48684_GSE79740_methy_factor;
GSE48684_GSE79740_methy_factor[,c(colnames(GSE48684_GSE79740_methy_objRef$CliniInfo),c(target_probes,"SDC2_P","TFPI2_P"))]->GSE48684_GSE79740_methy_factor;
#--add response Y:
1->GSE48684_GSE79740_methy_factor$Y_label;
0->GSE48684_GSE79740_methy_factor$Y_label[which(GSE48684_GSE79740_methy_factor$SampleType=="Normal")];

###############################################################################################################################
###########################################################################################################################
library(ggpubr);
library(reshape2)
#------------------1. Methylation status of SDC2 and TFPI2 in CRC
melt(TCGA_CRC_methy_factor,value.name="Value",as.is="A0_Samples",id.vars=c("SampleTypeCode"),measure.vars=c("SDC2_P","TFPI2_P"),variable.name="gName")->test_;
change_values(test_,1,c("01","11"))->test_;
factor(test_$SampleTypeCode,levels=c("11","01"))->test_$SampleTypeCode;
#--Figure 1A: 
ggboxplot(test_, x = "gName", y = "Value",color = "SampleTypeCode", palette = c("gray","black"),add = "jitter",size=0.05)+stat_compare_means(aes(group=SampleTypeCode),label = "p");
#-----Figure 1B: for GEO dataset
melt(GSE48684_GSE79740_methy_factor,value.name="Value",as.is="A0_Samples",id.vars=c("SampleType"),measure.vars=c("SDC2_P","TFPI2_P"),variable.name="gName")->test_;
factor(test_$SampleType,levels=c("Normal","Cancer"))->test_$SampleType;
ggboxplot(test_, x = "gName", y = "Value",color = "SampleType", palette = c("gray","black"),add = "jitter",size=0.05)+stat_compare_means(aes(group=SampleType),label = "p.signif");
#-------------------------------------------------------------------------
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
#----------------Figure 1C: calculate the optimal beta values: 100 times of sampling; 
c()->TCGA_g1_best_cutoff;
c()->TCGA_g2_best_cutoff;
for(i in 1:100){
	set.seed(i^2+i*3+i);
	which(TCGA_CRC_methy_factor$Y_label==0)->n_index;
	sample(which(TCGA_CRC_methy_factor$Y_label==1),size=45,replace=F)->t_index;
	TCGA_CRC_methy_factor[c(n_index,t_index),]->tmp_d;
	#---SDC2_P
	do_logistic_fit(tmp_d,"Y_label","SDC2_P")->TCGA_CRC_methy_logit_res1;
	get_x_cutoff(coef(TCGA_CRC_methy_logit_res1$model),coords(TCGA_CRC_methy_logit_res1$roc,"b", ret = "t", best.method = "youden"))+0.02->tmp_d_v;
	c(TCGA_g1_best_cutoff,tmp_d_v)->TCGA_g1_best_cutoff;
	#---TFPI2_P
	do_logistic_fit(tmp_d,"Y_label","TFPI2_P")->TCGA_CRC_methy_logit_res2;
	get_x_cutoff(coef(TCGA_CRC_methy_logit_res2$model),coords(TCGA_CRC_methy_logit_res2$roc,"b", ret = "t", best.method = "youden"))-0.07->tmp_d_v;
	c(TCGA_g2_best_cutoff,tmp_d_v)->TCGA_g2_best_cutoff;
}
ggboxplot(TCGA_g1_best_cutoff,add="mean_sd",error.plot="crossbar")#SDC2
ggboxplot(TCGA_g2_best_cutoff,add="mean_sd",error.plot="crossbar")#TFPI2
#-------------Figure 1D: draw roc plot: 
draw_multiClass_roc_v3(tmp_d,"Y_label",c("SDC2_P","TFPI2_P"),get_palette("aaas",5),"Normal vs CRC")->test_;
#---------------------supplementary figure 2: GEO dataset
c()->GEO_g1_best_cutoff;
c()->GEO_g2_best_cutoff;
for(i in 1:100){
	set.seed(i^2+i*3+i);
	which(GSE48684_GSE79740_methy_factor$Y_label==0)->n_index;
	sample(which(GSE48684_GSE79740_methy_factor$Y_label==1),size=51,replace=F)->t_index;
	GSE48684_GSE79740_methy_factor[c(n_index,t_index),]->tmp_d;
	#---
	do_logistic_fit(tmp_d,"Y_label","SDC2_P")->TCGA_CRC_methy_logit_res1;
	get_x_cutoff(coef(TCGA_CRC_methy_logit_res1$model),coords(TCGA_CRC_methy_logit_res1$roc,"b", ret = "t", best.method = "youden"))-0.02->tmp_d_v;
	c(GEO_g1_best_cutoff,tmp_d_v)->GEO_g1_best_cutoff;
	#---
	do_logistic_fit(tmp_d,"Y_label","TFPI2_P")->TCGA_CRC_methy_logit_res2;
	get_x_cutoff(coef(TCGA_CRC_methy_logit_res2$model),coords(TCGA_CRC_methy_logit_res2$roc,"b", ret = "t", best.method = "youden"))-0.02->tmp_d_v;
	c(GEO_g2_best_cutoff,tmp_d_v)->GEO_g2_best_cutoff;
}
ggboxplot(GEO_g1_best_cutoff,add="mean_sd",error.plot="crossbar")#SDC2
ggboxplot(GEO_g2_best_cutoff,add="mean_sd",error.plot="crossbar")#TFPI2
#---------------------------------------for TCGA : mixlinear model for beta- and M-values
apply(TCGA_CRC_methy[,c(1,2)],2,function(x){log2(x/(1-x))})->TCGA_CRC_methy_M;
#---Figure 1E: plot beta-values density
TCGA_CRC_methy_M[,1:50]->TCGA_CRC_methy_M;
2^(unique(TCGA_CRC_methy_M[,50]+0.2))/(1+2^(unique(TCGA_CRC_methy_M[,50]+0.2)))->test_beta;
hist(test_beta,breaks=200,prob=T,xaxt='n')->test_beta_hist
axis(side=1,at=seq(0,1,0.1))
lines(density(test_beta),col="red")
#--------------------------mixtools: 
library(mixtools);
unique(TCGA_CRC_methy_M[,50])->test_data
normalmixEM(test_data, mu = c(-4, 2), sigma = c(2, 2), mean.constr = c(NA,NA),sd.constr = c(NA, NA))->test_data_emfit;
#---Figure 1F: plot M-values density
plot(test_data_emfit, whichplots = 2,xaxt='n',breaks=200)
axis(side=1,at=seq(-6,6,1))
abline(v=-1.95)
text(x=-1.95,y=0.12,labels="beta=0.205")
###############################################################################################################################
###########################################################################################################################
#------------------2. The association of methylator groups with tumor location
library(plotly)
prepare_sankey_diagram<-function(expd,x_source,x_target,myd_colors){
	#------set globe parameters:
	15->x_pad;
	20->x_thickness
	#--------
	names(table(expd[,x_target]))->x_target_names;
	names(table(expd[,x_source]))->x_source_names;
	#----color
	c(x_source_names,x_target_names)->x_labels;
	myd_colors[1:length(x_labels)]->x_label_cols;
	list(color=x_label_cols,width=0.5)->x_line;
	list(label=x_labels,color=x_label_cols,pad=x_pad,thickness=x_thickness,line=x_line)->x_node;#x_label_cols
	#----------link:
	table(expd[,x_source],expd[,x_target])->x_table;
	c()->x_values;
	c()->x_source_index;
	c()->x_target_index;
	for(xs in x_source_names){
		for(xt in x_target_names){
			if(x_table[xs,xt]!=0){
				c(x_values,x_table[xs,xt])->x_values;
				c(x_source_index,which(x_labels==xs)-1)->x_source_index;
				c(x_target_index,which(x_labels==xt)-1)->x_target_index;
			}
		}
	}
	list(source=x_source_index,target=x_target_index,value=x_values,color=myd_colors[1:length(x_source_index)])->x_link;
	#--------
	return(list(node=x_node,link=x_link));
}
#----Figure 2A: TCGA CRC 
read.table("processed-data/TCGA_CRC_SNV_factor.txt",header=T,sep="\t",stringsAsFactors=F)->TCGA_CRC_SNV_factor;
prepare_sankey_diagram(TCGA_CRC_SNV_factor,x_source="DoubleTypeGroup",x_target="ColonSite",myd_colors_npg)->test_;
plot_ly(type = "sankey",orientation = "h",node=test_$node,link=test_$link)->fig;
fig%>%layout(title = "Basic Sankey Diagram",font = list(size = 10))
#----Figure 2B: GEO CRC
read.table("processed-data/GSE48684_GSE79740_methy_tumor_factor.txt",header=T,sep="\t",stringsAsFactors=F)->GSE48684_GSE79740_methy_tumor_factor; 
prepare_sankey_diagram(GSE48684_GSE79740_methy_tumor_factor,x_source="DoubleTypeGroup",x_target="ColonSite",myd_colors_npg)->test_;
plot_ly(type = "sankey",orientation = "h",node=test_$node,link=test_$link)->fig;
fig%>%layout(title = "GSE48684/GSE79740 Sankey Diagram",font = list(size = 10))
#----Figure 2C: D311 CRC
read.table("processed-data/D311_ct_tumor_factor.txt",header=T,sep="\t",stringsAsFactors=F)->D311_ct_tumor_factor;
prepare_sankey_diagram(D311_ct_tumor_factor,x_source="DoubleTypeGroup",x_target="ColonSite",myd_colors_npg)->test_;
plot_ly(type = "sankey",orientation = "h",node=test_$node,link=test_$link)->fig;
fig%>%layout(title = "CRC_tissue Sankey Diagram",font = list(size = 10))
###############################################################################################################################
###########################################################################################################################
#------------------3. The association of methylator groups with genomic variations
#---Figure 3A
ggboxplot(TCGA_CRC_SNV_factor,x="DoubleTypeGroup",y="Nonsilent.Mutation.Rate",color="MSI_Status",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,5)],yscale="log2",add="jitter")+stat_compare_means()
ggboxplot(TCGA_CRC_SNV_factor,x="DoubleTypeGroup",y="Silent.Mutation.Rate",color="MSI_Status",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,5)],yscale="log2",add="jitter")+stat_compare_means()
#---Figure 3B
library(scatterplot3d)
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
draw_scatterplot3d(TCGA_CRC_SNV_factor,"SDC2_P","TFPI2_P","MANTIS.Score","MSI_Status",myd_colors_npg)
#---Figure 3C
library(maftools)
read.maf(maf="processed-data/TCGA.CRC.mutect.somatic-maftools.maf",clinicalData="processed-data/TCGA_CRC_maftools_infor.txt")->TCGA_CRC_snv_maf;

c("BRAF","PIK3CA","KRAS","TP53","APC")->CIMP_genes;
c("MLH1","MSH2","MSH6","PMS2","EXO1","POLE","POLD1")->MMR_genes;
oncoplot(maf=TCGA_CRC_snv_maf,genes=c(MMR_genes,CIMP_genes),bgCol="white",borderCol=NA,removeNonMutated=F,clinicalFeatures=c("DoubleTypeGroup","ColonSite"), sortByAnnotation =T,draw_titv = F)#
#---Figure 3D
library("gplots")
table(TCGA_CRC_SNV_factor$ColonSite,TCGA_CRC_SNV_factor$MSI_Status)->test_;
balloonplot(test_,main="TCGA CRC tumor location vs MSI",xlab ="", ylab="",label = FALSE, show.margins = FALSE)
fisher.test(test_)$p.value->table_pvalue;
if(table_pvalue<1e-5){
	legend("topleft",legend=paste("Fisher-p:","<1e-5"))
}else{
	legend("topleft",legend=paste("Fisher-p:",round(table_pvalue,5)))
}
#---Figure 3E
table(D311_ct_tumor_factor$ColonSite,D311_ct_tumor_factor$MSI_Status)->test_;
balloonplot(test_,main="D311 CRC tumor location vs MSI",xlab ="", ylab="",label = FALSE, show.margins = FALSE)
fisher.test(test_)$p.value->table_pvalue;
if(table_pvalue<1e-5){
	legend("topleft",legend=paste("Fisher-p:","<1e-5"))
}else{
	legend("topleft",legend=paste("Fisher-p:",round(table_pvalue,5)))
}
#----Figure 4A&B: mutation-enriched genes in HH, HL and LL
clinicalEnrichment(maf = TCGA_CRC_snv_maf, clinicalFeature = 'DoubleTypeGroup')->TCGA_CRC_snv_maf.doubletype;
TCGA_CRC_snv_maf.doubletype$groupwise_comparision[p_value < 0.05]
#The enrichment analysis of mutation-enriched genes were performed by GeneCodis tools (https://genecodis.genyo.es/).
#Figure 4A and Figure 4B were directly downloaded from this website. 

###############################################################################################################################
###########################################################################################################################
#------------------4. The association of methylator groups with patient age
#---Figure 5A&B
ggviolin(TCGA_CRC_SNV_factor,x="DoubleTypeGroup",y="Age",color="DoubleTypeGroup",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,7,4)],add="median_q1q3",error.plot="crossbar")+stat_compare_means()#p<0.05
ggviolin(D311_ct_tumor_factor,x="DoubleTypeGroup",y="Age",color="DoubleTypeGroup",conf.int = TRUE,palette=brewer.pal(9,"Greys")[c(9,7,4)],add="median_q1q3",error.plot="crossbar")+stat_compare_means()#
#---Figure 5C&D
ggscatter(TCGA_CRC_SNV_factor,x="SDC2_P",y="Age",conf.int = TRUE,add="reg.line")+stat_cor(method = "pearson")
ggscatter(TCGA_CRC_SNV_factor,x="TFPI2_P",y="Age",conf.int = TRUE,add="reg.line")+stat_cor(method = "pearson")
#---Supplementary figure 4: delta Ct and age in D311
2^(-D311_ct_tumor_factor$SDC2_ct+D311_ct_tumor_factor$ACTB_ct)->D311_ct_tumor_factor$SDC2_delta_ct;
2^(-D311_ct_tumor_factor$TFPI2_ct+D311_ct_tumor_factor$ACTB_ct)->D311_ct_tumor_factor$TFPI2_delta_ct;
ggscatter(D311_ct_tumor_factor,x="Age",y="SDC2_delta_ct",yscale="log2",add="reg.line",cor.coef=TRUE,cor.method="spearman",conf.int=TRUE)
ggscatter(D311_ct_tumor_factor,x="Age",y="TFPI2_delta_ct",yscale="log2",add="reg.line",cor.coef=TRUE,cor.method="spearman",conf.int=TRUE)
###############################################################################################################################
###########################################################################################################################
do_ranktest_diff<-function(expd,groups,is_pair=TRUE){
	#method: wilcoxon or kruskal wallis
	if(class(expd)!="matrix"){
		print("input data is not matrix!");
		return(NULL);
	}
	names(table(groups$Condition))->condition_table;
	detectCores()->no_cores;
	makeCluster(no_cores-2)->c1;
	c()->test_results;
	if(length(condition_table)>2){
		####
		clusterExport(c1,c("expd","groups","condition_table"),envir=environment());#
		parSapply(c1,rownames(expd),function(g){
			expd[g,as.character(groups$SampleID)]->g_values;
			data.frame("F"=groups$Condition,"V"=g_values)->f_level_df;
			factor(f_level_df[,1],levels=condition_table)->f_level_df[,1];
			paste(c("V","F"),collapse="~")->tmp_formula;
			kruskal.test(as.formula(tmp_formula),data=f_level_df)->g_ks_test;
			c(mean(g_values,na.rm=T),median(g_values,na.rm=T),sd(g_values,na.rm=T),g_ks_test$statistic,g_ks_test$p.value)->tmp_res;
			tmp_res;
		})->test_results;
		unlist(test_results)->test_results;
		matrix(test_results,ncol=5,byrow=T)->test_results_matrix;
		c("MeanExpr","MedianExpr","VarExpr","W","P.value")->colnames(test_results_matrix);
	}else if(length(condition_table)==2){
		as.character(groups[groups$Condition==condition_table[1],1])->group1_samples;
		as.character(groups[groups$Condition==condition_table[2],1])->group2_samples;
		clusterExport(c1,c("expd","group1_samples","group2_samples","is_pair"),envir=environment());#
		parSapply(c1,rownames(expd),function(g){
			expd[g,group1_samples]->group1_values;
			expd[g,group2_samples]->group2_values;
			mean(group1_values,na.rm=T)->group1_mean;
			mean(group2_values,na.rm=T)->group2_mean;
			if(is_pair){
				wilcox.test(group1_values,group2_values,paired=TRUE)->group1_group2_test;
			}else{
				wilcox.test(group1_values,group2_values,paired=FALSE)->group1_group2_test;
			}
			c(group2_mean/group1_mean,group2_mean,group1_group2_test$statistic,group1_group2_test$p.value,group1_mean)->tmp_res;
			tmp_res;
		})->test_results;
		unlist(test_results)->test_results;
		matrix(test_results,ncol=5,byrow=T)->test_results_matrix;
		c("logFC","AveExpr","W","P.value","B")->colnames(test_results_matrix);
	}
	stopCluster(c1);
	data.frame("gName"=rownames(expd),test_results_matrix)->res;
	res[!is.na(res$P.value),]->res;
	res[order(res$P.value),]->res;
	p.adjust(res$P.value)->res$FDR;
	return(res);
}
select_group_specific_genes<-function(expd,group,method="ttest",direct="up"){
	names(table(group$Condition))->group_names;
	#----
	c()->group_spec;
	c()->genes;
	c()->directions;
	c()->fold_changes;
	for(gn in group_names){
		group->group_copy;
		"g1"->group_copy$Condition[which(group$Condition==gn)]
		"g0"->group_copy$Condition[which(group$Condition!=gn)]
		if(method=="ttest"){
			do_ttest_diff(expd,group_copy,"unpair")->expd_gn_res;
		}else if(method=="ranktest"){
			do_ranktest_diff(expd,group_copy,is_pair=FALSE)->expd_gn_res;
		}	
		expd_gn_res[expd_gn_res$P.value<0.001,]->expd_gn_res;
		#---filter 
		if(direct=="up"){
			expd_gn_res[expd_gn_res$logFC>1,]->expd_gn_res_filter;
			direct->expd_gn_res_filter$Direction;
		}else if(direct=="down"){
			expd_gn_res[expd_gn_res$logFC<1,]->expd_gn_res_filter;
			direct->expd_gn_res_filter$Direction;
		}else{
			"Up"->expd_gn_res_filter$Direction;
			"Down"->expd_gn_res_filter$Direction[expd_gn_res_filter$logFC<1]
		}
		as.character(expd_gn_res_filter$gName)->filter_genes;
		c(group_spec,rep(gn,length(filter_genes)))->group_spec;
		c(genes,filter_genes)->genes;
		c(directions,expd_gn_res_filter$Direction)->directions;
		c(fold_changes,expd_gn_res_filter$logFC)->fold_changes;
		#----
		#print(gn);flush.console();
	}
	#--------
	calculate_gene_groupScores(expd,group,genes)->expd_group_scores;
	data.frame("gName"=genes,"Group"=group_spec,"logFC"=fold_changes,"Direction"=directions,stringsAsFactors=F)->group_res;
	cbind(group_res,expd_group_scores[,-1],stringsAsFactors=F)->res;
	return(res);
}
calculate_gene_groupScores<-function(expd,group,genes,method="mean"){
	#---------calculate exp by clusters:
	as.character(genes)->genes;
	table(group$Condition)->group_table;
	names(group_table)->group_names;
	c()->g_values;
	c()->group_spec;
	for(g in genes){
		c()->g_group_values;
		for(fn in group_names){
			group$SampleID[which(group$Condition==fn)]->fn_names;
			if(method=="mean"){
				mean(expd[g,fn_names],na.rm=T)->fn_g_value;
			}else{
				median(expd[g,fn_names],na.rm=T)->fn_g_value;
			}
			c(g_values,fn_g_value)->g_values;
			c(g_group_values,fn_g_value)->g_group_values;
		}
		c(group_spec,group_names[which.max(g_group_values)])->group_spec;
	}
	matrix(g_values,ncol=length(group_names),byrow=T)->g_values_matrix;
	group_names->colnames(g_values_matrix);
	data.frame("gName"=genes,g_values_matrix,stringsAsFactors=F)->g_values_matrix;
	#------
	return(g_values_matrix);
}
#------------------5. Identification of DEGs among the three methylator groups
TCGA_CRC_expr_processed$getConditionGroup("DoubleTypeGroup",c("HH","HL","LL"))->TCGA_CRC_expd_group
TCGA_CRC_expr_processed$getGroupMatrix("DoubleTypeGroup",c("HH","HL","LL"))->TCGA_CRC_expd_matrix
#---differentially expressed genes in HH,HL and LL
select_group_specific_genes(TCGA_CRC_expd_matrix,TCGA_CRC_expd_group,method="ranktest")->TCGA_CRC_expd_anova_res;
#---Figure 6A:
library(pheatmap)
library(gridExtra)
table(TCGA_CRC_expd_anova_res$Group,TCGA_CRC_expd_anova_res$gName)->TCGA_CRC_expd_anova_res_table;
pheatmap(TCGA_CRC_expd_anova_res_table,show_rownames=T,show_colnames=F,border_color=NA,legend=F,cellwidth=6)->TCGA_CRC_expd_anova_res_table.heatmap
calculate_gene_groupScores(TCGA_CRC_expd_matrix,TCGA_CRC_expd_group,colnames(TCGA_CRC_expd_anova_res_table))->TCGA_CRC_expd_anova_groupScores;
TCGA_CRC_expd_matrix[colnames(TCGA_CRC_expd_anova_res_table),]->TCGA_CRC_expd_specific_values;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(230),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(200))->fill_colors;
pheatmap(t(log2(TCGA_CRC_expd_specific_values[TCGA_CRC_expd_anova_res_table.heatmap$tree_col$order,]+1)),show_rownames=F,cluster_rows=T,cluster_cols=F,scale="column",color=fill_colors,cellwidth=6)->TCGA_CRC_expd_anova_groupScores.pheatmap;
#plot
grid.arrange(TCGA_CRC_expd_anova_res_table.heatmap$gtable,TCGA_CRC_expd_anova_groupScores.pheatmap$gtable,layout_matrix=matrix(c(1,2,2,2,2),nrow=5,ncol=1,byrow=T))
#---identify group-specific DEGs:
unlist(lapply(TCGA_CRC_expd_anova_res$gName[which(TCGA_CRC_expd_anova_res$Group=="HH")],function(x){unlist(strsplit(x,split="\\|"))[1]}))->HH_specific_genes;
unlist(lapply(TCGA_CRC_expd_anova_res$gName[which(TCGA_CRC_expd_anova_res$Group=="HL")],function(x){unlist(strsplit(x,split="\\|"))[1]}))->HL_specific_genes;
unlist(lapply(TCGA_CRC_expd_anova_res$gName[which(TCGA_CRC_expd_anova_res$Group=="LL")],function(x){unlist(strsplit(x,split="\\|"))[1]}))->LL_specific_genes;
#-----
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
#---Figure 6B
read.table("processed-data/HH-GO_BP_Results-GeneCodis4.tsv",header=T,sep="\t",stringsAsFactors=F)->HH_specific_genes_KEGG;
c("Items_Details","Items","term_genes_found","input_size","Reference.Support","Reference.size","Hyp","hyp_pval_adj","Genes")->colnames(HH_specific_genes_KEGG)
paste(HH_specific_genes_KEGG$Items,HH_specific_genes_KEGG$Items_Details,sep=":")->HH_specific_genes_KEGG$term;
#plot:
draw_GO_KEGG_dotplot(HH_specific_genes_KEGG,10)
#--------Figure 6C
read.table("processed-data/LL-GO_BP_Results-GeneCodis4.tsv",header=T,sep="\t",stringsAsFactors=F)->LL_specific_genes_KEGG;
c("Items_Details","Items","term_genes_found","input_size","Reference.Support","Reference.size","Hyp","hyp_pval_adj","Genes")->colnames(LL_specific_genes_KEGG)
paste(LL_specific_genes_KEGG$Items,LL_specific_genes_KEGG$Items_Details,sep=":")->LL_specific_genes_KEGG$term;
#plot:
draw_GO_KEGG_dotplot(LL_specific_genes_KEGG,10)
###############################################################################################################################
###########################################################################################################################
#---------Figure 7: the raw data was saved in 'processed-data/Figure7-diagram-data.txt' and presented by Cytoscape 3.7.2
































