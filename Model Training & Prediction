#FUNCTIONS ONLY
#load packages
library(matrixStats)
library(randomForest)
library(pROC)
library(ggplot2)
library(pheatmap)
#generate data function
generate_data<-function(data,target_snp){
  duplicated_rows<-duplicated(data[,2:ncol(data)])
  data<-data[!duplicated_rows,]
  x_data<-data[,2:ncol(data)]
  y_data<-data[,1]
  x_data_standard<-x_data
  for(i in 1:ncol(x_data_standard)){
    average<-mean(x_data[,i])
    stdev<-sd(x_data[,i])
    for(j in x_data_standard[,i]){
      x_data_standard[j,i]<-(x_data[j,i]-average)/stdev
    }
  }
  if(target_snp==ncol(data)){
    target_snp = target_snp-1
  }
  column<-target_snp
  trial_column<-x_data[,target_snp]
  trial_data<-x_data_standard
  trial_data[,target_snp]<-trial_column
  colnames(trial_data)[column]<-"target"
  trial_indexes<-setdiff(1:ncol(trial_data),column)
  new_data_package<-list(
    new_ind=trial_indexes,
    new_data=trial_data,
    target_col = target_snp,
    output_data = y_data
  )
  return(new_data_package)
}

#build model original
build_model_orig<-function(new_data_package){
  trial_data<-new_data_package$new_data
  column<-new_data_package$target_col
  raw_data_indexes<-1:nrow(trial_data)
  sample_size<-floor(floor(length(raw_data_indexes)/5)/2)
  predicted_differences<-c()
  correlations<-c()
  for(i in 1:5){
    training_indexs<-c()
    predicted_differences<-c()
    start_pos<-((i-1)*sample_size)+1
    end_pos<-start_pos+sample_size
    train_test_indexes<-raw_data_indexes[start_pos:end_pos]
    training_indexes_1<-setdiff(raw_data_indexes,c(train_test_indexes,train_test_indexes+floor(length(raw_data_indexes)/2)))
    for(i in seq(1,length(training_indexes_1)/2)){
      training_indexs<-c(training_indexs,sample(c(training_indexes_1[i],training_indexes_1[i+(length(training_indexes_1)/2)]),1))
    }
    training_data<-trial_data[training_indexs,]
    rf_model<-randomForest(x=training_data,y=combined_ys[training_indexs])
    
    for(i in train_test_indexes){
      test_data<-rbind(trial_data[i,],trial_data[i+floor(length(raw_data_indexes)/2),])
      prediction<-predict(rf_model,newdata=test_data)
      predicted_difference<-prediction[1]-prediction[2]
      predicted_differences<-c(predicted_differences,predicted_difference)
    }
    correlations<-c(correlations,cor(predicted_differences,actual_y_diffs[train_test_indexes]))
  }
  orig_model_package<-list(
    final_correlation = mean(correlations),
    #orig_predicted_diff = predicted_differences
  )
  return(orig_model_package)
}

    #pca model building function
    build_model_pca<-function(new_data_package,enlarged_target_x,small_target_x){
      samplings<-c(5,10,20,30,40,50,60)
      #unloaded stuff from the input data package
      trial_data<-new_data_package$new_data
      column<-new_data_package$target_col
      trial_indexes<-new_data_package$new_ind
      raw_data_indexes<-1:nrow(trial_data)
      correlations_pca<-c()
      predictions<-rep(0,length(raw_data_indexes))
      combined_elements<-c(enlarged_target_x,small_target_x)
      combined_ys<-new_data_package$output_data
      size<-floor(length(raw_data_indexes)/5)
      size<-size/2
      #START OF TRAINING THE MODEL TO SELECT THE MOST OPTIMAL MODEL
      for(j in samplings){
        correlations_pca_separate<-c()
        transformed_x_controls<-prcomp(trial_data[,trial_indexes],rank=j,center=TRUE,scale=TRUE)
        transformed_data<-cbind(combined_ys,trial_data[,column],transformed_x_controls$x)
        colnames(transformed_data)[2]="target"
        for(n in 1:5){
          training_indexs<-c()
          start_pos<-((n-1)*size)+1
          end_pos<-start_pos+size-1
          train_test_indexes<-raw_data_indexes[start_pos:end_pos]
          training_indexes_1<-setdiff(raw_data_indexes,c(train_test_indexes,train_test_indexes+floor(nrow(transformed_data)/2)))
          training_data<-transformed_data[training_indexes_1,2:ncol(transformed_data)]
          rf_model<-randomForest(x=training_data,y=combined_ys[training_indexes_1])
          test_data<-rbind(transformed_data[train_test_indexes,],transformed_data[train_test_indexes+floor(nrow(transformed_data)/2),])
          
          prediction<-predict(rf_model,newdata=test_data[,2:ncol(test_data)])
          predicted_differences<-prediction[seq(length(train_test_indexes)+1,length(prediction))]-prediction[1:length(train_test_indexes)]
          actual_differences<-transformed_data[nrow(transformed_data)/2+train_test_indexes,1]-transformed_data[1:length(train_test_indexes),1]
          correlations_pca_separate<-c(correlations_pca_separate,cor(predicted_differences,actual_differences))
        }
        correlations_pca<-c(correlations_pca,mean(correlations_pca_separate))
      } 
      #END OF TRAINING THE MODEL TO FIND THE MOST OPTIMAL MODEL
      final_pca_controls<-prcomp(trial_data[,trial_indexes],rank=samplings[which.max(correlations_pca)],center=TRUE,scale=TRUE)
      final_controlled_x<-cbind(trial_data[,column],final_pca_controls$x)
      final_data<-cbind(combined_ys,final_controlled_x)
      colnames(final_data)[2]="target"
      final_predicted_differences<-c()
      predicted_gene_expression<-rep(0,length(raw_data_indexes))
      rf_model<-randomForest(x=final_data[,2:ncol(final_data)],y=final_data[,1])
      for(i in 1:(nrow(final_data)/2)){
        tissue_1<-final_data[i,]
        tissue_1_rep<-tissue_1
        tissue_2<-final_data[i+nrow(final_data)/2,]
        original_predictions<-predict(rf_model,rbind(tissue_1[2:length(tissue_1)],tissue_2[2:length(tissue_2)]))
        tissue_1[2]<-tissue_2[2]
        tissue_2[2]<-tissue_1_rep[2]
        individual_predicts<-predict(rf_model,rbind(tissue_1[2:length(tissue_1)],tissue_2[2:length(tissue_2)]))
        predictions[c(i,i+nrow(final_data)/2)]<-individual_predicts-original_predictions
      }
      tissue_3<-final_data[which(final_data[,2]==enlarged_target_x[length(enlarged_target_x)]),]
      tissue_4<-final_data[which(final_data[,2]==small_target_x[1]),]
      tissue_3_rep<-tissue_3
      if(length(tissue_3)!=ncol(final_data)){
        original_predictions<-predict(rf_model,tissue_3[1,2:ncol(tissue_3)])
        tissue_3[2]<-tissue_4[2]
        tissue_4[2]<-tissue_3_rep[2]
        individual_predicts<-predict(rf_model,tissue_3[1,2:ncol(tissue_3)])
      }else{
        original_predictions<-predict(rf_model,tissue_3[2:length(tissue_3)])
        tissue_3[2]<-tissue_4[2]
        tissue_4[2]<-tissue_3_rep[2]
        individual_predicts<-predict(rf_model,tissue_3[2:length(tissue_3)])
      }
      predictions[c(which(final_data[,2]==enlarged_target_x[length(enlarged_target_x)]))]<-individual_predicts-original_predictions
      predictions<-cbind(rownames(final_data),predictions)
      return(predictions)
    }
  
#generate output function
generate_output<-function(data,snp){
  setwd("C:/Users/kelly/R/GenePrediction/training_data_gwas_genes/training_data_gwas_genes")
  
    generated_data<-generate_data(data,snp)
    target_x<-generated_data$target_col
    generated_data$new_data<-generated_data$new_data[order(generated_data$new_data[,target_x],decreasing=TRUE),]
    generated_data_indexes<-order(generated_data$new_data[,target_x],decreasing=TRUE)
    generated_data$output_data<-generated_data$output_data[generated_data_indexes]
    pairs<-c()
    small_target_x<-c()
    enlarged_target_x<-c()
    new_ys<-rep(0,nrow(generated_data$new_data))
    measured_differences<-rep(0,nrow(generated_data$new_data))
    
    for (i in 1:floor(nrow(generated_data$new_data)/2)){
      pair <- c(generated_data$new_data[i,target_x], generated_data$new_data[nrow(generated_data$new_data)-i+1,target_x])
      pairs <- c(pairs, list(pair))
      small_target_x<-c(small_target_x,pair[2])
      enlarged_target_x<-c(enlarged_target_x,pair[1])
    }
    last_pair <- c(generated_data$new_data[(nrow(generated_data$new_data))/2+1,target_x],
                   generated_data$new_data[nrow(generated_data$new_data),target_x])
    pairs <- c(pairs, list(last_pair))
    enlarged_target_x<-c(enlarged_target_x,last_pair[1])
    pca_model_package<-build_model_pca(generated_data,enlarged_target_x,small_target_x)
    pca_model_package[,2]<-as.numeric(pca_model_package[,2])
    pca_model_package<-pca_model_package[order(pca_model_package[,2],decreasing=TRUE),]
  
  
  #analyze results

  violin_tissues<-c()
  violin_nums<-c()
  mean_values<-c()
  
  for(i in unique(bulk_tissue$X.2)){
    tissue_indexes<-which(bulk_tissue$X.2==i)
    tissues<-bulk_tissue$Accession_RNA[tissue_indexes]
    violin_nums<-c(violin_nums,pca_model_package[which(pca_model_package[,1]%in%tissues),2])
    violin_tissues<-c(violin_tissues,rep(i,length(pca_model_package[which(pca_model_package[,1]%in%tissues),2])))
  }
  
  violin_data_hypertension<-cbind(violin_tissues,violin_nums)
  violin_data<-data.frame(
    group=c(violin_data_hypertension[,1]),
    values=c(as.numeric(violin_data_hypertension[,2]))
  )
  
  violin_data$group <- reorder(violin_data$group, -violin_data$values)
  #for(i in unique(violin_tissues)){
  #  values<-violin_nums[which(violin_tissues==i)]
  #  mean_values<-c(mean_values,mean(as.numeric(values)))
  #}
  #mean_data<-cbind(unique(violin_tissues),mean_values)
    #violin_data_package<-list(
    #violin_data=violin_data,
    #identifier=full_list_genes$Identifier[gene],
    #gene_access=full_list_genes$Genes[gene],
    #trait=disease_gwas$DISEASE.TRAIT[which(disease_gwas$CHR_POS==full_list_genes$SNP.Position[gene]&disease_gwas$CHR_ID==full_list_genes$Chr.Number[gene])[1]]
    #)
  
    return(violin_data)
}

#Main:

master_snp_data<-read.csv("C:/Users/kelly/R/GenePrediction/matched_snp_ccres/qualified_genes1.csv")
setwd("C:/Users/kelly/R/GenePrediction/training_data_gwas_genes/training_data_gwas_genes")
bulk_tissue<-read.csv("C:/Users/kelly/R/GenePrediction/training_data_gwas_genes/training_data_gwas_genes/meta/bulk_tissue.csv")
snp_map<-c()
for(i in 1:nrow(master_snp_data)){
  data<-readRDS(paste0(master_snp_data$Genes[i],"_bulk.rds"))
  snp<-which(colnames(data)==master_snp_data$Identifier[i])
  generated_output<-generate_output(data,snp)
  if(i==1){
    snp_map<-rbind(snp_map,generated_output[,1])
    snp_map<-rbind(snp_map,generated_output[,2])
  }else{
    snp_map<-rbind(snp_map,generated_output[,2])
  }
}
