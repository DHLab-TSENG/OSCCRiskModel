
modelList<-c("XGBoost","RF","SVM","LR","Stacking")
AgeList<-c("O","T","B")
MFIList<-c("O","B","L","S")
HyperparameterAll<-NULL
HyperparameterBest<-NULL
for(model in modelList){
  for(MFI in MFIList){
    for(Age in AgeList){
      filename_all<-paste0("Output/",model,"allPar_",MFI,Age,"_new.csv")
      filename_best<-paste0("Output/",model,"bestPar_",MFI,Age,"_new.csv")
      if(file.exists(filename_all)){
        para_all<-read.csv(filename_all)
        para_unique<-para_all %>% dplyr::select(-Class,-Algorithm) %>% unique()
        for (col in colnames(para_unique)){
          r<-para_unique %>% dplyr::select(col) %>% unique() %>% arrange() %>%
            rename(value=col)
          r$par<-col
          if(MFI=="O"){
            r$MFI<-"Original"
          }else if(MFI=="B"){
            r$MFI<-"Binary"
          }else if(MFI=="L"){
            r$MFI<-"Logarithmic"
          }else if(MFI=="S"){
            r$MFI<-"Standardized"
          }
          if(Age=="O"){
            r$Age<-"Original"
          }else if(Age=="B"){
            r$Age<-"Binary"
          }else if(Age=="T"){
            r$Age<-"Ternary"
          }else if(Age=="S"){
            r$Age<-"Standardized"
          }
          r$Algorithm<-model
          if(model=="STK"){
            r$Algorithm<-"Stacking"
          }
          HyperparameterAll<-rbind(HyperparameterAll,r)
        }
      }else{
        print(paste(filename_all,"not exist"))
      }
      if(file.exists(filename_best)){
        para_best<-read.csv(filename_best)
        parList<-colnames(para_best)[c((-length(colnames(para_best))+1):(-length(colnames(para_best))))]
        for (col in parList){
          r<-para_best %>% dplyr::select(col) %>% group_by_all() %>% 
            summarise(Count=n()) %>% arrange(desc(Count)) %>% rename(value=col)
          r$par<-col
          if(MFI=="O"){
            r$MFI<-"Original"
          }else if(MFI=="B"){
            r$MFI<-"Binary"
          }else if(MFI=="L"){
            r$MFI<-"Logarithmic"
          }else if(MFI=="S"){
            r$MFI<-"Standardized"
          }
          if(Age=="O"){
            r$Age<-"Original"
          }else if(Age=="B"){
            r$Age<-"Binary"
          }else if(Age=="T"){
            r$Age<-"Ternary"
          }else if(Age=="S"){
            r$Age<-"Standardized"
          }
          r$Algorithm<-model
          if(model=="STK"){
            r$Algorithm<-"Stacking"
          }
          HyperparameterBest<-rbind(HyperparameterBest,r)
        }
      }else{
        print(paste(filename_best,"not exist"))
      }
      
    }
  }
}


write.csv(HyperparameterAll,"HyperparameterAll.csv",row.names = F)
HyperparameterAll %>% group_by(Algorithm,par,value) %>%
  summarise(Count=n()) %>% arrange(Algorithm,par,as.numeric(value)) %>% 
  mutate(all = paste0(value, collapse = ", ")) %>%
  mutate(all=ifelse(par=="sigma",paste0(round(range(as.numeric(value)),3), collapse = ", "),all)) %>%
  write.csv("HyperparameterAll_sum.csv",row.names = F)
write.csv(HyperparameterBest,"HyperparameterBest.csv",row.names = F)
HyperparameterBest %>% group_by(Algorithm,par,value) %>%
  summarise(CountSum=sum(Count)) %>%
  arrange(Algorithm,par,desc(CountSum)) %>%
  dplyr::slice(1)%>% 
  write.csv("HyperparameterBest_sum.csv",row.names = F)
  