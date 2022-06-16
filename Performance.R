
modelList<-c("XGBoost","RF","SVM","LR","Stacking")
AgeList<-c("O","T","B") #,"OPatient","BPatient"
MFIList<-c("O","B","L","S")
PerformanceAll<-NULL
for(model in modelList){
  for(MFI in MFIList){
    for(Age in AgeList){
      filename<-paste0("Output/",model,"5f_",MFI,Age,"_new.csv")
      if(!file.exists(filename)){
        print(paste(filename,"not exist"))
        next
      }
      perf<-read.csv(filename)
      #perf$Class<-paste0(MFI,Age)
      if(MFI=="O"){
        perf$MFI<-"Original"
      }else if(MFI=="B"){
        perf$MFI<-"Binary"
      }else if(MFI=="L"){
        perf$MFI<-"Logarithmic"
      }else if(MFI=="S"){
        perf$MFI<-"Standardized"
      }
      if(Age=="O"){
        perf$Age<-"Original"
      }else if(Age=="B"){
        perf$Age<-"Binary"
      }else if(Age=="T"){
        perf$Age<-"Ternary"
      }else if(Age=="S"){
        perf$Age<-"Standardized"
      }
      perf$Algorithm<-model
      if(model=="STK"){
        perf$Algorithm<-"Stacking"
      }
      PerformanceAll<-rbind(PerformanceAll,perf)
    }
  }
}

write.csv(PerformanceAll,"Performance_all.csv",row.names = F)
