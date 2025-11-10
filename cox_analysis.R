library('survival')
library('forestplot')
library(ggplot2)
library(survminer)

colon <-  read.csv('/Users/yuw/Documents/Development/Colon_Data/data/read.csv')

###########Select effective rows#############
colon$til_status <- factor(colon$til_status, levels=c('low', 'high'))
colon$MSI <- factor(colon$MSI, levels=c('MSS/MSI-L', 'MSI-H'))
colon$AJCC_stage <- factor(colon$AJCC_stage, levels=c('Stages I-II', 'Stages III-IV'))
colon$tumorSize_status <-factor(colon$tumorSize_status, levels=c('low', 'high'),
                                labels = c("≤ 4.9cm", "> 4.9cm"))
colon$histoGrade <- factor(colon$histoGrade, levels=c('Well-Mod', 'Poorly'))
colon$immuneSubtype <- factor(colon$immuneSubtype, levels=c('C1', 'C2'))
#colon$risk_status <- factor(colon$risk_status, levels=c('low', 'high'), labels=c('≤ 0', '> 0'))
#colon$buddingScore_status <- factor(colon$buddingScore_status, levels=c('low', 'high'),
#                                    labels = c('Bd1/Bd2', 'Bd3'))

###########Modify the column name#############
names(colon)[names(colon)=='til_status'] <- 'TILs%'
names(colon)[names(colon)=='MSI'] <- 'MSI status'
names(colon)[names(colon)=='AJCC_stage'] <- 'AJCC stages'
names(colon)[names(colon)=='tumorSize_status'] <- 'Tumor size'
names(colon)[names(colon)=='histoGrade'] <- 'Histological grades'
names(colon)[names(colon)=='immuneSubtype'] <- 'Immune subtypes'
#names(colon)[names(colon)=='risk_status'] <- 'Risk scores'
#names(colon)[names(colon)=='buddingScore_status'] <- 'Budding scores'

##########Generate cox forest plot############
png(filename = '/Users/yuw/Desktop/multi-cox-READ-pfi.png',
    width = 13, height = 7, units = "in", res = 300)
model <- coxph(Surv(PFI.month, PFI)~ `TILs%` + `AJCC stages` + 
                 `Tumor size` + `MSI status`+ `Histological grades` + `Immune subtypes`, data=colon)
ggforest(model, colon, main='Progression Free Interval: TCGA READ',fontsize = 1.4)
dev.off()









png(filename = "/Users/yuw/Desktop/os-multiCox-UK.png",
    width = 10, height = 5, units = "in", res = 300)
model <- coxph(Surv(OS.month, OS)~ `TILs%` + `AJCC stages` +
                 `Tumor size` + `Budding scores`, data=colon)
ggforest(model, colon, main='Overall Survival:UKentucky')
dev.off()



png(filename = "/Users/yuw/Desktop/pfi-Cox-TCGA.png",
    width = 10, height = 2, units = "in", res = 300)
names(colon)[names(colon)=='risk_score'] <- 'Risk score'
model <- coxph(Surv(PFI.month, PFI)~ `Risk score`, data=colon)
ggforest(model, colon, main='Progression Free Interval: TCGA')
dev.off()

