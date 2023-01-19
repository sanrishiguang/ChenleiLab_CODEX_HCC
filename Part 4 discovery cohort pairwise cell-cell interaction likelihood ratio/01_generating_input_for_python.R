
###R 4.1.0

########spit table : for python input######
load("All_CN_TP_0424.Rdata")

unique(All_CN_TP$SP)

data_SP1 <- dplyr::filter(All_CN_TP, SP == "SP_LII")
data_SP2 <- dplyr::filter(All_CN_TP, SP == "SP_HII" )
data_SP3 <- dplyr::filter(All_CN_TP, SP == "SP_Pf")
data_SP4 <- dplyr::filter(All_CN_TP, SP == "SP_MII_Tumor")
data_SP5 <- dplyr::filter(All_CN_TP, SP == "SP_MII_Fibro")

write.csv(data_SP1, file='E:\\myPythonProject\\TrainingData\\AllCN_r50_CN40_0331\\spatialanalysis\\data_SP_LII.csv')
write.csv(data_SP2, file='E:\\myPythonProject\\TrainingData\\AllCN_r50_CN40_0331\\spatialanalysis\\data_SP_HII.csv')
write.csv(data_SP3, file='E:\\myPythonProject\\TrainingData\\AllCN_r50_CN40_0331\\spatialanalysis\\data_SP_Pf.csv')
write.csv(data_SP4, file='E:\\myPythonProject\\TrainingData\\AllCN_r50_CN40_0331\\spatialanalysis\\data_SP_MII_Tumor.csv')
write.csv(data_SP5,file = 'E:\\myPythonProject\\TrainingData\\AllCN_r50_CN40_0331\\spatialanalysis\\data_SP_MII_Fibro.csv')

