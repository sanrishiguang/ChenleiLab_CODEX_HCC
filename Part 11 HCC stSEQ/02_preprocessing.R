
####R 4.1.0


library(Seurat)  #‘4.1.1’
library(SeuratData) #‘0.2.2’
library(ggplot2)   #‘3.3.6’
library(patchwork) #‘1.1.1’
library(dplyr)   #‘1.0.9’
library(SPATA2)   #‘0.1.0’


filename <- dir('./Spata Object')
pathname <- sapply(filename, function(x){strsplit(x, '[.]')[[1]][1]})
pathname
length(filename)



for (i in c(21:29)) {
  dir.create(pathname[i])
  spata_obj <- readRDS(file = paste0('./Spata Object/', filename[i]))
  
  #denoising
  spata_obj <-
    runAutoencoderAssessment(
      object = spata_obj,
      activations = c("relu", "selu", "sigmoid"),
      #activations = "selu",
      bottlenecks = c(32, 40, 48, 56, 64),
      epochs = 20,
      layers = c(128, 64, 32),
      dropout = 0.1
    )
  
  assdf <- getAutoencoderAssessment(spata_obj)$df
  activation_name <- assdf$activation[which(assdf$total_var == max(assdf$total_var) )]
  bottlenecks_name <- assdf$bottleneck[which(assdf$total_var == max(assdf$total_var) )] %>% as.character()
  
  spata_obj <-
    runAutoencoderDenoising(
      object = spata_obj,
      activation = activation_name,
      bottleneck = bottlenecks_name,
      epochs = 20,
      layers = c(128, 64, 32),
      dropout = 0.1
    )
  
  saveRDS(spata_obj, file =paste0('./', pathname[i], '/', pathname[i], '_denoising.RDS' ))
  
  
  #DEA
  spata_obj <- 
    runDeAnalysis(object = spata_obj,
                  across = "seurat_clusters", # across which identity groups
                  method_de = c("wilcox") # with which methods
    )
  
  saveRDS(spata_obj, file =paste0('./', pathname[i], '/', pathname[i], '_DEA.RDS' ))
  # 
  #   
  #   #GSEA
  #   spata_obj <- 
  #     runGsea(
  #       object = spata_obj, 
  #       across = "seurat_clusters",
  #       methods_de = "wilcox" 
  #     )
  #   
  #   saveRDS(spata_obj, file =paste0('./', pathname[i], '/', pathname[i], '_GSEA.RDS' ))
  #   
  print(paste0("have processed", i, " files", ", filename[i]))
  rm(spata_obj)
  gc()
  
}

