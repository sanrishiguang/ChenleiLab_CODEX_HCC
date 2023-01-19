

####R 4.1.0
library(Seurat)  #‘4.1.1’
library(SeuratData) #‘0.2.2’
library(ggplot2)   #‘3.3.6’
library(patchwork) #‘1.1.1’
library(dplyr)   #‘1.0.9’
library(SPATA2)   #‘0.1.0’


###SA indicated Data used in Science Advances articles

####SA HCC-1 225847

spata_obj_225847_A <- # nolint
  initiateSpataObject_10X(
    directory_10X = 'D:\\Spatial Transcriptomic-20201223\\225847\\225847-ST\\225847-ST-outs\\225847_A4\\225847_A4', # the directory from which to load the data # nolint
    sample_name = "225847_A"
  )

saveRDS(spata_obj_225847_A, file='./Spata Object/225847_A_T_orig.RDS')





####SA HCC-4 226766 


spata_obj_226766_C <- # nolint
  initiateSpataObject_10X(
    directory_10X = 'D:\\Spatial Transcriptomic-20201223\\226766\\226766-ST\\226766_ST_outs\\226766_C\\226766_C', # the directory from which to load the data # nolint
    sample_name = "226766_C"
  )

saveRDS(spata_obj_226766_C, file='./Spata Object/226766_C_L_orig.RDS')




####SA HCC-5 084717
spata_obj_084717_A <- # nolint
  initiateSpataObject_10X(
    directory_10X = 'D:\\Spatial Transcriptomic-20201223\\084717\\084717-ST\\Project_chenlei_SP_20200513N_20200526N_4s_Visium_084717_outs\\084717_A\\084717_A', # the directory from which to load the data # nolint
    sample_name = "084717_A"
  )

saveRDS(spata_obj_084717_A, file='./Spata Object/084717_A_orig.RDS')





##YangFu  HCC2

counts <- Read10X('./stSEQ data from YangFu/HCC2/raw_feature_bc_matrix') # nolint
images <- Read10X_Image('./stSEQ data from YangFu/HCC2/spatial')
coords_df <- select(images@coordinates, imagerow,imagecol )
coords_df$barcodes <- rownames(coords_df)
coords_df <- coords_df %>% rename(x=imagecol, y=imagerow)
library(EBImage)
xximage <- readImage('./stSEQ data from YangFu/HCC2/spatial/aligned_fiducials.jpg', type = 'jpg', names='YangFu_HCC2_image')

##
spata_obj_HCC2 <- 
  initiateSpataObject_CountMtr(
    coords_df = coords_df, # your coordinate data.frame
    count_mtr = counts, # a matrix with unprocessed count values
    sample_name = "YangFu_HCC2",
    image = xximage
  )


saveRDS(spata_obj_HCC2, file = './stSEQ data from YangFu/YangFu_HCC2.RDS')

plotSurface(object = spata_obj_HCC2, color_by = "seurat_clusters", pt_size = 2, pt_alpha = 1, display_image = T)
plotSurface(object = spata_obj_HCC2, color_by = "nCount_RNA", smooth_span = 0.1, pt_size = 2)




###YangFu  HCC3

counts <- Read10X('./stSEQ data from YangFu/HCC3/raw_feature_bc_matrix') # nolint
images <- Read10X_Image('./stSEQ data from YangFu/HCC3/spatial')
coords_df <- select(images@coordinates, imagerow,imagecol )
coords_df$barcodes <- rownames(coords_df)
coords_df <- coords_df %>% rename(x=imagecol, y=imagerow)
library(EBImage)
xximage <- readImage('./stSEQ data from YangFu/HCC3/spatial/aligned_fiducials.jpg', type = 'jpg', names='YangFu_HCC3')

##
spata_obj_HCC3 <- 
  initiateSpataObject_CountMtr(
    coords_df = coords_df, # your coordinate data.frame
    count_mtr = counts, # a matrix with unprocessed count values
    sample_name = "YangFu_HCC3",
    image = xximage
  )


saveRDS(spata_obj_HCC3, file = './stSEQ data from YangFu/YangFu_HCC3.RDS')

plotSurface(object = spata_obj_HCC3, color_by = "seurat_clusters", pt_size = 2, pt_alpha = 0, display_image = T)
plotSurface(object = spata_obj_HCC3, color_by = "nCount_RNA", smooth_span = 0.1, pt_size = 2)




