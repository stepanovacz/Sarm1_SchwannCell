#load the libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(SeuratObject)

#Read in and merge all sham data for 1 dpi
Sham.combined_1dpi <- merge(Sham.singletB1, y = c(Sham.singletB2, Sham.singletB3, Sham.singletB4), 
                            add.cell.ids = c("Sham B1", "Sham B2", "Sham B3", "Sham B4"), 
                            project = "Sham")
#Read in and merge all 1dpi data
Injury.combined_1dpi <- merge(Injury.singletB1, y = c(Injury.singletB2, Injury.singletB3, Injury.singletB4),
                              add.cell.ids = c("Injury B1","Injury B2","Injury B4", "Injury B5"), 
                              project = "1dpi")

#Read in the Naive data and the 2hours post injury data
Cond.combined <- merge(x = Sham.combined_1dpi,
                       y = c(Injury.combined_1dpi,
                             Sham_for_TwoHpi_singlet.clean,
                             Two_hpi_hashtag_singlets.clean,
                             Naive.singlet),
                       add.cell.ids = c("Sham_1dpi",
                                        "Injury_1dpi",
                                        "Sham_2hpi",
                                        "Injury_2hpi",
                                        "Naive"),
                       project = "Combined")

#Count total number of cells present
length(unique(WhichCells(Cond.combined))) #60802

#Now select only WT and Sarm1KOs (other data is present here and not needed for this project)
unique(Cond.combined@meta.data$MULTI_ID)
Idents(Cond.combined)<- "MULTI_ID"
Cond.combined <- subset(Cond.combined, 
                        subset = MULTI_ID %in% 
                          unique(MULTI_ID)[!grepl("DR6", unique(MULTI_ID))])
#Count total number of cells
length(unique(WhichCells(Cond.combined)))#verify that have less cells #45176
unique(Cond.combined@meta.data$MULTI_ID)#verify that no DR6 data is here


# Create a unified layer
DefaultAssay(Cond.combined) <- "RNA"
Cond.combined <- JoinLayers(Cond.combined)

# Verify the names of each 10x lane
unique(Cond.combined@meta.data$orig.ident)

######
DefaultAssay(Cond.combined) <- "RNA"

Cond.combined[["RNA"]] <- split(Cond.combined[["RNA"]], f = Cond.combined$orig.ident)
Cond.combined
Cond.combined <- NormalizeData(Cond.combined)
Cond.combined <- FindVariableFeatures(Cond.combined)
Cond.combined <- ScaleData(Cond.combined)
Cond.combined <- RunPCA(Cond.combined)

# CCA integration
Cond.combineds_integrated <- IntegrateLayers(object = Cond.combined, method = CCAIntegration, 
                                             orig.reduction = "pca", new.reduction = "integrated.cca",
                                             verbose = FALSE,k.weight=30)

# re-join layers after integration
Cond.combineds_integrated[["RNA"]] <- JoinLayers(Cond.combineds_integrated[["RNA"]])
ElbowPlot(Cond.combineds_integrated, n=50)

x<-25
Cond.combineds_integrated <- FindNeighbors(Cond.combineds_integrated, reduction = "integrated.cca", dims = 1:x)
Cond.combineds_integrated <- FindClusters(Cond.combineds_integrated, resolution = 0.5)#0.4
Cond.combineds_integrated <- RunUMAP(Cond.combineds_integrated, dims = 1:x, reduction = "integrated.cca")


##properly name each sample, current names are difficult to understand
##Do not change this as only the authors know the nomenclature
unique(Cond.combineds_integrated@meta.data$MULTI_ID)

#properly name each sample
mapping <- c("Sham-WT-Male"= 'Sham-WT-Male',
             'Sham-Sarm1KO-male'="Sham-Sarm1-KO-Male",
             "Sham-Sarm1KO-Male"="Sham-Sarm1-KO-Male",
             "Sham-Sarm-Male" ="Sham-Sarm1-KO-Male",
             "Sham-WT-Female" = 'Sham-WT-Female',
             "Sham-Sarm-Female"= "Sham-Sarm1-KO-Female",
             "Sham-WT-male" = "Sham-WT-Male",
             "Sham-Sarm1KO-male" = "Sham-Sarm1-KO-Male", 
             "Sham-WT-female" = "Sham-WT-Female",
             "Sham-Sarm1KO-female" = "Sham-Sarm1-KO-Female",
             "Female-10-min" = "Naive-Female",
             "Male-10-min" = "Naive-Male",
             "Male-30-min"='Naive-Male',
             "Female-30-min"="Naive-Female",
             "Sham-Sarm1KO-Male" = "Sham-Sarm1-KO-Male",
             "Sham-Sarm1KO-Female"="Sham-Sarm1-KO-Female",
             "Injury-Sarm-Female" = "1dpi-Sarm1-KO-Female",
             "Injury-Sarm-Male" = "1dpi-Sarm1-KO-Male",
             "Injury-WT-Female" = "1dpi-WT-Female",
             "Injury-WT-Male"= "1dpi-WT-Male",
             "Injury-Sarm1KO-Female" = "1dpi-Sarm1-KO-Female",
             "Injury-Sarm1KO-Male" = "1dpi-Sarm1-KO-Male",
             "3dpi-Sarm1KO-Male-Batch3" = "3dpi-Sarm1-KO-Male",
             "Sham-WT-Female-Batch3" = "Sham-for-3dpi-WT-Female",
             "3dpi-WT-Male-Batch3" = "3dpi-WT-Male",
             "3dpi-Sarm1KO-Female-Batch3" ="3dpi-Sarm1-KO-Female",
             "Sham-WT-Male-Batch3" = "Sham-for-3dpi-WT-Male",
             "3dpi-WT-Female-Female-Batch3" = "3dpi-WT-Female",
             "Sham-male-for-2hpi-R3" = "Sham-for-2hpi-WT-Male",
             "Sham-female-for-2-hpi-R2" = "Sham-for-2hpi-WT-Female",
             "Sham-male-for-2hpi-R4" = "Sham-for-2hpi-WT-Male",
             "Sham-female-for-2hpi-R1" = "Sham-for-2hpi-WT-Female",
             "2hpi-male-R4" = "2hpi-WT-Male",
             "2hpi-female-R2" = "2hpi-WT-Female",
             "2hpi-female-R1" ="2hpi-WT-Female",
             "2hpi-male-R3" = "2hpi-WT-Male",
             "Injury-Sarm1KO-female"="1dpi-Sarm1-KO-Female",
             "Injury-Sarm1KO-male"="1dpi-Sarm1-KO-Male",
             "Injury-WT-female"="1dpi-WT-Female",
             "Injury-WT-male"="1dpi-WT-Male")

Cond.combineds_integrated@meta.data$Hashtags <- mapping[Cond.combineds_integrated@meta.data$MULTI_ID]
unique(Cond.combineds_integrated@meta.data$Hashtags)
tail(x = Cond.combineds_integrated[[]])


##properly name each sample, current names are difficult to understand
##Do not change this as only the authors know the nomenclature
##add a new column "label" that assigns the name of sample
unique(Cond.combineds_integrated@meta.data$Hashtags)
mapping<-c("Sham-WT-Male"="Sham-for-1dpi-WT",
           "Sham-Sarm1-KO-Male"="Sham-for-1dpi-Sarm1-KO",
           "Sham-WT-Female"="Sham-for-1dpi-WT",
           "Sham-Sarm1-KO-Female"="Sham-for-1dpi-Sarm1-KO",
           "1dpi-Sarm1-KO-Female"="1dpi-Sarm1-KO",
           "1dpi-Sarm1-KO-Male"="1dpi-Sarm1-KO",
           "1dpi-WT-Female"="1dpi-WT",
           "1dpi-WT-Male"="1dpi-WT",
           "3dpi-Sarm1-KO-Male"="3dpi-Sarm1-KO",
           "Sham-for-3dpi-WT-Female"="Sham-for-3dpi-WT",
           "3dpi-WT-Male"="3dpi-WT",
           "3dpi-Sarm1-KO-Female"="3dpi-Sarm1-KO",
           "Sham-for-3dpi-WT-Male"="Sham-for-3dpi-WT",
           "3dpi-WT-Female"="3dpi-WT",
           "Sham-for-2hpi-WT-Male"="Sham-for-2hpi-WT",
           "Sham-for-2hpi-WT-Female"="Sham-for-2hpi-WT",
           "2hpi-WT-Male"="2hpi-WT",
           "2hpi-WT-Female"="2hpi-WT",
           "Naive-Male" = "Naive",
           "Naive-Female" = "Naive")
Cond.combineds_integrated@meta.data$label <- mapping[Cond.combineds_integrated@meta.data$Hashtags]
unique(Cond.combineds_integrated@meta.data$label)
tail(x = Cond.combineds_integrated[[]])

##properly name each sample, current names are difficult to understand
##Do not change this as only the authors know the nomenclature
####Shams and 2hpi, 1dpi, 3dpi information
unique(Cond.combineds_integrated@meta.data$Hashtags)
mapping<-c("Sham-WT-Male"="Sham",
           "Sham-Sarm1-KO-Male"="Sham",
           "Sham-WT-Female"="Sham",
           "Sham-Sarm1-KO-Female"="Sham",
           "1dpi-Sarm1-KO-Female"="1dpi",
           "1dpi-Sarm1-KO-Male"="1dpi",
           "1dpi-WT-Female"="1dpi",
           "1dpi-WT-Male"="1dpi",
           "3dpi-Sarm1-KO-Male"="3dpi",
           "Sham-for-3dpi-WT-Female"="Sham",
           "3dpi-WT-Male"="3dpi",
           "3dpi-Sarm1-KO-Female"="3dpi",
           "Sham-for-3dpi-WT-Male"="Sham",
           "3dpi-WT-Female"="3dpi",
           "Sham-for-2hpi-WT-Male"="Sham",
           "Sham-for-2hpi-WT-Female"="Sham",
           "2hpi-WT-Male"="2hpi",
           "2hpi-WT-Female"="2hpi",
           "Naive-Male" = "Naive",
           "Naive-Female" = "Naive")
Cond.combineds_integrated@meta.data$condition <- mapping[Cond.combineds_integrated@meta.data$Hashtags]
unique(Cond.combineds_integrated@meta.data$condition)
# Check the number of cells for each condition
#length(WhichCells(Cond.combineds_integrated, expression = condition == "Sham"))#19115
#length(WhichCells(Cond.combineds_integrated, expression = condition == "Naive"))#9123
#length(WhichCells(Cond.combineds_integrated, expression = condition == "1dpi"))#14678
#length(WhichCells(Cond.combineds_integrated, expression = condition == "2hpi"))#2260
#length(WhichCells(Cond.combineds_integrated))

# Reorder the levels so that order is logical (progression of injury)
Cond.combineds_integrated$condition <- factor(Cond.combineds_integrated$condition, 
                                              levels = c("Naive","Sham", "2hpi",
                                                         "1dpi"))

##properly name each sample, current names are difficult to understand
##Do not change this as only the authors know the nomenclature
####add batch information
#add a new column batch.orig that assigns what batch a certain sample was from
unique(Cond.combineds_integrated@meta.data$orig.ident)

current.idsB <- c('Sham.singletB1', 'Sham.singletB2', 'Sham.singletB3', 'Sham.singletB4',
                  'Injury.singletB1', 'Injury.singletB2', 'Injury.singletB3','Injury.singlet.B4',
                  "Sham_for2hpi.singlet", "2hpi.singlet","Naive")
new.idsB <- c('Batch1','Batch2',"Batch3",'Batch4', 
              "Batch1","Batch2","Batch3",'Batch4',
              "Batch7","Batch7","Batch5")
batch.orig <- plyr::mapvalues(x = Cond.combineds_integrated@meta.data$orig.ident, from = current.idsB, to = new.idsB)


names(batch.orig) <- colnames(x = Cond.combineds_integrated)
Cond.combineds_integrated <- AddMetaData(
  object = Cond.combineds_integrated,
  metadata = batch.orig,
  col.name = 'batch.orig'
)
unique(Cond.combineds_integrated@meta.data$batch.orig)


##properly name each sample, current names are difficult to understand
##Do not change this as only the authors know the nomenclature

####Add replicate information without splitting by sex
unique(Cond.combineds_integrated@meta.data$MULTI_ID)
unique(Cond.combineds_integrated@meta.data$batch.orig)
# Create new column by pasting MULTI_ID and batch.orig together
Cond.combineds_integrated@meta.data$HTO_classification_batch.orig <- paste(
  Cond.combineds_integrated@meta.data$MULTI_ID,
  Cond.combineds_integrated@meta.data$batch.orig,
  sep = "_")
unique(Cond.combineds_integrated@meta.data$HTO_classification_batch.orig)

mapping <- c(
  "Sham-WT-male_Batch1" = "Rep1",
  "Sham-WT-female_Batch1" = "Rep2", 
  "Sham-Sarm1KO-male_Batch1" = "Rep1",
  "Sham-Sarm1KO-female_Batch1" = "Rep2",
  "Sham-WT-male_Batch2" = "Rep3",
  "Sham-WT-female_Batch2" = "Rep4",
  "Sham-Sarm1KO-male_Batch2" = "Rep4",
  "Sham-Sarm1KO-female_Batch2" = "Rep3",
  
  "Sham-WT-Male_Batch1"="Rep1",
  "Sham-WT-Female_Batch1"="Rep2",
  "Sham-Sarm-Male_Batch1"="Rep1",
  "Sham-Sarm-Female_Batch1"="Rep2",
  "Sham-Sarm-Female_Batch2"="Rep3",
  "Sham-Sarm-Male_Batch2"="Rep4",
  "Sham-WT-Male_Batch2"="Rep3",
  "Sham-WT-Female_Batch2"="Rep4",
  "Sham-Sarm1KO-Male_Batch3"="Rep5",
  "Sham-WT-Female_Batch3"="Rep5",
  "Sham-WT-Male_Batch3"="Rep6",
  "Sham-Sarm1KO-Female_Batch3"="Rep6",
  "Sham-Sarm1KO-Female_Batch4"="Rep7",
  "Sham-WT-Male_Batch4"="Rep7",
  "Sham-WT-Female_Batch4"="Rep8",
  "Sham-Sarm1KO-Male_Batch4" = "Rep8",  
  
  "Injury-WT-female_Batch1" = "Rep1",       
  "Injury-WT-male_Batch1" = "Rep2",          
  "Injury-Sarm1KO-female_Batch1" = "Rep1", 
  "Injury-Sarm1KO-male_Batch1" = "Rep2",    
  "Injury-WT-female_Batch2"="Rep3",
  "Injury-Sarm1KO-female_Batch2"="Rep3",
  "Injury-Sarm1KO-male_Batch2"="Rep4",
  "Injury-WT-male_Batch2"="Rep4",
  "Injury-WT-Male_Batch3"="Rep5",
  "Injury-Sarm1KO-Female_Batch3"="Rep5",
  "Injury-Sarm1KO-Male_Batch3"="Rep6",
  "Injury-WT-Female_Batch3"="Rep6",
  "Injury-Sarm1KO-Male_Batch4"="Rep7",
  "Injury-WT-Male_Batch4"="Rep7",
  "Injury-Sarm1KO-Female_Batch4"="Rep8",
  "Injury-WT-Female_Batch4"="Rep8",
  
  
  "Sham-male-for-2hpi-R3_Batch7"="Rep1",
  "Sham-male-for-2hpi-R4_Batch7"="Rep2",
  "Sham-female-for-2-hpi-R2_Batch7"="Rep3",
  "Sham-female-for-2hpi-R1_Batch7"="Rep4",
  "2hpi-female-R1_Batch7"="Rep1",
  "2hpi-male-R4_Batch7"="Rep2",
  "2hpi-female-R2_Batch7"="Rep3",
  "2hpi-male-R3_Batch7"="Rep4",
  
  "Male-30-min_Batch5"="Rep1",
  "Male-10-min_Batch5"="Rep2",
  "Female-10-min_Batch5"="Rep3",
  "Female-30-min_Batch5"="Rep4"
  
)

Cond.combineds_integrated@meta.data$replicate<- mapping[Cond.combineds_integrated@meta.data$HTO_classification_batch.orig]
unique(Cond.combineds_integrated@meta.data$replicate)


#Save the Seurat file so you do not have to rerun all the lines of code every time
#saveRDS(Cond.combineds_integrated, "Data_Files/Files/Cond.combineds_integrated.rds")


