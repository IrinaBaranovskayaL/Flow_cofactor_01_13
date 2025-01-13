#####
# transforming the raw data using channel-specific cofactor
#####
# make sure you have CD45+, live and singlets
#####
getwd() #check the working directory -it should be in the project folder
wd <- "/Users/ibaranovskaya/Projects/flow_cofactors"
setwd(wd)
getwd()

#renv:: activate() #-will work but the working directory should be the same whre activate.R is located

#?file.path #to write the path, will work for all system WindrowsLinux
#source- run r-script that are located in a directory

source(file.path(getwd(), "renv", 'activate.R'))

renv:: status()

#when I create a new project on new computer. FIle-create project
#renv::init()
#renv::restore() # only if I copied renv.lock file and it did not download all packages automatically

#renv::snapshot() #record all packages that were used in the project



# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("flowVS")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("flowCore")




library(flowVS)
library(flowCore)

# install.packages("ggplot2")
library(ggplot2)

getwd()

#fcs_files_folder<- "~/Data/flowdata/Flow_singlets_cof/PA_mel_12-11-24"
fcs_files_folder<- "~/Data/flowdata/Flow_singlets_cof/BATF3_WT_12-24_singlets_CD45"
setwd(fcs_files_folder)  

# Make a "panel" file we use for CyTOF analysis and mark "type" for every channel you want to transform.
files <- list.files(fcs_files_folder, pattern="fcs")

fs <- read.flowSet(files = files , path = fcs_files_folder , truncate_max_range = FALSE)
colnames(fs)    # The variable name that corresponds to the column names in the flowSet object

x <- list(colnames(fs))
x

my_panel <- data.frame(fcs_colname =  row.names(data.frame(markers =  markernames(fs))), 
                       markers =  markernames(fs), 
                       data.frame(antigen = sub("^[^_]*_", "", markernames(fs))), 
                       marker_class = "")

#my_panel$antigen[my_panel$markers == my_panel$antigen] <- ""
my_panel$antigen[my_panel$antigen %in% c("EQ Bead", "Intercalator-Ir", "Cisplatin")] <- ""
my_panel$marker_class[my_panel$antigen == ""] <- "none"
my_panel$antigen[my_panel$antigen == ""] <- "x"

write.csv(my_panel, "panel.csv")
getwd()
#####
# In the "panel" file, mark "type" as marker_class for every channel you want to transform.
#####
rm(list = ls())

panel_file_name <- "panel_to_use.csv"

getwd()

#?fcs_files_folder<- "~/Data/flowdata/Flow_singlets_cof/PA_mel_12-11-24"
#?setwd(fcs_files_folder)  

fcs_files_folder<- "~/Data/flowdata/Flow_singlets_cof/BATF3_WT_12-24_singlets_CD45"


# read the panel
panel <- read.csv(paste0(fcs_files_folder, "/",  panel_file_name), header = TRUE)
panel


# read fcs files as flowSet
files <- list.files(fcs_files_folder, pattern="fcs")
fs <- read.flowSet(files = files , path = fcs_files_folder , truncate_max_range = FALSE)

## identify optimum cofactor for your marker
markers_to_be_transformed <- panel$fcs_colname[panel$marker_class == "type"]
as.data.frame(markernames(fs))[markers_to_be_transformed,]


# estimate the cofactors
#sink(paste0(fcs_files_folder, "/estParamFlowVS_out.txt"))
sink(file.path(fcs_files_folder, "estParamFlowVS_out.txt"))


#markers_to_be_transformed2 <- markers_to_be_transformed[-7]
#cofactors = estParamFlowVS(fs,channels=markers_to_be_transformed2)



cofactors = estParamFlowVS(fs,channels=markers_to_be_transformed)



sink() #return output into console

# save optimum cofactors
write.csv(merge(panel[panel$fcs_colname %in% markers_to_be_transformed, ], 
                data.frame(fcs_colname = markers_to_be_transformed, opt_cofactors = round(cofactors, 2))),
          paste0(fcs_files_folder, "/opt.cofactors.csv"))

# transform selected channels
fs.VS = transFlowVS(fs, channels=markers_to_be_transformed, cofactors)


# if you want to save the flowSet
# saveRDS(fs.VS, file = paste0(fcs_files_folder, "/fs.VS_transformed.rds"), refhook = NULL)

# move the Bartlett's test plot to your folder
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
dir.create(paste0(fcs_files_folder, "/", "Bartlett"))
file.copy(from=plots.png.paths, to=paste0(fcs_files_folder, "/", "Bartlett"))

# write the transformed data to fcs file
dir.create(paste0(fcs_files_folder, "/", "transformed_fcs"))
flowset_list <- flowSet_to_list(fs.VS)
flowframe_names <- names(flowset_list)
for (flowframe_names_temp in flowframe_names) {
  write.FCS(fs.VS[[flowframe_names_temp]],paste0(fcs_files_folder, "/", "transformed_fcs/", "transformed_",  flowframe_names_temp))
}
