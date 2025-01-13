temp_marker <- "Yb171Di"

# I tested every individual FCS file and found that "WT_3_IgG.fcs" and "OLFR2_KO_8_PD1.fcs" have problems with this marker. Probably low-expressed cells. 
# So, I did estParamFlowVS once with other markers with all FCS files. And another time with only "Yb171Di" and without those two FCS files. Then, I combined the two tables. 
fs
View(fs[c(1:7, 9:11, 13:18)])

sink(paste0(fcs_files_folder, "/estParamFlowVS_out.txt"))
cofactors = estParamFlowVS(fs[c(1:7, 9:11, 13:18)],channels=temp_marker)
sink()


# save optimum cofactors
write.csv(merge(panel[panel$fcs_colname %in% temp_marker, ], 
                data.frame(fcs_colname = temp_marker, opt_cofactors = round(cofactors, 2))),
          paste0(fcs_files_folder, "/opt.cofactors_Yb171Di.csv"))



setwd("/Users/hdivakaran/Desktop/CyTOF_data_analysis/Layne_WT_Olfr2KO_IgG_PD1_11152024/fcs_files_folder")
other_markers <- read_csv("fcs_files_folder/opt.cofactors_others.csv")
Yb171Di <- read.csv("fcs_files_folder/opt.cofactors_Yb171Di.csv")

all_markers <- rbind(other_markers, Yb171Di)

markers_to_be_transformed <- panel$fcs_colname[panel$marker_class == "type"]

# make sure the order of markers are the same as "markers_to_be_transformed"
all_markers <- all_markers[match(markers_to_be_transformed, all_markers$fcs_colname), ]
all_markers

all(all_markers$fcs_colname == markers_to_be_transformed)

# save the cofactors 
write.csv(all_markers, "opt.cofactors_all.csv")

# cofactors to transform the data
cofactors <- all_markers$opt_cofactors

# continue the rest of the pipeline