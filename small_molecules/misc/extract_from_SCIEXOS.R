# Script to extract relevant columns from SCIEX OS-exported data

# find .txt file with data
# NOTE: if there is more than one .txt file in working directory, then the correct file will need to be manually selected from the list
data_file <- list.files(pattern = ".txt", full.names = T)

# read in data from SCIEXOS .txt file
data <- read.delim(data_file)

# select for the relevant columns
#data_select <- data[,c(2,4,17,21,38,43,56,64,66,86)]
data_select <- data[,c(3,2,24,30,46,49,58,68,70,91)]

# rename columns to typical column names used in data quality check, data analysis, etc.
colnames(data_select) <- c("sample_no","lcms_name","component_name","mrm","area","height","rt","fwhm","sn","pts_across_pk")

# save file to then edit in spreadsheet program to add sample/project-specific information
write.table(data_select, "exported_data.txt", row.names = F, quote = F, sep = "\t")
