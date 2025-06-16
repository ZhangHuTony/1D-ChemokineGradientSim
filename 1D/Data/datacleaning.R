# Script to take raw data and split into a dataframe of just CD4 cells, YFP cells and RFP cells.

library(readxl) #package for reading excel files
raw_df = read_excel('3mm_tumor/6235_3mm X coordinates for cross sections.xlsx', sheet = 1)

##function to split dataset into cancer and CD4 dataframes (removes any rows with multiple markers)
library(dplyr)
clean_raw_data <- function(df) {
  # Clean the data: Trim spaces and convert to lowercase
  df <- df %>%
    rename_with( ~ "pos_x", .cols = 7) %>%
    mutate(
      `AF647_clean` = trimws(tolower(`AF647?`)),
      `YFP_clean` = trimws(tolower(`YFP?`)),
      `RFP_clean` = trimws(tolower(`RFP?`)),
      pos_x = round(pos_x) #round to nearest integer
    ) 
    
  
  # Filter CD4 cells (AF647)
  CD4_Cells <- df %>%
    filter( 
      AF647_clean == "af647" & is.na(YFP_clean) & is.na(RFP_clean)#keep rows with only AF647
    )%>% 
    select(pos_x)
  # Filter YFP cells
  YFP_Cells <- df %>%
    filter(
      (YFP_clean == "yfp" & is.na(RFP_clean) & is.na(AF647_clean) ) #keep rows with only YFP
    )%>%
    
    select(pos_x)
  
  
  # Filter RFP cells
  RFP_Cells <- df %>%
    filter(
      (RFP_clean == "rfp" & is.na(YFP_clean) &  is.na(AF647_clean) ) #keep rows with only RFP
    )%>%
    select(pos_x)
  
  
  return(list(CD4_df = CD4_Cells, YFP_df = YFP_Cells, RFP_df = RFP_Cells))
}

output = clean_raw_data(raw_df)
CD4_Output = output$CD4_df
YFP_Output = output$YFP_df
RFP_Output = output$RFP_df

#Output files to desired location
write.csv(CD4_Output, "3mm_tumor/CrossSection1/CD4_XPos.csv",row.names=FALSE)
write.csv(YFP_Output, "3mm_tumor/CrossSection1/YFP_XPos.csv",row.names=FALSE)
write.csv(RFP_Output, "3mm_tumor/CrossSection1/RFP_XPos.csv",row.names=FALSE)
