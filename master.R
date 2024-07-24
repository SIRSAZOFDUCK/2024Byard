### Interrupted time series analysis of primary care prescribing post-lockdown
### Byard et al. 2024

## SET PATHS AND LOAD PACKAGES -------------

# Clear environment

rm(list=ls())

# Define project directories

path.project <- "C:/Users/sirsa/OneDrive/Documents/2024Byard"
path.data.prescribing <- "E:/EPD_DATA"

# Identify uninstalled packages

list.of.packages <- c("data.table","dplyr","plyr","RColorBrewer","ggplot2","corrplot","rgdal","rgeos","raster","maptools","scales","plotly","latticeExtra","maps","classInt","grid","pals","reshape2","cowplot","sandwich","msm","Cairo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages

library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(corrplot)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)
library(scales)
library(plotly)
library(latticeExtra)
library(maps)
library(classInt)
library(grid)
library(pals)
library(reshape2)
library(cowplot)
library(sandwich)
library(msm)
library(Cairo)
library(lme4)
library(RCurl)
library(arm)
library(stringr)
library(MASS)
library(jtools)
library(Epi) 
library(tsModel)


## LOAD BNF CODES -------------------------

# Set working directory

setwd(path.project)

# Read in master BNF codes spreadsheet

bnf.codes <- read.csv("bnf_codes_table.csv")

# Identify individual drug classes 

bnf.codes.classes <- bnf.codes %>%
  # Ensure chapter code numeric
  mutate(BNF.Chapter.Code = as.numeric(BNF.Chapter.Code)) %>%
  # Remove chapter 15 onwards (anaesthesia and devices)
  filter(BNF.Chapter.Code < 15) %>%
  # Arrange in numerical order of BNF code
  arrange(BNF.Subparagraph.Code) %>%
  # Select code and name columns
  dplyr::select(BNF.Subparagraph, BNF.Subparagraph.Code) %>%
  # Keep unique entires
  distinct() %>%
  # Ensure all entries are 7 digits by adding a leading zero if needed
  mutate(BNF.Subparagraph.Code = ifelse(nchar(BNF.Subparagraph.Code) == 6, str_pad(BNF.Subparagraph.Code, width = 7, pad = "0"), BNF.Subparagraph.Code))
  
# Identify individual drugs

bnf.codes.drugs <- bnf.codes %>%
  # Ensure chapter code numeric
  mutate(BNF.Chapter.Code = as.numeric(BNF.Chapter.Code)) %>%
  # Remove chapter 15 onwards (anaesthesia and devices)
  filter(BNF.Chapter.Code < 15) %>%
  # Arrange in numerical order of BNF code
  arrange(BNF.Chemical.Substance.Code) %>%
  # Select code and name columns
  dplyr::select(BNF.Chemical.Substance, BNF.Chemical.Substance.Code) %>%
  # Keep unique entires
  distinct() %>%
  # Ensure all entries are 9 digits by adding a leading zero if needed
  mutate(BNF.Chemical.Substance.Code = ifelse(nchar(BNF.Chemical.Substance.Code) == 8, str_pad(BNF.Chemical.Substance.Code, width = 9, pad = "0"), BNF.Chemical.Substance.Code))

# Create lists of bnf drug and drug class codes

bnf.codes.drugs.list <- bnf.codes.drugs %>% pull(BNF.Chemical.Substance.Code)
  
bnf.codes.classes.list <- bnf.codes.classes %>% pull(BNF.Subparagraph.Code) 


## LOAD PRESCRIBING DATA ---------------

# Define data months to include

included.months <- list("201501","201502","201503", "201504", "201505", "201506", "201507", "201508", "201509", "201510", "201511", "201512",
                        "201601","201602","201603", "201604", "201605", "201606", "201607", "201608", "201609", "201610", "201611", "201612",
                        "201701","201702","201703", "201704", "201705", "201706", "201707", "201708", "201709", "201710", "201711", "201712",
                        "201801","201802","201803", "201804", "201805", "201806", "201807", "201808", "201809", "201810", "201811", "201812",
                        "201901","201902","201903", "201904", "201905", "201906", "201907", "201908", "201909", "201910", "201911", "201912",
                        "202001","202002","202003", "202004", "202005")

# Create directory to save precribing outputs and results

dir.create(paste0(path.project,"/aggregates"), recursive = TRUE)
dir.create(paste0(path.project,"/results_class"), recursive = TRUE)
dir.create(paste0(path.project,"/results_drugs"), recursive = TRUE)

# Read in data and aggregate by drug class

for(i in included.months){
  setwd(path.data.prescribing) # set to prescribing data folder
  file = paste0("EPD_",i,".csv") # define file naming convention
  data <- fread(file = file, header = T, sep = ",") # read in data file
  data$BNF_CODE <- str_sub(data$BNF_CODE,1,9) # Keep first 9 digits of BNF code (drug) 
  data1 <- setDT(data)[,.(ITEMS = sum(ITEMS)), by = .(BNF_CODE)] # Aggregate
  assign(paste0("data_drug_",i), data1) # keep as dataframe object
  setwd(paste0(path.project,"/aggregates")) # set to project folder
  write.csv(data1, paste0("data_drug_",i,".csv"), row.names = F) # save to folder
  data1$BNF_CODE <- str_sub(data1$BNF_CODE,1,7) # Keep first 7 digits of BNF code (drug class)
  data2 <- setDT(data1)[,.(ITEMS = sum(ITEMS)), by = .(BNF_CODE)] # Aggregate
  assign(paste0("data_class_",i), data2) # keep as dataframe object
  write.csv(data2, paste0("data_class_",i,".csv"), row.names = F) # save to folder
  rm(data, data1, data2)
}


## ADD LIST SIZE DATA -------------

## Add to drug class data

# Get a list of all files that start with "data_class"
files <- list.files(path = paste0(path.project,"/aggregates"), pattern = "^data_class", full.names = TRUE)

# Function to read each file
read_file <- function(file) {
  read.csv(file)
}

# Read all the files into a list
data_list <- lapply(files, read_file)

# Name each element of the list with the corresponding file name (without the path and extension)
names(data_list) <- basename(files)

# Get the names of the list elements
list_names <- names(data_list)

# Add year, month and list size to each file
for(i in 1:length(data_list)){
  
  # Get dataframe
  data1 <- data_list[[i]]
  
  # Identify year and month of data
  year <- str_sub(list_names[i], -10, -7)
  month <- str_sub(list_names[i], -6,-5)
  time <- paste0(year,month) %>% as.numeric()
  
  # Add year and month to prescribing data
  data2 <- data1 %>% 
    mutate(year = year) %>%
    mutate(month = month) %>%
    setNames(c("bnf.code", "items", "year", "month"))
  
  # Get list size for that month
  setwd(paste0(path.project,"/listsize"))
  
  data.listsize <- read.csv(paste0(time,"_listsize.csv"))
  
  # Add list size of all practices to prescribing data - note different formats of list size in different years
  
  if (time <201507){
    
    # Sum list size of each practice
    total.listsize <- sum(data.listsize[6])
    
    # Add to prescribing data
    data3 <- data2 %>%
      mutate(listsize = total.listsize)
    
  } else if (201506 < time  & time < 201707) {
    
    # Sum list size of each practice
    total.listsize <- sum(data.listsize[9])
    
    # Add to prescribing data
    data3 <- data2 %>%
      mutate(listsize = total.listsize)
    
  } else if (201706 < time) {
    
    # Get "all" ages rows only
    data.listsize <- data.listsize %>% 
      filter(AGE == "ALL")
    
    # Sum list size of each practice
    total.listsize <- sum(data.listsize[8])*2
    
    # Add to prescribing data
    data3 <- data2 %>%
      mutate(listsize = total.listsize)
    
  }                         
  
  # Assign to object
  assign(paste0("data_",time), data3) # keep as dataframe object
  
  # remove extraneous data
  rm(data.listsize, data1, data2, data3)
  
}

# Bind all prescribing dataframes  

# List all objects
all_objects <- ls()

# Keep all those starting with "data_2"
df_names <- all_objects[grep("^data_2", all_objects)]

# Get names of all dataframes
df_list <- lapply(df_names, get)

# Bind
alldata <- bind_rows(df_list)

# Save
setwd(path.project)
write.csv(alldata,"aggregate_data_by_class.csv", row.names = F)


## Add to individual drug data

# Get a list of all files that start with "data_class"
files <- list.files(path = paste0(path.project,"/aggregates"), pattern = "^data_drug", full.names = TRUE)

# Function to read each file
read_file <- function(file) {
  read.csv(file)
}

# Read all the files into a list
data_list <- lapply(files, read_file)

# Name each element of the list with the corresponding file name (without the path and extension)
names(data_list) <- basename(files)

# Get the names of the list elements
list_names <- names(data_list)

# Add year, month and list size to each file
for(i in 1:length(data_list)){
  
  # Get dataframe
  data1 <- data_list[[i]]
  
  # Identify year and month of data
  year <- str_sub(list_names[i], -10, -7)
  month <- str_sub(list_names[i], -6,-5)
  time <- paste0(year,month) %>% as.numeric()
  
  # Add year and month to prescribing data
  data2 <- data1 %>% 
    mutate(year = year) %>%
    mutate(month = month) %>%
    setNames(c("bnf.code", "items", "year", "month"))
  
  # Get list size for that month
  setwd(paste0(path.project,"/listsize"))
  
  data.listsize <- read.csv(paste0(time,"_listsize.csv"))
  
  # Add list size of all practices to prescribing data - note different formats of list size in different years
  
  if (time <201507){
    
    # Sum list size of each practice
    total.listsize <- sum(data.listsize[6])
    
    # Add to prescribing data
    data3 <- data2 %>%
      mutate(listsize = total.listsize)
    
  } else if (201506 < time  & time < 201707) {
    
    # Sum list size of each practice
    total.listsize <- sum(data.listsize[9])
    
    # Add to prescribing data
    data3 <- data2 %>%
      mutate(listsize = total.listsize)
    
  } else if (201706 < time) {
    
    # Get "all" ages rows only
    data.listsize <- data.listsize %>% 
      filter(AGE == "ALL")
    
    # Sum list size of each practice
    total.listsize <- sum(data.listsize[8])*2
    
    # Add to prescribing data
    data3 <- data2 %>%
      mutate(listsize = total.listsize)
    
  }                         
  
  # Assign to object
  assign(paste0("data_",time), data3) # keep as dataframe object
  
  # remove extraneous data
  rm(data.listsize, data1, data2, data3)
  
}

## Bind all prescribing dataframes  

# List all objects
all_objects <- ls()

# Keep all those starting with "data_2"
df_names <- all_objects[grep("^data_2", all_objects)]

# Get names of all dataframes
df_list <- lapply(df_names, get)

# Bind
alldata <- bind_rows(df_list)

# Save
setwd(path.project)
write.csv(alldata,"aggregate_data_by_drug.csv", row.names = F)


## DRUG CLASS: PROCESS DATA -----------

# Clear environment (except project paths)
keep <-c("path.project", "path.data.prescribing")
rm(list = setdiff(ls(), keep))


# Read in drug class data
setwd(path.project)
all.data <- read.csv("aggregate_data_by_class.csv") %>%
  filter(bnf.code <1500000) # get rid of chapters 15 onwards

# Create bnf code-to-name lookup
bnf.lookup <- read.csv("bnf_codes_table.csv") %>%
  mutate(bnf.code = BNF.Subparagraph.Code) %>%
  mutate(bnf.name = BNF.Subparagraph) %>%
  dplyr::select(15,16) %>%
  unique()

# Identify unique drug classes
drug.list <- all.data %>% pull(bnf.code) %>% unique()

# Prime a dataframe for the results
column_names <- c("bnf.code", "drug.class", "irr", "lcl", "ucl", "pval")
results <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(results) <- column_names

# Set to results directory
setwd(paste0(path.project,"/results_class"))


## DRUG CLASS: ANALYSE DATA AND SAVE OUTPUTS: LOOP ------------

for (i in 1:length(drug.list)){

data <- all.data %>% 
  filter(bnf.code == drug.list[i]) 

# Check at least 1000 items prescribed each month
low.rx <- data %>%
  filter(items <1000)

if (nrow(low.rx) == 0){

# Calculate prescribing rates per 1000 people per year
data <- data %>%
  mutate(items.per.1000 = 1000*items/listsize)

# Allocate to pre or post-lockdown timeframes
data$intervention <- 0 #starts column of zeros
data$intervention[c(64:65)] <- 1 #add 1 to intervention rows which are after intervention

# Add time column
data$time <- 1:nrow(data) #add "time" starting at 1 from month 1

# Poisson ITS (to include time trend and seasonal adjustment)
model <- glm(items ~ offset(log(listsize)) + intervention + time + harmonic(month,2,12), family=quasipoisson, data)
summary <- summary(model)
cis <- ci.exp(model)
#round(ci.lin(model,Exp=T),3)

# Add key results to results dataframe

  # BNF code
  results[i,1] <- drug.list[i]
  
  # BNF name
  bnf.name <- data.frame(drug.list[i]) %>%
    rename("bnf.code" = 1) %>%
    left_join(bnf.lookup, by = "bnf.code")
  
  results[i,2] <- bnf.name[1,2]
  
  # IRR
  results[i,3] <- round(cis[2,1],3)
  
  # Lower CI
  results[i,4] <- round(cis[2,2],3)
  
  # Upper CI
  results[i,5] <- round(cis[2,3],3)
  
  # P value
  results[i,6] <- round(summary$coefficients["intervention", "Pr(>|t|)"],4)
  

# Reform data frame with 0.1 time units to improve plotting
datanew <- data.frame(listsize=mean(data$listsize),intervention=rep(c(0,1),c(630,90)),
                      time= 1:720/10,month=rep(1:120/10,6))
datanew <- datanew[1:650,]

# Save plot

Cairo(file=paste0("Drug class ",drug.list[i]," ",bnf.name[1,2],".png"), 
      type="png",
      units="in", 
      width=5, 
      height=4, 
      pointsize=8, 
      dpi=1200)

# Identify max precribing rate
ymax <- data %>%
  pull(items.per.1000) %>%
  max()

# Define y-axis upper limit
ymax <- ceiling(ymax/10)*10

# Identify max precribing rate
ymin <- data %>%
  pull(items.per.1000) %>%
  min()

# Define y-axis upper limit
ymin <- floor(ymin/10)*10

# Plot model... ADD: flexible y axis
pred <- predict(model,type="response",datanew)/mean(data$listsize)*10^3
plot(data$items.per.1000,type="n",ylim=c(ymin,ymax),xlab="Year",ylab="Items per 1000 registered patients",
     bty="l",xaxt="n")
rect(66,ymin,63,ymax,col=grey(0.9),border=F)
points(data$items.per.1000,cex=0.7, pch=19)
axis(1,at=0:5*12,labels=F)
axis(1,at=0:5*12,tick=F,labels=2015:2020)
lines(1:650/10,pred,col=2)
title(paste0("Drug class: ", bnf.name[1,2],"\nBNF code: ",bnf.name[1,1]), cex.main = 1)

# Predict and plot the deseasonalised trend
pred2b <- predict(model,type="response",transform(datanew,month=6))/
  mean(data$listsize)*10^3
lines(1:650/10,pred2b,col=3,lty=2)


dev.off()

}  else {
  
  # Do nothing
  
}

}


## DRUG CLASS: SAVE ALL RESULTS ---------------

# Clean results output 

results <- results %>%
  # Remove empty rows
  filter(rowSums(is.na(.)) != ncol(results)) %>%
  # Order by lowest p value
  arrange(pval) %>%
  # Create single column with 95% CIs
  mutate(ci = paste0("(",format(round(lcl,3), nsmall = 3)," - ",format(round(ucl,3), nsmall = 3),")")) %>%
  # Remove defunct columns
  dplyr::select(-c(ucl, lcl)) %>%
  # Reorder columns
  dplyr::select(bnf.code, drug.class, irr, ci, pval) %>%
  rename(
    "BNF code" = bnf.code,
    "Drug class" = drug.class,
    "IRR change" = irr,
    "95% CI" = ci,
    "p value" = pval
  )

# Say if p value is small
results$`p value` <- ifelse(results$`p value` < 0.0001, "< 0.0001", results$`p value`)

# Save
write.csv(results, "Results_drug_class.csv", row.names = F)  
  



