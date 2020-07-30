#source("ROMS Data Extraction Functions.R")


devtools::source_url("https://raw.github.com/HaydenSchilling/ROMS-Functions/master/ROMS_Data_Extraction_Functions.R")

mydata <- read.csv("Test points.csv")
FULL_DATA <- Extract_Profile_Data_ROMS(data = mydata, variable = "temp")
head(FULL_DATA)



mydata <- read.csv("Test points.csv")

FULL_DATA <- Extract_2D_Point_Data_ROMS(data = mydata, variable = "Hsbl")
FULL_DATA

