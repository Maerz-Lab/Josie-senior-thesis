DipnetData <- DipnetData[!DipnetData$Year %in% c("View(DipnetData)", 
                                                 "DipnetData <- read.csv(ARWMA DipnetSurveyData_Aug-24-2023.csv)" ),]
DipnetData$Year <- str_replace_all(DipnetData$Year, "2021 June", "2021")
##Species Updates
replacements <- c("Acris gryllys"="Acris gryllus", 
                  "Pseudacris crucifer "="Pseudacris crucifer", 
                  "Acris gryllus "="Acris gryllus", 
                  "Pseudacris ocularis "="Pseudacris ocularis", 
                  "Pseudocris ocularis"="Pseudacris ocularis", 
                  "Toad spp."="Anaxyrus spp.", 
                  "Acris crepitans "="Acris crepitans", 
                  "Rana catebesbeiana"="Rana catesbeiana", 
                  "Psuedacris ocularis"="Pseudacris ocularis", 
                  "Notophtalmus viridescens "="Notophthalmus viridescens", 
                  "Scaphiopus holbrooki"="Scaphiopus holbrookii", 
                  "acris crepitans "="Acris crepitans", 
                  "Pseudacris nigrita "="Pseudacris nigrita",
                  "Scaphiopus holbrookiii" = "Scaphiopus holbrookii")
DipnetData$Species <- str_replace_all(DipnetData$Species, replacements)

##DateUpdates
replacement <- c("13/02/2021" = "2/12/2021",
                 "13/03/2021" = "3/13/2021",
                 "31-Jan-19" = "1/31/2019",
                 "1-Mar-19" = "3/1/2019",
                 "5-Apr-19" = "4/5/2019",
                 "3-May-19" = "5/3/2019",
                 "3-Jun-19" = "6/3/2019",
                 "12-Feb-21" = "2/12/2021",
                 "12-Jan-22" = "1/12/2022",
                 "10-Feb-22" = "2/10/2022",
                 "14-Mar-22" = "3/14/2022",
                 "7-Apr-22" = "4/7/2022",
                 "13/02/2023" = "2/13/2023",
                 "14/04/2023" = "4/14/2023",
                 "15/02/2020" = "2/15/2020",
                 "17/03/2023" = "3/17/2023",
                 "18/05/2021" = "5/18/2021",
                 "24/01/2020" = "1/24/2020",
                 "9-Feb-22" = "2/9/2022",
                 "1-Feb-19" = "2/1/2019",
                 "2-Mar-19" = "3/2/2019",
                 "16-Nov-22" = "11/16/2022",
                 "15/12/2021" = "12/15/2021",
                 "18/01/2023" = "1/18/2023",
                 "19/03/2022" = "3/19/2022",
                 "17/5/2021" = "5/17/2021",
                 "17/02/2023" = "2/17/2023",
                 "16/12/2021" = "12/16/2021",
                 "12/7/2022" = "7/12/2022",
                 "9/5/2023" = "5/9/2023",
                 "9/1/2021" = "1/9/2021",
                 "11/12/2020" = "12/11/2020",
                 "9/1/2022" = "1/9/2022",
                 "10/3/2023" = "3/10/2023",
                 "12-Feb-2021" = "2/12/2021",
                 "16-Nov-2022" = "11/16/2022",
                 "7-Apr-2022" = "4/7/2022",
                 "2020-11-12" = "12/11/2020",
                 "2021-07-04" = "4/7/2021",
                 "2021-08-04" = "4/8/2021",
                 "2021-08-01" = "1/8/2021",
                 "2021-09-01" = "1/9/2021",
                 "2021-07-13" = "7/13/2021",
                 "2021-07-14" = "7/14/2021",
                 "2021-12-03" = "3/12/2021",
                 "16/11/2021" = "11/16/2021",
                 "2022-04-13" = "4/13/2022",
                 "2022-09-01" = "1/9/2022",
                 "2022-12-07" = "7/12/2022",
                 "2022-02-06" = "6/2/2022",
                 "2022-03-05" = "5/3/2022",
                 "2022-04-05" = "5/4/2022",
                 "2023-09-06"= "6/9/2023",
                 "2023-10-03" = "3/10/2023",
                 "2023-09-05" = "5/9/2023",
                 "2021-11-06" = "6/11/2021"
)
DipnetData$Date <- str_replace_all(DipnetData$Date, replacement)

# Create the named vector of replacements
sitereplacements <- c("Plum Orchard" = "4",
                      "4 (4)" = "4",
                      "4- 4" = "4",
                      "snot bonnett" = "8", 
                      "4 (Plum Orchard)" = "4", 
                      "Cutgrass" = "10", 
                      "10- Cutgrass" = "10", 
                      "10- 10" = "10",
                      "8- Snot Bonnet" = "8", 
                      "4- Plum Orchard" = "4",
                      "Little mill" = "6",
                      "Little Mill" = "6",
                      "Millfoil" = "6",
                      "Deadwood" = "9",
                      "22/29" = "29",
                      "4 " = "4",
                      "8 " = "8")
##22/29?
##removing ceylon sites
DipnetData <- DipnetData[!DipnetData$WetlandID %in% c("CB4", "CB3", "CB2", "CB1", "CB5", "CB6", "CB7"), ]
# Apply the replacements to the WetlandID column
for (pattern in names(sitereplacements)) {
  replacement <- sitereplacements[pattern]
  DipnetData$WetlandID <- str_replace_all(DipnetData$WetlandID, fixed(pattern), replacement)
}

##Site Updates
replacement1 <- c("Townsend WMA" = "TWMA")
DipnetData$Site <- str_replace_all(DipnetData$Site, replacement1)

#correct m/a in egg data
DipnetData$Egg <- str_replace_all(DipnetData$Egg, "m/a", "0")

DipnetData$Total <- NA


##2024 Data Set
ARWMA_2024 <- read.csv("ARWMA_2024_Dipnet_Data.csv")

##Pond Updates
pondreplacements <- c("X"="3", 
                      "428"="42",
                      "Pond 9" = "8")
ARWMA_2024$Wetland.ID <- str_replace_all(ARWMA_2024$Wetland.ID, pondreplacements)
#removing ceylon data
ARWMA_2024 <- ARWMA_2024[!ARWMA_2024$Wetland.ID %in% c("999", "49", "21", "55", "54", "127", "41B", "109", "69", "59"), ]

unique(ARWMA_2024$Wetland.ID)

datetimeARWMA <- unique(ARWMA_2024$Date.and.Time)




# Convert Date.and.Time to Date and time columns
ARWMA_2024 <- ARWMA_2024 %>%
  mutate(
    Date.and.Time = mdy_hms(Date.and.Time),  # Convert to POSIXct
    Date = format(Date.and.Time, "%Y-%m-%d"), # Reformat to YYYY-MM-DD
    Time = format(Date.and.Time, "%H:%M:%S")    # Reformat to HH:MM:SS (24-hour format)
  )

write.csv(DipnetData, "CleanedDipnetData.csv", row.names = FALSE)
write.csv(ARWMA_2024, "CleanedDipnet2024.csv", row.names = FALSE)
unique(ARWMA_2024$Wetland.ID)
unique(ARWMA_2024)