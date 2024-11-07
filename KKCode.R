weather.dat$RH<-str_replace_all(weather.dat$RH,c("37.10%"="37.1"))
weather.dat$RH<-str_replace_all(weather.dat$RH,c("40.00%"="40"))
weather.dat$RH<-str_replace_all(weather.dat$RH,c("58.50%"="58.5"))
weather.dat$RH<-str_replace_all(weather.dat$RH,c("88.20%"="88.2"))
weather.dat$RH<-str_replace_all(weather.dat$RH,c("43.20%"="43.2"))
weather.dat$RH<-str_replace_all(weather.dat$RH,c("n/a"="NA"))
weather.dat$RH<-str_replace_all(weather.dat$RH,c("na"="NA"))
weather.dat$RH<-str_replace_all(weather.dat$RH,c("6236"="62.36"))
##convert weather data RH to numeric
weather.dat$RH<-as.numeric(weather.dat$RH)
##convert weather data RH to percentage
weather.dat$RH <- sprintf("%.2f%%", weather.dat$RH)
##fix errors in weather data SkyCode
unique(weather.dat$SkyCode)
weather.dat$SkyCode<-str_replace_all(weather.dat$SkyCode,c("n/a"="NA"))
weather.dat$SkyCode<-str_replace_all(weather.dat$SkyCode,c("na"="NA"))
##fix errors in weather data WindCode
unique(weather.dat$WindCode)
weather.dat$WindCode<-str_replace_all(weather.dat$WindCode,c("n/a"="NA"))
weather.dat$WindCode<-str_replace_all(weather.dat$WindCode,c("N/A"="NA"))
weather.dat$WindCode<-str_replace_all(weather.dat$WindCode,c("na"="NA"))
##fix errors in weather data WindSpeed
unique(weather.dat$WindSpeed)
weather.dat$WindSpeed<-str_replace_all(weather.dat$WindSpeed,c("n/a"="NA"))
weather.dat$WindSpeed<-str_replace_all(weather.dat$WindSpeed,c("na"="NA"))
weather.dat$WindSpeed<-str_replace_all(weather.dat$WindSpeed,c("N/A"="NA"))
##convert weather data WindSpeed to numeric
weather.dat$WindSpeed<-as.numeric(weather.dat$WindSpeed)
##fix errors in weather data AmountofStandingWater
unique(weather.dat$AmountofStandingWater)
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("n/a"="NA"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("<1"=".5"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("<5"="4"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c(">6"="7"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("<10"="9"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("na"="NA"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("5 to 10"="7.5"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("~5"="5"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c(">5"="6"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("35-40"="37.5"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("60 to 70"="65"))
weather.dat$AmountofStandingWater<-str_replace_all(weather.dat$AmountofStandingWater,c("NA"="0.0"))
##convert weather data AmountofStandingWater to numeric
weather.dat$AmountofStandingWater<-as.numeric(weather.dat$AmountofStandingWater)
##fix errors in weather data WaterTemp
unique(weather.dat$WaterTemp)
weather.dat$WaterTemp<-str_replace_all(weather.dat$WaterTemp,c("n/a"="NA"))
weather.dat$WaterTemp<-str_replace_all(weather.dat$WaterTemp,c("N/A"="NA"))
weather.dat$WaterTemp<-str_replace_all(weather.dat$WaterTemp,c("na"="NA"))
##convert weather data WaterTemp to numeric
weather.dat$WaterTemp<-as.numeric(weather.dat$WaterTemp)
##fix errors in weather data DissolvedOxygen (no data)
unique(weather.dat$DissolvedOxygen)
weather.dat$DissolvedOxygen<-str_replace_all(weather.dat$DissolvedOxygen,c("n/a"="NA"))
weather.dat$DissolvedOxygen<-str_replace_all(weather.dat$DissolvedOxygen,c("N/A"="NA"))
weather.dat$DissolvedOxygen<-str_replace_all(weather.dat$DissolvedOxygen,c("na"="NA"))
##fix errors in weather data WaterpH (no data)
unique(weather.dat$WaterpH)
weather.dat$WaterpH<-str_replace_all(weather.dat$WaterpH,c("n/a"="NA"))
weather.dat$WaterpH<-str_replace_all(weather.dat$WaterpH,c("N/A"="NA"))
weather.dat$WaterpH<-str_replace_all(weather.dat$WaterpH,c("na"="NA"))
##fix errors in weather data MaxWaterDepth
unique(weather.dat$MaxWaterDepth)
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("n/a"="NA"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("50-70"="70"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("N/A"="NA"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("30 cm"="30"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("60 cm"="60"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("<5-10"="10"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("18-20"="20"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("na"="NA"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("~35"="35"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("~30"="30"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("~15"="15"))
weather.dat$MaxWaterDepth<-str_replace_all(weather.dat$MaxWaterDepth,c("~18"="18"))
##convert weater data MaxWaterDepth to numeric
weather.dat$MaxWaterDepth<-as.numeric(weather.dat$MaxWaterDepth)
##fix errors in weather data CanopyClosure (need to see how to fix, how to define range)
unique(weather.dat$CanopyClosure)
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("6-25%"="6-25"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("51-75%"="51-75"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("N/A"="NA"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("n/a"="NA"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("<5%"="<5"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c(">75%"=">75"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("26-50%"="26-50"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("na"="NA"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c("6-25, 51-75"="51-75"))
weather.dat$CanopyClosure<-str_replace_all(weather.dat$CanopyClosure,c(">75"="76-100"))
##fix errors in weather data LeafOut
unique(weather.dat$LeafOut)
weather.dat$LeafOut<-str_replace_all(weather.dat$LeafOut,c("y"="Y"))
weather.dat$LeafOut<-str_replace_all(weather.dat$LeafOut,c("n"="N"))
weather.dat$LeafOut<-str_replace_all(weather.dat$LeafOut,c("N/A"="NA"))
weather.dat$LeafOut<-str_replace_all(weather.dat$LeafOut,c("n/a"="NA"))
weather.dat$LeafOut<-str_replace_all(weather.dat$LeafOut,c("na"="NA"))
##fix errors in weather data EmergentVegCover (need to see how to fix)
unique(weather.dat$EmergentVegCover)
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("6-25%"="6-25"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("n/a"="NA"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("51-75%"="51-75"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("N/A"="NA"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("6-25%, >75"=">75"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("<5, 6-25%"="6-25"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("51-75, >75"=">75"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c(">75%"=">75"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("<5%"="<5"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("26-50%"="26-50"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("na"="NA"))
weather.dat$EmergentVegCover<-str_replace_all(weather.dat$EmergentVegCover,c("None"="0"))
##convert weather data EmergentVegCover to numeric
weather.dat$EmergentVegCover<-as.numeric(weather.dat$EmergentVegCover)
##creating waterdepth.weather.dat (water depth)
waterdepth.weather.dat<-weather.dat %>%
  group_by(WetlandID, Date, Site, SampleWindow)%>%
  summarize(MaxWaterDepth=mean(MaxWaterDepth))
##filter to only ARWMA
waterdepth.weather.dat <- waterdepth.weather.dat %>%
  filter(Site =="ARWMA")
##replacing NA with zero in waterdepth.weather.dat
waterdepth.weather.dat[is.na(waterdepth.weather.dat)] <- 0



##extract unique dates, wetland IDs, and species from amph.surv.dat data frame
dates <- unique(amph.surv.dat$Date)
wetland_ids <- unique(amph.surv.dat$WetlandID)
species <- unique(amph.surv.dat$Species)
##creating combinations
combinations <- expand.grid(Date = dates, WetlandID = wetland_ids, Species = species)
##convert to data frame
combinations_df <- as.data.frame(combinations)
##view the resulting data frame
print(combinations_df)
##create combinations data frame
combinations_df<-left_join(combinations_df,amph.surv.dat %>%
                             dplyr::select(Date, WetlandID, Species, Egg, Larvae, JuveAdult, Call, Notes, Total),
                           by=c("Date","WetlandID","Species"))
##replace NA with zero
combinations_df[is.na(combinations_df)] <- 0
##creating If/else statement
combinations_df$Detected<- ifelse(combinations_df$Total==0, "0", "1")
##filtered other combinations_df to only Pseudacris ornata
orn.dat<- combinations_df %>%
  filter(Species =="Pseudacris ornata")
##create hydroperiod from waterdepth.weather.dat
waterdepth.weather.dat$SampleWindow<-str_replace_all(waterdepth.weather.dat$SampleWindow,c("March "="March"))
waterdepth.weather.dat$SampleWindow<-str_replace_all(waterdepth.weather.dat$SampleWindow,c("July "="July"))
unique(waterdepth.weather.dat$SampleWindow)
unique(waterdepth.weather.dat$Date)
waterdepth.weather.dat$WetlandID<-str_replace_all(waterdepth.weather.dat$WetlandID,c("22/29"="29"))
hydroperiod_data <- waterdepth.weather.dat %>%
  group_by(WetlandID) %>%
  summarise(water.present = sum(MaxWaterDepth > 0, na.rm = TRUE),total.obs=n())
hydroperiod_data$prop.water<-hydroperiod_data$water.present/hydroperiod_data$total.obs
View(manage.dat)
View(hydroperiod_data)
View(orn.dat)
library(unmarked)
install.packages("unmarked")
library(unmarked)
##convert long to wide format orn.dat
model.dat<-orn.dat%>%pivot_wider(names_from=Date,values_from=Detected)
View(model.dat)
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
View(model.dat)
model.dat<-model.dat%>%pivot_wider(names_from=Date,values_from=Detected)
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date)
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
model.dat<-model.dat%>%pivot_wider(names_from=Date,values_from=Detected_Date)
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
model.dat<-model.dat%>%pivot_wider(names_from=Vist,values_from=Detected_Date)
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
model.dat<-model.dat%>%pivot_wider(names_from=Visit,values_from=Detected_Date)
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
model.dat<-model.dat[,c("WetlandID","Detected_Date","Visit")]
model.dat<-model.dat%>%pivot_wider(names_from=Visit,values_from=Detected_Date)
##adding covariates
model.dat<-left_join(model.dat,manage.dat,by=WetlandID)
##adding covariates
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
##adding covariates
manage.dat$WetlandID<-as.character(manage.dat$WetlandID)
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
model.dat<-left_join(model.dat,hydroperiod_data,by="WetlandID")
##adding covariates
manage.dat$WetlandID<-as.character(manage.dat$WetlandID)
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
model.dat<-left_join(model.dat,hydroperiod_data[,c"WetlandID","prop.water"],by="WetlandID")
##adding covariates
manage.dat$WetlandID<-as.character(manage.dat$WetlandID)
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
model.dat<-left_join(model.dat,hydroperiod_data[,c("WetlandID","prop.water")],by="WetlandID")
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
model.dat<-model.dat[,c("WetlandID","Detected_Date","Visit")]
model.dat<-model.dat%>%pivot_wider(names_from=Visit,values_from=Detected_Date)
##adding covariates
manage.dat$WetlandID<-as.character(manage.dat$WetlandID)
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
model.dat<-left_join(model.dat,hydroperiod_data[,c("WetlandID","prop.water")],by="WetlandID")
View(amph.surv.dat)
##checking for errors in amphibian survey data SampleWindow
unique(amph.surv.dat$SampleWindow)
##
amph.surv.dat%>%
  distinct(Year,SampleWindow)
##
amph.surv.dat%>%
  arrange(Date)%>%
  distinct(Year,SampleWindow)
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  arrange(WetlandID,Date)%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
##extract unique dates, wetland IDs, and species from amph.surv.dat data frame
dates <- unique(amph.surv.dat$Date)
wetland_ids <- unique(amph.surv.dat$WetlandID)
species <- unique(amph.surv.dat$Species)
##creating combinations
combinations <- expand.grid(Date = dates, WetlandID = wetland_ids, Species = species)
##convert to data frame
combinations_df <- as.data.frame(combinations)
##view the resulting data frame
print(combinations_df)
##create combinations data frame
combinations_df<-left_join(combinations_df,amph.surv.dat %>%
                             dplyr::select(Date, SampleWindow, WetlandID, Species, Egg, Larvae, JuveAdult, Call, Notes, Total),
                           by=c("Date","WetlandID","Species"))
##replace NA with zero
combinations_df[is.na(combinations_df)] <- 0
##creating If/else statement
combinations_df$Detected<- ifelse(combinations_df$Total==0, "0", "1")
##filtered other combinations_df to only Pseudacris ornata
orn.dat<- combinations_df %>%
  filter(Species =="Pseudacris ornata")
##create hydroperiod from waterdepth.weather.dat
waterdepth.weather.dat$SampleWindow<-str_replace_all(waterdepth.weather.dat$SampleWindow,c("March "="March"))
waterdepth.weather.dat$SampleWindow<-str_replace_all(waterdepth.weather.dat$SampleWindow,c("July "="July"))
unique(waterdepth.weather.dat$SampleWindow)
unique(waterdepth.weather.dat$Date)
waterdepth.weather.dat$WetlandID<-str_replace_all(waterdepth.weather.dat$WetlandID,c("22/29"="29"))
hydroperiod_data <- waterdepth.weather.dat %>%
  group_by(WetlandID) %>%
  summarise(water.present = sum(MaxWaterDepth > 0, na.rm = TRUE),total.obs=n())
hydroperiod_data$prop.water<-hydroperiod_data$water.present/hydroperiod_data$total.obs
##convert long to wide format orn.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat<-model.dat%>%
  arrange(WetlandID,Date)%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
##extract unique dates, wetland IDs, and species from amph.surv.dat data frame
dates <- unique(amph.surv.dat$Date)
wetland_ids <- unique(amph.surv.dat$WetlandID)
species <- unique(amph.surv.dat$Species)
samplewindow <- unique(amph.surv.dat$SampleWindow)
##creating combinations
combinations <- expand.grid(Date = dates, WetlandID = wetland_ids, Species = species, SampleWindow = samplewindow)
##convert to data frame
combinations_df <- as.data.frame(combinations)
##view the resulting data frame
print(combinations_df)
##create combinations data frame
combinations_df<-left_join(combinations_df,amph.surv.dat %>%
                             dplyr::select(Date, SampleWindow, WetlandID, Species, Egg, Larvae, JuveAdult, Call, Notes, Total),
                           by=c("Date","WetlandID","Species"))
##replace NA with zero
combinations_df[is.na(combinations_df)] <- 0
##creating If/else statement
combinations_df$Detected<- ifelse(combinations_df$Total==0, "0", "1")
##filtered other combinations_df to only Pseudacris ornata
orn.dat<- combinations_df %>%
  filter(Species =="Pseudacris ornata")
##extract unique dates, wetland IDs, and species from amph.surv.dat data frame
dates <- unique(amph.surv.dat$Date)
wetland_ids <- unique(amph.surv.dat$WetlandID)
species <- unique(amph.surv.dat$Species)
samplewindow <- unique(amph.surv.dat$SampleWindow)
##creating combinations
combinations <- expand.grid(Date = dates, WetlandID = wetland_ids, Species = species)
##convert to data frame
combinations_df <- as.data.frame(combinations)
##view the resulting data frame
print(combinations_df)
##create combinations data frame
combinations_df<-left_join(combinations_df,amph.surv.dat %>%
                             dplyr::select(Date, SampleWindow, WetlandID, Species, Egg, Larvae, JuveAdult, Call, Notes, Total),
                           by=c("Date","WetlandID","Species"))
##replace NA with zero
combinations_df[is.na(combinations_df)] <- 0
##creating If/else statement
combinations_df$Detected<- ifelse(combinations_df$Total==0, "0", "1")
##filtered other combinations_df to only Pseudacris ornata
orn.dat<- combinations_df %>%
  filter(Species =="Pseudacris ornata")
##extract unique dates, wetland IDs, and species from amph.surv.dat data frame
dates <- unique(amph.surv.dat$Date)
wetland_ids <- unique(amph.surv.dat$WetlandID)
species <- unique(amph.surv.dat$Species)
##creating combinations
combinations <- expand.grid(Date = dates, WetlandID = wetland_ids, Species = species)
##convert to data frame
combinations_df <- as.data.frame(combinations)
##view the resulting data frame
print(combinations_df)
##create combinations data frame
combinations_df<-left_join(combinations_df,amph.surv.dat %>%
                             dplyr::select(Date, WetlandID, Species, Egg, Larvae, JuveAdult, Call, Notes, Total),
                           by=c("Date","WetlandID","Species"))
##replace NA with zero
combinations_df[is.na(combinations_df)] <- 0
##creating If/else statement
combinations_df$Detected<- ifelse(combinations_df$Total==0, "0", "1")
##filtered other combinations_df to only Pseudacris ornata
orn.dat<- combinations_df %>%
  filter(Species =="Pseudacris ornata")
##create hydroperiod from waterdepth.weather.dat
waterdepth.weather.dat$SampleWindow<-str_replace_all(waterdepth.weather.dat$SampleWindow,c("March "="March"))
waterdepth.weather.dat$SampleWindow<-str_replace_all(waterdepth.weather.dat$SampleWindow,c("July "="July"))
unique(waterdepth.weather.dat$SampleWindow)
unique(waterdepth.weather.dat$Date)
waterdepth.weather.dat$WetlandID<-str_replace_all(waterdepth.weather.dat$WetlandID,c("22/29"="29"))
hydroperiod_data <- waterdepth.weather.dat %>%
  group_by(WetlandID) %>%
  summarise(water.present = sum(MaxWaterDepth > 0, na.rm = TRUE),total.obs=n())
hydroperiod_data$prop.water<-hydroperiod_data$water.present/hydroperiod_data$total.obs
##convert long to wide format from orn.dat to model.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat$month<-month(model.dat$Date)
install.packages("lubridate")
library(lubridate)
##convert long to wide format from orn.dat to model.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat$month<-month(model.dat$Date)
model.dat<-model.dat%>%
  arrange(WetlandID,Date)%>%
  group_by(WetlandID,Date)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Date,Detected_Date)
##convert long to wide format from orn.dat to model.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat$Month<-month(model.dat$Date)
model.dat$yYear<-year(model.dat$Date)
##convert long to wide format from orn.dat to model.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat$Month<-month(model.dat$Date)
model.dat$Year<-year(model.dat$Date)
model.dat<-model.dat%>%
  arrange(WetlandID,Year,Month)%>%
  group_by(WetlandID,Year,Month)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Detected_Date,Year,Month)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
model.dat<-model.dat[,c("WetlandID","Detected_Date","Visit")]
model.dat<-model.dat%>%pivot_wider(names_from=Visit,values_from=Detected_Date)
##adding covariates
manage.dat$WetlandID<-as.character(manage.dat$WetlandID)
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
model.dat<-left_join(model.dat,hydroperiod_data[,c("WetlandID","prop.water")],by="WetlandID")
##creating dynamic occupancy model
Occ.dat<-model.dat[,2:32]
View(Occ.dat)
##creating dynamic occupancy model
occ.dat<-model.dat[,2:32]
View(occ.dat)
##adding covariates
manage.dat$WetlandID<-as.character(manage.dat$WetlandID)
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
model.dat<-left_join(model.dat,hydroperiod_data[,c("WetlandID","prop.water")],by="WetlandID")
model.dat$New1<-NA
model.dat$New2<-NA
model.dat$New3<-NA
model.dat$New4<-NA
model.dat$New5<-NA
model.dat$New6<-NA
model.dat$New7<-NA
model.dat$New8<-NA
model.dat$New9<-NA
##convert long to wide format from orn.dat to model.dat
model.dat<-orn.dat[,c("Date","WetlandID","Detected")]
model.dat$Month<-month(model.dat$Date)
model.dat$Year<-year(model.dat$Date)
model.dat<-model.dat%>%
  arrange(WetlandID,Year,Month)%>%
  group_by(WetlandID,Year,Month)%>%
  mutate(Detected_Date=max(Detected))%>%
  distinct(WetlandID,Detected_Date,Year,Month)
model.dat<-model.dat%>%
  group_by(WetlandID)%>%
  mutate(Visit=row_number())
model.dat<-model.dat[,c("WetlandID","Detected_Date","Visit")]
model.dat<-model.dat%>%pivot_wider(names_from=Visit,values_from=Detected_Date)
##adding covariates
manage.dat$WetlandID<-as.character(manage.dat$WetlandID)
model.dat<-left_join(model.dat,manage.dat,by="WetlandID")
model.dat<-left_join(model.dat,hydroperiod_data[,c("WetlandID","prop.water")],by="WetlandID")
model.dat$New1<-NA
model.dat$New2<-NA
model.dat$New3<-NA
model.dat$New4<-NA
model.dat$New5<-NA
model.dat$New6<-NA
model.dat$New7<-NA
model.dat$New8<-NA
model.dat$New9<-NA
model.dat<-relocate(model.dat,New1,.after=6)
model.dat<-relocate(model.dat,New1,.after=7)
model.dat<-relocate(model.dat,New1,.after="6")
model.dat<-relocate(model.dat,New2,.after="6")
model.dat<-relocate(model.dat,c(New3,New4,New5,New6,New7,New8).after="8")
model.dat<-relocate(model.dat,c(New3,New4,New5,New6,New7,New8),.after="8")
model.dat<-relocate(model.dat,New9,.after="31")
summary(umfMS)
##creating dynamic occupancy model
umfMS <- unmarkedMultFrame(y=model.dat[,2:41], siteCovs=model.dat[,42:46], numPrimary=5)
summary(umfMS)
ModelMS <- colext(psiformula = ~1, # First-year occupancy
                  gammaformula = ~1, # Colonization
                  epsilonformula = ~1, # Extinction
                  pformula = ~1, # Detection
                  data = umfMS)
summary(ModelMS)
backTransform(ModelMS, type="psi")
backTransform(ModelMS, type="col")
backTransform(ModelMS, type="ext")
backTransform(nullModelMS, type="det")
backTransform(ModelMS, type="det")
