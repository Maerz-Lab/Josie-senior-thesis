library(rmarkdown)

library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jagsUI)
library(rjags)
library(coda)

DipnetData <- read.csv("ARWMA DipnetSurveyData_Aug-24-2023.csv")
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


##Converting Egg, Larvae, and Adult to numeric values
DipnetData$Egg <- gsub("[^0-9.-]", "", DipnetData$Egg)
DipnetData$Larvae <- gsub("[^0-9.-]", "", DipnetData$Larvae)
DipnetData$JuveAdult <- gsub("[^0-9.-]", "", DipnetData$JuveAdult)

# Convert to numeric
DipnetData$Egg <- as.numeric(DipnetData$Egg)
DipnetData$Larvae <- as.numeric(DipnetData$Larvae)
DipnetData$JuveAdult <- as.numeric(DipnetData$JuveAdult)

# Check the result
str(DipnetData)
# Ensure every WetlandID-Year combination is present (even if no data exists)

##*********************************************************************************************************
  
  
  
  
# Create a full grid of WetlandID and Year combinations, and then join with presence/absence matrix
wetland_year_grid <- expand.grid(WetlandID = unique(DipnetData$WetlandID),
                                 Year = unique(DipnetData$Year))

# Check the unique values of WetlandID
unique(DipnetData$WetlandID)

#ChatGPT experiement

# Extract unique values for Year, WetlandID, and Species from the data
dates <- unique(DipnetData$Year)
wetland_ids <- unique(DipnetData$WetlandID)
species <- unique(DipnetData$Species)
sum(is.na(DipnetData$Species))
#Summarize by group and year
# Summarize by Year and WetlandID, and include the count of unique species per combination
# Summarize by Year and WetlandID, and include the count of unique species per combination
total_summary_by_wetland <- DipnetData %>%
  group_by(Year, WetlandID) %>%
  summarise(
    Egg = sum(Egg, na.rm = TRUE),           # Sum of Egg for each Year-Wetland combination
    Larvae = sum(Larvae, na.rm = TRUE),      # Sum of Larvae for each Year-Wetland combination
    JuveAdult = sum(JuveAdult, na.rm = TRUE), # Sum of JuveAdult for each Year-Wetland combination
    Total = sum(Total, na.rm = TRUE),         # Sum of Total for each Year-Wetland combination
    NumSpecies = n_distinct(Species),         # Count of unique species per Year-Wetland combination
    SpeciesList = paste(unique(Species), collapse = ", ")  # List of species for each Year-Wetland combination
  ) %>%
  # Add a new column for presence of Ambystoma tigrinum
  mutate(
    Ambystoma_tigrinum_Present = ifelse(grepl("Ambystoma tigrinum", SpeciesList), 1, 0)
  )

# Print the summarized table with the new column
print(total_summary_by_wetland)
##populating data frame with those ones and zeros
# Define the years and WetlandIDs
years <- 2019:2023
wetland_ids <- unique(DipnetData$WetlandID)


#dataframe without eggs, juveadult, larvae
total_summary_without <- total_summary_by_wetland %>%
  select(Year, WetlandID, NumSpecies, SpeciesList, Ambystoma_tigrinum_Present)

# Print the second data frame (without Egg, Larvae, and JuveAdult)
print(total_summary_without)


total_summary_without <- total_summary_without[-1,]

# Step 2: Pivot wider to flip the headers (Years) to columns and fill with the WetlandID as row identifiers
flipped_df <- total_summary_without %>%
  pivot_wider(names_from = Year, values_from = Ambystoma_tigrinum_Present)


print(flipped_df)
flipped_df <-flipped_df[,-c(2,3)]


# Merging rows by wetlandID and keeping the correct status for each year
merged_df <- flipped_df %>%
  group_by(WetlandID) %>%
  summarise(
    `2019` = ifelse(any(`2019` == "1", na.rm = TRUE), "1", 
                    ifelse(any(`2019` == "0", na.rm = TRUE), "0", NA)),
    `2020` = ifelse(any(`2020` == "1", na.rm = TRUE), "1", 
                    ifelse(any(`2020` == "0", na.rm = TRUE), "0", NA)),
    `2021` = ifelse(any(`2021` == "1", na.rm = TRUE), "1", 
                    ifelse(any(`2021` == "0", na.rm = TRUE), "0", NA)),
    `2022` = ifelse(any(`2022` == "1", na.rm = TRUE), "1", 
                    ifelse(any(`2022` == "0", na.rm = TRUE), "0", NA)),
    `2023` = ifelse(any(`2023` == "1", na.rm = TRUE), "1", 
                    ifelse(any(`2023` == "0", na.rm = TRUE), "0", NA))
  )

print(merged_df)

#SelectSpecies, dates, wetland ID
tig.2024.arwma <- ARWMA_2024 %>%
  select(Date, Wetland.ID, Species.1.observed, Species.2.observed, Species.3.observed, Species.4.observed, Species.5.observed, Species.6.observed, Species.7.observed, Species.8.observed, Species.9.observed, Species.10.observed, Species.11.observed, Species.12.observed, Species.13.observed, Species.14.observed, Species.15.observed)
##Replace NA with 0, ambystoma with 1
tig.2024.arwma[is.na(tig.2024.arwma)] <- 0
tig.2024.arwma <- tig.2024.arwma %>%
  mutate_at(vars(Species.1.observed, Species.2.observed, Species.3.observed, Species.4.observed, Species.5.observed, Species.6.observed, Species.7.observed, Species.8.observed, Species.9.observed, Species.10.observed, Species.11.observed, Species.12.observed, Species.13.observed, Species.14.observed, Species.15.observed), ~ ifelse(. == "Ambystoma tigrinum", 1, 0))
#Combine species columns
tig.2024.arwma <- tig.2024.arwma %>%
  mutate(Date = year(Date))  # Extracts the year
tig.2024.arwma <- tig.2024.arwma %>%
  mutate(Presence = ifelse(rowSums(select(., starts_with("Species"))) > 0, 1, 0))
tig.2024.arwma <- tig.2024.arwma %>% select(Date, Wetland.ID, Presence)
#merge wetland ID
tig.2024.arwma <- tig.2024.arwma %>%
  group_by(Wetland.ID) %>%
  summarise(
    Presence = as.integer(any(Presence == 1))  # 1 if any presence, 0 if no presence
  )

##merged data frame
colnames(tig.2024.arwma)[which(colnames(tig.2024.arwma) == "Wetland.ID")] <- "WetlandID"
colnames(tig.2024.arwma)[which(colnames(tig.2024.arwma) == "Presence")] <- "2024"
final.dat <- merge(merged_df, tig.2024.arwma, by = "WetlandID")

#plotattempt
final.dat <- final.dat %>%
  mutate_at(vars(`2019`, `2020`, `2021`, `2022`, `2023`), as.integer)

final.dat$distance <- c(314, 1519, 775, 259, 270, 478, 557, 83, 83, 643, 500, 121, 557, 489, 489, 121, 259)

#hydroperiod data calculation
hydroperiod.dat <- ARWMA_2024 %>%
  select(Date, Max.water.depth..cm., Wetland.ID)
hydroperiod.dat <- hydroperiod.dat %>%
  mutate(Date = month(Date))  # Extracts the year

colnames(hydroperiod.dat)[which(colnames(hydroperiod.dat) == "Wetland.ID")] <- "WetlandID"
hydroperiod.dat <- hydroperiod.dat %>%
  group_by(WetlandID) %>% 
  summarise(hydroperiod.24 = n_distinct(Date[Max.water.depth..cm. > 0])) # View the resulting dataframe print(hydroperiod_data)

#add hydroperiod to data frame
salamander <- merge(final.dat, hydroperiod.dat, by = "WetlandID")

#DipnetData Sample period to Months
month_names <- c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')
DipnetData$Month <- match(DipnetData$SampleWindow, month_names)

##subset data from other years
DipnetData.hydro <- DipnetData %>%
  select(Month, WetlandID, Year)

##summarize number of months
hydro.data.dip <- DipnetData.hydro %>%
  group_by(WetlandID, Year) %>%
  summarize(unique_months_count = n_distinct(Month), .groups = 'drop')
#data frame
hydro.data.dip <- hydro.data.dip %>%
  pivot_wider(
    names_from = Year,
    values_from = unique_months_count,
    names_prefix = "year_"
  )
hydro.data.dip <- hydro.data.dip[-1,]

salamander1 <- salamander %>%
  left_join(hydro.data.dip, by = "WetlandID")
salamander1 <- salamander1 %>%
  select(-year_)




##adding in management as variable

treatment.data <- read.csv("Wetland_Level_Data.csv")

new_row <- data.frame(WetlandID = 4, CanopyCover = "NA", CanopyThinning = 0, Fire = 0, DitchFilling = 0)
treatment.data2 <- rbind(treatment.data, new_row)

salamander2 <- merge(salamander1, treatment.data2, by = "WetlandID")


##Occupancy Model
nSites <- length(unique(salamander1$WetlandID))
nVisits <- 6


##Occupancy Matrix 

salamander.matrix <- salamander1[,c("2019","2020", "2021","2022","2023","2024")]


salamander.matrix <- as.matrix(salamander.matrix)
effort <- as.matrix(salamander.matrix)
effort[!is.na(effort)] <- -50
effort[is.na(effort)] <- 0
effort[effort == -50] <- 1
effort
#site covariates 
distance <- salamander1$distance

Fire <- as.numeric(salamander2$Fire)
DitchFilling <- as.numeric(salamander2$DitchFilling)
CanopyThinning <- as.numeric(salamander2$CanopyThinning)


#observation covariates
hydro.matrix <- as.matrix(salamander2)
hydro.matrix <- apply(hydro.matrix, 2, as.numeric)

hydro.mat <- c("year_2019", "year_2020", "year_2021", "year_2022", "year_2023", "hydroperiod.24")
hydro <- hydro.matrix[,hydro.mat]

hydro[hydro == ""] <-0
hydro[is.na(hydro)] <-0


#standardize 
# Standardize column-wise, handling NAs
hydro.scaled <- apply(hydro, 2, function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))  # Return NA if all values are NA
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
})


distance_scaled <- (distance-mean(distance))/sd(distance)





cat(file = "tiger.model.txt",
    "
    model {
    #priors
    #priors for occupancy/site coefficients
    beta0 ~ dnorm(0, 0.5)
    beta1 ~ dnorm(0, 0.5)
    beta2 ~ dlogis(0, 1)
    beta3 ~ dlogis(0, 1)
    beta4 ~ dlogis(0, 1)
    
    #detection/observation covariate
    alpha0 ~ dnorm(0, 0.5)
    alpha1 ~ dnorm(0, 0.5)
    
    for(i in 1:nSites) {
      logit(psi[i]) <- beta0 + beta1*dist[i]+ beta2*plug[i]+ beta3*thin[i]+ beta4*fire[i]
      z[i] ~ dbern(psi[i])
    
    for(j in 1:nOccasions) {
      logit(p[i, j]) <- alpha0 + alpha1*hydro[i,j]
      y[i,j] ~ dbern(z[i]*pf[i,j])
      pf[i,j] <- p[i,j]*effort[i,j]
      }
    }
    sitesOccupied <- sum(z)
   }
    "
)


#data into named list
jags.data <- list(y=salamander.matrix, hydro=hydro.scaled, dist=distance_scaled, nSites=nSites, nOccasions=nVisits, effort=effort, plug = DitchFilling, thin = CanopyThinning, fire = Fire)

#initial values
jags.inits <- function(){list(beta0=rnorm(1), alpha0=rnorm(1), z=rep(1, nSites))
}

#parameters to be monitored
jags.pars <- c("beta0", "beta1", "beta2", "beta3", "beta4", "alpha0", "alpha1", "sitesOccupied")

#MCMC

jags.post.samples <- jags.basic(data=jags.data, inits=jags.inits,
                                parameters.to.save = jags.pars,
                                model.file = "tiger.model.txt",
                                n.chains=3, n.adapt=100, n.burnin=500,
                                n.iter=7000, parallel=TRUE)
summary(jags.post.samples)


#### Test for convergence ####
quartz()
plot(jags.post.samples[,1:7])

# summarizing posteriors into a table#
(jagssum<- summary(jags.post.samples))



#means
Dist.Mean.JAGS <- jagssum$statistics["beta1", "Mean"]
Ditch.Mean.JAGS <- jagssum$statistics["beta2", "Mean"]
Canopy.Mean.JAGS <- jagssum$statistics["beta3", "Mean"]
Fire.Mean.JAGS <- jagssum$statistics["beta4", "Mean"]
Hydro.Mean.JAGS <- jagssum$statistics["alpha1", "Mean"]

#lower CI
LCI.dist.JAGS <- jagssum$quantiles["beta1", "2.5%" ]
LCI.ditch.JAGS <- jagssum$quantiles["beta2", "2.5%"]
LCI.canopy.JAGS <- jagssum$quantiles["beta3", "2.5%"]
LCI.fire.JAGS <- jagssum$quantiles["beta4", "2.5%"]
LCI.hydro.JAGS <- jagssum$quantiles["alpha1", "2.5%" ]

#upper CI
UCI.dist.JAGS <- jagssum$quantiles["beta1", "97.5%" ]
UCI.ditch.JAGS <- jagssum$quantiles["beta2", "97.5%"]
UCI.canopy.JAGS <- jagssum$quantiles["beta3", "97.5%"]
UCI.fire.JAGS <- jagssum$quantiles["beta4", "97.5%"]
UCI.hydro.JAGS <- jagssum$quantiles["alpha1", "97.5%" ]

#Bayesian p value


comparison_table <- data.frame(
  "Species" = "Ambystoma tigrinum",
  "Parameter" = c("Distance","Hydroperiod", "Ditch Filling", "Canopy Thinning", "Fire"),
  "Mean" = c(Dist.Mean.JAGS, Hydro.Mean.JAGS, Ditch.Mean.JAGS, Canopy.Mean.JAGS, Fire.Mean.JAGS), 
  "Lower CI" = c(LCI.dist.JAGS, LCI.hydro.JAGS, LCI.ditch.JAGS, LCI.canopy.JAGS, LCI.fire.JAGS),
  "Upper CI" = c(UCI.dist.JAGS, UCI.hydro.JAGS, UCI.ditch.JAGS, UCI.canopy.JAGS, UCI.fire.JAGS)
)

# compare
print(comparison_table)


#  Shows parameters on x-axis and mean estimate as a dot and the 95% credible interval as a line.


########## OCCUPANCY PROBABILITY
#### Plot predicted occupancy probability: Distance to last occupied ####
# Check if psi.post.pred is correctly filled
# Set up a sequence for prediction data
##is this right

pred.data <- data.frame(distance_scaled = seq(from = min(distance_scaled), to = max(distance_scaled), length = 100), hydro.scaled =0)
pred.data$distance <- pred.data$distance_scaled*(sd(distance))+(mean(distance))
pred.data$hydro <- rep(0, nrow(pred.data))  # holding hydro constant at zero if that's the intention
psi.coef.post <- as.matrix(jags.post.samples[,c("beta0", "beta1", "beta2", "beta3", "beta4", "alpha0", "alpha1")])
head(psi.coef.post)

# Create a matrix to store predictions across MCMC samples
n.iter <- nrow(psi.coef.post)  # number of iterations in posterior samples
psi.post.pred <- matrix(NA, nrow = n.iter, ncol = nrow(pred.data))

# Calculate predictions
for(i in 1:n.iter) {
  psi.post.pred[i, ] <- plogis(psi.coef.post[i, "beta0"] +
                                 psi.coef.post[i, "beta1"] * pred.data$distance_scaled
                                )
}

# Check the structure of psi.post.pred to ensure it's filled
head(psi.post.pred)

# Calculate mean and credible intervals for the predictions
pred.post.mean <- colMeans(psi.post.pred)
pred.post.lower <- apply(psi.post.pred, 2, quantile, prob = 0.025)
pred.post.upper <- apply(psi.post.pred, 2, quantile, prob = 0.975)

# Plot the predicted occupancy probability against distance
quartz()
plot(pred.data$distance, psi.post.pred[1,], type="l", xlab="Distance to nearest occupied wetland (meters)",  
     ylab="Occurrence probability", ylim=c(0, 1), col=gray(0.8)) #prediction line for first posterior samples
for(i in 1:n.iter) {
  lines(pred.data$distance, psi.post.pred[i,], col=gray(0.8))
} # posterior predictive distribution
lines(pred.data$distance, pred.post.mean, col="blue") #mean
lines(pred.data$distance, pred.post.lower, col="blue", lty=2) #lower CI
lines(pred.data$distance, pred.post.upper, col="blue", lty=2) #upper CI







####### DETECTION PROBABILITY
pred.data2 <- data.frame(hydro.scaled = seq(from = min(hydro.scaled), to =5, length = 100), distance_scaled =0)
pred.data2$hydro <- pred.data2$hydro.scaled*(sd(hydro))+(mean(hydro))
pred.data2$distance <- rep(0, nrow(pred.data2))  # holding hydro constant at zero if that's the intention


# Create a matrix to store predictions across MCMC samples
n.iter <- nrow(psi.coef.post)  # number of iterations in posterior samples
psi.post.pred2 <- matrix(NA, nrow = n.iter, ncol = nrow(pred.data2))

# Calculate predictions
for(i in 1:n.iter) {
  psi.post.pred2[i, ] <- plogis(psi.coef.post[i, "alpha0"] +
                                 psi.coef.post[i, "alpha1"] * pred.data2$hydro.scaled
  )
}

# Check the structure of psi.post.pred to ensure it's filled
head(psi.post.pred2)

# Calculate mean and credible intervals for the predictions
pred.post.mean <- colMeans(psi.post.pred2)
pred.post.lower <- apply(psi.post.pred2, 2, quantile, prob = 0.025)
pred.post.upper <- apply(psi.post.pred2, 2, quantile, prob = 0.975)

# Plot the predicted occupancy probability against distance
quartz()
plot(pred.data2$hydro, psi.post.pred2[1,], type="l", xlab="Hydroperiod (months)", xlim = c(0,12),  
     ylab="Detection probability", ylim=c(0, 1), col=gray(0.8)) #prediction line for first posterior samples
for(i in 1:n.iter) {
  lines(pred.data2$hydro, psi.post.pred2[i,], col = rgb(0,1,0, alpha = 0.01)) } # posterior predictive distribution
lines(pred.data2$hydro, pred.post.mean, col="black") #mean
lines(pred.data2$hydro, pred.post.lower, col="black", lty=2) #lower CI
lines(pred.data2$hydro, pred.post.upper, col="black", lty=2) #upper CI
                       # Use minimal theme


#potentially cool graph ???????
# Check the structure of jags.post.samples

samples_matrix <- as.matrix(jags.post.samples)
beta2_samples <- samples_matrix[, "beta2"]  # DitchFilling
beta3_samples <- samples_matrix[, "beta3"]  # Canopy Thinning
beta4_samples <- samples_matrix[, "beta4"]  # Fire
management_effects <- data.frame(
  Parameter = rep(c("Ditch Filling", "Canopy Thinning", "Fire"), each = length(beta2_samples)),
  Value = c(beta2_samples, beta3_samples, beta4_samples)
)


library(ggplot2)
library(ggridges)


# Plot the ridge density plot
quartz()
ridges <- ggplot(data = management_effects, 
       aes(y = factor(Parameter), x = Value)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Ditch Filling", "Canopy Thinning", "Fire")) + 
  ylab("") + 
  xlab("Estimated Impact on Occupancy Probability") +  # Adjust the label based on your parameter
  labs(title = "Posterior Distribution of Management Strategies Impact") +
  theme_minimal()
print(ridges)


ggplot(data = obsprob.female.jags, 
       aes(y = factor(x), x = pp)) + 
  stat_density_ridges(quantile_lines = TRUE, 
                      quantiles = c(0.025, 0.5, 0.975), vline_color = "white") + 
  scale_y_discrete(labels = c("Male", "Female")) + 
  ylab("") + 
  xlab("Estimated probability of volunteering") + 
  labs(title = "Probability based on observed-case approach") +
  theme_minimal()

library(ggplot2)

bayes.p <- as.matrix(jags.data)
T_obs <- sum(bayes.p$y, na.rm = TRUE)  # Sum of observed detections
# Initialize storage for simulated test statistics
T_rep <- numeric(length = length(beta0_samples))

for (i in 1:length(jags.pars)) {
  # Sample posterior parameters
  beta0_i <- beta0_samples[i]
  beta1_i <- beta1_samples[i]
  beta2_i <- beta2_samples[i]
  beta3_i <- beta3_samples[i]
  beta4_i <- beta4_samples[i]
  alpha0_i <- alpha0_samples[i]
  alpha1_i <- alpha1_samples[i]
  
  # Simulate occupancy probabilities
  psi_sim <- plogis(beta0_i + beta1_i * jags.data$dist +
                      beta2_i * jags.data$plug + 
                      beta3_i * jags.data$thin +
                      beta4_i * jags.data$fire)
  
  # Simulate detection probabilities
  p_sim <- plogis(alpha0_i + alpha1_i * jags.data$hydro)
  
  # Simulate occupancy and detection
  z_sim <- rbinom(jags.data$nSites, 1, psi_sim)
}               
bayesian_p_value <- mean(T_rep >= T_obs)
print(bayesian_p_value)



library(ggplot2)

# Convert to a data frame for ggplot
parameters_df <- as.data.frame(samples_matrix)
quartz()
# Example: Plot for beta1 (Distance)
parameterplot <- ggplot(parameters_df, aes(x = beta2)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = mean(parameters_df$beta1), color = "red", linetype = "dashed") +
  labs(title = "Posterior Distribution of Ditch Filling Parameter (beta2)",
       x = "Value",
       y = "Density") +
  theme_minimal()
print(parameterplot)
# Repeat for other parameters (beta2, beta3, etc.)

# Convert to a data frame for ggplot
parameters_df <- as.data.frame(samples_matrix)
quartz()
# Example: Plot for beta1 (Distance)
parameterplot2 <- ggplot(parameters_df, aes(x = beta3)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = mean(parameters_df$beta1), color = "red", linetype = "dashed") +
  labs(title = "Posterior Distribution of Canopy Thinning Parameter (beta3)",
       x = "Value",
       y = "Density") +
  theme_minimal()
print(parameterplot2)

# Convert to a data frame for ggplot
parameters_df <- as.data.frame(samples_matrix)
quartz()
# Example: Plot for beta1 (Distance)
parameterplot3 <- ggplot(parameters_df, aes(x = beta4)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = mean(parameters_df$beta1), color = "red", linetype = "dashed") +
  labs(title = "Posterior Distribution of Fire Parameter (beta4)",
       x = "Value",
       y = "Density") +
  theme_minimal()
print(parameterplot3)
# Repeat for other para
# Repeat for other para



library(tidyr)
parameters_long <- parameters_df %>%
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Value")

# Plot all parameter distributions in one graph
library(ggplot2)
quartz()
combined_plot <- ggplot(parameters_long, aes(x = Value, fill = Parameter)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Parameter, scales = "free") +
  labs(title = "Posterior Distributions of Parameters",
       x = "Value",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none") # Hide legend since we use faceting
print(combined_plot)



parameters_long <- parameters_df %>%
  pivot_longer(cols = everything(),
               names_to = "Parameter",
               values_to = "Value")

# Plot all parameter distributions on one graph
library(ggplot2)
quartz()
combined_plot <- ggplot(parameters_long, aes(x = Value, color = Parameter, fill = Parameter)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior Distributions of Parameters",
       x = "Value",
       y = "Density") +
  theme_minimal()
print(combined_plot)


# Convert samples_matrix to a data frame
parameters_df <- as.data.frame(samples_matrix)

# Reshape to long format and filter for specific parameters
library(tidyr)
parameters_long <- parameters_df %>%
  pivot_longer(cols = c(beta2, beta3, beta4),  # Include only beta2, beta3, and beta4
               names_to = "Parameter",
               values_to = "Value")

# Plot the selected parameter distributions on one graph
library(ggplot2)
quartz()
combined_plot <- ggplot(parameters_long, aes(x = Value, color = Parameter, fill = Parameter)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior Distributions of Selected Parameters",
       x = "Value",
       y = "Density") +
  theme_minimal()
print(combined_plot)

library(tidyr)
parameters_long2 <- parameters_df %>%
  pivot_longer(cols = c(alpha1, beta1, beta2, beta3, beta4),  # Include only beta2, beta3, and beta4
               names_to = "Parameter",
               values_to = "Value")

quartz()
combined_plot2 <- ggplot(parameters_long2, aes(x = Value, color = Parameter, fill = Parameter)) +
  geom_density(alpha = 0.3) +
  labs(title = "Posterior Distributions of Selected Parameters",
       x = "Value",
       y = "Density") +
  theme_minimal()
print(combined_plot2)



# Create the initial plot to set up the window
plot(distance_scaled, rep(NA, length(distance_scaled)), type = "n", 
     xlab = "Distance", ylab = "Predicted Occupancy",
     xlim = range(distance_scaled), ylim = c(0, 1))

# Loop through each management strategy
management_strategies <- c("fire", "ditch", "canopy")
colors <- c("red", "blue", "green")  # Color for each management strategy

for (i in 1:length(management_strategies)) {
  strategy <- management_strategies[i]
  pred.data[[strategy]] <- rep(1, nrow(pred.data))  # Set the management strategy to "active"
  
  # Recalculate predictions holding each strategy constant
  psi.coef.post <- as.matrix(jags.post.samples[, c("beta0", "beta1", "beta2", "beta3", "beta4", "alpha0", "alpha1")])
  psi.post.pred <- matrix(NA, nrow = n.iter, ncol = nrow(pred.data))
  
  for (j in 1:n.iter) {
    # Adjust for the current management strategy using the correct beta value
    if (strategy == "fire") {
      psi.post.pred[j, ] <- plogis(psi.coef.post[j, "beta0"] + 
                                     psi.coef.post[j, "beta1"] * pred.data$distance_scaled + 
                                     psi.coef.post[j, "beta4"] * pred.data$fire)  # Use beta4 for fire
    } else if (strategy == "ditch") {
      psi.post.pred[j, ] <- plogis(psi.coef.post[j, "beta0"] + 
                                     psi.coef.post[j, "beta1"] * pred.data$distance_scaled + 
                                     psi.coef.post[j, "beta2"] * pred.data$ditch)  # Use beta2 for ditch
    } else if (strategy == "canopy") {
      psi.post.pred[j, ] <- plogis(psi.coef.post[j, "beta0"] + 
                                     psi.coef.post[j, "beta1"] * pred.data$distance_scaled + 
                                     psi.coef.post[j, "beta3"] * pred.data$canopy)  # Use beta3 for canopy
    }
  }
  
  # Plot for each management strategy
  pred.post.mean <- colMeans(psi.post.pred)
  pred.post.lower <- apply(psi.post.pred, 2, quantile, prob = 0.025)
  pred.post.upper <- apply(psi.post.pred, 2, quantile, prob = 0.975)
  
  lines(pred.data$distance, pred.post.mean, col = colors[i], lwd = 2)
  lines(pred.data$distance, pred.post.lower, col = colors[i], lty = 2)
  lines(pred.data$distance, pred.post.upper, col = colors[i], lty = 2)
}

# Add a legend to describe each strategy
legend("topright", legend = management_strategies, col = colors, lty = 1:2, cex = 0.8)

is.na(salamander)
newsalamander <- salamander[, -c(1, 8, 9)]
newsalamander <- colSums(newsalamander, na.rm = TRUE) / nrow(newsalamander) * 100

percentage_occupied <- colSums(newsalamander) / nrow(newsalamander) * 100
years <- 2019:2024  # Adjust based on your actual year range
percentage_df <- data.frame(year = years, percentage_occupied = percentage_occupied)
quartz()
percentagegraph <- ggplot(percentage_df, aes(x = year, y = percentage_occupied)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 3) +
  labs(title = "Percentage of Wetlands Occupied Over Time",
       x = "Year",
       y = "Percentage of Wetlands Occupied (%)") +
  theme_minimal()
print(percentagegraph)
# View the percentage occupied for each year
percentage_occupied


hydro.sal <- hydro.data.dip[, -c(1,2)]






# Ensure that your matrix is correctly structured with years as columns
# Assuming each column represents a year and rows represent wetlands.

# Calculate the percentage of wetlands occupied for each year
percentage_occupied <- colSums(salamander.matrix != 0) / nrow(salamander.matrix) * 100

# Adjust the years range to match the number of columns in your matrix
years <- 2019:(2019 + ncol(salamander.matrix) - 1)

# Create a data frame for plotting
percentage_df <- data.frame(year = years, percentage_occupied = percentage_occupied)

# Plot the data
quartz()
percentagegraph <- ggplot(percentage_df, aes(x = year, y = percentage_occupied)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 3) +
  labs(title = "Percentage of Wetlands Occupied Over Time",
       x = "Year",
       y = "Percentage of Wetlands Occupied (%)") +
  theme_minimal()

print(percentagegraph)




quartz()
distancegraph <- ggplot(final.dat, aes(x = year, y = distance)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 3) +
  labs(title = "Percentage of Wetlands Occupied Over Time",
       x = "Year",
       y = "Percentage of Wetlands Occupied (%)") +
  theme_minimal()

print(distancegraph)
