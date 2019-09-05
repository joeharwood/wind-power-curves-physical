library(RODBC)
require(lubridate)
require(ggplot2)
#require(hexbin)
require(tidyr)
library(splines)
library(nlme)
library(dplyr)
library("RColorBrewer")
library(grid)
library(gridExtra)

source("S:\\OandT\\OptRisk\\Energy_Requirements\\09 - Demand forecasting (Forecasting team)\\Renewable Generation Capacity\\R function\\running 32 bit query function.R")

# check for 32 bit set up
ifelse(.Machine$sizeof.pointer == 4, "Congrats you are running 32bit", "Caution - Its a trap 64bit")

#obtain GEN_IDs of current BMU's
#set up connection to DEAF
EFS <- odbcConnect("DEAFP", uid = "tde", pwd = "control1")

root_dir <- "S:/OandT/OptRisk/Energy_Requirements/19 - Wind Energy/WPFS/Wind farm models"
dev_dir <- "S:/OandT/OptRisk/Energy_Requirements/19 - Wind Energy/WPFS/Wind farm models/Physical models/development"

files_name_only <- list.files(path=paste0(root_dir, "/Historic_data"), pattern="*.RDS", recursive=FALSE)
files_name_only <- gsub("wind_id_", "", files_name_only)
files_name_only <- as.numeric(gsub(".RDS", "", files_name_only))

move_plots <- list.files(path=paste0(root_dir, "/Physical models/plots"), pattern="*.png", recursive=FALSE)

for (i in 1 : length(move_plots)){
  file.rename(from=paste0(root_dir, "/Physical models/plots/", move_plots[i]),
              to=paste0(root_dir, "/Physical models/plots/archive/", move_plots[i]))
}

#################################################################################################
##### Add any wind farms which need to have there curves update manually). These are excluded from update)
#################################################################################################

wind_curve_updated_manually <- c(15,91536)

######
## loop through list of wind farms

for(i in 1:length(files_name_only)){
  tryCatch({
## assign gen_id to i
Gen_id <- files_name_only[i]

## obtain historic data
metered_BOA_removed <- readRDS(paste0(root_dir, "/Historic_data/wind_id_", Gen_id, ".RDS"))

##assign capacity
capacity <- unique(metered_BOA_removed$capacity)

##create load factor
metered_BOA_removed <- metered_BOA_removed %>% mutate(LF = GEN_MW/capacity)

## Adding a Weighting column with higher values given to more recent dates
metered_BOA_removed <- metered_BOA_removed %>% 
  mutate(Weighting = (seq(1, nrow(metered_BOA_removed)))/nrow(metered_BOA_removed))

#set up connection args to DEAF
connection_args <-  list("DEAFP", "tde", "control1")

## Acquire existing physical model parameters
## create sqlquery
query_PM_parameters <- "select 
GTM.GEN_ID, 
wg.gen_name,
wg.gen_full_name,
GTM.TRBNE_ID, 
TT.SCALE, 
TT.MULT, 
TT.SENS, 
TT.PWROFF
from VTURBINE_TYPES TT
LEFT JOIN VGEN_TRBNE_MAP GTM on GTM.TRBNE_ID=TT.TRBNE_ID
LEFT JOIN WIND_GENERATOR wg on GTM.GEN_ID=wg.gen_id
where GTM.GEN_ID =Gen_id;"

query_PM_parameters <- gsub("Gen_id", Gen_id, query_PM_parameters)
query_PM_parameters <- gsub("\n", " ", query_PM_parameters)

#obtain existing parameters for a GEN_ID
existing_PM_parameters <- tbl_df(RODBC_query(sqlQuery(EFS, query_PM_parameters)))

s <- NA
s[1] <- existing_PM_parameters$SCALE
s[2] <- existing_PM_parameters$MULT
s[3] <- existing_PM_parameters$SENS
s[4] <- existing_PM_parameters$PWROFF

###########################################################################################################################
#############################               DATA CLEANING               ###################################################
###########################################################################################################################

#create functions for lines with data that is to be removed on one side and kept on the other
# the first two are combine in the filter in an OR statement so both coloured blue in the graph
# the production of equations and functions is explained in the document

right_removal_side = function(x) {(-.35*x + 9)}
right_removal_down = function(x) {(.6)}
left_removal = function(x) {(.18*x)}
negative_removal = function(x) {(0)}

cut_off_line = function(x) {(1*x)}

##  Plot data with selected removal lines
ggplot(data = metered_BOA_removed, aes(x = WS, y = LF)) + 
  geom_point(size = 0.1) +
  stat_function(fun = function(x) 1*s[1]/(1+s[2]*exp(-s[3]*x)), colour='red', size = 1.5)+
  stat_function(fun = negative_removal, colour='purple', size = 1.5)+
  stat_function(fun = right_removal_side, colour='blue', size = 1.5)+
  stat_function(fun = right_removal_down, colour='blue', size = 1.5)+
  stat_function(fun = left_removal, colour='red', size = 1.5) + coord_cartesian(ylim = c(-0.1, 1.1))

##  Use removal lines to clean data set
metered_BOA_removed_clean <- metered_BOA_removed %>% filter(LF < right_removal_side(WS) | LF > right_removal_down(WS)) 
metered_BOA_removed_clean <- metered_BOA_removed_clean %>% filter(LF < left_removal(WS))
metered_BOA_removed_clean <- metered_BOA_removed_clean %>% filter(LF > 0 | WS < 8)



#############################################
#############################################

##  Use graph to choose cut off value

cut_off <- 24

#############################################
#############################################

## perform linear regression
#mod1 <- nls(LF ~ 1*SCALE/(1 + MULT*exp(-SENS*WS)), 
#                              data = metered_BOA_removed, 
#                              start = list(SCALE = 0.9,MULT = 250, SENS = 0.7))

## perform linear regression with weightings
mod_w <- nls(LF ~ 1*SCALE/(1 + MULT*exp(-SENS*WS)), 
            data = metered_BOA_removed, 
            start = list(SCALE = 0.9,MULT = 250, SENS = 0.7), weights=Weighting^3)


# define esp equation
power_curve <- function(p, w) {s[1]/(1+s[2]*exp(-s[3]*w)) }

## use 2 models and defined equation(with existing parameters) to make dataframe
new_d1 <- data_frame(WS = seq(0, 28, by = 0.1))
new_d1 <- new_d1 %>%  mutate(nls_w = predict(mod_w, new_d1)) %>%
                      mutate(exist = power_curve(s, new_d1$WS))

#new_d1 <- new_d1 %>% mutate(nls = predict(mod1, new_d1)) %>% 
#  mutate(nls_w = predict(mod_w, new_d1)) %>%
#  mutate(exist = power_curve(s, new_d1$WS))

## remodel data for use with ggplot
long_new_d1 <- gather(new_d1, key = model, value = lf, nls_w, exist)

##  save results
gen_model <- data_frame(GEN_ID = Gen_id, SCALE = coefficients(mod_w)[["SCALE"]], 
                                          MULT = coefficients(mod_w)[["MULT"]], 
                                          SENS = coefficients(mod_w)[["SENS"]], 
                                          CUT_OFF = cut_off)

physical_models <- tbl_df(read.csv(paste0(root_dir, "/Physical models/physical_model.csv")))
physical_models <- physical_models %>% filter(GEN_ID != Gen_id)
physical_models <- physical_models %>% filter(!GEN_ID %in% wind_curve_updated_manually)
physical_models <- rbind(physical_models, gen_model) %>% arrange(GEN_ID)
write.csv(physical_models, paste0(root_dir, "/Physical models/physical_model.csv"), row.names = FALSE)

n <- NA
n[1] <- filter(physical_models, GEN_ID == Gen_id) %>% pull(SCALE)
n[2] <- filter(physical_models, GEN_ID == Gen_id) %>% pull(MULT)
n[3] <- filter(physical_models, GEN_ID == Gen_id) %>% pull(SENS)



##  Plot cleaned data set to check removal correctly performed
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1))

g <-   ggplot(metered_BOA_removed_clean, aes(x = WS, y = LF)) + 
  geom_point(size = 2, aes(colour=Weighting)) + sc

g <-   g + stat_function(fun = function(x) 1*n[1]/(1+n[2]*exp(-n[3]*x)), size = 1.5, aes(linetype="new")) +  
  stat_function(fun = function(x) 1*s[1]/(1+s[2]*exp(-s[3]*x)), size = 1.5, aes(linetype="existing")) 


h <-   ggplot(metered_BOA_removed, aes(x=GEN_MW))+ geom_histogram(binwidth=1, colour="white") +
geom_vline(aes(xintercept=capacity), color="red", linetype="dashed", size=1) +
  ggtitle(paste0("MW Values, Red line is Capacity"))

p6 <- grid.arrange(g, h, top = paste0("Gen_id = ", Gen_id, ",  Gen_Name = ", Gen_name, ",  Gen_Name = ", gen_full_name ,existing_PM_parameters$GEN_NAME), widths = c(2, 1))



## save plot

png(filename=paste0(root_dir, "/Physical models/plots/GEN_ID_",Gen_id,".png"),
    width = 1000, height = 750, units = "px")
grid.arrange(g, h, top = paste0("Gen_id = ", Gen_id, ",  Gen_Name = ", existing_PM_parameters$GEN_NAME), widths = c(2, 1))

dev.off()


  }, error = function(e){write(paste0(dev_dir, "  ", Sys.Date(), "/Gen_ID ", i , " not extracted ", conditionMessage(e), "\n"),file = "log", append = TRUE)}
)
}


