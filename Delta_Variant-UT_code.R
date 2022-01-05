
library(ggplot2)
library(dplyr)
library(lubridate)
library(cowplot)
library(zoo)
library(boot)


# Load Data ---------------------------------------------------------------
setwd("~/Downloads/Utah_COVID19_data_Oct15/")

## Read in the case data + vaccination data
vaccine_dat<-read.csv("Risk_Daily Case Counts by Vaccination Status_2021-10-15.csv")
vaccine_pop<-read.csv("Risk_Weekly Case Rates in Vaccinated and Unvaccinated People2021-10-15.csv")




## Deal with dates 
#vaccine_dat$Case.Date <- lubridate::ymd(vaccine_dat$Case.Date)
vaccine_pop$Case.Date <- lubridate::ymd(vaccine_pop$Week.Date..Sunday.)


## Read in the lineage data
lineage_data <- read.csv("Trends_Sequencing Results by Week_2021-10-15.csv")
lineage_data<-lineage_data[!(lineage_data$Percentage.of.all.Sequenced.Results.per.Week=="Sample size too small"),]
lineage_data$Percentage.of.all.Sequenced.Results.per.Week<-as.numeric(as.character(
  lineage_data$Percentage.of.all.Sequenced.Results.per.Week))
## Deal with dates
lineage_data$Week.Sample.was.Collected <- lubridate::ymd(lineage_data$Week.Sample.was.Collected)


## Filter the vaccination data to match the dates of the case data-- 
#### vaccine_pop starts in 2020, vaccine_dat (on breakthrough infections) starts 2021
#vaccine_pop <- vaccine_pop %>% filter(Case.Date >= min(vaccine_dat$Case.Date)) %>% filter(Case.Date <= max(vaccine_dat$Case.Date))

## Merge the two datasets which are now on the same time frame
vaccine_data <- vaccine_pop %>% select("Week.Date..Sunday.",  "Case.Date", "Breakthrough.Cases", "Unvaccinated.Cases", 
                                       "Unvaccinated.Population", "Fully.vaccinated.Population") %>% filter(Fully.vaccinated.Population >= 1)

vaccine_data$Total.Cases <-vaccine_data$Unvaccinated.Cases + vaccine_data$Breakthrough.Cases
  #cbind(vaccine_dat, vaccine_pop$Unvaccinated.Population, vaccine_pop$X14.Days.Post.Full.Vaccination.Population)

## Rename the colums of the data being merged in
# colnames(vaccine_data) <- c("Case.Date", "Total.Cases", "Breakthroughs..Vaccinated.Cases.", 
#                             "Unvaccinated.Cases", "Unvaccinated.Population", "X14.Days.Post.Full.Vaccination.Population")



# Calculate VE from attack rate -------------------------------------------

## VE = ARU - ARV / ARV *100

ARU <- vaccine_data$Unvaccinated.Cases/vaccine_data$Unvaccinated.Population
ARV <- vaccine_data$Breakthrough.Cases/vaccine_data$Fully.vaccinated.Population

vaccine_data$VE <- ((ARU-ARV)/ARU)*100


## Calculate the 14-day rolling average for VE 
ARU_avg<- ARU #rollmean(vaccine_data$Unvaccinated.Cases, 14, na.pad=TRUE)/rollmean(vaccine_data$Unvaccinated.Population,  7, na.pad=TRUE)
ARV_avg<- ARV #rollmean(vaccine_data$Breakthroughs..Vaccinated.Cases., 14, na.pad=TRUE)/rollmean(vaccine_data$X14.Days.Post.Full.Vaccination.Population, 7, na.pad=TRUE)

vaccine_data$VE_avg <- ((ARU_avg-ARV_avg)/ARU_avg)*100


# Estimating VE_delta -----------------------------------------------------

## VE against Delta Varaint

## Split the dataset into pre and post delta
Vaccine_notdelta <- vaccine_data %>% filter(Case.Date >= "2021-04-04")  %>% filter(Case.Date <= "2021-05-02")
Vaccine_delta <- vaccine_data %>% filter(Case.Date > "2021-05-23")  

## Calculate the mean VE during the not delta 
mean_notdelta <- mean(Vaccine_notdelta$VE)

## Define the weeks that are in the delta period to make an index to loop over
delta_weeks<-distinct(lineage_data %>% filter(Week.Sample.was.Collected >= "2021-05-23"), Week.Sample.was.Collected)
# week <- seq(1, dim(delta_weeks)[1])

## Loop to estimate VE of delta and take the mean over that period
VE_d <- NULL

for(i in 1:dim(delta_weeks)[1]){
  pct_delta <- as.numeric(lineage_data %>% filter(Lineage == "B.1.617.2", Week.Sample.was.Collected ==delta_weeks[1,]) %>% 
                            select(Percentage.of.all.Sequenced.Results.per.Week))/100
  VE_t <- Vaccine_delta%>% filter(Case.Date >= delta_weeks[i,]) %>% filter(Case.Date < delta_weeks[i,]+7) %>% select(VE)
  VE_est <- (VE_t - (mean_notdelta*(1-pct_delta)))/pct_delta
  VE_d <- rbind(VE_d, VE_est)
}
mean(VE_d$VE)





# BOOTSTRAP FOR CIs -------------------------------------------------------

N_draws = 1000


n_iter = 1000

vacc_data_all <- NULL
ve_d <- ve_d_sens <- NULL
for (j in 1:n_iter){
  
  vaccine_data_boot <- vaccine_data %>% mutate(pop = (Fully.vaccinated.Population + Unvaccinated.Population),
                                               p_vacc_case = Breakthrough.Cases / Total.Cases)
  vaccine_data_boot <- vaccine_data_boot %>% mutate(
    Breakthrough.Cases = round(rbinom(rep(1, length(p_vacc_case)), Total.Cases, prob=p_vacc_case)),
    Unvaccinated.Cases = Total.Cases - Breakthrough.Cases,
    break_throughs_sens = round(rbinom(rep(1, length(p_vacc_case)), Total.Cases, prob=(p_vacc_case*.9))),
    unvacc_case_sens = Total.Cases - break_throughs_sens)
  
  
  ## Calculate VE from unvaccinated and vaccinated attack rate
  
  ARU <- vaccine_data_boot$Unvaccinated.Cases/vaccine_data_boot$Unvaccinated.Population
  ARV <- vaccine_data_boot$Breakthrough.Cases/vaccine_data_boot$Fully.vaccinated.Population
  vaccine_data_boot$VE <- ((ARU-ARV)/ARU)*100
  
  ARU_ <- vaccine_data_boot$unvacc_case_sens / vaccine_data_boot$Unvaccinated.Population
  ARV_ <- vaccine_data_boot$break_throughs_sens / vaccine_data_boot$Fully.vaccinated.Population
  vaccine_data_boot$VE_sens <- ((ARU_-ARV_)/ARU_)*100
  
  
  # Calculate the 14-day rolling average
  ARU_avg<- rollmean(vaccine_data_boot$Unvaccinated.Cases, 14, na.pad=TRUE)/rollmean(vaccine_data_boot$Unvaccinated.Population,
                                                                                7, na.pad=TRUE)
  ARv_avg<- rollmean(vaccine_data_boot$Breakthrough.Cases, 14,
                     na.pad=TRUE)/rollmean(vaccine_data_boot$Fully.vaccinated.Population, 7, na.pad=TRUE)
  vaccine_data_boot$VE_avg <- ((ARU_avg-ARv_avg)/ARU_avg)*100


  ARU_avg_ <- rollmean(vaccine_data_boot$unvacc_case_sens, 14, na.pad=TRUE)/rollmean(vaccine_data_boot$Unvaccinated.Population,
                                                                                7, na.pad=TRUE)
  ARv_avg_ <- rollmean(vaccine_data_boot$break_throughs_sens, 14,
                     na.pad=TRUE)/rollmean(vaccine_data_boot$Fully.vaccinated.Population, 7, na.pad=TRUE)
  vaccine_data_boot$VE_avg_sens <- ((ARU_avg_-ARv_avg)/ARU_avg_)*100
  
  
  
  ## VE against Delta Varaint
  
  # Estimating VE_delta -----------------------------------------------------
  ## VE against Delta Varaint
  ## Split the dataset into pre and post delta
  Vaccine_notdelta <- vaccine_data_boot %>% filter(Case.Date >= "2021-04-04")  %>% filter(Case.Date <= "2021-05-02")
  Vaccine_delta <- vaccine_data_boot %>% filter(Case.Date > "2021-05-23")  
  
  ## Calculate the mean VE during the not delta 
  mean_notdelta <- mean(Vaccine_notdelta$VE)
  
  ## Define the weeks that are in the delta period to make an index to loop over
  delta_weeks<-distinct(lineage_data %>% filter(Week.Sample.was.Collected >= "2021-05-23"), Week.Sample.was.Collected)

  ## Loop to estimate VE of delta and take the mean over that period
  VE_d <- VE_d_sens <- NULL
  for(i in 1:dim(delta_weeks)[1]){
    pct_delta <- as.numeric(lineage_data %>% filter(Lineage == "B.1.617.2", Week.Sample.was.Collected ==delta_weeks[i,]) %>% 
                              select(Percentage.of.all.Sequenced.Results.per.Week))/100
    VE_t <- Vaccine_delta%>% filter(Case.Date >= delta_weeks[i,]) %>% filter(Case.Date < delta_weeks[i,]+7) %>% select(VE)
    VE_est <- (VE_t - (mean_notdelta*(1-pct_delta)))/pct_delta
    VE_d <- rbind(VE_d, VE_est)
    
    VE_t_sens <- Vaccine_delta%>% filter(Case.Date >= delta_weeks[i,]) %>% filter(Case.Date < delta_weeks[i,]+7) %>% select(VE_sens)
    VE_est_sens <- (VE_t_sens - (mean_notdelta*(1-pct_delta)))/pct_delta
    VE_d_sens <- rbind(VE_d_sens, VE_est_sens)
  }
  
  
  # Put it into one boot form
  
  vacc_data_all <- vacc_data_all %>% 
    bind_rows(vaccine_data_boot %>% mutate(boot = j) %>% select(Case.Date, boot, VE, VE_avg, VE_sens, VE_avg_sens))
  
  ve_d <- c(ve_d, mean(VE_d$VE))
  ve_d_sens <- c(ve_d_sens, mean(VE_d_sens$VE))
  
}


vacc_data_all_cis <- vacc_data_all %>% 
  group_by(Case.Date) %>%
  summarise(VE_mean = mean(VE, na.rm=TRUE),
            VE_low = quantile(VE, probs=0.025, na.rm=TRUE),
            VE_high = quantile(VE, probs=0.975, na.rm=TRUE),
            VE_avg_mean = mean(VE_avg, na.rm=TRUE),
            VE_avg_low = quantile(VE_avg, probs=0.025, na.rm=TRUE),
            VE_avg_high = quantile(VE_avg, probs=0.975, na.rm=TRUE),
            VE_sens_mean = mean(VE_sens, na.rm=TRUE),
            VE_sens_low = quantile(VE_sens, probs=0.025, na.rm=TRUE),
            VE_sens_high = quantile(VE_sens, probs=0.975, na.rm=TRUE),
            VE_avg_sens_mean = mean(VE_avg_sens, na.rm=TRUE),
            VE_avg_sens_low = quantile(VE_avg_sens, probs=0.025, na.rm=TRUE),
            VE_avg_sens_high = quantile(VE_avg_sens, probs=0.975, na.rm=TRUE))

mean(ve_d)
quantile(ve_d, probs=c(0.025, 0.975))

mean(ve_d_sens)
quantile(ve_d_sens, probs=c(0.025, 0.975))

# Plots -------------------------------------------------------------------

## PLOT 1: 14 day rolling average of VE

p1<- ggplot(vacc_data_all_cis) + 
  geom_ribbon(aes(x =Case.Date, ymin=VE_low, ymax=VE_high), fill="black", alpha=0.25) +
  geom_line(aes(x=Case.Date, y = VE_mean)) +
  #geom_line(aes(x=Case.Date, y = VE_avg_low), alpha = 0.7, linetype = "dashed") + 
  #geom_line(aes(x=Case.Date, y = VE_avg_high), alpha = 0.7, linetype = "dashed") + 
  ylab("Vaccine Effectiveness") + 
  xlab("Date") + 
  theme_bw() + 
  theme(legend.position = c(0.25, 0.125)) + 
  xlim(as_date("2021-01-21"), as_date("2021-10-09"))
  



## PLOT 2: % of each lineage identified each week 
p2<- ggplot(lineage_data %>% filter(Week.Sample.was.Collected >= as.Date("2021-01-17")), 
            aes(fill=Lineage, y=Percentage.of.all.Sequenced.Results.per.Week, x=Week.Sample.was.Collected)) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("Week Sample Collected") + 
  ylab("% Identified ") + 
  theme_bw() + 
  scale_fill_discrete( 
    labels = c("B.1.1.7" = "Alpha", "B.1.351" = "Beta", "B.1.427" = "Epsilon1",
               "B.1.429" = "Epsilon2", "B.1.617.2" = "Delta", 
               "Other Lineage" = "Other Lineage", "P.1" = "Gamma")) +
  theme(legend.position = "bottom") + 
  xlim(as_date("2021-01-21"), as_date("2021-10-09"))



## PLOT 3: % of casses vaccinated by % of population vully vaccinated, for 3 VEs
A<- ggplot() +  
  theme_bw() + 
  geom_line(data = Expected, aes(x = PPV*100, y = PCV_95*100), color = "black", alpha=0.65, linetype = "dashed") + 
  geom_line(data = Expected, aes(x = PPV*100, y = PCV_90*100), color = "black",alpha=0.65, linetype = "dashed") + 
  geom_line(data = Expected, aes(x = PPV*100, y = PCV_85*100), color = "black",alpha=0.65, linetype = "dashed") + 
  geom_line(data = Expected, aes(x = PPV*100, y = PCV_80*100), color = "black", alpha=0.65,linetype = "dashed") + 
  annotate("text", x = 43, y = 12, label = "VE = 80%", angle = 42,alpha=0.65, color = "black") + 
  annotate("text", x = 46, y = 10.3, label = "VE = 85%", angle = 39,alpha=0.65, color = "black") + 
  annotate("text", x = 54, y = 4.9, label = "VE = 95%", angle = 22, alpha=0.65,color = "black") + 
  annotate("text", x = 50, y = 8, label = "VE = 90%", angle = 30, alpha=0.65,color = "black") + 
  geom_line(data = vaccine_data, aes(
    x=(Fully.vaccinated.Population/(Unvaccinated.Population + Fully.vaccinated.Population)*100),
    y=rollmean(Frac_Breakthrough,1, na.pad=TRUE)), color = "Navy", size = 1.1) + 
  xlim(0, 70) + 
  ylim(0,30) + 
  xlab("% Population Fully Vaccinated") + 
  ylab("% of Cases Vaccinated") + 
  annotate("segment", x = 25.78836, xend = 25.78836, y = 0.3, yend = 2.3, size=0.8,  arrow=arrow(length = unit(3, "mm"), type = "closed")) + 
  annotate("text", x = 25.78836, y = 0.05, label = "Delta = 4.9%") 


## Plot of the RHS of figure 
BC<-plot_grid(p1, p2, nrow=2, labels = c('A', 'B'), rel_heights = c(0.7,1))


## Make the whole figure
plot_grid(A, BC, ncol=2, labels = c('A', ' '))

