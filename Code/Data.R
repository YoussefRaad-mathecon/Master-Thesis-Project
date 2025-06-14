
#------------------------------------- Packages and WD ---------------------------------------#

library(tidyverse) ### Data manipulations
library(dtplyr) ### Data manipulations - merge datasets
library(ggplot2) ### Plots
library(gridExtra) ### Plots
library("fHMM") ### Model fitting
library(Quandl) ### Data
library(dplyr) ### Data manipulations
library(lubridate) ### Time-variables
library(stats) ### ACF plots
library(matrixcalc) ### Matrix calculations
library("RColorBrewer") ### Colors
library(latex2exp) ### Text for plots
library(matrixStats) ### ColSds
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis")
load("yields.RData")

#------------------------------------- Data --------------------------------------------------#



### Create data

#codes <- c("FRED/DGS3MO",
           #"FRED/DGS6MO",
           #"FRED/DGS1",
           #"FRED/DGS2",
           #"FRED/DGS3",
           #"FRED/DGS5",
           #"FRED/DGS7",
           #"FRED/DGS10",
           #"FRED/DGS20",
           #"FRED/DGS30")

#yields <- Quandl(codes, collapse = "daily", type = "zoo", start_date = "1984-01-01", end_date = "2023-12-31")
#colnames(yields) <- c("3M", "6M", "1Y", "2Y", "3Y", "5Y", "7Y", "10Y", "20Y", "30Y")
#yields$Date <- rownames(as.data.frame(yields))
head(yields)



### Load data FRED (Federal Reserve Economic Data) treasury bill


yields$"3M" <- as.numeric(yields$"3M")/100
yields$"6M" <- as.numeric(yields$"6M")/100
yields$"1Y" <- as.numeric(yields$"1Y")/100
yields$"2Y" <- as.numeric(yields$"2Y")/100
yields$"3Y" <- as.numeric(yields$"3Y")/100
yields$"5Y" <- as.numeric(yields$"5Y")/100
yields$"7Y" <- as.numeric(yields$"7Y")/100
yields$"10Y" <- as.numeric(yields$"10Y")/100
yields$"20Y" <- as.numeric(yields$"20Y")/100
yields$"30Y" <- as.numeric(yields$"30Y")/100

yields_df <- as.data.frame(yields)
yields_df$DateCont <- year(yields_df$Date) + (month(yields_df$Date) - 1)/12 + (day(yields_df$Date) - 1)/365





yields_df <- yields_df[yields_df$`3M` != 0, ]
hist(as.numeric(yields_df$"3M"))

#------------------------------------- Timeseries plot -------------------------------------------------------------#

yields_df %>%
  ggplot(aes(x = DateCont, y = as.numeric(yields_df$"3M"))) +
  geom_line(color = "#901a1E") +
  theme_bw() +
  xlab("Year") +
  ylab("3M Rates") +
  theme(plot.title = element_text(size=17, hjust= 0.5)) +
  theme(axis.title = element_text(size=13)) +
  theme(axis.text.x = element_text(size=13, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size=13))
+
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11))

#------------------------------------- ACF plot -----------------------------------------------------------------------#

# Compute the autocorrelation
acf_data <- acf(as.numeric(yields_df$"3M"), lag.max = 100, plot = FALSE)

# Extract values
acf_values <- data.frame(
  Lag = acf_data$lag[, 1, 1],
  ACF = acf_data$acf[, 1, 1]
)

# Define confidence limits  
conf_limit <- qnorm((1 + 0.95) / 2) / sqrt(length(yields_df$"3M"))

ggplot(acf_values, aes(x = Lag, y = ACF)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_segment(aes(xend = Lag, yend = 0), color = "#901a1E") +
  geom_point(color = "#901a1E", size = 2) +
  geom_hline(yintercept = c(-conf_limit, conf_limit), color = "#901a1E", linetype = "dashed") +
  labs(
    x = "Lag",
    y = "ACF"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  )
