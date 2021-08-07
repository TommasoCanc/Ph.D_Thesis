#######################
# MOdels of Chapter 6 #
#######################

# Loading Libraries ------------------------------------------------------------
library(MASS) # <- Stat
library(ggplot2) # <- Plot
library(gridExtra) # <- Plot
library(psych) # <- Stat 
library(sjstats) # <- Stat
library(osfr) # <- Connect to OSF

#####################
# Connection to OSF #
#####################

main.path <- "/Users/tommasocancellario/Desktop/OSF_prove"
# Create folder where you can download the files. 
# For more information: ?dir.create
dir.create(main.path)
setwd(main.path)

# Connect to Open Science Framework
cr_project <- osf_retrieve_node("4p63t")
osf_files <- osf_ls_files(cr_project, type = "folder")

file.number <- which(osf_ls_files(osf_files[4, ])[ ,1] == "Variables_GeoDef.csv")
osf_download(osf_ls_files(osf_files[4, ])[file.number, ], path = main.path)
rm(file.number)

file.number <- which(osf_ls_files(osf_files[4, ])[ ,1] == "Variables_TaxDef.csv")
osf_download(osf_ls_files(osf_files[4, ])[file.number, ], path = main.path)
rm(file.number)

######################### 
# Geographic deficiency #
#########################

# Loading database -------------------------------------------------------------
db <- read.csv("./Variables_GeoDef.csv", sep = ";")

# Rename variables -------------------------------------------------------------
colnames(db) <- c("Country_code", "Country", "Area", "Population", "GDP", "x",
                      "y", "Elev_min", "Elev_max", "gap", "gap2", "richness",
                      "N_univ", "LIFE_fund")

# Creating useful new variables ------------------------------------------------
db$GDP           <- db$GDP/db$Population
db$N_univ        <- db$N_univ/(db$Population/100000) 
db$elev_range    <- db$Elev_max - db$Elev_min

# Data exploration -------------------------------------------------------------

####### Data exploration ##############
# A. Look at collinearity
# B. Look at outliers 
# C. Look at relationships
# D. Other stuffs

### Checking Collinearity ###
MyVar <- c(3 :9,12:14)
psych::pairs.panels(db[,MyVar])

# y 
# GDP
# Richness 
# N_univ
# LIFE
# Log_area

MyVar <- c(3,5,7,12:14,11)

### Checking outliers ###
colnames(db[,MyVar])
par(mfrow=c(3,3))
for(i in MyVar)
  dotchart(db[,i], main= as.character(colnames(db)[i]) ) ; rm(i)

# outlier in N_univ & Area & LIFE
db$Area   <- log(db$Area +1) # <- log-transform Area
par(mfrow=c(1,1)) ; dotchart(db$Area) # fine now

plot(db$N_univ)
#identify(db$N_univ) #<---- outlier at line 6
# manually removing the outlier
db <- db[-6,] 

### Cheking relationship between Y and X 
psych::pairs.panels(db[,c(10,3,5,7,12:14)])

### Checking Zero inflation ###
100 * sum(db$gap == 0, na.rm = TRUE) / nrow(db) #0% of zero

# Fitting the model -------------------------------------------------------

#Poisson Model
M0 <- glm(gap ~ Area + N_univ + y + GDP, data = db, family = "poisson")

performance::check_overdispersion(M0) # <- huge overdispersion

#Negative binomial model
M1 <- MASS::glm.nb(gap ~ scale(Area) + scale(N_univ) + scale(y) + scale(GDP), data = db) #scale to facilitate convergence

performance::check_model(M1) 

sink("./performance_scale.txt")
print("---------------------------VARIANCE---------------------------")
print(performance::r2(M1))  # <- Model explain 0.56 Variance
print("---------------------------VIF---------------------------")
car::vif(M1) # VIF score is <3. It's a good indicator. Low score of VIF indicates there is no collinearity in the variables.
print("---------------------------SUMMARY---------------------------")
summary(M1) 
print("---------------------------ANOVA---------------------------")
anova(M1)
sink()

# Visualizing significant relationships

# For simplicity, I re-do a model just with significant var
M_plot <- MASS::glm.nb(gap ~ N_univ + Area + y + GDP, data = db) #removing scale to facilitate plotting

# Generating data to predict
MyData_Univ <- data.frame(N_univ=seq(from = range(db$N_univ,na.rm=T)[1], to = range(db$N_univ,na.rm=T)[2], length = 200),
                          Area = mean(db$Area),
                          y = mean(db$y), 
                          GDP = mean(db$GDP))

MyData_Area <- data.frame(N_univ = mean(db$N_univ),
                          Area=seq(from = range(db$Area,na.rm=T)[1], to = range(db$Area,na.rm=T)[2], length = 200),
                          y = mean(db$y), 
                          GDP = mean(db$GDP))

# Prediction
P1 <- predict(M_plot, newdata = MyData_Univ, se = TRUE, type="response")
head(P1)
MyData_Univ$mu     <- P1$fit  # Fitted values
MyData_Univ$selow  <- P1$fit - 1.9 * P1$se.fit  # lower bound 
MyData_Univ$seup   <- P1$fit + 1.9 * P1$se.fit  # upper bound

P2 <- predict(M_plot, newdata = MyData_Area, se = TRUE,type="response")

MyData_Area$mu       <- P2$fit  # Fitted values
MyData_Area$selow    <- P2$fit - 1.9 * P2$se.fit  # lower bound
MyData_Area$seup     <- P2$fit + 1.9 * P2$se.fit  # upper bound

# Plot
plot1 <- ggplot() + 
    xlab("Number of university / population") +
    ylab("Gap")+
    geom_line(data = MyData_Univ, aes(x = N_univ, y = mu),size = 1,col = "blue") +
    geom_ribbon(data = MyData_Univ, aes(x = N_univ, ymax = seup,ymin = selow),
                fill = "lightblue",alpha = 0.4)+
    geom_point(data = db, aes(y = gap, x = N_univ), shape = 16, size = 2, 
               col = "grey20", alpha = 1) +
    theme_classic()

plot2 <- ggplot() + 
    xlab("Country area (log-transformed)") +
    ylab("Gap")+
    geom_line(data = MyData_Area, aes(x = Area, y = mu),size = 1,col = "purple") +
    geom_ribbon(data = MyData_Area, aes(x = Area, ymax = seup, ymin = selow), 
                fill = "plum", alpha = 0.4)+
    geom_point(data = db, aes(y = gap, x = Area), shape = 16, size = 2, 
               col = "grey20", alpha = 1) +
    theme_classic()

gridExtra::grid.arrange(plot1, plot2, ncol = 2, nrow = 1)

rm(list=setdiff(ls(), c("main.path")))


########################
# Taxonomic deficiency #
########################

# Load data --------------------------------------------------------------------
db.g <- read.csv("./Variables_TaxDef.csv", sep = ";")
head(db.g); str(db.g)

# Longitudinal
# Body size
# Endemic

# Select variables
MyVar <- c(4:6)

# Checking outliers ------------------------------------------------------------

# Deficiency
dotchart(db.g$gap)
db.g <- db.g[db.g$gap <1000,]
dotchart(db.g$gap) # <- Better

# Endemic
dotchart(db.g[,6]+1, main= as.character(colnames(db.g)[6]) )
# Log transformation
db.g$Number.of.Endemic  <- log(db.g$Number.of.Endemic +1)
dotchart(db.g[,6], main= as.character(colnames(db.g)[6]) ) # <- Better

# Body size
table(db.g$Body.size) # Unbalanced levels!
db.g$Body.size <- as.factor(db.g$Body.size)

levels(db.g$Body.size) <- c("0-2cm",">2cm",">2cm","0-2cm") # Let's try this
table(db.g$Body.size) # <- Better

# Longitudinal band
table(db.g$Longitudinal.distribution) # Unbalance

### Cheking relationship between Y and X ###
psych::pairs.panels(db.g[,c(2,6)])
boxplot(db.g$gap~db.g$Body.size) # maybe an effect
boxplot(db.g$gap~db.g$Longitudinal.distribution) # no effect
boxplot(db.g$Number.of.Endemic~db.g$Body.size) # no effect

### Checking Zero inflation ###
100 * sum(db.g$gap == 0, na.rm = TRUE) / nrow(db.g) # 8% of zero --- Very low

# Fitting the model ------------------------------------------------------------

#Poisson Model
M0 <- glm(gap ~ Body.size + Number.of.Endemic, data = db.g, family = "poisson")

performance::check_overdispersion(M0) # <- Overdispersion detected.

#Negative binomial model -------------------------------------------------------
M1 <- MASS::glm.nb(gap ~ as.factor(Body.size) + scale(Number.of.Endemic), data = db.g) 

# Checking
performance::r2(M1) 
car::vif(M1) 
performance::check_model(M1) 

summary(M1)

# Generating data to predict ---------------------------------------------------
MyData_endemic1 <- data.frame(
  Body.size = "0-2cm",
  Number.of.Endemic =seq(from = range(db.g$Number.of.Endemic,na.rm = TRUE)[1], 
                         to = range(db.g$Number.of.Endemic,na.rm = TRUE)[2], 
                         length = 200))

MyData_endemic2 <- data.frame(
  Body.size = ">2cm",
  Number.of.Endemic =seq(from = range(db.g$Number.of.Endemic,na.rm = TRUE)[1], 
                         to = range(db.g$Number.of.Endemic,na.rm = TRUE)[2], 
                         length = 200))

P1 <- predict(M1, newdata = MyData_endemic1, se = TRUE, type="response")

MyData_endemic1$mu     <- P1$fit  # Fitted values
MyData_endemic1$selow  <- P1$fit - 1.9 * P1$se.fit  # lower bound 
MyData_endemic1$seup   <- P1$fit + 1.9 * P1$se.fit  # upper bound

P2 <- predict(M1, newdata = MyData_endemic2, se = TRUE, type="response")

MyData_endemic2$mu     <- P2$fit  #F itted values
MyData_endemic2$selow  <- P2$fit - 1.9 * P2$se.fit  # lower bound
MyData_endemic2$seup   <- P2$fit + 1.9 * P2$se.fit  # upper bound

MyData_endemic3 <- rbind(MyData_endemic1, MyData_endemic2)

plot1 <- ggplot() + 
  xlab("Number of Endemic species") +
  ylab("gap")+
  facet_grid(~ Body.size)+
  geom_line(data = MyData_endemic3,aes(x = Number.of.Endemic, y = mu)) +
  geom_ribbon(data = MyData_endemic3, aes(x = Number.of.Endemic, ymax = seup, ymin = selow ),alpha = 0.2)+
  geom_point(data = db.g, aes(y = gap, x = Number.of.Endemic),shape = 16, size = 2, col="black") + 
  geom_point(data = db.g, aes(y = gap, x = Number.of.Endemic),shape = 16, size = 1, col="darkslategray2") + 
  theme_bw()

rm(list = ls())