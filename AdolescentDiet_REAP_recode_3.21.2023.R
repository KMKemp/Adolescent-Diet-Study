

###########################################################
#################### Recode Reap Data #####################
###########################################################
library("readxl")
library("tidyverse")

############# Set working directory #############
setwd("~/boxdrive/_Mrug/Diet")

############# Import metadata #############
reap<-read_excel("~/boxdrive/_Mrug/Diet/REAP/_reap_before_recode/W1reap.xlsx")

#identify all character columns
chars <- sapply(reap, is.character)

#convert all character columns to numeric
reap[ , chars] <- as.data.frame(apply(reap[ , chars], 2, as.numeric))

#check
str(reap)

############## Current Codes ##############
# usually/often 0
# sometimes 1   
# rarely/never 2 

############################################
################## Recode ##################
############################################


### c1reap1: In an average week, how often do you skip breakfast?
reap <- mutate(reap,
                c1reap1 = case_when(
                  c1reap1 == 0 ~ 3, 
                  c1reap1 == 1 ~ 2,
                  c1reap1 == 2 ~ 1, 
                    TRUE ~ NA_real_ # This is for all other values 
                  ))                # not covered by the above.

### c1reap2: In an average week, how often do eat 4 or more meals from a 
### sit-down or take out restaurant?
reap <- mutate(reap,
                c1reap2 = case_when(
                  c1reap2 == 0 ~ 3, 
                  c1reap2 == 1 ~ 2, 
                  c1reap2 == 2 ~ 1,
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap3: In an average week, how often do you eat less than 2 servings of 
### whole grain products or high fiber starches a day?
reap <- mutate(reap,
                c1reap3 = case_when(
                  c1reap3 == 0 ~ 3, 
                  c1reap3 == 1 ~ 2, 
                  c1reap3 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap4: In an average week, how often do you eat less than 2 servings of fruit a day? 
### fruit a day? 
reap <- mutate(reap,
                c1reap4 = case_when(
                  c1reap4 == 0 ~ 3, 
                  c1reap4 == 1 ~ 2, 
                  c1reap4 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.


### c1reap4a: In an average week, how often do you eat less than 2 servings of vegetables a day?
### This has already been recoded by Catherine and Anna- see the excel file


### c1reap5: In an average week, how often do you eat less than 2 servings of milk, yogurt, 
### or cheese a day?
reap <- mutate(reap,
                c1reap5 = case_when(
                  c1reap5 == 0 ~ 3, 
                  c1reap5 == 1 ~ 2, 
                  c1reap5 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap6: In an average week, how often do you eat more than 8 ounces (see sizes below) of
### meat, chicken, turkey, or fish per day?
reap <- mutate(reap,
                c1reap6 = case_when(
                  c1reap6 == 0 ~ 3, 
                  c1reap6 == 1 ~ 2, 
                  c1reap6 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap7: In an average week, how often do you use regular processed meats (like bologna,
### salami, corned beef, hotdogs, sausage, or bacon) instead of low fat processed meats (like roast
### beef, turkey, lean ham, low-fat cold cuts/hotdogs)?
reap <- mutate(reap,
                c1reap7 = case_when(
                  c1reap7 == 0 ~ 3, 
                  c1reap7 == 1 ~ 2, 
                  c1reap7 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap8: In an average week, how often do you eat fried foods such as fried chicken, fried fish,
### French fries, fried plantains, tostones or fried yucca?
reap <- mutate(reap,
                c1reap8 = case_when(
                  c1reap8 == 0 ~ 3, 
                  c1reap8 == 1 ~ 2, 
                  c1reap8 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap9: In an average week, how often do you eat eat regular potato chips, nacho chips, corn
### chips, crackers, regular popcorn, nuts, instead of pretzels, low-fat chips, or low-fat crackers, 
### air popped popcorn?
reap <- mutate(reap,
                c1reap9 = case_when(
                  c1reap9 == 0 ~ 3, 
                  c1reap9 == 1 ~ 2, 
                  c1reap9 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap10: In an average week, how often do you add butter, margarine, or oil to bread, potatoes,
### rice, or vegetables at the table?
reap <- mutate(reap,
                c1reap10 = case_when(
                  c1reap10 == 0 ~ 3, 
                  c1reap10 == 1 ~ 2, 
                  c1reap10 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.


### c1reap11: In an average week, how often do you eat sweets like cake, cookies, pastries, donuts,
### muffins, chocolate and candies, more than 2 times per day.
reap <- mutate(reap,
                c1reap11 = case_when(
                  c1reap11 == 0 ~ 3, 
                  c1reap11 == 1 ~ 2, 
                  c1reap11 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

### c1reap12: In an average week, how often do you drink 16 ounces or more of non-diet soda, fruit
### drink/punch or Kool-Aid a day? Note: 1 can of soda = 12 ounces
reap <- mutate(reap,
               c1reap12 = case_when(
                 c1reap12 == 0 ~ 3, 
                 c1reap12 == 1 ~ 2, 
                 c1reap12 == 2 ~ 1, 
                  TRUE ~ NA_real_ # This is for all other values 
                ))                # not covered by the above.

############# Write Files #############

write.table(reap,"./REAP/W1reap.recode.txt", sep="\t", col.names = TRUE)



