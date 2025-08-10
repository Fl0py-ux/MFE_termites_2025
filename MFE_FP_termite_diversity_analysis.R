### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#-------- EvALUATION OF TERMITE DIVERSITY IN FMNP, TOGO -------------

#Florence Pivetta
#MFE 2025

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# ___________________________________________________________________

# --------------------------1. Loading data -------------------------

#____________________________________________________________________

library(fossil)
library(readxl)
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(iNEXT)
library(ggtext)

#Dataset has one column with an ID, one column with the family name and one column with the speices naem.
#Each line of the dataset is one occurence.

termites = read_xlsx(" ")#file path


termites <-termites %>%  #pipe
  mutate(famille_espece = paste(Genus, species, sep = " "))


termites <- termites %>% 
  mutate(presence = 1) %>%
  mutate(across(c(Genus, species, famille_espece), tolower))
termites$`IDcolumn`<-as.numeric(as.character(termites$`IDcolumn`)) #mettre les ID des échantillons en numérique



# ___________________________________________________________________

#---------------------2. Accumulation curve ------------------------

# ___________________________________________________________________


matrix_pres <- termites %>%
  pivot_wider(names_from = famille_espece, values_from = presence, values_fill = 0, values_fn = max )
matrix_pres <- matrix_pres[, -c(2:4)]  #Creating a presence-absence matrix

#creating the accumulation curve

accum = specaccum(matrix_pres%>% select(where(is.numeric)), method = "random", permutations = 1000)
accum_df = data.frame(
  sites = accum$sites,
  richness = accum$richness,
  sd = accum$sd)


ggplot(accum_df, aes(x = sites, y = richness)) +
  geom_line(color = "#117a65", linewidth = 1) +
  
  geom_ribbon(
    aes(ymin = richness - sd, ymax = richness + sd),
    fill = "#117a65", alpha = 0.2  
  ) +
  labs(
    x = "Number of samples",
    y = "Species richness",
    title = "Species accumulation curve"
  ) +
  theme_minimal(base_size = 13) +  # police plus lisible
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# ___________________________________________________________________

#-------------------------3. Species sampled ------------------------

# ___________________________________________________________________


#especes les plus communes
#faire plot 

matrix_work <- matrix_pres[, -1] # I remove the ID column from the matrix
head(matrix_work)

species_counts <- colSums(matrix_work) # total number of occurences per species


freq <- colSums(matrix_work)
df_freq <- data.frame(
  Espece = names(freq),
  Frequence = as.numeric(freq)) %>%
  arrange(desc(Frequence))  #Sorting the species from the more to the less numerous


df_freq$Espece <- tolower(df_freq$Espece)  
df_freq$Espece <- sapply(strsplit(df_freq$Espece, " "), function(x) {
  if (length(x) == 2) {
    genre <- paste0(toupper(substring(x[1], 1, 1)), substring(x[1], 2))
    espece <- x[2]
    paste0("<i>", genre, " ", espece, "</i>")
  } else {
    paste0("<i>", paste(x, collapse = " "), "</i>")
  }
})   #writing genus and species in italics



#making a barplot

ggplot(df_freq, aes(x = reorder(Espece, Frequence), y = Frequence)) +
  geom_bar(stat = "identity", fill = "#117a65", color = "black") +
  geom_text(aes(label = Frequence), hjust = -0.2, size = 3.5) + 
  coord_flip() +
  labs(title = "Number of samples per species", 
       x = "Species", 
       y = "Occurence frequency") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_markdown(),
    axis.text = element_text(color = "black") 
  )



# ___________________________________________________________________

#----------------------4. Diversity indices -------------------------

# ___________________________________________________________________

## Simpson Index ----

#I group all of the data to compute Simpson index
total_abundances <- colSums(matrix_pres%>% select(-1))
total_pool <- matrix(total_abundances, nrow = 1)

simpson_global <- diversity(total_pool, index = "simpson") 
print(simpson_global)


#bootstrapping to estimate a confidence interval
set.seed(123)
n_boot <- 1000

boot_simpson <- numeric(n_boot)

for (i in 1:n_boot) {
  # Resampling the rows of matrix_pres with replacement
  resample <- matrix_pres[sample(1:nrow(matrix_pres), replace = TRUE), ]
  
  # Sum of columns (presences per species) after removing the ID column
  total_abund <- colSums(resample[, -1])
  total_pool <- matrix(total_abund, nrow = 1)
  
  # Simpson calculation
  boot_simpson[i] <- diversity(total_pool, index = "simpson")
}

# Confidence interval 95 %
quantile(boot_simpson, c(0.025, 0.975))



## Shannon Index ----

#I group all of the data to compute Shannon index
shannon_global <- diversity(colSums(matrix_pres[, -1]), index = "shannon")
print(shannon_global)


#bootstrapping to estimate a confidence interval
set.seed(123)
n_boot <- 1000

boot_shannon <- numeric(n_boot)

for (i in 1:n_boot) {
  # Resampling the rows of matrix_pres with replacement
  resample <- matrix_pres[sample(1:nrow(matrix_pres), replace = TRUE), ]
  
  # Sum of columns (presences per species) after removing the ID column
  total_abund <- colSums(resample[, -1])
  
  # Shannon calculation
  boot_shannon[i] <- diversity(total_abund, index = "shannon")
}


# Confidence Interval 95%
quantile(boot_shannon, c(0.025, 0.975))



##Pielou's Evenness ----

#Pielou index computation - Shannon/log(species richness)
S <- specnumber(total_abundances)
print(S)
Pielou <- (shannon_global / log(S))
print(Pielou)


set.seed(123)
n_boot <- 1000

boot_pielou <- numeric(n_boot)


#bootstrapping to estimate a confidence interval
for (i in 1:n_boot) {
  # Resampling the rows of matrix_pres with replacement
  resample <- matrix_pres[sample(1:nrow(matrix_pres), replace = TRUE), ]
  
  # Sum of columns (presences per species) after removing the ID column
  total_abund <- colSums(resample[, -1])
  
  # Shannon
  H <- diversity(total_abund, index = "shannon")
  
  # Species richness
  S <- specnumber(total_abund)
  
  # avoiding log(0)
  if (S > 1) {
    boot_pielou[i] <- H / log(S)
  } else {
    boot_pielou[i] <- NA
  }
}

boot_pielou <- na.omit(boot_pielou)

# Confidence Interval 95\%
quantile(boot_pielou, c(0.025, 0.975))


# ___________________________________________________________________

#----------------5. Species richness estimators ---------------------

# ___________________________________________________________________

## Chao Estimator ----

chao2(matrix_work, taxa.row = FALSE)

## Jacknife Estimators ----

jack1(matrix_work, taxa.row = FALSE, abund = FALSE)
jack2(matrix_work, taxa.row = FALSE, abund = FALSE)

##Incidence-based Coverage Estimator ----

ICE(matrix_work, taxa.row = FALSE) 



## calculation of the estimators(via 'specpool') ----
specpool_result <- specpool(matrix_work)



# Backup: calulation with package 'fossil' as a cross-check

# Chao2 estimator
chao2(matrix_work, taxa.row = FALSE)

# Jackknife Estimators

jack1(matrix_work, taxa.row = FALSE, abund = FALSE)
jack2(matrix_work, taxa.row = FALSE, abund = FALSE)

# Incidence-based Coverage Estimator (ICE)

ICE(matrix_work, taxa.row = FALSE) 

#building a plot

#I manually write the dataframe to later compute confidence itervals
estimates <- data.frame(
  method = c("Chao2", "ICE", "Jackknife1", "Jackknife2", "Bootstrap"),
  estimate = c(37.71012,37.15957, 41.96889, 40.02649, 38.91457),
  se = c(2.779068, NA, 2.633992, NA, 1.636109)  
)


# Confidence intervals 95%
estimates <- estimates %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )
head(estimates$lower)
head(estimates$upper)

estimates$method <- factor(estimates$method,
                           levels = c("Chao2", "ICE", "Jackknife1", "Jackknife2", "Bootstrap"))


#summary plot for the species richness estimators

ggplot(estimates, aes(x = method, y = estimate)) +
  geom_point(size = 4, color = "#117a65") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, color = "#117a65") +
  geom_hline(yintercept = 35, linetype = "dashed", color = "indianred") +
  labs(
    title = "Termite species richness estimators",
    x = "Estimation method",
    y = "Estimated species richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black") 
  )


