# Jaccard distance-based permutation test between presence/absence gene profiles
# within households, with randomization of a SINGLE compartment of choice 
# (Patient, Animal, Environment or Food).

# This script is identical to the other three compartment-based randomization scripts (patient/animal/environment/food), with the exception of the randomization section (lines 117-131), which randomizes only the Environment compartment

# Load the necessary packages
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)

# Open the data file (csv with sample_ID | patient_ID | sample_type | gene1 | gene2 | ... | geneX)
table <- as.data.frame(read.csv2("data_file"))

# Build the dataframe of combinations
ID <- expand.grid(table$sample_ID,table$sample_ID)
HH <- expand.grid(table$patient_ID,table$patient_ID)
ST <- expand.grid(table$sample_type,table$sample_type)
expanded_df <- cbind(ID,HH,ST)
colnames(expanded_df) <- c("ID1", "ID2", "HH1", "HH2", "sample_type1", "sample_type2")      # rename the columns 

# Extract only the gene matrix and set the sample_ID as rownames
gene_matrix <- (table[,4:ncol(table)])
rownames(gene_matrix) <- table$sample_ID

# Calculation of the Jaccard distance matrix
distance_jaccard<- vegdist(gene_matrix, method = "jaccard", binary = TRUE)
distance_matrix_jaccard <- as.matrix(distance_jaccard)                          # Convert to matrix
distance_table_jaccard <- as.data.frame(as.table(distance_matrix_jaccard))      # Put in long table format: sample_x | sample_y | distance
colnames(distance_table_jaccard) <- c("ID1", "ID2", "jaccard_distance")
expanded_df <- merge(expanded_df, distance_table_jaccard, by = c("ID1","ID2"))  # Merging distances into the combinations dataframe

# Delete the rows with self-comparisons (ID1=ID2), duplicates (keep only ID1 vs ID2 or ID2 vs ID1), and interhousehold comparisons (HH1 =/= HH2)
expanded_df_no_self_comparison <- expanded_df %>%
  filter(ID1 != ID2)                                                            # delete self-comparisons
expanded_df_no_self_comparison_unique <- expanded_df_no_self_comparison %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(ID1, ID2)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)                                                             # Delete duplicates
expanded_df_no_self_comparison_unique_intra_HH <- expanded_df_no_self_comparison_unique %>%
  filter(HH1 == HH2)                                                            # keep only intra household comparisons



# Compute the observed mean of distances (Mtot) for each comparison category (HH - HA - HE - HF - AA - AE - AF - EE - EF - FF)

## List of unique households
foyers <- unique(expanded_df_no_self_comparison_unique_intra_HH$HH1)

## Creating the empty results tables
summary_table_jaccard <- data.frame(
  HH = foyers,
  HuHu = NA, HuAni = NA, HuEnv = NA, HuFood = NA,
  AniAni = NA, AniEnv = NA, AniFood = NA,
  EnvEnv = NA, EnvFood = NA,
  FoodFood = NA, 
  stringsAsFactors = FALSE
)

## Filter comparisons by comparison type (including both directions, for example A-E and E-A)
get_subdf <- function(df, type1, type2) {
  df %>%
    filter((sample_type1 == type1 & sample_type2 == type2) |
           (sample_type1 == type2 & sample_type2 == type1))
}

## Loop the distance computation per comparison on each household
for (i in seq_along(foyers)) {
  hh <- foyers[i]
  df_hh <- expanded_df_no_self_comparison_unique_intra_HH %>% filter(HH1 == hh)          # Subset of the household
  
  # Calculation of means for each type of comparison
  summary_table_jaccard$HuHu[i]   <- mean(get_subdf(df_hh, "Patient", "Patient")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$HuAni[i]  <- mean(get_subdf(df_hh, "Patient", "Animal")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$HuEnv[i]  <- mean(get_subdf(df_hh, "Patient", "Environment")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$HuFood[i] <- mean(get_subdf(df_hh, "Patient", "Food")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$AniAni[i] <- mean(get_subdf(df_hh, "Animal", "Animal")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$AniEnv[i] <- mean(get_subdf(df_hh, "Animal", "Environment")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$AniFood[i]<- mean(get_subdf(df_hh, "Animal", "Food")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$EnvEnv[i] <- mean(get_subdf(df_hh, "Environment", "Environment")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$EnvFood[i]<- mean(get_subdf(df_hh, "Environment", "Food")$jaccard_distance, na.rm = TRUE)
  summary_table_jaccard$FoodFood[i]<- mean(get_subdf(df_hh, "Food", "Food")$jaccard_distance, na.rm = TRUE)
}

## Calculate the mean of each comparison category
Mtot_jaccard <- colMeans(summary_table_jaccard[ , -1], na.rm = TRUE)
names(Mtot_jaccard) <- paste0("Mtot_", names(Mtot_jaccard))                     ## Rename with the prefix "Mtot_"
Mtot_jaccard <- as.data.frame(t(Mtot_jaccard))                                  ## Convert to data.frame if necessary
write.csv(Mtot_jaccard, "Mtot_jaccard.csv")                                     ## Save distance table 



# Randomisation : Compute the null distribution of Mtot for each comparison category

n_iter <- 10000                                                                 ## Loop over 10000 rounds

## Preallocate a matrix 
n_var <- 10
mat_jaccard <- matrix(NA, nrow = n_iter, ncol = n_var)
colnames(mat_jaccard) <- names(Mtot_jaccard)

## Loop the Mtot computation on randomised dataframes for n_iter 
for (i in 1:n_iter) {
  
  ### Randomise the dataframe
  #### Separation of "table" by sample type
  df_patient <- table %>% filter(sample_type == "Patient")
  df_animal  <- table %>% filter(sample_type == "Animal")
  df_env     <- table %>% filter(sample_type == "Environment")
  df_food    <- table %>% filter(sample_type == "Food")

  #### Randomize by sample type:
  identifiants_env <- df_env %>% select(sample_ID, patient_ID, sample_type)     # Keep the metadata columns fixed
  genes_env_random <- df_env %>%
    select(-sample_ID, -patient_ID, -sample_type) %>%
    slice_sample(prop = 1)                                                      # Randomize the rows of the 0/1 gene matrix
  df_env_random <- bind_cols(identifiants_env, genes_env_random)       # Reconnect identifiers and randomized data

  ### Reassemble the table while maintaining the structure of each household in terms of number of H, A, E, and F:  
  table_metadata <- table %>% select(sample_ID, patient_ID, sample_type)        # Keep the metadata columns fixed
  genes_random <- rbind(df_patient,df_animal, df_env_random, df_food)           # Link the randomized tables of the different sample types together
  table_random <- merge(
    table_metadata,
    genes_random,
    by = c("sample_ID", "patient_ID", "sample_type"),
    all.x = TRUE
  )                                                                             # merge them with the first 3 columns, matching the sample IDs


  ### Compute the Jaccard distances for this new, randomized, table.

  #### Build the dataframe of combinations
  ID_r <- expand.grid(table_random$sample_ID,table_random$sample_ID)
  HH_r <- expand.grid(table_random$patient_ID,table_random$patient_ID)
  ST_r <- expand.grid(table_random$sample_type,table_random$sample_type)
  expanded_df_r <- cbind(ID_r,HH_r,ST_r)
  colnames(expanded_df_r) <- c("ID1", "ID2", "HH1", "HH2", "sample_type1", "sample_type2")       

  #### Extract only the gene matrix and set the sample_ID as rownames.
  gene_matrix_r <- (table_random[,4:ncol(table_random)])
  rownames(gene_matrix_r) <- table_random$sample_ID

  #### Calculation of the Jaccard distance matrix
  distance_jaccard_r <- vegdist(gene_matrix_r, method = "jaccard", binary = TRUE)
  distance_matrix_jaccard_r <- as.matrix(distance_jaccard_r)                    # Convert to matrix
  distance_table_jaccard_r <- as.data.frame(as.table(distance_matrix_jaccard_r))       # Put in long table format: sample_x | sample_y | distance
  colnames(distance_table_jaccard_r) <- c("ID1", "ID2", "jaccard_distance")
  expanded_df_r <- merge(expanded_df_r, distance_table_jaccard_r, by = c("ID1","ID2"))      # Merging distances into the combinations dataframe

  #### Delete the rows with self-comparisons (ID1=ID2), duplicates (keep only ID1 vs ID2 or ID2 vs ID1), and interhousehold comparisons (HH1 =/= HH2)
  expanded_df_no_self_comparison_r <- expanded_df_r %>%
    filter(ID1 != ID2)                                                          # delete self-comparisons
  expanded_df_no_self_comparison_unique_r <- expanded_df_no_self_comparison_r %>%
    rowwise() %>%
    mutate(pair_id = paste(sort(c(ID1, ID2)), collapse = "_")) %>%
    ungroup() %>%
    distinct(pair_id, .keep_all = TRUE) %>%
    select(-pair_id)                                                            # Delete duplicates
  expanded_df_no_self_comparison_unique_intra_HH_r <- expanded_df_no_self_comparison_unique_r %>%
    filter(HH1 == HH2)                                                          # keep only intra household comparisons

  #### Compute the observed mean of distances (Mtot) for each comparison category (HH - HA - HE - HF - AA - AE - AF - EE - EF - FF)
  ##### List of unique households
  foyers_r <- unique(expanded_df_no_self_comparison_unique_intra_HH_r$HH1)

  ##### Create the empty results table
  summary_table_jaccard_r <- data.frame(
    HH = foyers_r,
    HuHu = NA, HuAni = NA, HuEnv = NA, HuFood = NA,
    AniAni = NA, AniEnv = NA, AniFood = NA,
    EnvEnv = NA, EnvFood = NA,
    FoodFood = NA, 
    stringsAsFactors = FALSE
  )
  
  ##### Loop the distance computation per comparison on each household
  for (j in seq_along(foyers_r)) {
    hh <- foyers_r[j]
    df_hh_r <- expanded_df_no_self_comparison_unique_intra_HH_r %>% filter(HH1 == hh)          # Subset of household
    
    # Calculation of means for each type of comparison
    summary_table_jaccard_r$HuHu[j]   <- mean(get_subdf(df_hh_r, "Patient", "Patient")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$HuAni[j]  <- mean(get_subdf(df_hh_r, "Patient", "Animal")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$HuEnv[j]  <- mean(get_subdf(df_hh_r, "Patient", "Environment")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$HuFood[j] <- mean(get_subdf(df_hh_r, "Patient", "Food")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$AniAni[j] <- mean(get_subdf(df_hh_r, "Animal", "Animal")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$AniEnv[j] <- mean(get_subdf(df_hh_r, "Animal", "Environment")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$AniFood[j]<- mean(get_subdf(df_hh_r, "Animal", "Food")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$EnvEnv[j] <- mean(get_subdf(df_hh_r, "Environment", "Environment")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$EnvFood[j]<- mean(get_subdf(df_hh_r, "Environment", "Food")$jaccard_distance, na.rm = TRUE)
    summary_table_jaccard_r$FoodFood[j]<- mean(get_subdf(df_hh_r, "Food", "Food")$jaccard_distance, na.rm = TRUE)
  }
  
  ##### Calculate the mean of each comparison category
  Mtot_jaccard_r <- colMeans(summary_table_jaccard_r[ , -1], na.rm = TRUE)
  names(Mtot_jaccard_r) <- paste0("Mtot_", names(Mtot_jaccard_r))               # Rename with the prefix "Mtot_"
  Mtot_jaccard_r <- as.data.frame(t(Mtot_jaccard_r))                            # Convert to data.frame if necessary
  mat_jaccard[i, ] <- as.numeric(Mtot_jaccard_r)                                # Store results
}

## Convert the matrix to df
Mtot_jaccard_r <- as.data.frame(mat_jaccard)


# Graphical comparison of the observed mean (Mtot) with Mtots' null distributions (randomised Mtots, Mtot_r) 
## Observe the null distribution (can be gaussian) and compare the "real" Mtot (in blue) to the 0,05 thresholds of the AUC (in red).
real_val_jaccard_Huhu <- as.numeric(Mtot_jaccard$Mtot_HuHu) 
real_val_jaccard_HuAni <- as.numeric(Mtot_jaccard$Mtot_HuAni) 
real_val_jaccard_HuEnv <- as.numeric(Mtot_jaccard$Mtot_HuEnv) 
real_val_jaccard_HuFood <- as.numeric(Mtot_jaccard$Mtot_HuFood)
real_val_jaccard_AniAni <- as.numeric(Mtot_jaccard$Mtot_AniAni) 
real_val_jaccard_AniEnv <- as.numeric(Mtot_jaccard$Mtot_AniEnv) 
real_val_jaccard_AniFood <- as.numeric(Mtot_jaccard$Mtot_AniFood) 
real_val_jaccard_EnvEnv <- as.numeric(Mtot_jaccard$Mtot_EnvEnv) 
real_val_jaccard_EnvFood <- as.numeric(Mtot_jaccard$Mtot_EnvFood) 
real_val_jaccard_FoodFood <- as.numeric(Mtot_jaccard$Mtot_FoodFood) 

# Human-Human comparison
tiff("jac_HuHu.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")
ggplot(Mtot_jaccard_r, aes(x = Mtot_HuHu)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_Huhu, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuHu, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuHu, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuHu, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuHu, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Hu-Hu (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_HuHu (in blue = real value of Mtot_jaccard_HuHu)")
dev.off()

# Human-Animal comparison
tiff("jac_HuAni.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_HuAni)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_HuAni, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuAni, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuAni, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuAni, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuAni, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Hu-Ani (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_HuAni (in blue = real value of Mtot_jaccard_HuAni)")
dev.off()

# Human-Environment comparison
tiff("jac_HuEnv.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_HuEnv)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_HuEnv, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuEnv, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuEnv, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuEnv, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuEnv, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Hu-Env (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_HuEnv (in blue = real value of Mtot_jaccard_HuEnv)")
dev.off()

# Human-Food comparison
tiff("jac_HuFood.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")
ggplot(Mtot_jaccard_r, aes(x = Mtot_HuFood)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_HuFood, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuFood, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuFood, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuFood, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_HuFood, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Hu-Food (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_HuFood (in blue = real value of Mtot_jaccard_HuFood)")
dev.off()

# Animal-Animal comparison
tiff("jac_AniAni.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_AniAni)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_AniAni, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniAni, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniAni, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniAni, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniAni, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Ani-Ani (randomisations)", y = "Density",
     title = "Null distribution of Mtot_jaccard_AniAni (in blue = real value of Mtot_jaccard_AniAni)")
dev.off()

# Animal-Environment comparison
tiff("jac_AniEnv.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_AniEnv)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_AniEnv, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniEnv, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniEnv, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniEnv, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniEnv, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Ani-Env (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_AniEnv (in blue = real value of Mtot_jaccard_AniEnv)")
dev.off()

# Animal-Food comparison
tiff("jac_AniFood.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_AniFood)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_AniFood, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniFood, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniFood, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniFood, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_AniFood, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Ani-Food (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_AniFood (in blue = real value of Mtot_jaccard_AniFood)")
dev.off()


# Environment-Environment comparison
tiff("jac_EnvEnv.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_EnvEnv)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_EnvEnv, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvEnv, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvEnv, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvEnv, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvEnv, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Env-Env (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_EnvEnv (in blue = real value of Mtot_jaccard_EnvEnv)")
dev.off()

# Environment-Food comparison
tiff("jac_EnvFood.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_EnvFood)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_EnvFood, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvFood, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvFood, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvFood, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_EnvFood, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Env-Food (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_EnvFood (in blue = real value of Mtot_jaccard_EnvFood)")
dev.off()

# Food-Food comparison
tiff("jac_FoodFood.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")

ggplot(Mtot_jaccard_r, aes(x = Mtot_FoodFood)) +
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "white") +
  geom_density() +
  geom_vline(xintercept = real_val_jaccard_FoodFood, color = "blue", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_FoodFood, probs = 0.05, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_FoodFood, probs = 0.95, na.rm = TRUE),
             color = "red", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_FoodFood, probs = 0.025, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  geom_vline(xintercept = quantile(Mtot_jaccard_r$Mtot_FoodFood, probs = 0.975, na.rm = TRUE),
             color = "orange", size = 1, linetype = "dashed") +
  labs(x = "Mtot Food-Food (randomisations)", y = "Density",
       title = "Null distribution of Mtot_jaccard_FoodFood (in blue = real value of Mtot_jaccard_FoodFood)")
dev.off()


# Calculation of the p-value : proportion of randomised values lower than the observed value
p_value_jaccard_HuHu <- mean(Mtot_jaccard_r$Mtot_HuHu <= real_val_jaccard_Huhu, na.rm = TRUE)
p_value_jaccard_HuAni <- mean(Mtot_jaccard_r$Mtot_HuAni <= real_val_jaccard_HuAni, na.rm = TRUE)
p_value_jaccard_HuEnv <- mean(Mtot_jaccard_r$Mtot_HuEnv <= real_val_jaccard_HuEnv, na.rm = TRUE)
p_value_jaccard_HuFood <- mean(Mtot_jaccard_r$Mtot_HuFood <= real_val_jaccard_HuFood, na.rm = TRUE)
p_value_jaccard_AniAni <- mean(Mtot_jaccard_r$Mtot_AniAni <= real_val_jaccard_AniAni, na.rm = TRUE)
p_value_jaccard_AniEnv <- mean(Mtot_jaccard_r$Mtot_AniEnv <= real_val_jaccard_AniEnv, na.rm = TRUE)
p_value_jaccard_AniFood <- mean(Mtot_jaccard_r$Mtot_AniFood <= real_val_jaccard_AniFood, na.rm = TRUE)
p_value_jaccard_EnvEnv <- mean(Mtot_jaccard_r$Mtot_EnvEnv <= real_val_jaccard_EnvEnv, na.rm = TRUE)
p_value_jaccard_EnvFood <- mean(Mtot_jaccard_r$Mtot_EnvFood <= real_val_jaccard_EnvFood, na.rm = TRUE)
p_value_jaccard_FoodFood <- mean(Mtot_jaccard_r$Mtot_FoodFood <= real_val_jaccard_FoodFood, na.rm = TRUE)

## Print p values
comparisons <- c("HuHu", "HuAni", "HuEnv", "HuFood",
                  "AniAni", "AniEnv", "AniFood",
                  "EnvEnv", "EnvFood", "FoodFood")                              # Define the types of comparisons
table_pvalues <- data.frame(
  Comparisons = comparisons,
  Jaccard = NA
)                                                                               # Create an empty dataframe

table_pvalues$Jaccard     <- c(p_value_jaccard_HuHu, p_value_jaccard_HuAni, p_value_jaccard_HuEnv, p_value_jaccard_HuFood, p_value_jaccard_AniAni, p_value_jaccard_AniEnv, p_value_jaccard_AniFood, p_value_jaccard_EnvEnv, p_value_jaccard_EnvFood, p_value_jaccard_FoodFood)                 # fill the table with p-values
write.csv(table_pvalues, "table_pvalues.csv")                                   # Save to .csv format


# Create a stabilization by randomization graph based on the p-value calculated at each round : to check that the number of randomisation iterations is sufficient
## Fonction : p-value + 95% IC
calc_pvalues <- function(M_real, M_rand) {
  n_iter <- nrow(M_rand)
  res_list <- list()
  for (col in colnames(M_rand)) {
    pvals <- numeric(n_iter)
    lower <- numeric(n_iter)
    upper <- numeric(n_iter)
    for (i in 1:n_iter) {
      p_hat <- mean(M_rand[[col]][1:i] < M_real[[col]], na.rm = TRUE)
      se <- ifelse(i > 1, sqrt(p_hat * (1 - p_hat) / i), NA)
      pvals[i] <- p_hat
      lower[i] <- max(p_hat - 1.96 * se, 0)
      upper[i] <- min(p_hat + 1.96 * se, 1)
    }
    res_list[[col]] <- data.frame(
      iteration = 1:n_iter,
      category = col,
      pvalue = pvals,
      lower = lower,
      upper = upper
    )
  }
  return(bind_rows(res_list))
}

## Calculate p-values
p_jac_long <- calc_pvalues(Mtot_jaccard, Mtot_jaccard_r)

## Automatic annotation : retrieve the latest iteration for each category.
last_points_jac <- p_jac_long %>%
  group_by(category) %>%
  slice_tail(n = 1) %>%
  mutate(label = sprintf("p = %.3f [%.3f–%.3f]", pvalue, lower, upper))

## Graph with CI and annotations
tiff("p-values_evolution_jac.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")
ggplot(p_jac_long, aes(x = iteration, y = pvalue, color = category, fill = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  geom_text(
    data = last_points_jac,
    aes(label = label),
    hjust = -0.05, vjust = 0.5,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  expand_limits(x = max(p_jac_long$iteration) * 1.5) + # space for text
  labs(title = "P-values Jaccard distance (with 95% CI)",
       x = "Iteration",
       y = "P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()

## “faceting” (one chart per category)
tiff("p-values_evolution_jac_faceting.tiff",
     width = 2000, height = 1500, res = 300, compression = "lzw", type = "cairo")
ggplot(p_jac_long, aes(x = iteration, y = pvalue, color = category, fill = category)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1) +
  geom_text(
    data = last_points_jac,
    aes(label = label),
    hjust = 1, vjust = 0,
    size = 3.5,
    show.legend = FALSE
  ) +
  facet_wrap(~ category, scales = "free_y") +
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "P-values Jaccard distance (with 95% CI)",
       x = "Iteration",
       y = "P-value") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()