library(tidyverse)


tooltip_css <- "background-color:black;color:white;padding:10px;border-radius:10px;font-family:sans-serif;font-weight:bold"

##########################

df <- read.csv('S2-SATTs for radar plot.csv', na.strings=c("",".","NA", "NaN", "#N/A", "#VALUE!"), stringsAsFactors = TRUE, header = F) 

new_header <- paste((df[1,]) %>% unlist(use.names=FALSE),(df[2,]) %>% unlist(use.names=FALSE), sep = '_')

df <- df[-1:-2, ]
colnames(df) <- new_header

# Remove duplicate GENE_NAME
dup <- df$GENE_NAME[duplicated(df$GENE_NAME)]

if (length(dup) == 0) dup <- 'none!'

print(paste("These Gene(s) were present multiple times and deleted from the analysis:", dup))

df <- df %>% filter(GENE_NAME != dup)

df_tidy <- df %>% gather(KEY, Value, 2:ncol(df))

df_tidy <- df_tidy %>% separate(KEY, into = c("Sumo", "metric"), sep = "_") %>% 
  spread(metric, Value)

# df_tidy <- df_tidy %>% filter(!is.na(Diff))
# Convert factors to numbers
df_tidy <- df_tidy  %>% 
  mutate_at(vars(Diff,minLogP,`SATT index`), as.numeric)

#Rename Z3 and ZATT
df_tidy <- df_tidy %>% mutate(Sumo =
                                case_when(Sumo == "Z3" ~ "LAZSUL",
                                          Sumo == "ZATT" ~"ZNF451",
                                          TRUE ~ Sumo)
)



# Minimal dataframe for processing is saved
df_tidy %>% write.csv('S2-SATTs_tidy.csv', row.names = F)



