library(tidyverse)
library(ggrepel)
library(ggiraph)
library(dplyr)
library(glue)
library(scales)

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

df_tidy <- df_tidy %>% group_by(GENE_NAME) %>% mutate(ID = row_number(GENE_NAME))

# auto
# max_Diff <- max(df_tidy$Diff, na.rm=T)
# Manual
max_Diff <- 14

max_P <- max(df_tidy$minLogP, na.rm=T)

df_tidy <- df_tidy %>% mutate(x = (x= ID +(Diff/max_Diff)))

#Select top hits
# df_top <- df_tidy %>% arrange(desc(`Diff`)) %>% top_n(5,`Diff`)

# top_hits <- c("RANGAP1", "SP100", "PML", "RANGAP1", "YEATS2")
# top_hits <- c("RANGAP1", "SP100")

# df_top <- df_tidy %>% filter(GENE_NAME %in% top_hits)

# df_top$Diff[is.na(df_top$Diff)] <- 0
# df_top$minLogP[is.na(df_top$minLogP)] <- 0
# df_top <- df_top %>% mutate(x = (x= ID +(Diff/max_Diff)))

df_sectors <- df_tidy %>% ungroup() %>% select(Sumo, ID) %>% distinct()
df_sectors <- df_sectors %>% mutate(x=ID, y=max_P/2, width=1, height=max_P*1.1)
  
  
#######

# p <- ggplot(df_tidy, aes(x=x, y=minLogP, color=`SATT index`)) + geom_point()  + theme_light() + scale_color_viridis_c(limits=c(0, 1))


p2 <- ggplot(df_tidy, aes(x=x, y=minLogP))+
  geom_tile(data=df_sectors, aes(x=x+0.5, y=y, width=width, height=height, fill=factor(x)), color='black', alpha=0.1)+
  
  annotate(geom='segment', x = 1, xend=1, y=max_P*1.1, yend=max_P*1.15, color='black', size=.25) +
  annotate(geom="text", x = 1, y=max_P*1.25, label = '0.0', size=3) +  
  annotate(geom="text", x = 1.95, y=max_P*1.4, label = 'Change', size=4) +  
  annotate(geom='segment', x = 1, xend=2, y=max_P*1.1, yend=max_P*1.1, color='black', size=.25) +
  annotate(geom='segment', x = 2, xend=2, y=max_P*1.1, yend=max_P*1.15, color='black', size=.25) +
  annotate(geom="text", x = 2, y=max_P*1.25, label = paste(round(max_Diff,1)), size=3) +

  geom_point_interactive(aes(color=`SATT index`, tooltip = glue("{GENE_NAME}\nDifference: {Diff}\n-Log[p]: {minLogP}"), data_id = GENE_NAME))+
  
    
  # geom_point(data=df_top, aes(x=x, y=minLogP), color='black', shape=21,size=5,fill=NA, alpha=1) + 

  # geom_text_repel(
  #   data = df_top,
  #   aes(label = GENE_NAME),
  #   size = 4,
  #   box.padding = unit(0.4, "lines"),
  #   # nudge_x = 0.2,
  #   # nudge_y       = 2,
  #   segment.size  = 0.2,
  #   segment.color = "black"
  # ) + 
  # geom_polygon(data=df_top, aes(x=x, y=minLogP, group=GENE_NAME),color='black', fill=NA)+
  coord_polar(start = 0)+

  # geom_path(aes(group=GENE_NAME)) +
  scale_color_viridis_c(limits=c(0, 1), oob=squish)+
  theme_light()

p2 <- p2 + geom_label(data=df_sectors, aes(x=x+0.5, y=height*1.25, label=Sumo, fill=factor(x)), alpha=0.2)

#Remove legend for sectors
p2 <- p2 + guides(fill = FALSE)

#Remove grid
p2 <- p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())

p2 <- p2 + labs(x = NULL, title ='Polar Volcano Plot')


# Dim non-identified
g <- girafe(ggobj = p2, options = list(
  opts_hover_inv(css = "opacity:0.5;"),
  opts_hover(css = "stroke:black;r:4pt"),
  opts_tooltip(css= tooltip_css,delay_mouseover = 0, offy = -90,opacity = .75),
  opts_zoom(min = 1, max = 4)
))

g

#Save as HTML widget
htmlwidgets::saveWidget(file = "PolarVolcano.html", g)

write.csv(df_tidy, 'tidy_S2-SATTs for radar plot.csv')

