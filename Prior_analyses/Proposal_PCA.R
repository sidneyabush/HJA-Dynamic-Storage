## Adapting HJA post-fire scripts for Jaime's stream chem analysis:
## Creating Figure 8 PCA with experimental watersheds mixed in with 13 month synoptic sampling:

# Clear environment
rm(list = ls())

librarian:: shelf(ggpmisc, ggpubr, dplyr, reshape, tidyr, ggalt, GGally, stringr, PCAtools, remotes, TigR, ggstream)

# Set working directory and import data
setwd("/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/HJA_LongTerm_Stream_Chem")
synoptic <- read.csv("./PCA_Full_Sidney.csv", header = TRUE, na.strings = "n.a.", stringsAsFactors = FALSE) %>%
  select(-Unique, -SO4) %>%
  filter(!Site %in% c("LO5", "MC")) %>%  # Exclude specific sites
  mutate(SamplingDate = as.Date(SamplingDate, format = "%m/%d/%Y"),
         Year = format(SamplingDate, "%Y"))

# Read and preprocess the CSV file
long_term <- read.csv("./CF00201_v6.csv", header = TRUE, na.strings = "n.a.", stringsAsFactors = FALSE) %>%
  mutate(DATE_TIME = as.POSIXct(DATE_TIME, format = "%Y-%m-%d %H:%M:%S"),  # Correct time format with minutes (%M)
         Year = format(DATE_TIME, "%Y")) %>%
  dplyr::rename(Si = SI, Na = NA., Ca = CA, Mg = MG, Cl = CL, Site = SITECODE, SamplingDate = DATE_TIME) %>%
  select(-contains("CODE"), -STCODE, -ENTITY, -WATERYEAR, -LABNO, -TYPE, -INTERVAL, -MEAN_LPS, -Q_AREA_CM, 
         -PH, -COND, -ALK, -ALKCODE, -SSED, -UTP, -TDP, -PARTP, -UTN, -TDN, -DON, -PARTN, -UTKN, -TKN, -PVOL, -ANCA,
         -PO4P, -NH3N, -NO3N, -SO4S, -DOC) %>%
  na.omit()  # Remove rows with NA values

# Step 1: Calculate the total number of unique sites
total_sites <- long_term %>% 
  select(Site) %>% 
  distinct() %>% 
  nrow()

# Step 2: Identify years with complete data for all sites
complete_years <- long_term %>%
  group_by(Year) %>%
  summarise(site_count = n_distinct(Site)) %>%
  filter(site_count == total_sites) %>%
  pull(Year)

# Step 3: Filter the original dataframe to keep only the complete years
long_term_filtered <- long_term %>%
  filter(Year %in% complete_years)

# Remove outliers of all tracer solutes
tracers <- c("Ca", "Mg", "Na", "K", "Si")
cols_to_consider <- tracers
sd_limit <- 2
remove_outliers <- function(data_to_filter, cols = cols_to_consider, limit = sd_limit){
  z_scores <- as.data.frame(sapply(data_to_filter[cols], function(data) (abs(data-mean(data, na.rm = TRUE))/sd(data, na.rm = TRUE))))    
  return(subset(data_to_filter, !rowSums(z_scores>limit, na.rm = TRUE)))}

synoptic <- remove_outliers(synoptic, tracers)
long_term_filtered <-remove_outliers(long_term_filtered)

# Combine the data sets and add a month column
all_chemistry <- bind_rows(synoptic, long_term_filtered) %>%
  mutate(month = format(as.Date(SamplingDate), "%m") %>% as.factor())

## Now compute averages by group: 
chem_averages <- all_chemistry %>% 
  group_by(Site)%>% 
  summarize(
    Ca = mean(Ca, na.rm = TRUE),
    Mg = mean(Mg, na.rm = TRUE),
    Na = mean(Na, na.rm = TRUE),
    K = mean(K, na.rm = TRUE),
    Si = mean(Si, na.rm = TRUE),
  )

## Need to match the types/sources of each sample to three categories: trib/ main stem for synoptic and longterm for longterm data:
types <- read.csv("./type_guide_Figure8.csv", header =TRUE, na.strings=c("n.a."),stringsAsFactors=FALSE)
chem_averages <-merge(chem_averages, types, by="Site")

chem_averages$unique<-paste0(chem_averages$Site, chem_averages$Type) #create a column of "unique" values - type of sample and date collected
chem_averages <-chem_averages[!duplicated(chem_averages$unique),] #remove all duplicated samples

#set up dataframe for centered and scaled data
center_all <- chem_averages

# create numerical columns variable
num_cols <- c(2:6)
center_all[num_cols] <- sapply(center_all[num_cols], as.numeric)

#center and scale data
center_all[num_cols]<-data.frame(sapply(center_all[num_cols], scale))

#create datafame of just data to use in PCA (pull out only numeric columns)
final_HJA_mat <- center_all[num_cols]

#transpose data to match required format of PCAtools
final_HJA_mat_t<-data.frame(t(final_HJA_mat))

# #rename columns in df as the "unique" values created above
colnames(final_HJA_mat_t)<-center_all$unique

#set all other columns in dataframe as metadata
metadata <-center_all[,c(1,7)]

#rename columns in df as the "unique" values created above
rownames(metadata) <- center_all$unique

##run PCA
pca<-pca(final_HJA_mat_t, metadata = metadata)

## order the sources, now precipitation will always appear first and gw will always appear last
pca$metadata$Site <- factor(pca$metadata$Site,
                            levels = c("CC", "LC", "NC", "MR", "GSWS06", "GSWS07", "GSWS08", "GSMACK", "GSWS01", "GSWS02", "GSWS09","GSWS10","GSLOOK"),
                            labels = c("Cold", "Longer", "Nostoc", "McRae", "WS6", "WS7", "WS8", "Mack", "WS1", "WS2", "WS9", "WS10", "Lookout"))

pca$metadata$Type <- factor(pca$metadata$Type,
                            levels = c("main", "trib", "longterm_up", "longterm_low"),
                            labels = c("Main Stem", "Tributary", "Long-term Upper", "Long-term Lower"))

Fig1C <- biplot(pca, 
               showLoadings = TRUE, 
               colby = "Site", 
               shape = "Type",
               pointSize = 6,
               shapekey = c(17, 15, 15),
               colkey = c(
                 "#facba6", 
                 "#f8b58b", 
                 "#f59e72", 
                 "#f2855d", 
                 "#c1e7ff",
                 "#abd2ec",
                 "#94bed9",
                 "#7faac6",
                 "#6996b3",
                 "#5383a1",
                 "#3d708f",
                 "#255e7e",
                 "#004c6d"
               ),
               legendPosition = NULL,
               legendLabSize = NULL,
               colLegendTitle = NULL,
               shapeLegendTitle = NULL,
               legendTitleSize = 16,
               legendIconSize = 5,
               sizeLoadingsNames = 5,
               boxedLoadingsNames = TRUE,
               fillBoxedLoadings = alpha("white", 1/4),
               lab = NULL, 
               ylim = c(-5.5, 3),
               xlim = c(-5, 4),
               max.overlaps = Inf, 
               axisLabSize = 18,
               titleLabSize = 20,
               borderWidth = 0.5,
               borderColour = "gray",
               subtitle = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title = element_text(size = 16, face = "bold", vjust = 1)
  ) +
  ggrepel::geom_label_repel(aes(label = pca$metadata$Site, color = pca$metadata$Site), 
                            size = 5,
                            fontface = "bold",  # Make text bold
                            box.padding = 0.9, 
                            point.padding = 0.8, 
                            max.overlaps = 30,  # Adjust overlap threshold
                            fill = alpha("white", 3),  # White background
                            segment.color = "darkgray") +
  geom_point(aes(color = pca$metadata$Site, shape = pca$metadata$Type), 
             size = 6) + # Connector line color
  scale_color_manual(values = c(
    "#facba6", 
    "#f8b58b", 
    "#f59e72", 
    "#f2855d", 
    "#c1e7ff",
    "#abd2ec",
    "#94bed9",
    "#7faac6",
    "#6996b3",
    "#5383a1",
    "#3d708f",
    "#255e7e",
    "#004c6d"
  ))

Fig1C

#ggsave("LongTerm_HJA_proposal.png", plot = Fig1C, width = 7, height = 7.5, units = "in")
