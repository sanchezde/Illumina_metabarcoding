##### R code for conducting read-minimum and read-threshold filtering
##### using a mock community sequenced in parallel to test samples.

### Needed are a TSV formatted biom file from Qiime2 and a list of
### ESV identifiers that do not align (97% ID) to a sequence file of expected sequences
### only in the mock communities. We also need a list of Sample names for the mock communities

### The script could also be adapted for ambiguous read percentage thresholds

library(tidyverse)

### read in entire feature table of all samples ###

ft <- read_tsv("feature-table_coi_arth.lenhostfilt.final.tsv",
               skip = 1)


#### 1. Read minimum filtering ####
# Keeping samples with reads greater or equal to 2000 (From a histogram)


# Making a column of samples to keep that have at least 2000 reads

samples_to_keep <- ft %>%
  gather(key = Sample, value = Reads, 2:ncol(ft)) %>% 
  filter(Reads > 0) %>% 
  group_by(Sample) %>%
  summarize(Read_num = sum(Reads)) %>% 
  ungroup() %>% 
  filter(Read_num >= 2000) %>% 
  select(Sample)


# Pivot the feature table to long-form for tidyverse manipulation
# Filter low-read samples from feature table

ft_f <- ft %>%
  gather(key = Sample, value = Reads, 2:ncol(ft)) %>% 
  filter(Sample %in% samples_to_keep$Sample) %>% 
  filter(Reads > 0)

### Reading in the mock community file for miss

da.miss <- read_tsv("mock.missseqs.txt", col_names = "HashID") %>% 
  mutate(MockAlign = "miss")

## Filter the mock_communities of interest

of_interest <- c("IM3Xa", "IM3Xb")

da.mock_basic <- ft_f %>%
  filter(Sample %in% of_interest)

## Calculate how many reads misses occur in

da.mock_basic %>% 
  group_by(Sample) %>%
  mutate(Prop.Read = Reads / sum(Reads) * 100) %>% 
  ungroup() %>% 
  filter(`#OTU ID` %in% da.miss$HashID) %>% 
  group_by(Sample) %>% 
  summarize(MaxReads = max(Reads),
            MaxPercent = max(Prop.Read))

# Output from above code chunk
#|Sample | MaxReads| MaxPercent|#
#|:------|--------:|----------:|#
#|IM3Xa  |     6828|       3.07|#
#|IM3Xb  |     6638|       3.58|#
   
### So will filter the bottom 3.6% of arthropod reads; seems a bit huge

## Filter out bottom 3.6% of reads per sample -- example
ft_f_p <- ft_f %>%
  group_by(Sample) %>% 
  mutate(p = Reads / sum(Reads) * 100) %>% 
  arrange(Sample, desc(Reads)) %>% 
  filter(p >= 3.6) %>% 
  select(-p)

## Back to biom tsv format

# back to wide form

ft_w <- ft_f_p %>%
  spread(Sample, Reads)

# replace na with 0

ft_w[is.na(ft_w)] <- 0

write_tsv(ft_w, "../final_final_feature.table.tsv")
