# Plot the region with a view to finding common variants with better signals than rs62625034

library(tidyverse)
library(ggrepel)
a <- read.csv("Downloads/results_combined.csv")
subset(a, snp == "rs62625034")

a %>% 
filter(binary == 1) %>%
ggplot(aes(x=eaf)) +
geom_histogram(aes(fill=binary)) +
facet_grid(. ~ binary)

a <- group_by(a, snp) %>%
do({
	.$eaf <- subset(., binary == 0)$eaf
	.
	})

a %>% 
filter(binary == 1) %>%
ggplot(aes(x=pos, y=-log10(a1_p))) +
geom_point(aes(size=snp=="rs62625034", colour=eaf>0.05)) +
geom_label_repel(
	data=subset(a, binary == 1 & eaf > 0.05 & a1_p < subset(a, binary == 1 & snp == "rs62625034")$a1_p), aes(label=snp))
ggsave("regional_plot.pdf")


