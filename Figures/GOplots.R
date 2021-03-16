#Rscript to create GO plots in Figure 3 and SI Appendix 

library("ggplot2")
library("ggrepel")

#Read in GO results file
mega_conv3B <- read.csv("MegALL_corr_2020-0827.csv", stringsAsFactors = FALSE)

#Rename columns
colnames(mega_conv3B) <- c("Mega", "Outgroup1", "Outgroup1$Conv.comp", "Outgroup2", "Outgroup2$Conv.path", "Count", "AA.conv.possible.GO", "AA.conv.genome", "AA.conv.possible", "Odds.ratio", "Pval", "Genes", "BH_corr")

##Subset colobine and horse as well as the different outgroups (micMur and Propithecus)##

mega_colobine3b <- subset(mega_conv3B, micMur1.bosTau4=="micMur1$colobine")
colnames(mega_colobine3b) <- c("Mega", "Outgroup1", "Conv.comp", "Out2", "Conv.path", "Count", "AA.conv.possible.GO", "AA.conv.genome", "AA.conv.possible", "Odds.ratio", "Pval", "Genes", "BH_corr")
mega_colobine3c <- subset(mega_conv3B, micMur1.bosTau4=="Propithecus$colobine")
colnames(mega_colobine3c) <- c("Mega", "Outgroup1", "Conv.comp", "Out2", "Conv.path", "Count", "AA.conv.possible.GO", "AA.conv.genome", "AA.conv.possible", "Odds.ratio", "Pval", "Genes", "BH_corr")

mega_horse3b <- subset(mega_conv3B, micMur1.bosTau4=="micMur1$equCab2")
colnames(mega_horse3b) <- c("Mega", "Outgroup1", "Conv.comp", "Out2", "Conv.path", "Count", "AA.conv.possible.GO", "AA.conv.genome", "AA.conv.possible", "Odds.ratio", "Pval", "Genes", "BH_corr")
mega_horse3c <- subset(mega_conv3B, micMur1.bosTau4=="Propithecus$equCab2")
colnames(mega_horse3c) <- c("Mega", "Outgroup1", "Conv.comp", "Out2", "Conv.path", "Count", "AA.conv.possible.GO", "AA.conv.genome", "AA.conv.possible", "Odds.ratio", "Pval", "Genes", "BH_corr")

##Characters as numberic##
mega_colobine3b$AA.conv.possible.GO = as.numeric(as.character(mega_colobine3b$AA.conv.possible.GO))
mega_colobine3b$Count = as.numeric(as.character(mega_colobine3b$Count))

mega_colobine3c$AA.conv.possible.GO = as.numeric(as.character(mega_colobine3c$AA.conv.possible.GO))
mega_colobine3c$Count = as.numeric(as.character(mega_colobine3c$Count))

mega_horse3b$Count = as.numeric(as.character(mega_horse3b$Count))
mega_horse3b$AA.conv.possible.GO = as.numeric(as.character(mega_horse3b$AA.conv.possible.GO))

mega_horse3c$Count = as.numeric(as.character(mega_horse3c$Count))
mega_horse3c$AA.conv.possible.GO = as.numeric(as.character(mega_horse3c$AA.conv.possible.GO))

##Colobine plots with micMur and Propithecus outgroups##

#micMur outgroup
mega_colobine3b_plot <- ggplot(mega_colobine3b, aes(x=AA.conv.possible.GO, y=Count, label=Conv.path, color = ifelse( Count >= 5, "5+ events", "<5 events"))) +
  theme_bw() +
  geom_point(size=6, alpha=0.45) + 
  stat_smooth(method="lm", se=TRUE, color="navy", size=0.5,) +
  scale_color_manual(name="Count", values = c("darkgray","blue"))+
  geom_text_repel(data = subset(mega_colobine3b, mega_colobine3b$Conv.path == "panTro2$hydrolase activity"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_colobine3b, mega_colobine3b$Conv.path == "panTro2$hydrolase activity"), color="darkred", size=6) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line = element_line())
mega_colobine3b_plot 

#propithecus outgroup
mega_colobine3c_plot <- ggplot(mega_colobine3c, aes(x=AA.conv.possible.GO, y=Count, label=Conv.path, color = ifelse( Count >= 5, "5+ events", "<5 events"))) +
  theme_bw() +
  geom_point(size=6, alpha=0.45) + 
  stat_smooth(method="lm", se=TRUE, color="navy", size=0.5,) +
  scale_color_manual(name="Count", values = c("darkgray","blue"))+
  geom_text_repel(data = subset(mega_colobine3c, mega_colobine3c$Conv.path == "panTro2$hydrolase activity"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_colobine3c, mega_colobine3c$Conv.path == "panTro2$hydrolase activity"), color="darkred", size=6) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line = element_line())
mega_colobine3c_plot

##Horse plots with micMur and Propithecus outgroups##

#micMur
mega_horse3b_plot <- ggplot(mega_horse3b, aes(x=AA.conv.possible.GO, y=Count, label=Conv.path, color = ifelse( Count >= 5,"5+ events", "<5 events"))) +
  theme_bw() +
  geom_point(size=6, alpha=0.45) + 
  stat_smooth(aes(group = 1), method="lm", se=TRUE, color="navy", size=0.5) +
  scale_color_manual(name="Count", values = c("darkgrey","blue"))+
  geom_text_repel(data = subset(mega_horse3b, mega_horse3b$Conv.path == "myoLuc1$brush border"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_horse3b, mega_horse3b$Conv.path == "myoLuc1$brush border"), color="darkred", size=6) +
  geom_text_repel(data = subset(mega_horse3b, mega_horse3b$Conv.path == "myoLuc1$innate immune response"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_horse3b, mega_horse3b$Conv.path == "myoLuc1$innate immune response"), color="red", size=6) +
  geom_text_repel(data = subset(mega_horse3b, mega_horse3b$Conv.path == "myoLuc1$serine-type endopeptidase activity"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_horse3b, mega_horse3b$Conv.path == "myoLuc1$serine-type endopeptidase activity"), color="pink", size=6) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line = element_line())+
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, by = 25)) +
  scale_x_continuous(limits=c(0,600000))
mega_horse3b_plot

#Propithecus
mega_horse3c_plot <- ggplot(mega_horse3c, aes(x=AA.conv.possible.GO, y=Count, label=Conv.path, color = ifelse( Count >= 5, "5+ events", "<5 events"))) +
  theme_bw() +
  geom_point(size=6, alpha=0.45) + 
  stat_smooth(aes(group = 1), method="lm", se=TRUE, color="navy", size=0.5) +
  scale_color_manual(name="Count", values = c("darkgray","blue"))+
  geom_text_repel(data = subset(mega_horse3c, mega_horse3c$Conv.path == "myoLuc1$brush border"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_horse3c, mega_horse3c$Conv.path == "myoLuc1$brush border"), color="darkred", size=6) +
  geom_text_repel(data = subset(mega_horse3c, mega_horse3b$Conv.path == "myoLuc1$innate immune response"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_horse3c, mega_horse3c$Conv.path == "myoLuc1$innate immune response"), color="red", size=6) +
  geom_text_repel(data = subset(mega_horse3c, mega_horse3c$Conv.path == "myoLuc1$serine-type endopeptidase activity"), nudge_x = 5, nudge_y=2, point.padding = .6) +
  geom_point(data = subset(mega_horse3c, mega_horse3c$Conv.path == "myoLuc1$serine-type endopeptidase activity"), color="magenta", size=6) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.line = element_line())+
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 20)) +
  scale_x_continuous(limits=c(0,600000))
mega_horse3c_plot
