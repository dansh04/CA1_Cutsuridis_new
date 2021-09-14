# INSTALL GGPLOT2 
library(ggplot2)
library(glue)
library(dplyr)
library(hrbrthemes)

# AVG CURRENT VS. MG CONCENTRATION
mgconcs <- c('0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6', '1.8', '2.0')
setwd('C:/Users/minse/OneDrive/Documents/GitHub/CA1_Cutsuridis')
nonpattern_avgcur <- c()
pattern_avgcur <- c()

for (mg in mgconcs) {
  nonpattern_dir <- 'pyresults/mgf_{mg}_celli_0.dat'
  nonpattern_data <- read.table(glue(nonpattern_dir))
  nonpattern_cot <- nonpattern_data$V2
  nonpattern_avgcur <- append(nonpattern_avgcur, mean(nonpattern_cot))
  pattern_dir <- 'pyresults/mgf_{mg}_celli_1.dat'
  pattern_data <- read.table(glue(pattern_dir))
  pattern_cot <- pattern_data$V2
  pattern_avgcur <- append(pattern_avgcur, mean(pattern_cot))
}

mgconcs_df <- as.data.frame(as.numeric(mgconcs))
nonpattern_avgcur_df <- as.data.frame(nonpattern_avgcur, col.names=)
nonpattern_avgcur_df <- cbind(mgconcs_df, nonpattern_avgcur_df)
colnames(nonpattern_avgcur_df) <- c('Mg.Conc', 'Avg.Cur')
pattern_avgcur_df <- as.data.frame(pattern_avgcur)
pattern_avgcur_df <- cbind(mgconcs_df, pattern_avgcur_df)
colnames(pattern_avgcur_df) <- c('Mg.Conc', 'Avg.Cur')

cols <- c('Nonpattern'='#f04546','Pattern'='#3591d1')
ggplot() + 
  ggtitle('Average Current vs. Magnesium Concentration') + 
  geom_point(data=nonpattern_avgcur_df, aes(x=Mg.Conc, y=Avg.Cur, color='Nonpattern'), shape=15, size=2.5) + 
  geom_line(data=nonpattern_avgcur_df, aes(x=Mg.Conc, y=Avg.Cur, color='Nonpattern'), size=0.7) + 
  geom_point(data=pattern_avgcur_df, aes(x=Mg.Conc, y=Avg.Cur, color='Pattern'), shape=16, size=2.5) + 
  geom_line(data=pattern_avgcur_df, aes(x=Mg.Conc, y=Avg.Cur, color='Pattern'), size=0.7) +
  labs(x='Extracellular Mg Concentration (mM)', y='Average Current (nA)') +
  scale_color_manual(name='Cell Type', values=cols) +
  scale_x_continuous(breaks=seq(from=0, to=2.0, by=0.2)) + 
  scale_y_continuous(breaks=seq(from=-0.02, to=0, by=0.0025)) + 
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5, size=14, margin=margin(t=0, r=0, b=15, l=0))) +
  theme(legend.title=element_text(size=12), legend.text=element_text(size=10)) + 
  theme(axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0))) + 
  theme(axis.title.y=element_text(margin=margin(t=0, r=15, b=0, l=0))) + 
  theme(panel.grid.minor=element_blank())

# MEMBRANE POTENTIAL VS. TIME


# DEAD CELL COUNT VS. MG CONCENTRATION
nonpattern_count <- read.table('pyresults/mgf_nonpattern_count.dat', header=TRUE)
nonpattern_count_df <- as.data.frame(nonpattern_count)
colnames(nonpattern_count_df) <- c('Idx', 'Dead.Count')
nonpattern_count_df$Mg.Conc = c(0.0, 0.5, 1.0, 2.0, 4.0)
pattern_count <- read.table('pyresults/mgf_pattern_count.dat', header=TRUE)
pattern_count_df <- as.data.frame(pattern_count)
colnames(pattern_count_df) <- c('Idx', 'Dead.Count')
pattern_count_df$Mg.Conc = c(0.0, 0.5, 1.0, 2.0, 4.0)

# AVG ALIVE TIME VS. MG CONCENTRATION
nonpattern_avgtime <- read.table('pyresults/mgf_nonpattern_avgtime.dat', header=TRUE)
nonpattern_avgtime_df <- as.data.frame(nonpattern_avgtime)
colnames(nonpattern_avgtime_df) <- c('Idx', 'Avg.Time')
nonpattern_avgtime_df$Mg.Conc = c(0.0, 0.5, 1.0, 2.0, 4.0)
pattern_avgtime <- read.table('pyresults/mgf_pattern_avgtime.dat', header=TRUE)
pattern_avgtime_df <- as.data.frame(pattern_avgtime)
colnames(pattern_avgtime_df) <- c('Idx', 'Avg.Time')
pattern_avgtime_df$Mg.Conc = c(0.0, 0.5, 1.0, 2.0, 4.0)

cols <- c('Nonpattern'='#f04546','Pattern'='#3591d1')
ggplot() + 
  ggtitle('Average Survival Time vs. Magnesium Concentration') + 
  geom_point(data=nonpattern_avgtime_df, aes(x=Idx, y=Avg.Time, color='Nonpattern'), shape=15, size=2.5) + 
  geom_line(data=nonpattern_avgtime_df, aes(x=Idx, y=Avg.Time, color='Nonpattern'), size=0.7) + 
  geom_point(data=pattern_avgtime_df, aes(x=Idx, y=Avg.Time, color='Pattern'), shape=16, size=2.5) + 
  geom_line(data=pattern_avgtime_df, aes(x=Idx, y=Avg.Time, color='Pattern'), size=0.7) +
  labs(x='Extracellular Mg Concentration (mM)', y='Average Survival Time (ms)') +
  scale_color_manual(name='Cell Type', values=cols) +
  # scale_x_discrete(labels=c('0.0', '0.5', '1.0', '2.0', '4.0')) +
  scale_y_continuous(breaks=seq(from=400, to=1700, by=200)) + 
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5, size=14, margin=margin(t=0, r=0, b=15, l=0))) +
  theme(legend.title=element_text(size=12), legend.text=element_text(size=10)) + 
  theme(axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0))) + 
  theme(axis.title.y=element_text(margin=margin(t=0, r=15, b=0, l=0))) + 
  theme(panel.grid.minor=element_blank())

