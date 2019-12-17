<h1 text-align="center">R code for box plot representation of Editing Index distributions</h1>

Data table should be tab delimited and formatted as <a href="https://github.com/BioinfoUNIBA/QEdit/blob/master/Figures_R_code/volcano_plot_recoding_cerebVSarttib2.txt"> AEI_GTEx_selected.txt</a>

<pre>library (ggplot2)
data <- read.table("AEI_GTEx_selected.txt", header=TRUE, check.names = FALSE, sep = "\t")
cols <- c("Brain - Frontal Cortex (BA9)" = "royalblue2", "Brain - Hippocampus" = "steelblue3", "Brain - Cerebellum" = "deepskyblue4", "Brain - Spinal cord (cervical c-1)" = "steelblue1", "Muscle - Skeletal" = "lightsalmon3", "Brain - Hypothalamus" = "slateblue1", "Lung" = "darkslategray2", "Artery - Tibial" = "lightcoral", "Artery - Aorta" = "indianred2", "Brain - Amygdala" = "lightblue3")
png("AEI_GTEx_selected.png", w=25, h=16, res = 300, units = 'in', pointsize=25)
ggplot(data, aes(x=Condition, y=Value, fill=Condition)) + geom_boxplot(alpha=0.7) + xlab(" ") + theme(axis.title=element_text(size="30"), axis.text.y=element_text(size="20"), legend.text=element_text(size="25"), legend.title=element_blank(), strip.text.x=element_text(size="25"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position = "none", plot.title=element_text(size="30", hjust = 0.5)) + scale_fill_manual(values=cols)+ ylab("AEI (%)") + ggtitle("Alu Editing Index")
dev.off()
</pre>
<h1 text-align="center">R code for Volcano plot representation of differentially edited sites</h1>

Data table should be tab delimited and formatted as <a href="https://github.com/BioinfoUNIBA/QEdit/blob/master/Figures_R_code/volcano_plot_recoding_cerebVSarttib2.txt"> volcano_plot_recoding_cerebVSarttib.txt</a>

In particular, the "Delta" value is the difference between the average level of editing between the two compared groups and the "log" value is the -log(10) of Benjamini Hochberg corrected pvalue calculated by Mann-Whitney two tailed test.

<pre>library(ggplot2)
library(ggrepel)
data <- read.table("volcano_plot_recoding_cerebVSarttib.txt", header=TRUE, check.names = FALSE)
names <- subset(data, log > 1.3 & Delta > 0.05 | log > 1.3 & Delta < -0.05)
cols <- c("UP" = "salmon", "DOWN" = "dodgerblue", "NS." = "darkgrey")
title <- expression("Cerebellum"~italic("vs.")~"Artery Tibial")
png("volcano_plot.png", w=24, h=20, res = 300, units = 'in', pointsize=25)
ggplot(data, aes(x=Delta, y=log, color=color)) + geom_point(size = 5,alpha = 1) + geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.1, colour="#990000", linetype="dashed") + ggtitle(label =title, subtitle = "n. sites = 85") + scale_colour_manual(values = cols)+ ylab(expression(paste("",-log[10],"(", italic("p"),"adj" ,")",sep="")))+ xlab(expression(Delta~editing)) + scale_x_continuous(limits=c(-0.8,0.8))+ theme(axis.title=element_text(size="30"), axis.text.x=element_text(size="20"), axis.text.y=element_text(size="20"), legend.text=element_text(size="20"), legend.title=element_blank(),plot.subtitle = element_text(size="30"),plot.title = element_text(size="35")) + geom_text_repel(data=names ,aes(label=names$Site),size=8,show.legend=F)
dev.off()
</pre>

<h1 text-align="center">R code for heatmap showing editing levels at recoding sites</h1>

Data table, containing editing levels (%) at recoding poisitions covered by at least 10 reads, should be tab delimited and formatted as <a href="https://github.com/BioinfoUNIBA/QEdit/blob/master/Figures_R_code/filtered_full_recoding_table2.txt"> filtered_full_recoding_table.txt</a>

<pre>library(ComplexHeatmap)
library(circlize)
library('RColorBrewer')
data <- read.table("filtered_full_recoding_table.txt", header=TRUE, row.names=1,check.names = FALSE, sep = "\t")
mat <- as.matrix(data)
cols <- colorRampPalette(brewer.pal(9, 'Blues'))(100)
ha_column = HeatmapAnnotation(df = data.frame(tissue = c(rep("Artery - Aorta", 14), rep("Artery - Tibial", 14), rep("Brain - Amygdala",14), rep("Brain - Cerebellum", 14), rep("Brain - Frontal Cortex (BA9)", 14), rep("Brain - Hippocampus", 14), rep("Brain - Hypothalamus", 14), rep("Brain - Spinal cord (cervical c-1)", 14), rep("Lung", 14), rep("Muscle - Skeletal", 14))), col = list(tissue = c("Brain - Frontal Cortex (BA9)" = "royalblue2", "Brain - Hippocampus" = "steelblue3", "Brain - Cerebellum" = "deepskyblue4", "Brain - Spinal cord (cervical c-1)" = "steelblue1", "Muscle - Skeletal" = "lightsalmon3", "Brain - Hypothalamus" = "slateblue1", "Lung" = "darkslategray2", "Artery - Tibial" = "lightcoral", "Artery - Aorta" = "indianred2", "Brain - Amygdala" = "lightblue3")), annotation_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 15) , title = " ", grid_height = unit(2, "cm"), grid_width = unit(1.5, "cm"), title_position = "leftcenter-rot"))
png("recoding_editing_heatmap.png", w=28, h=22, res = 300, units = 'in', pointsize=25)
Heatmap(mat, cluster_rows = TRUE, top_annotation =ha_column, show_row_dend = FALSE, cluster_columns = FALSE, col = cols, na_col = "gray95", show_column_names = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20) ,col_fun = cols, title = "editing level (%)", at = c(0, 50, 100), legend_height = unit(10, "cm"), grid_width = unit(1.5, "cm"), title_position = "leftcenter-rot"))
dev.off()
</pre>





