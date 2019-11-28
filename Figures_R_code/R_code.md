<h1 text-align="center">R code for Volcano plot representation of differentially edited sites</h1>

Data table should be tab delimited and formatted as <a href="https://github.com/BioinfoUNIBA/QEdit/blob/master/Figures_R_code/volcano_plot_recoding_cerebVSarttib.txt"> volcano_plot_recoding_cerebVSarttib.txt</a>

In particular, the "Delta" value is the difference between the average level of editing between the two compared groups and the "log" value is the -log(10) of Benjamini Hochberg corrected pvalue calculated by Mann-Whitney two tailed test.

<pre>library(ggplot2)
library(ggrepel)
data <- read.table("volcano_plot_recoding_cerebVSarttib.txt", header=TRUE, check.names = FALSE)
test <- subset(countdata, log > 1.3 & Delta > 0.05 | log > 1.3 & Delta < -0.05)
cols <- c("UP" = "salmon", "DOWN" = "dodgerblue", "NS." = "darkgrey")
t <- expression("The"~italic("alpha")-"male percentage")
png("volcano_plots.png", w=24, h=20, res = 300, units = 'in', pointsize=25)
ggplot(data, aes(x=Delta, y=log, color=color)) + geom_point(size = 5,alpha = 1) + geom_hline(yintercept = 1.3, colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.1, colour="#990000", linetype="dashed") + ggtitle(label =t, subtitle = "n. sites = 85") + scale_colour_manual(values = cols)+ ylab(expression(paste("",-log[10],"(", italic("p"),"adj" ,")",sep="")))+ xlab(expression(Delta~editing)) + scale_x_continuous(limits=c(-0.8,0.8))+ theme(axis.title=element_text(size="30"), axis.text.x=element_text(size="20"), axis.text.y=element_text(size="20"), legend.text=element_text(size="20"), legend.title=element_blank(),plot.subtitle = element_text(size="30"),plot.title = element_text(size="35")) + geom_text_repel(data=test ,aes(label=test$Site),size=8,show.legend=F) + guides(color = guide_legend(reverse = TRUE))
dev.off()
</pre>
