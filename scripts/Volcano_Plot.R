library(plotly)
my_file <- read.csv("condition_results.csv", header = TRUE)
p <- ggplot(data = my_file[my_file$padj >= 0, ], aes(x = log2FoldChange, y = -log10(pvalue)))+
  geom_point(aes(color = my_file$padj < 0.05), size = 2.5, alpha = 0.4)+
  guides(color=FALSE)+
  xlab("log2 Fold Change")+
  ylab("-log10 (p-value)")+
  ggtitle("Volcano Plot showing Up and Down regulated genes")+
  theme_minimal()


x <- list(
  title = "log2 Fold Change"
)
y <- list(
  title = "-log10 (p-value)"
)
plot_ly(data = my_file[my_file$padj >= 0, ], x=~log2FoldChange, y=~-log10(pvalue), color = my_file$padj < 0.05, alpha = 0.7) %>%
  add_markers(text=my_file$X) %>%
  layout(
    title = "Upregulated and Downregulated transcripts in a Volcano Plot",xaxis = x, yaxis = y,showlegend = FALSE)
