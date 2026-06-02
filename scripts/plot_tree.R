suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
output_image <- args[2] 

tree <- read.tree(tree_file)

num_tips <- length(tree$tip.label)

plot_height <- max(10, num_tips * 0.15) 

p <- ggtree(tree) + 
 
  geom_tiplab(size=1.5, color="#333333", offset=0.01) + 

  geom_nodelab(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 70), 
               size=1.5, color="blue", vjust=-0.5) +    
  theme_tree2() + 

  xlim(0, 3.5) 

ggsave(output_image, plot=p, width=12, height=plot_height, dpi=150, bg="white", limitsize=FALSE)

pdf_output <- sub("\\.png$", ".pdf", output_image)
ggsave(pdf_output, plot=p, width=12, height=plot_height, bg="white", limitsize=FALSE)