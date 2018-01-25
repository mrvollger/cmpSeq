#!/usr/bin/env Rscript
library(ggplot2)
require(scales)
library(grid)
path = "~/Desktop/data/genomeWide/Mitchell_CHM1/LocalAssemblies/clones/cmps/group.0_collapse.000145F.235600.249399_id.1/perID/slidePerID.tsv"
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
#parser$add_argument("-n", "--add_numbers", action="store_true", default=FALSE, help="Print line number at the beginning of each line [default]")
parser$add_argument("input", nargs='?', help="File to be displayed", default=path)
parser$add_argument("out",   nargs='?', help="output png file", default="~/Desktop/slideAln.svg")
args <- parser$parse_args()

path = args$input
mydata <- read.table(path, header=TRUE)
dim(mydata)
selfAlign = (as.character(mydata$qname) == as.character(mydata$rname))
mydata = mydata[!selfAlign, ]
dim(mydata)
#
# ggplot stuff 
#
ymin=90
textSize = 54
theme_set(theme_minimal(base_size = textSize))
th = theme(axis.title.y = element_text(margin = ggplot2::margin(t = 50, r = 50, b = 50, l = 50)))
size = 2
alpha = 0.1
h=2000
w=3500

#
# this is to make text boxes
#
qpre = ""
rpre=""
changes = c()
for(idx in 1:nrow(mydata)){
  row = mydata[idx, ]
  if(row$qname != qpre | row$rname != rpre){
    changes = c(changes, idx-1, idx)
  }
  qpre = row$qname
  rpre = row$rname
}
changes = changes[2:length(changes)]
changes = c(changes, nrow(mydata))
textdf = mydata[changes,]
textdf

#mydata = mydata[mydata$rname == "AC270130.1", ]
p1 = ggplot(mydata, aes(Reference_Coordinate, perID_by_matches, color=qname, fill=rname) )  +
  geom_point(size = size, stroke = size)  + 
  geom_area(alpha = alpha) +
  geom_text(data=textdf, aes(Reference_Coordinate, perID_by_matches, label=Query_Coordinate), color="black", size=textSize/3) +
  scale_x_continuous(labels = comma) + th + facet_grid(.~rname, scales = "free_x") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  ylab("% Identity") + coord_cartesian(ylim = c(ymin, 100)) + theme(legend.position = "none") 

p2 = ggplot(mydata, aes(Reference_Coordinate, Deletions, color=qname, fill=rname))  + 
  geom_area(alpha = alpha)+
  geom_point(size = size, stroke = size)  + 
  scale_x_continuous(labels = comma) + th + facet_grid(.~rname, scales = "free_x") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  ylab("# Deletions") + theme(legend.position = "none")


p3 = ggplot(mydata, aes(Reference_Coordinate, Insertions, color=qname, fill=rname))  +  
  geom_area(alpha = alpha)+
  geom_point(size = size, stroke = size)  + 
  scale_x_continuous(labels = comma) + th + facet_grid(.~rname, scales = "free_x") +
  ylab("# Insertions") + xlab("Alignment position (bp)")  + theme(legend.position = "bottom")



grid.newpage()
svg(args$out, height = h, width = w)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))
## Stop writing to the  file
dev.off()



