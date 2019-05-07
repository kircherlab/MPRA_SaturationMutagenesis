library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

standard_SatMut_region_style <- function() {
  ## styles
  size.text = 10;
  size.title = 12;
  size.line = 1.2;
  size.geom_line = 1;
  size.geom_point = 2
  standard_style <- theme_bw() + theme(plot.title = element_text(size = size.title, face="bold",hjust = 0.5),
                                       panel.grid.major = element_blank() , panel.grid.minor = element_blank(), panel.border = element_blank(),
                                       axis.text = element_text(colour = "black",size=size.text), axis.title = element_text(colour = "black",size=size.title), axis.ticks = element_line(colour = "black", size=1), axis.line.y = element_line(color="black", size = size.line), axis.line = element_line(colour = "black", size=size.line),
                                       legend.key =  element_blank(), legend.text = element_text(size=size.text),
                                        legend.position="top", legend.box.just = "left",  legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.margin = margin(0, 0, 0, 0),
                                       legend.key.size = unit(2, 'lines'), legend.title=element_text(size=size.text))+
    theme(axis.line.x = element_line(color="black", size = size.line))
}

modify.filterdata <- function(data,barcodes=10, threshold=1e-5, deletions=TRUE, range=NULL) {
  data <- data %>% select(Chrom,Pos,Ref,Alt,Barcodes, Coefficient, pValue) %>%
    filter(Barcodes >= barcodes) %>% 
    mutate(significance=if_else(pValue<threshold,"Significant", "Not significant")) %>%
    mutate(printpos=if_else(Alt=="A",as.double(Pos)-0.4,if_else(Alt=="T",as.double(Pos)-0.2, if_else(Alt=="G",as.double(Pos)+0.0,if_else(Alt=="C",as.double(Pos)+0.2,as.double(Pos)+0.4)))))
  if (!deletions) {
    data <- data %>% filter(Alt != "-")
  }
  if (!is.null(range)) {
    data <- data %>% filter(Pos >= range[1] & Pos <= range[2])
  }
  return(data)
}


defaultColours=c("A"="#0f9447","C"="#235c99","T"="#d42638","G"="#f5b328","-"="#cccccc", 
                         "Significant"="#005500","Not significant"="red")
colorblindColors=c("A"="#1B9E77","C"="#7570B3","T"="#D95F02","G"="#E6AB02","-"="#A6761D", 
                   "Significant"="#66A61E","Not significant"="#E7298A")
getPlot <- function(data,name, release, colourPalette="default") {
  colours <- defaultColours
  
  if (colourPalette == "colorblind") {
    print(colourPalette)
    colours <- colorblindColors
  }
  
  refs <- data$Ref %>% unique()
  aesRefsValues <- c(if_else("A" %in% refs,15,c()),if_else("C" %in% refs,16,c()),if_else("G" %in% refs,17,c()),if_else("T" %in% refs,18,c()))
  aesRefsShape <- c(if_else("A" %in% refs,0,c()),if_else("C" %in% refs,1,c()),if_else("G" %in% refs,2,c()),if_else("T" %in% refs,5,c()))
  aesRefsValues <- aesRefsValues[!is.na(aesRefsValues)]
  aesRefsShape <- aesRefsShape[!is.na(aesRefsShape)]
  
  alts <- data$Alt %>% unique()
  sigs <- data$significance %>% unique()
  aesSize<-c(rep(7,length(alts)), rep(3,length(sigs)))
  aesLine<-c(rep(0,length(alts)), rep(1,length(sigs)))
  aesShape<-c(if_else("A" %in% alts,15,c()),if_else("C" %in% alts,16,c()),if_else("G" %in% alts,17,c()),if_else("T" %in% alts,18,c()),if_else("-" %in% alts,19,c()), rep(32,length(sigs)))
  aesShape <- aesShape[!is.na(aesShape)]
  altBreaks<-c(as.character(alts), sigs)
  
  chr <- data$Chr %>% unique()
  data <- data %>% select(printpos,Coefficient,significance,Alt,Ref) %>% dplyr::rename(Position=printpos)
  p <- ggplot() +
    geom_segment(data = data, aes(x=Position, xend=Position,y=0,yend=Coefficient, colour=significance), size=0.3, show.legend = TRUE) +
    geom_point(data= data, aes(x=Position,y=Coefficient,colour=Alt, shape=Ref), size=1, show.legend = TRUE) +
    scale_shape_manual("",values=aesRefsValues, guide=guide_legend(override.aes = list(size=1, linetype=0, shape=aesRefsShape),nrow=2)) +
    scale_colour_manual("", values = colours, breaks=altBreaks, labels=altBreaks,
                        guide= guide_legend(override.aes = list(size=aesSize, linetype=aesLine, shape=aesShape))) +
    ggtitle(name) + labs(x = paste0("Chromosome ",chr," (",release,")"), y= "Log2 variant effect") + standard_SatMut_region_style()
  return(p)
}
