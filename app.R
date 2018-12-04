#
# This is a Shiny web application to access the MPRA saturation mutagenisis data.
# 
# It should be hosted on https://mpra.gs.washington.edu/satMutagenesis and has teh following functionality
#
# 1. Inform about the project/data and the people/labs behind it (About website)
# 2. Browse the different elements (separated by Enhancers and Promoters)
# 2.1 Filter the variants according to Reference genome, p-value, number of Tags, Region (Sidebarpanel)
# 2.2 Show the variants in a table, filteded on the input settings (Table View)
# 2.3 Plot the Variants of the region. use the barcode, and p-value as significance theshhold.


# Required libraries
library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(plotly)
# import plots.R for visualisations
source("plots.R")


# load the elements seperated by enhacners and promoters
loadElementList <- function(path) {
  elements <- read.delim(path, header=FALSE, quote="") 
  colnames(elements) <- c("Element")
  elements <- elements %>% arrange(Element)
  return(elements)
}

enhancers <- loadElementList("data/enhancers.tsv") 
promoters <- loadElementList("data/promoters.tsv")

# Combine the elements (for download)
elements <- bind_rows(enhancers,promoters)  %>% arrange(Element)

# Read the variant values!
data_table <- read.delim("data/elements.tsv.gz", quote="") %>%  na.omit() %>% arrange(Element,Chrom,Pos,Ref,Alt) %>% mutate(Coefficient=round(Coefficient,2),pValue=round(pValue,5))

# Separate the Varants by elements. This is to speedup the process later because this step has to be done only once.
data_tables <- list()
for (name in elements$Element) {
  data_tables[[name]] <- data_table %>% filter(Element==name)
}


# The list of different experiments is dynamic. 
# This function creates different tabs for an element with asidebar panel for filtering
# and a main panel with Tabset table view and Plot. It is debendend on the name element
experimentPage <- function(name){
  tabPanel(name,
           fluidPage(
             
             # Application title
             headerPanel(name),
             mainPanel(
               sidebarLayout(
                 # User inputs
                 sidebarPanel(
                   # Genome release
                   selectInput(paste0(name,'_reference'), "Genome Release:",c("GRCh37","GRCh38")),
                   # Filtering
                   h3("Filter"),
                   checkboxInput(paste0(name,'_deletions'), "With 1 bp deletions", value = TRUE),
                   # Seperate panel to select significance level (for plots)
                   wellPanel(
                    h4("Significance level"),
                   numericInput(paste0(name,'_pvalue'), "P-Value:",1e-5,min=1e-5,max=1.0),
                   numericInput(paste0(name,'_barcodes'), "Minimum Tags:",10,min=0,step=1)
                   ),
                   sliderInput(paste0(name,"_region"), "Region", min=0, max=100, value=c(0,100), step = 1)#,
                   #actionButton(paste0(name,"_filter"),"Filter")
                 ),
                 # Main panel for tab view and plot view
                 mainPanel(
                   tabsetPanel(type = "tabs",
                               # table View using DT table
                               tabPanel("Variant table",
                                        fluidRow(
                                          column(12,
                                            DT::DTOutput(paste0(name,'_table'))
                                          ),
                                          column(12,id="firefox_warning",
                                                 tags$script(HTML('
                                            //Javascript
                                            var FIREFOX = /Firefox/i.test(navigator.userAgent);
                                            
                                            if (FIREFOX) {
                                              document.getElementById("firefox_warning").style.color="red";
                                            }
                                            ')),
                                                 "We have some problems visualizing the variant table in Firefox or Edge correctly. 
                                                 If you cannot see the variant table please use the Chrome browser."
                                           )
                                        )
                              ),
                               # Plot View using ggplot
                               tabPanel("Variant plot", plotlyOutput(paste0(name,'_plot'))
                    )
                  )
                )
                          
              )
             )
           )
  )
}

# Define UI for application that draws a histogram
ui <- fluidRow(id="canvas",
        column(12,
               navbarPage("MPRA-Data",
                 theme = "mpra.css",
                 tabPanel("About",
                          h1("MPRA data access portal"),
                          hr(),
                          includeMarkdown("mrkdown/about.md")
                          ),
                 tabPanel("Promoter",
                          tabsetPanel(type = "tabs",id="promoterNavigation")
                 ),
                 tabPanel("Enhancer",
                          tabsetPanel(type = "tabs",id="enhancerNavigation")),
                 tabPanel("Download", 
                          fluidPage(title="Download elements",
                            h2("Download elements"),
                    
                            # Application title
                            hr(),
                            fluidRow(
                              column(4,
                                     id="download_panel",
                                     selectInput("download_reference_all", "Genome release:",c("GRCh37","GRCh38")),
                                     selectInput("download_promoter", "Promoters:", promoters, multiple=TRUE, selectize=TRUE),
                                     selectInput("download_enhancer", "Enhancers:", enhancers, multiple=TRUE, selectize=TRUE),
                                     selectInput("download_format_all", "Format:",c(".tsv",".csv")),
                                     downloadButton("downloadData_selected", "Download Selected Elements"),
                                     downloadButton("downloadData_all", "Download All Elements")
                                     
                              ),
                              column(8,
                                     id="file_format",
                                     includeMarkdown("mrkdown/file_format.md")
                                     )
                            )
                          )
                 )
               )
         ),
        column(12,
               id="footer",
                HTML('<a href="https://www.washington.edu/online/terms">Terms and Conditions</a> and the <a href="https://www.washington.edu/online/privacy">Online Privacy Statement</a> of the University of Washington apply.')
        )
                 
)


# function to render the table/plot for a specific selected element (by name)
# is dependend on the release, barcodes, p-value threshsold, true/false deletions, 
# range of the plot and needs output and session element
renderElement <- function(name,release,barcodes, threshold, deletions, range, output,session) {
  # check first if elements are not null
  if(!is.null(release) & !is.null(barcodes) & !is.null(threshold)) {
    # values cannot be NA (possible if out of range is selected, like <1 barcoded)
    if (!is.na(threshold) & !is.na(barcodes) & threshold > 0 & barcodes > 0) {
      # get name depeneded data and filter on release
      data_general <- data_tables[[name]] %>% filter(Release==release)
      # filter on barcodes as well as pvalue for the table, remove not needed rows
      data <- data_general %>% filter(Barcodes >= barcodes) %>% filter(pValue < threshold) %>%  
        select(-Element,-Release)
      # remove deletions if unselected
      if (!deletions) {
        data <- data %>% filter(Alt != "-")
      }
      # filter down to range if selected
      if (!is.null(range)) {
        data <- data %>% filter(Pos >= range[1] & Pos <= range[2])
      }
      # render the table with the data
      output[[paste0(name,"_table")]] <- DT::renderDT(data, rownames = FALSE, filter="top",
                                                  options = list(lengthMenu = list(c(25,50, 100, 500, 1000, -1), list('25', '50', '100','500','1000', 'All'))),
                                                  colnames = c('Chromosome', 'Position', 'Ref', 'Alt', 'Tags','DNA','RNA','Value','P-Value'))
      
      # plot the data. here the helper   modify.filterdata is used to filter the data (from plots.R
      # and the function getPlot from the same source to plot the ggplot
      plotData <- modify.filterdata(data_general, barcodes,threshold,deletions,range)
      if (plotData %>% nrow > 0) {
        output[[paste0(name,"_plot")]] <- renderPlotly({
          p <- getPlot(plotData,name,release)
            ggplotly(p) %>% 
              layout(autosize=TRUE) %>% layout(xaxis=list(fixedrange=TRUE)) %>% layout(yaxis=list(fixedrange=TRUE)) %>% 
              config(displayModeBar = F)
          })
      }
    }
  }
}


# update the slider. We have to set the min/max of a slider depending on the element and on the selected genome release
updateRegionSlider <- function(session, id, name, release) {
  # always check if release is not NULL
  if(!is.null(release)) {
    # filter data on release and get min/max
    data <- data_tables[[name]] %>% filter(Release==release)
    min <- min(data$Pos)
    max <- max(data$Pos)
    
    updateSliderInput(session, id, value = c(min,max), min=min, max=max, step=1)
  }
}

# Define server logic required for this website
server <- function(input, output, session) {
  
  # need to store the selected genome release here to update min/max if needed
  session$userData$actualSelectedRelease <- list()
  
  # initialize the tabs with the elements. promoter and enhancer wise. redenr element and set slider first
  for (name in promoters$Element){
    session$userData$actualSelectedRelease[[name]] <- "GRCh37"
    isolate({
      appendTab("promoterNavigation", experimentPage(name))
      updateRegionSlider(session, paste0(name,"_region"), name, "GRCh37")
      renderElement(name,"GRCh37",10,1e-5,TRUE,NULL,output,session)
    })
  }
  for (name in enhancers$Element){
    session$userData$actualSelectedRelease[[name]] <- "GRCh37"
    isolate({
      appendTab("enhancerNavigation", experimentPage(name))
      updateRegionSlider(session, paste0(name,"_region"), name, "GRCh37")
      renderElement(name,"GRCh37",10,1e-5,TRUE,NULL,output,session)
    })
  }
  
  
  # observer to observe filter and other selections
  # FIXME: Right now it observes everything and acts on that.
  # It would be better to act on the direct klick on a specifi slider 
  # or other selection element in the filter settings not on all.
  # But the names are dynamically. So it will be difficult
  
  observe({
      for (name in elements$Element){
        
        # get the inputs
        release <- input[[paste0(name,"_reference")]]
        barcodes <- input[[paste0(name,"_barcodes")]]
        pValue <- input[[paste0(name,"_pvalue")]]
        deletions <- input[[paste0(name,"_deletions")]]
        
        
        # update the slider if the release changed
        if (!is.null(release) && session$userData$actualSelectedRelease[[name]] != release) {
          session$userData$actualSelectedRelease[[name]] <- release
          updateRegionSlider(session, paste0(name,"_region"), name, release)
        }
        # update data table and plot
        range <- input[[paste0(name,"_region")]]
        isolate({
          renderElement(name,release,barcodes,pValue,deletions,range,output,session)
        })
      }
    })
  
  ## Download button. Observe this specific event
  # selected
  output$downloadData_selected <- downloadHandler(
    # function for filename
    filename = function() {
      selectedElements <- c(input$download_promoter,input$download_enhancer)
      return(paste0(paste(input$download_reference_all,paste(selectedElements,sep="_",collapse = "_"), sep="_"), input$download_format_all))
    },
    content = function(file) {
      release=input$download_reference_all
      selectedElements <- c(input$download_promoter,input$download_enhancer)
      data <- data_table %>% filter(Release==release) %>% filter(Element %in% selectedElements) %>% select(-Release)
      colnames(data) <- c('Chromosome', 'Position', 'Ref', 'Alt', 'Tags','DNA','RNA','Value','P-Value','Element')
      sep="\t"
      if (input$download_format_all == ".csv") {
        sep=","
      }
      # download table
      write.table(data, file, row.names = FALSE, sep=sep, quote=FALSE)
    }
  )
  # all
  output$downloadData_all <- downloadHandler(
    # function for filename
    filename = function() {
      return(paste0(paste(input$download_reference_all,"ALL", sep="_"), input$download_format_all))
    },
    content = function(file) {
      release=input$download_reference_all
      data <- data_table %>% filter(Release==release) %>% select(-Release)
      colnames(data) <- c('Chromosome', 'Position', 'Ref', 'Alt', 'Tags','DNA','RNA','Value','P-Value','Element')
      sep="\t"
      if (input$download_format_all == ".csv") {
        sep=","
      }
      # download table
      write.table(data, file, row.names = FALSE, sep=sep, quote=FALSE)
    }
  )
  
  ## ABOUT page right now simple HTML
#  output$about <- renderUI({
#    includeMarkdown("mrkdown/about.md")
#  })
}

# Run the application 
shinyApp(ui = ui, server = server)

