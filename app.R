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
                   checkboxInput(paste0(name,'_deletions'), "Deletions", value = TRUE),
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
                               tabPanel("Variant table",DTOutput(paste0(name,'_table'))),
                               # Plot View using ggplot
                               tabPanel("Variant plot", plotOutput(paste0(name,'_plot'))
                    )
                  )
                )
                          
              )
             )
           )
  )
}

# Define UI for application that draws a histogram
ui <- navbarPage("MPRA-Data",
                 tabPanel("About",
                          h1("MPRA data access portal"),
                          hr(),
                          htmlOutput("about")
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
                            mainPanel(
                                     selectInput("download_reference_all", "Genome release:",c("GRCh37","GRCh38")),
                                     selectInput("download_promoter", "Promoters:", promoters, multiple=TRUE, selectize=TRUE),
                                     selectInput("download_enhancer", "Enhancers:", enhancers, multiple=TRUE, selectize=TRUE),
                                     selectInput("download_format_all", "Format:",c(".tsv",".csv")),
                                     downloadButton("downloadData_selected", "Download Selected Elements"),
                                     downloadButton("downloadData_all", "Download All Elements")
                            )
                          )
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
      output[[paste0(name,"_table")]] <- renderDT(data, rownames = FALSE, 
                                                  colnames = c('Chromosome', 'Posistion', 'Ref', 'Alt', 'Tags','DNA','RNA','Value','P-Value'))
      
      # plot the data. here the helper   modify.filterdata is used to filter the data (from plots.R
      # and the function getPlot from the same source to plot the ggplot
      plotData <- modify.filterdata(data_general, barcodes,threshold,deletions,range)
      if (plotData %>% nrow > 0) {
        output[[paste0(name,"_plot")]] <- renderPlot(getPlot(plotData,name,release))
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
  output$about <- renderUI({
    HTML('We selected 21 regulatory elements, including 20 commonly studied, 
         disease-relevant promoter and enhancer sequences from the literature, and one ultraconserved enhancer (UC88). 
         For the former, we focused primarily on regulatory sequences in which specific mutations are known to cause disease, 
         both for their clinical relevance and to provide for positive control variants. 
         Selected elements were limited to 600 base pairs (bp) for technical reasons related to the mapping of variants to 
         barcodes by subassembly. In addition, we selected only sequences where cell line-based reporter assays were 
         previously established. <br>
        <h2>Promoter</h2>
<style type="text/css">
	table.tableizer-table {
         font-size: 12px;
         border: 1px solid #CCC; 
         font-family: Arial, Helvetica, sans-serif;
  } 
         .tableizer-table td {
         padding: 4px;
         margin: 3px;
         border: 1px solid #CCC;
         }
         .tableizer-table th {
         background-color: #104E8B; 
         color: #FFF;
         font-weight: bold;
         }
         </style>
         <table class="tableizer-table">
<thead><tr class="tableizer-firstrow"><th>Name</th><th>Genomic coordinates (GRCh37/GRCh38)</th><th>Transcript</th><th>Associated Phenotype</th><th>Luciferase vector</th><th>MPRA vector</th><th>Cell line</th><th>Transf. time (hr)</th><th>Fold Ch.</th><th>Construct size (bp)</th></tr></thead><tbody>
 <tr><td>F9</td><td>X:138,612,622-138,612,924</td><td>NM_000133.3</td><td>Hemophilia B</td><td>pGL4.11b</td><td>pGL4.11c</td><td>HepG2</td><td>24</td><td>2.6</td><td>303</td></tr>
 <tr><td>&nbsp;</td><td>X:139,530,463-139,530,765</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>FOXE1</td><td>9:100,615,537-100,616,136</td><td>NM_004473.3</td><td>Thyroid cancer</td><td>pGL4.11b</td><td>pGL4.11c</td><td>HeLa</td><td>24</td><td>6.6</td><td>600</td></tr>
 <tr><td>&nbsp;</td><td>9:97,853,255-97,853,854</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>GP1BB</td><td>22:19,710,789-19,711,173</td><td>NM_000407.4</td><td>Bernard-Soulier Syndrome</td><td>pGL4.11b</td><td>pGL4.11c</td><td>HEL 92.1.7</td><td>24</td><td>22.1</td><td>385</td></tr>
 <tr><td>&nbsp;</td><td>22:19,723,266-19,723,650</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>HBB</td><td>11:5,248,252-5,248,438</td><td>NM_000518.4</td><td>Thalassemia</td><td>pGL4.11b</td><td>pGL4.11c</td><td>HEL 92.1.7</td><td>24</td><td>14.3</td><td>187</td></tr>
 <tr><td>&nbsp;</td><td>11:5,227,022-5,227,208</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>HBG1</td><td>11:5,271,035-5,271,308</td><td>NM_000559.2</td><td>Hereditary persistence of fetal hemoglobin</td><td>pGL4.11b</td><td>pGL4.11c</td><td>HEL 92.1.7</td><td>24</td><td>118.1</td><td>274</td></tr>
 <tr><td>&nbsp;</td><td>11:5,249,805-5,250,078</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>HNF4A (P2)</td><td>20:42,984,160-42,984,444</td><td>NM_175914.4</td><td>Maturity-onset diabetes of the young (MODY)</td><td>pGL4.11b</td><td>pGL4.11c</td><td>HEK293T</td><td>24</td><td>2.8</td><td>285</td></tr>
 <tr><td>&nbsp;</td><td>20:44,355,520-44,355,804</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>LDLR</td><td>19:11,199,907-11,200,224</td><td>NM_000527.4</td><td>Familial hypercholesterolemia</td><td>pGL4.11b</td><td>pGL4.11b</td><td>HepG2</td><td>24</td><td>110.7</td><td>318</td></tr>
 <tr><td>&nbsp;</td><td>19:11,089,231-11,089,548 </td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>MSMB</td><td>10:51,548,988-51,549,578</td><td>NM_002443.3</td><td>Prostate cancer</td><td>pGL4.11b</td><td>pGL4.11c</td><td>HEK293T</td><td>24</td><td>8.4</td><td>593</td></tr>
 <tr><td>&nbsp;</td><td>10:46,046,244-46,046,834</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>PKLR</td><td>1:155,271,186-155,271,655</td><td>NM_000298.5</td><td>Pyruvate kinase deficiency</td><td>pGL4.11b</td><td>pGL4.11c</td><td>K562</td><td>48</td><td>29.4</td><td>470</td></tr>
 <tr><td>&nbsp;</td><td>1:155,301,395-155,301,864</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>TERT</td><td>5:1,295,104-1,295,362</td><td>NM_198253.2</td><td>Various types of cancer</td><td>pGL4.11b</td><td>pGL4.11b</td><td>HEK293T, GBM</td><td>24</td><td>231.8</td><td>259</td></tr>
 <tr><td>&nbsp;</td><td>5:1,294,989-1,295,247</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td></td></tr>
</tbody></table>
        <h2>Enhancer</h2>
        <table class="tableizer-table">
<thead><tr class="tableizer-firstrow"><th>Name</th><th>Genomic coordinates (GRCh37/GRCh38)</th><th>Associated Phenotype</th><th>Luciferase vector</th><th>MPRA vector</th><th>Cell line</th><th>Transf. time (hr)</th><th>Fold Ch.</th><th>Construct (bp)</th></tr></thead><tbody>
 <tr><td>BCL11A</td><td>2:60,722,075-60,722,674</td><td>Sickle cell disease </td><td>pGL4.23</td><td>pGL4.23d</td><td>HEL 92.1.7</td><td>24</td><td>2.5</td><td>600</td></tr>
 <tr><td>58</td><td>2:60,494,940-60,495,539</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>IRF4</td><td>6:396,143-396,593</td><td>Human pigmentation</td><td>pGL4.23</td><td>pGL4.23d</td><td>SK-MEL-28</td><td>24</td><td>44.5</td><td>451</td></tr>
 <tr><td>&nbsp;</td><td>6:396,143-396,593 </td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>IRF6</td><td>1:209,989,135-209,989,735 </td><td>Cleft lip</td><td>pGL4.23</td><td>pGL4.23c</td><td>HaCaT</td><td>24</td><td>17</td><td>600</td></tr>
 <tr><td>&nbsp;</td><td>1:209,815,790-209,816,390</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>MYC (rs 6983267)</td><td>8:128,413,074-128,413,673</td><td>Various types of cancer </td><td>pGL4.23</td><td>pGL4.23c</td><td>HEK293T</td><td>32, 20nM LiCl added after 24hr</td><td>0.8</td><td>600</td></tr>
 <tr><td>&nbsp;</td><td>8:127,400,829-127,401,428</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>MYC (rs 11986220)</td><td>8:128,531,515-128,531,977</td><td>Various types of cancer </td><td>pGL4.23</td><td>pGL4.23d</td><td>LNCaP + 100nM DHT</td><td>24</td><td>5.5</td><td>464</td></tr>
 <tr><td>&nbsp;</td><td>8:127,519,270-127,519,732</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>RET </td><td>10:43,581,927-43,582,526</td><td>Hirschsprung</td><td>pGL3</td><td>pGL3c</td><td>Neuro-2a</td><td>24</td><td>2</td><td>600</td></tr>
 <tr><td>&nbsp;</td><td>10:43,086,479-43,087,078 </td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>SORT1</td><td>1:109,817,274-109,817,873</td><td>Plasma low-density lipoprotein cholesterol & myocardial infraction</td><td>pGL4.23</td><td>pGL4.23</td><td>HepG2</td><td>24</td><td>235.3</td><td>600</td></tr>
 <tr><td>&nbsp;</td><td>1:109,274,652-109,275,251</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>TCF7L2 </td><td>10:114,757,999-114,758,598</td><td>Type 2 diabetes</td><td>pGL4.23</td><td>pGL4.23d</td><td>MIN6</td><td>24</td><td>9</td><td>600</td></tr>
 <tr><td>&nbsp;</td><td>10:112,998,240-112,998,839</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>UC88</td><td>2:162,094,919-162,095,508</td><td>-</td><td>pGL4.23</td><td>pGL4.23c</td><td>Neuro-2a</td><td>24</td><td>9.3</td><td>590</td></tr>
 <tr><td>&nbsp;</td><td>2:161,238,408-161,238,997 </td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>ZFAND3</td><td>6:37,775,275-37,775,853</td><td>Type 2 diabetes</td><td>pGL4.23</td><td>pGL4.23c</td><td>MIN6</td><td>24</td><td>14.3</td><td>579</td></tr>
 <tr><td>&nbsp;</td><td>6:37,807,499-37,808,077 </td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td></tr>
 <tr><td>ZRS</td><td>7:156,583,813-156,584,297</td><td>Limb malformations</td><td>TATA-pGL4m(EV087)*</td><td>pGL4Zc</td><td>NIH/3T3</td><td>24</td><td>4.2</td><td>485</td></tr>
 <tr><td>&nbsp;</td><td>7:156,791,119-156,791,603</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td></td></tr>
</tbody></table>
         ')
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

