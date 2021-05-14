# To run the Shiny app, from the terminal in the app folder:
# R -e "shiny::runApp('./app/app.R', launch.browser = TRUE)"

library(shiny)
library(shinyWidgets)
library(tidyverse)
library(DT)
library(shiny.router)

# Paths
figsFolder <- "../BBA/figs"
#homoSapiensDB <- "../dbTF/Homo_sapiens.tsv"
homoSapiensDB <- "Homo_sapiens.tsv"
tfIDFTables <- "../BBA/filtered"

# Read data
genesTable <- read.table(
  file = homoSapiensDB,
  header = TRUE,
  sep = "\t",
  na.strings = c("", "NA")
)

# Esternal links to other databases
resources_links <- c(
  "http://jaspar.genereg.net/matrix/",
  "https://hocomoco11.autosome.ru/motif/",
  "http://chip-atlas.org/view?id=",
  "#",
  "#",
  "https://www.ncbi.nlm.nih.gov/sra/",
  "https://www.ncbi.nlm.nih.gov/sra/",
  "http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/PWMs/Files/",
  "http://thebrain.bwh.harvard.edu/uniprobe/details34.php?id=",
  "https://www.ncbi.nlm.nih.gov/sra/",
  "http://remap.univ-amu.fr/target_page/"
)

names(resources_links) <- colnames(genesTable)[11:21]

# Select columns that will be displayed in home page table
genesTableSelectedCols <- genesTable %>% select(Gene.Name, Gene.ID, Species, UniProt.Accession, UniProt.Entry)

homePage <- div(
  class = "container-fluid wrapper-content-margin",
  DT::dataTableOutput("table")
)

genePage <- div(
  # Return menu
  tags$div(class = "return-menu",
    tags$a(class = "return-link", href = "/","Return")
  ),
  tags$div(class = "container-fluid wrapper-content",
    # Header container
    tags$div(
     class = "header-container",
     htmlOutput("gene_title"),
     tags$div(
       htmlOutput("gene_meta_data")
     )
    ),
    # Section: Word cloud
    tags$div(
      class = "section",
      tags$h2(
        class = "secondary-title",
        "Word cloud"
      ),
      fluidRow(
        # Column for table
        column(
          width = 5,
          DT::dataTableOutput("words_table")
        ),
        # Column for image
        column(
          width = 7,
          imageOutput(outputId = "wordcloud")
        )
      )
    ),
    # Section: Additional information
    tags$div(
      class = "section",
      tags$h2(
        class = "secondary-title",
        "Additional information"
      ),
      tabsetPanel(
        # Tab: Sequence
        tabPanel("Sequence",
                 htmlOutput("protein_seq")
        ),
        # Tab: Resources
        tabPanel("Resources",
                 tags$ul(
                   htmlOutput("html_list")
                 )
        )
      )
    )
  )
)

router <- make_router(
  route("/", homePage),
  route("gene", genePage)
)

ui <- fluidPage(
  includeCSS("www/app.css"),
  # Navbar
  HTML('
<nav class="navbar navbar-default navbar-static-top">
  <div class="container-fluid">
	<!-- Brand and toggle get grouped for better mobile display -->
	<div class="navbar-header">
	  <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#bs-example-navbar-collapse-1" aria-expanded="false">
		<span class="sr-only">Toggle navigation</span>
		<span class="icon-bar"></span>
		<span class="icon-bar"></span>
		<span class="icon-bar"></span>
	  </button>
	  <a class="navbar-brand" href="/">TF Word Clouds</a>
	</div>

	<!-- Collect the nav links, forms, and other content for toggling -->
	<div class="collapse navbar-collapse" id="bs-example-navbar-collapse-1">
	  <ul class="nav navbar-nav navbar-right">
	  </ul>
	</div><!-- /.navbar-collapse -->
  </div><!-- /.container-fluid -->
</nav>
  '),
  # Router
  router$ui
)

server <- function(input, output, session) {

  router$server(input, output, session)
  
  # Reactive values
  rv <- reactiveValues(
    gene = "ADNP",
    entrezid = "Blank",
    uniacc = "Blank",
    unientr = "Blank",
    gene_family = "Blank",
    gene_cluster = "Blank",
    gene_evidence = "Blank",
    wordsTable = NULL
  )
  
  ## Update all reactive values based on rv$gene
  observeEvent(rv$gene, {
    
    # Filter the row for the gene in the URL
    geneRow <- genesTable %>%
      filter(Gene.Name == rv$gene)
    
    #### Entrez ID
    rv$entrezid <- as.character(
      geneRow %>% select(Gene.ID)
    )
    
    ### Uniprot accession
    rv$uniacc <- as.character(
      geneRow %>% select(UniProt.Accession)
    )
    
    ### Uniprot entry
    rv$unientr <- as.character(
      geneRow %>% select(UniProt.Entry)
    )
    
    ### Family
    rv$gene_family <- as.character(
      geneRow %>% select(Family)
    )
    
    ### Cluster
    rv$gene_cluster <- as.character(
      geneRow %>% select(Cluster)
    )
    
    ### Evidence
    rv$gene_evidence <- as.character(
      geneRow %>% select(Evidence)
    )
    
  })
  
  # Observe gene parameter in URL
  observe({
    if (!is.null(get_query_param()$gene_id)) {
      rv$gene <- get_query_param()$gene_id
    }
  })
  
  # Home Page ------------------
  
  # Home page: Render data table
  output$table <- DT::renderDataTable(
    genesTableSelectedCols,
    options = list(scrollX = TRUE,
                   dom = '<"search-input"f><tip>',
                   pageLength = 50,
                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#56AB7F', 'color': '#fff'});",
                     "}"),
                   searchPanes = list(
                     layout = "columns-12"
                   ),
                   columnDefs = list(
                     list(
                       targets = 1,
                       render = JS(
                         "function ( data, type, row, meta ) {",
                         "return '<a href=\"/#!/gene?gene_id=' + data + '#!/\">' + data + '</a>'",
                         "}"
                       )
                     ),
                     list(
                       targets = 2,
                       orderable = JS("false")
                     )
                   )
    ),
    style = "bootstrap",
    # class = "custom"
  )
  
  
  # Gene page ---------
  
  ## Gene page: Render header
  output$gene_title <- renderUI({
    HTML(paste('<h1 class="title">', rv$gene, '</h1>'))
  })
  
  
  ## Gene page: Render metadata in header
  output$gene_meta_data <- renderUI({
    
    ncbiHTML <- paste('<p>NCBI Gene: ',
                      '<a target="_blank" href="https://www.ncbi.nlm.nih.gov/gene/',
                      rv$entrezid,
                      '">',
                      rv$entrezid,
                      '</a></p>',
                      sep = ""
    )
    
    uniprotAccHtml <- paste('<p>Uniprot Accession: ',
                         '<a target="_blank" href="https://www.uniprot.org/uniprot/',
                         rv$uniacc,
                         '">',
                         rv$uniacc,
                         '</a></p>',
                         sep = ""
    )
    
    uniprotEntryHtml <- paste('<p>Uniprot Entry: ', rv$unientr, '</p>', sep = "")
    
    familyHtml <- paste('<p>Family: ', rv$gene_family, '</p>', sep = "")
    
    clusterHtml <- paste('<p>Cluster: ', rv$gene_cluster, '</p>', sep = "")
    
    evidenceHtml <- paste('<p>Evidence: ', rv$gene_evidence, '</p>', sep = "")
    
    return(
      HTML(
        paste(
          ncbiHTML,
          uniprotAccHtml,
          uniprotEntryHtml,
          familyHtml,
          clusterHtml,
          evidenceHtml
        )
      )
    )
    
  })
  
  ## Gene page: Render list of resources from "Additional information"
  output$html_list <- renderUI({
    
    htmlContent <- ""
    
    for (resource_name in colnames(genesTable)[11:21]) {

      resources_with_link <- names(resources_links)[resources_links != "#"]
      htmlResourceContent <- ""
      
      cellValue <- genesTable %>% filter(Gene.Name == rv$gene) %>% select({{ resource_name }})
      
      ### If there is an accession
      if(!is.na(cellValue) & cellValue != "") {
        cellValue <- as.character(cellValue)
        accessions <- strsplit(cellValue, ";")
        
        # If resource has a link, create an a tag
        if(resource_name %in% resources_with_link) {
          
          for (accession in accessions[[1]]) {
            
            base_link <- resources_links[resource_name]
            
            # Modify link if necessary
            full_link <- 
              switch (resource_name,
                      "CIS.BP" = paste(base_link, accession, ".txt", sep = ""),
                      "ReMap" = paste(base_link, accession, ":9606", sep = ""),
                      "UniPROBE" = paste(base_link, str_remove(accession, "UP0"), sep = ""),
                      paste(base_link, accession, sep = "")
              )
            
            
            
            htmlResourceContent <- paste(
              htmlResourceContent,
              "<li><a target='_blank' href='",
              full_link,
              "'>",
              accession,
              "</a></li>",
              sep = ""
            )
          }
          
        } else {
          # If resource doesn't have a link, create just an li tag
          for (accession in accessions[[1]]) {
            htmlResourceContent <- paste(
              htmlResourceContent,
              "<li>",
              accession,
              "</li>",
              sep = ""
            )
          }
        }
        
        ### Wrap list of accessions into an li for the specific resource (eg. ChipSeq-Atlas)
        htmlResourceContent <- paste(
          "<li>",
          resource_name,
          "<ul>",
          htmlResourceContent,
          "</ul>",
          "</li>",
          sep = ""
        )
        
        ### Add the new resource item (would be an li wrapping it) to htmlContent
        htmlContent <- paste(
          htmlContent,
          htmlResourceContent,
          sep = ""
        )
      }
    }
    
    ### Return
    HTML(htmlContent)
    
  })


  ## Gene page: Read TF-IDFs table  
  observeEvent(rv$entrezid, {
    rv$wordsTable <- read.table(
      gzfile(
        paste(
          tfIDFTables, "/",
          rv$entrezid, ".tsv.gz",
          sep = ""
        )
      ),
      quote = "", header = FALSE, sep = "\t"
    )
    
    colnames(rv$wordsTable) <- c("Word", "TF-IDF")
  })
  
  ## Gene page: Render TF-IDFs table
  output$words_table <- DT::renderDataTable(
    rv$wordsTable,
    style = "bootstrap",
    options = list(
      scrollX = TRUE,
      dom = '<"search-input"f><"#words-table" t><ip>',
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#56AB7F', 'color': '#fff', 'text-align': 'left'});",
        "}"),
      columnDefs = list(
        list(
          targets = c(1,2),
          className = "dt-body-left"
        )
      )
    )
  )
  
  ## Gene page: Render word cloud image
  output$wordcloud <- renderImage({
    
    # Filename
    filename <- normalizePath(
      file.path(
        figsFolder,
        paste(
          rv$entrezid, ".png", sep = ""
        )
      )
    )
    
    list(src = filename,
         alt = "No image")
    
  }, deleteFile = FALSE)
  
  ## Gene page: Render protein sequence
  output$protein_seq <- renderUI({
    protein_seq <- as.character(
      genesTable %>%
        filter(Gene.Name == rv$gene) %>%
        select(Sequence)
    )
    
    HTML(
      paste(
        '<p class="sequence-text">',
        protein_seq,
        '</p>'
      )
    )
    
  })

}

shinyApp(ui, server)
