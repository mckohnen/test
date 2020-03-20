library(leaflet)
library(RColorBrewer)
library(scales)
library(ggplot2)
#library(lattice)
library(dplyr)
library(shinyjs)
#library(ggsci)

# library(reticulate)
# conda_list()
# use_condaenv("blast")


function(input, output, session) {
  #----------------------------------------------------------------------------#
  #                               read data                                    #
  #----------------------------------------------------------------------------#
  
  #-----------------------------------------------------------------#
  #                            read BLAST result                    #
  #-----------------------------------------------------------------#
  clearstatus <- reactiveValues(loaded = TRUE, clear = TRUE)
  observeEvent(input$file1, {
    clearstatus$clear <- FALSE
  }, priority = 1000)
  
  observeEvent(input$clear_stringSequence, {
    reset("stringSequence")
  }, priority = 1000)
  
  fasta_checker <- function() {
    
  }
  
  blaster_react <- eventReactive(input$searchSequence, {
    
    
    # if true will run demo or text from box
    if( clearstatus$loaded & clearstatus$clear ){
      #-----------------------------------------------------------------#
      #             Run from fasta sequence in box                      #
      #-----------------------------------------------------------------#
      if (input$stringSequence != "> my_corona") {
        my_path <- "data/inline.fasta"
        writeLines(input$stringSequence, my_path)
        fast_val <- validate_fasta(my_path, input$seq_type)
        if (!fast_val) {
          shiny::showNotification("Invalid sequence\n 
                                please check if the sequence contains header\n 
                                and is of the correct type ",
                                  duration = 10, closeButton = TRUE,
                                  type = "error")
          return(data.frame())
        }
        
      } else {
      #-----------------------------------------------------------------#
      #                   Run from examples                             #
      #-----------------------------------------------------------------#
        shiny::showNotification("Using default fasta sequence",
                                duration = 5, closeButton = TRUE,
                                type = "warning")
        if (input$seq_type == "nucleotide") {
          my_path <- "testquery/SARScov2_query_nucleotide.fasta"
        } else {
          my_path <- "testquery/SARScov2_query_protein.fasta"
        }
      }
      #-----------------------------------------------------------------#
      #                Run from uploaded fasta                          #
      #-----------------------------------------------------------------#
    } else { # will use fasta file
      clearstatus$clear <- TRUE
      clearstatus$loaded <- TRUE
      
      if (input$stringSequence != "> my_corona") {
        shiny::showNotification("Text box is not empty, but will use uploaded fasta file",
                                duration = 5, closeButton = TRUE,
                                type = "error")
      } 
      
      my_path <- input$file1$datapath
      if (input$seq_type == "nucleotide") {
        fast_val <- validate_fasta(my_path, "nucleotide")
      } else {
        fast_val <- validate_fasta(my_path, "protein")
      }
      if (!fast_val) {
        
        shiny::showNotification("Invalid sequence\n 
                                please check if the sequence contains header\n 
                                and is of the correct type ",
                                duration = 10, closeButton = TRUE,
                                type = "error")
        return(data.frame())
      }
      
      
      
    }
    reset("file1")
    reset("stringSequence")
    
    
    
    #--------------------------------------------------------------------------#
    #                                     Run BLAST                            #
    #--------------------------------------------------------------------------#
    print(my_path)
    
    #my_path <- "testquery/SARScov2_query_nucleotide.fasta"
    # Run Blast
    blast_strp <- paste0("blastp -query ", my_path, " -task 'blastp' -db db/protein/covid19 ")
    blast_strn <- paste0("blastn -query ", my_path, " -task 'megablast' -db db/nucleotide/covid19 ")
    blast_strr <- "-outfmt 6 -num_threads 1 > data/SARScov2_alignment.tsv"
    
    if (input$seq_type == "nucleotide") {
      system(paste0(blast_strn, blast_strr), wait = TRUE)
      anno_query <- read.csv("data/SARScov2_nucleotide_metadata.csv",
                             header=TRUE, stringsAsFactors = FALSE)
    } else {
      system(paste0(blast_strp, blast_strr), wait = TRUE)
      anno_query <- read.csv("data/SARScov2_protein_metadata.csv",
                             header=TRUE, stringsAsFactors = FALSE)
    }
    
    # Check if results are empty
    if (length(readLines("data/SARScov2_alignment.tsv")) == 0) {
      shiny::showNotification("No hits found, 
                                please check if sequence is of the correct type",
                              duration = 10, closeButton = TRUE,
                              type = "error")
      return(data.frame())
    } 
    
    align <- read.delim("data/SARScov2_alignment.tsv", 
                        header=FALSE, stringsAsFactors = FALSE)
    colnames(align) <- c("qaccver", "Accession", "pident", "length", 
                         "mismatch", "gapopen", "qstart", "qend", "sstart", 
                         "send", "evalue", "bitscore")
    
    
    align <- merge(align, anno_query, by = "Accession")
    align <- align[order(align$pident, - align$evalue, align$bitscore, 
                         decreasing = TRUE),]
    
    
    align %>% 
      filter(!Geo_Location == "") %>% 
      mutate(collection_months = substr(Collection_Date, 1, 7)) %>% 
      mutate(evalue_log10 = -log10(evalue + 0.000001))
    
  })
  
  
  #-----------------------------------------------------------------#
  #                        format BLAST result                      #
  #-----------------------------------------------------------------#
  blaster_form_react <- eventReactive({
    input$score_id
    input$searchSequence
  }, {
    
    if (nrow(blaster_react()) > 0 ) {
      
      score_ids <- c(pident = "Percent identity",
                     evalue_log10 = "evalue",
                     bitscore = "bitscore")
      score_id <- names(score_ids)[score_ids %in% input$score_id]
      
      #print(blaster_react())
      
      x <- blaster_react() %>% 
        filter(!Geo_Location == "") %>% 
        mutate(radius = cut(!! sym(score_id), 4)) 
      radius_levels <- setNames(seq(2, 8, 2), levels(x$radius))
      
      
      #radius_levels <- setNames(seq(3, 12, 3), levels(blaster$radius))
      x <- x %>% 
        mutate(radius = recode(radius, !!!radius_levels)) %>% 
        group_by(Geo_Location) %>% 
        mutate(idx_location = 1:n()) %>% 
        mutate(radiusfix = if_else(idx_location == 1 & sum(idx_location) > 1, 12, radius )) %>% 
        ungroup() %>% 
        mutate(radiusfix = factor(radiusfix, levels = sort(unique(radiusfix))))
      
      dots_pal <- colorFactor(c("grey20", "grey40", "grey60", "Tomato"), domain = levels(x$radius))
      #------------------------------------#
      #         Country mapper to sp       #
      #------------------------------------#
      #countries
      simp_id <- function(x){
        gsub(" ", "", tolower(x))
      }
      
      # Country mapper
      country_mapper <- data.frame(blast_id_orig = unique(x$Geo_Location),
                                   blast_id = simp_id(unique(x$Geo_Location)), 
                                   stringsAsFactors = FALSE) %>% 
        mutate(ADMIN  = match(blast_id, simp_id(countries$ADMIN))) %>% 
        mutate(ISO_A3 = match(blast_id, simp_id(countries$ISO_A3))) %>% 
        mutate(mapper = if_else(!is.na(ADMIN), ADMIN, ISO_A3)) %>% 
        filter(!is.na(mapper))
      
      
      
      # keep only countries in the blast results
      my_countries <- countries[country_mapper$mapper,]
      my_countries$blast_id <- country_mapper$blast_id_orig
      
      #print(country_mapper)
      
      
      # Add color by date
      
      x <- x %>% 
        # Fix collection date color
        mutate(Collection_Date2 = if_else(nchar(Collection_Date) == 7,
                                          paste0(Collection_Date, "-32"),
                                          Collection_Date)) %>% 
        mutate(fixdate = as.Date(gsub("32$", "15", Collection_Date2)))%>%
        mutate(fixdate2 = cut.Date(fixdate, 6, labels = FALSE)) %>%
        mutate(fixdate3 = cut.Date(fixdate, 6)) %>%
        mutate(col_collect = sort(unique(as.character(fixdate3)))[fixdate2]) %>% 
        # Fix Release date color
        mutate(fixdate = as.Date(Release_Date))%>%
        mutate(fixdate2 = cut.Date(fixdate, 6, labels = FALSE)) %>%
        mutate(fixdate3 = cut.Date(fixdate, 6)) %>%
        mutate(col_release = sort(unique(as.character(fixdate3)))[fixdate2]) 
      
      #------------------------------------#
      #     Ad coordinates to blaster      #
      #------------------------------------#
      #countries$blast_id
      idx <- match(x$Geo_Location, my_countries$blast_id)
      x$longitude <- my_countries$longitude[idx]
      x$latitude <- my_countries$latitude[idx]
      
      #print(head(x))
      x$dots_lab <- paste(sep = "<br/>",
                          paste0("<b><a href='https://www.ncbi.nlm.nih.gov/nuccore/", 
                                 x$Accession, 
                                 "'>", 
                                 x$Accession, 
                                 "</a></b>"),
                          paste0("Host : ", x$Host),
                          paste0("Geo Location : ", x$Geo_Location),
                          paste0("Collection Date : ", x$Collection_Date),
                          paste0("Release Date : ", x$Release_Date),
                          paste0("percent identity = ", x$pident),
                          paste0("evalue = ", x$evalue),
                          paste0("bitscore = ", x$bitscore)
      ) 
      #print(my_countries@data)
      list(df           = x,
           my_countries = my_countries,
           dots_pal     = dots_pal)
    } else {
      NULL
    }
    
  })
  
  
  
  #----------------------------------------------------------------------------#
  #                           Reactive widgets                                 #
  #----------------------------------------------------------------------------#
  # Country selector
  output$sel_country <- renderUI({
    selectInput(
      inputId = "sel_country",
      label = "Choose countries (default all):",
      choices = unique(na.omit(blaster_react()$Geo_Location)),
      multiple = TRUE
    )
  })
  
  # Date selector
  output$date_range <- renderUI({
    
    collection_months <- sort(unique(blaster_react()$collection_months))
    if (length(collection_months) == 0 ) {
      #collection_months <- c(0,0)
      collection_months <- NULL
      sliderTextInput(
        inputId = "date_range",
        label = "Date Range:", 
        choices = NULL,
        selected = NULL
      )
    } else {
      sliderTextInput(
        inputId = "date_range",
        label = "Date Range:", 
        choices = collection_months,
        selected = collection_months[c(1, length(collection_months))]
      )
    }
    
    
  })
  
  
  
  #----------------------------------------------------------------------------#
  #                                       Map                                  #
  #----------------------------------------------------------------------------#
  
  ## Interactive Map ###########################################

  # Create the map
  output$map <- renderLeaflet({
    leaflet() %>%
      setView(42, 16, 2) %>% 
      addProviderTiles(providers$CartoDB.DarkMatter)
  })
  
  fil_by_location <- function() {
    if (length(input$sel_country) == 0) {
      #input$sel_country
      unique(blaster_react()$Geo_Location)
    } else {
      input$sel_country
    }
  }
  
  fil_by_collection_date <- function() {
    if (length(input$date_range) == 0) {
      #input$sel_country
      c(min(blaster_react()$collection_months), max(blaster_react()$collection_months))
    } else {
      input$date_range
    }
  }
  #----------------------------------------------------------------------------#
  #                              Color areas                                   #
  #----------------------------------------------------------------------------#

  observe({
    #print(input$date_range)
    if (nrow(blaster_react()) > 0 & 
        !is.null(blaster_form_react()) #&
        #input$date_range[1] != 0  & input$date_range[2] != 0 
        ) {
      # print(blaster_form_react()$df)
      #print(blaster_form_react()$my_countries@data)
      blaster_summ <- blaster_form_react()$df %>%
        filter(Geo_Location %in% fil_by_location()) %>%
        filter(collection_months >= fil_by_collection_date()[1] &
                 collection_months <= fil_by_collection_date()[2]) %>%
        mutate(Collection_Date2 = if_else(nchar(Collection_Date) == 7,
                                          paste0(Collection_Date, "-32"),
                                          Collection_Date)) %>%
        # for each country keep only the top hit
        group_by(Geo_Location) %>%
        top_n(n = 1, pident) %>%
        top_n(n = 1, -evalue) %>%
        top_n(n = 1, bitscore) %>%
        top_n(n = -1, Collection_Date2) %>%
        top_n(n = -1, Release_Date) %>%
        # If there's more than one entry keep only the first one
        mutate(cums = 1) %>%
        mutate(cums = cumsum(cums)) %>%
        filter(cums == 1)

      
      
      countries_react <- blaster_form_react()$my_countries[blaster_form_react()$my_countries$blast_id %in% blaster_summ$Geo_Location,]
      colID <- names(color_area_IDs)[color_area_IDs %in% input$sel_area_col]

      #print(countries_react@data)
      
      if (colID == "none") {
        blaster_summ <- blaster_summ %>%
          mutate(area_col = "")

        leafletProxy("map", data = blaster_summ) %>%
          clearControls() %>%
          clearShapes()
      } else {
        
        blaster_summ <- blaster_summ %>%
          mutate(area_col = !! sym(colID))
        
        
        #print(countries_react$blast_id)
        #print(blaster_summ)
        countries_react$density <- blaster_summ$area_col[match(countries_react$blast_id, blaster_summ$Geo_Location)]
        
        
        #print(countries_react@data)
        pal <- colorFactor("YlOrRd", domain = sort(unique(countries_react$density)))

        
        leafletProxy("map", data = blaster_summ) %>%
          clearControls() %>%
          clearShapes() %>%
          addPolygons(data = countries_react,
                      fillColor = ~pal(countries_react$density),
                      weight = 0.5,
                      opacity = 1,
                      color = "white",
                      dashArray = "3",
                      fillOpacity = 0.7
          ) %>%
          addLegend(pal = pal, values = ~countries_react$density, opacity = 0.7, title = NULL,
                    position = "topright")
      }
    } else {
      leafletProxy("map") %>%
        clearControls() %>%
        clearShapes() %>%
        clearMarkerClusters() %>%
        clearMarkers()
    }

  })

  #----------------------------------------------------------------------------#
  #                              Add clusters                                  #
  #----------------------------------------------------------------------------#

  observe({
    if (nrow(blaster_react()) > 0 & !is.null(blaster_form_react())) {
      blaster_map <- blaster_form_react()$df %>%
        filter(Geo_Location %in% fil_by_location()) %>%
        filter(collection_months >= fil_by_collection_date()[1] &
                 collection_months <= fil_by_collection_date()[2])

      dots_pal <- blaster_form_react()$dots_pal
      #dots_pal <- colorFactor(c("grey20", "grey40", "grey60", "Tomato"), domain = levels(blaster_map$radius))
      #print(as.data.frame(blaster_map[is.na(blaster_map$latitude),]))

      leafletProxy("map", data = blaster_map) %>%
        clearMarkerClusters() %>%
        clearMarkers() %>%
        addCircleMarkers(lng = blaster_map$longitude,
                         lat = blaster_map$latitude,
                         radius = blaster_map$radiusfix,
                         color = ~dots_pal(blaster_map$radius),
                         clusterOptions = markerClusterOptions(
                           spiderfyDistanceMultiplier=1.2
                         ),
                         popup = blaster_map$dots_lab,
                         fillOpacity = 1)
    }

  })

  
  
  # 
  # ## Data Explorer ###########################################
  ##--------------------------------------------------------------------------##
  ##              Get main table and fitler according to options              ##
  ##--------------------------------------------------------------------------##
  
  output$blaster_ui <- DT::renderDataTable({
    
    df <- blaster_form_react()$df %>% 
      mutate(Release_Date = as.Date.character(Release_Date)) %>%
      filter(collection_months >= fil_by_collection_date()[1] &
               collection_months <= fil_by_collection_date()[2]) %>% 
      select(Accession, pident, evalue, bitscore, Geo_Location, Host,
             Release_Date, Collection_Date, length, mismatch, gapopen,
             qstart, qend, sstart, send, Length, Isolation_Source, Species) %>% 
      filter(Geo_Location %in% fil_by_location()) 
      
      
    action <- DT::dataTableAjax(session, df, outputId = "blaster_ui")
    
    DT::datatable(df, options = list(ajax = list(url = action)), escape = FALSE, class = "nowrap display")
  })
  
  
  
  
  observe(
    output$gg_data_months <- renderPlot({
      
      # Input data is filtered 
      blaster_form_react()$df %>% 
        filter(Geo_Location %in% fil_by_location()) %>% 
        filter(collection_months >= fil_by_collection_date()[1] &
                 collection_months <= fil_by_collection_date()[2]) %>% 
        mutate(collection_months = factor(collection_months, 
                                          levels = sort(unique(collection_months)))) %>% 
        ggplot(aes(x = collection_months, fill = collection_months)) +
        geom_bar(stat = "count") +
        theme_dark() + 
        #scale_fill_tron()  +
        theme(#legend.position = "none",
          text = element_text(colour = "white"),
          axis.text.x = element_text(colour = "white"),
          axis.text.y = element_text(colour = "white"),
          plot.background =element_rect(fill = "#2D2D2D"),
          panel.background = element_rect(fill = "#2D2D2D"),
          axis.line=element_blank(),
          
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          #axis.title.y=element_blank(),
          panel.grid = element_blank(),
          legend.position="none") 
      
        
    },
    #width  = 100, 
    height = 80
    )
  )
  
  
  
  
  
}
