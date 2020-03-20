library(leaflet)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)


dashboardPage(

  dashboardHeader(title = "Map my Corona"),
  
  #----------------------------------------------------------------------------#
  #                                Sidebar                                     #
  #----------------------------------------------------------------------------#
  dashboardSidebar(

    sidebarMenu(

      menuItem(
        "Map",
        tabName = "maps",
        icon = icon("globe")
      ),

      menuItem(
        "Data explorer",
        tabName = "dataexplore",
        icon = icon("bar-chart")
      ),

      # Horizontal line ----
      tags$hr(),
      h4("Search sequence", align = "center"),
      h6("from Fasta or from search box", align = "center"),

      # Input: Select a file ----
      fileInput("file1", "Choose Fasta File",
                multiple = FALSE,
                accept = c(".RDS", ".fa", ".fasta")),

      # Search fasta sequence
      textAreaInput("stringSequence", "Search sequence...", "> my_corona"),
      actionBttn(
        inputId = "clear_stringSequence",
        label = "Clear", 
        style = "minimal",
        size = "xs",
        color = "default"
      ),
      
      awesomeRadio(
        inputId = "seq_type",
        label = "Sequence Type",
        choices = c("nucleotide", "protein"),
        selected = "nucleotide",
      ),

      actionBttn(
        inputId = "searchSequence",
        label = "Submit",
        style = "jelly",
        size  = "sm",
        color = "danger"
      ),
      h6("Click Submit without any input", align = "center"),
      h6("to display results for MT188340 (nucleotide)", align = "center"),
      h6("or QHN73805 (protein)", align = "center"),


      br(),
      tags$hr(),
      h4("Options", align = "center"),
      #h6("from Fasta or from search box", align = "center"),

      awesomeRadio(
        inputId = "score_id",
        label = "Select score metric",
        choices = c("Percent identity", "evalue", "bitscore"),
        selected = "Percent identity",
      ),


      selectInput(
        inputId  = "sel_area_col",
        label    = "Color Area By:",
        choices  = unname(color_area_IDs),
        selected = unname(color_area_IDs)[1],
        multiple = FALSE
      )


    )
  ),

  #----------------------------------------------------------------------------#
  #                              dashboardBody                                 #
  #----------------------------------------------------------------------------#
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(
        HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             }
             "
        )
      )
    ),
    
    tabItems(

      tabItem(
        tabName = "maps",
        fluidRow(
          box(
            height = 100,
            width = 4,
            background = "black",
            plotOutput(outputId = "gg_data_months")
          ),

          box(
            height = 100,
            width = 4,
            background = "black",
            uiOutput("date_range")
          ),

          box(
            height = 100,
            width = 4,
            background = "black",
            # Input: Selector for choosing Regulatory Variance metric ----
            uiOutput("sel_country")
          )

        ),

        fluidRow(
          tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
          box(
            #height = 300,
            width = 12,
            leafletOutput("map")
          )
        )
      ),

      tabItem(
        tabName = "dataexplore",
        div(style = 'overflow-x: scroll',
            #DT::DTOutput("blaster"))
            #checkboxInput("dt_sel", "sel/desel all")
            DT::dataTableOutput('blaster_ui'))
        #DT::dataTableOutput("blaster")
      )

    )
  ),
  skin = "black"
)



# 
# dashboardPage(
#   
#   dashboardHeader(title = "Map my Corona"),
#   #----------------------------------------------------------------------------#
#   #                                Sidebar                                     #
#   #----------------------------------------------------------------------------#
#   dashboardSidebar(
#     
#     sidebarMenu(
#       
#       menuItem(
#         "Map", 
#         tabName = "maps", 
#         icon = icon("globe")
#       ),
#       
#       menuItem(
#         "Data explorer", 
#         tabName = "dataexplore", 
#         icon = icon("bar-chart")
#       ),
#       
#       # Horizontal line ----
#       tags$hr(),
#       h4("Search sequence", align = "center"),
#       h6("from Fasta or from search box", align = "center"),
#       
#       # Input: Select a file ----
#       fileInput("file1", "Choose Fasta File",
#                 multiple = FALSE,
#                 accept = c(".RDS", ".fa", ".fasta")),
#       
#       # Search fasta sequence
#       textAreaInput("stringSequence", "Search sequence...", "> my_corona\natatataaccddf"),
#       
#       actionBttn(
#         inputId = "searchSequence",
#         label = "Submit",
#         style = "jelly", 
#         size  = "sm",
#         color = "danger"
#       ),
#       h6("Click Submit without any input", align = "center"),
#       h6("to display results for MT188340", align = "center"),
#       
#       
#       br(),
#       tags$hr(),
#       h4("Options", align = "center"),
#       #h6("from Fasta or from search box", align = "center"),
#       
#       awesomeRadio(
#         inputId = "score_id",
#         label = "Select score metric", 
#         choices = c("Percent identity", "evalue", "bitscore"),
#         selected = "Percent identity",
#       ),
#       
#       
#       selectInput(
#         inputId  = "sel_area_col",
#         label    = "Color Area By:",
#         choices  = unname(color_area_IDs),
#         selected = unname(color_area_IDs)[1],
#         multiple = FALSE
#       )
#       
#       
#     )
#   ),
#   
#   #----------------------------------------------------------------------------#
#   #                              dashboardBody                                 #
#   #----------------------------------------------------------------------------#
#   dashboardBody(
#     useShinyjs(),
#     tabItems(
#       
#       
#       tabItem(
#         tabName = "maps",
#         fluidRow(
#           box(
#             height = 100, 
#             width = 4, 
#             background = "black",
#             plotOutput(outputId = "gg_data_months")
#           ),
#           
#           box(
#             height = 100, 
#             width = 4, 
#             background = "black",
#             uiOutput("date_range")
#           ),
#           
#           box(
#             height = 100,
#             width = 4, 
#             background = "black",
#             # Input: Selector for choosing Regulatory Variance metric ----
#             uiOutput("sel_country")
#           )
#           
#         ),
#         
#         fluidRow(
#           tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
#           box(
#             #height = 300,
#             width = 12, 
#             leafletOutput("map")
#           )
#         )
#       ),
#       
#       tabItem(
#         tabName = "dataexplore",
#         
#         fluidRow(
#           box(
#             height = 100, 
#             width = 4, 
#             background = "black",
#             plotOutput(outputId = "gg_data_months")
#           ),
#           
#           box(
#             height = 100, 
#             width = 4, 
#             background = "black",
#             uiOutput("date_range")
#           ),
#           
#           box(
#             height = 100,
#             width = 4, 
#             background = "black",
#             # Input: Selector for choosing Regulatory Variance metric ----
#             uiOutput("sel_country")
#           )
#           
#         ),
#         fluidRow(
#           div(style = 'overflow-x: scroll', 
#               #DT::DTOutput("blaster"))
#               #checkboxInput("dt_sel", "sel/desel all")
#               DT::dataTableOutput('blaster_ui'))
#         )
#         
#       )
#       
#     )
#   ), 
#   skin = "black"
# )


