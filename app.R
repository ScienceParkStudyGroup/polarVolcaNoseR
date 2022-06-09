# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PolaRVolcaNoseR - A Shiny app for multiple volcano plots
# Created by Joachim Goedhart (@joachimgoedhart) and Román Gonzalez-Prieto (@HombreRRo)
# first version 2021
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joachim Goedhart (C) 2021
# electronic mail address: j #dot# goedhart #at# uva #dot# nl
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# ToDo:
# warning when dots are located outside of sector (or clip!)
# Export to PDF
# Show table with hits

library(shiny)
library(ggplot2)
library(magrittr)
library(htmlwidgets)
library(ggrepel)
library(ggiraph)
library(dplyr)
library(glue)
library(scales)

source("themes.R")

#Read data
df_1 <- read.csv('S1-SATTs_tidy.csv') %>% mutate_at(3:5, round, 2)
df_2 <- read.csv('S2-SATTs_tidy.csv') %>% mutate_at(3:5, round, 2)

#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Create a reactive object here that we can share between all the sessions.
vals <- reactiveValues(count=0)

# Define UI for application that draws a histogram
ui <- fluidPage(
    #Required for proper formatting of the font of the girafeOutput
    tags$head(tags$style(type="text/css", "text {font-family: sans-serif}")),
    
    # Application title
    titlePanel("Polar Volcano plots"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width=3,
                     
                     radioButtons("dataset", "Data set:", choices = 
                                                     list(
                                                       "SUMO1-SATTs" = 1,
                                                       "SUMO2/3-SATTs" = 2),
                                                   selected =  1),
                     # 
                     # sliderInput("SATT_filter",
                     #             "SATT-index higher than:",
                     #             min = 0,
                     #             max = 1,
                     #             value = 0),
                     textInput("SATT_Range", "Range SATT values (min,max)", value = "0,1"),
                     
                     selectizeInput(inputId = 'user_gene_list',
                                    label = "User selected hits:",
                                    choices = "-",
                                    selected = "-",
                                    multiple = TRUE, # allow for multiple inputs
                                    options = list(create = TRUE)), 
            sliderInput("rotation",
                        "Rotate:",
                        min = 0,
                        max = 360,
                        value = 0),
            
            
            numericInput("pointSize", "Size of the datapoints", value=2),  

        checkboxInput(inputId = "polar",
                      label = "Polar plot",
                      value = TRUE),

        numericInput("scale_Diff", "Limit of the effect size (FC)", value = 14),
        numericInput("scale_P", "Limit of the significance (p-value)", value = 9),
        checkboxInput(inputId = "clip",
                      label = "Hide values that exceed the limits",
                      value = TRUE),
        # numericInput("plot_height", "Plot size (# pixels):", value = 600),
        # numericInput("plot_width", "Plot width (# pixels):", value = 800),
        checkboxInput(inputId = "dark", label = "Dark Theme", value = FALSE),
        NULL),

        # Show a plot of the generated distribution
        mainPanel(
            # plotOutput("coolplot", height = 'auto'),
          div(downloadButton("downloadHTML", "Download html file"),
              downloadButton("downloadPlotPDF", "Download pdf-file"),
              girafeOutput("iplot", height="800", width="800"), style = "float:left")
            # girafeOutput("iplot", height="600", width="800")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
  
  ##### Get Variables from the input ##############
  
  genelist.selected <- ""
  
  observe({
    updateSelectizeInput(session, "user_gene_list", selected = genelist.selected)
    
    df <- df_upload() %>% filter(!is.na(Diff))
    updateSelectizeInput(session, "user_gene_list", choices = df$GENE_NAME, selected = genelist.selected)
    })
  
  #Get the data
    df_upload <- reactive({
      
      if (input$dataset==1) {return(df_1)}
      else {return(df_2)}

      
    })
    
    plot_data <- reactive({
      
      if (input$dark) {line_color="white"} else {line_color="gray20"}
      
    df_tidy <- df_upload() %>% group_by(GENE_NAME) %>% mutate(ID = row_number(GENE_NAME))

    ##### SATT range #####
    rng_y <- as.numeric(strsplit(input$SATT_Range,",")[[1]])
    
    df_tidy <- df_tidy %>% filter(SATT.index > rng_y[1] & SATT.index < rng_y[2])
    
    df_minus <- df_tidy %>% filter(SATT.index < 0)

    
    # auto
    # scale_Diff <- max(df_tidy$Diff, na.rm=T)
    # Manual
    scale_Diff <- input$scale_Diff
    
    scale_P <- input$scale_P
    
    
    max_P <- max(df_tidy$minLogP, na.rm=T)
    max_Diff <- max(df_tidy$minLogP, na.rm=T)
    
    # Hide values that exceed those set by axis limits
    if (input$clip) {
      df_tidy <- df_tidy %>% filter(Diff <= scale_Diff) %>% filter(minLogP <= scale_P)
    }
    
    df_tidy <- df_tidy %>% mutate(x = (x= ID +(Diff/scale_Diff)))
    
    
    df_sectors <- df_tidy %>% ungroup() %>% select(Sumo, ID) %>% distinct()
    df_sectors <- df_sectors %>% mutate(x=ID, y=scale_P/2, width=1, height=scale_P*1.0)
    



    
    #######
    
    p <- ggplot(df_tidy, aes(x=x, y=minLogP))+
        geom_tile(data=df_sectors, aes(x=x+0.5, y=y, width=width, height=height, fill=factor(x)), color=line_color, alpha=0.1)
    
    
    #Generate x-axis scale when for polar Volcanos
    
    if (input$polar) {
      p <- p+
        ggplot2::annotate(geom='segment', x = 1, xend=1, y=scale_P*1.05, yend=scale_P*1.1, color='black', size=.25) +
        ggplot2::annotate(geom="text", x = 1, y=scale_P*1.15, label = '0.0', size=4) +
        ggplot2::annotate(geom="text", x = 1.95, y=scale_P*1.35, label = bquote(Log[2](FC)), size=5) +
        ggplot2::annotate(geom='segment', x = 1, xend=2, y=scale_P*1.05, yend=scale_P*1.05, color='black', size=.25) +
        ggplot2::annotate(geom='segment', x = 2, xend=2, y=scale_P*1.05, yend=scale_P*1.1, color='black', size=.25) +
        ggplot2::annotate(geom="text", x = 2, y=scale_P*1.15, label = paste(round(scale_Diff,1)), size=4)
        
        
    } else if (!input$polar) {
      
      #Generate x-axis scale when not presented as polar Volcanos
      p <- p+
        ggplot2::annotate(geom='segment', x = 1, xend=1, y=scale_P*1.01, yend=scale_P*1.03, color='black', size=.25) +
        ggplot2::annotate(geom="text", x = 1, y=scale_P*1.04, label = '0.0', size=4) +
        ggplot2::annotate(geom="text", x = 1.5, y=scale_P*1.08, label = bquote(Log[2](FC)), size=5) +
        ggplot2::annotate(geom='segment', x = 1, xend=2, y=scale_P*1.01, yend=scale_P*1.01, color='black', size=.25) +
        ggplot2::annotate(geom='segment', x = 2, xend=2, y=scale_P*1.01, yend=scale_P*1.03, color='black', size=.25) +
        ggplot2::annotate(geom="text", x = 2, y=scale_P*1.04, label = paste(round(scale_Diff,1)), size=4)
    }
    
    
    p <- p+ geom_point_interactive(aes(color=`SATT.index`, tooltip = glue("{GENE_NAME}\nDifference: {Diff}\n-Log[p]: {minLogP}\nSATT index: {SATT.index}"), data_id = GENE_NAME), size=input$pointSize)+
        scale_color_viridis_c(limits=c(0, 1), oob=squish)+
        theme_light(base_size = 16, base_family = "sans")
    
    
    if (input$dark) {p <- p+ theme_light_dark_bg(base_size = 16, base_family = "sans")}
    
    if(input$polar) {p <- p + coord_polar(start = input$rotation/360*2*pi)}
    
    df_top <- df_tidy %>% filter(GENE_NAME %in% input$user_gene_list) %>% filter(!is.na(Diff))
    observe({print(head(df_top))})
    
    # Highlight and label user selected hits
    p <-  p +
      geom_point(data=df_top, aes(x=`x`,y=`minLogP`), shape=1,color=line_color, size=input$pointSize+2)
    p <- p+  geom_text_repel(
        data = df_top,
        aes(label = GENE_NAME),
        # family="mono",
        # size = input$fnt_sz_cand,
        color=line_color,
        nudge_x = 0.2,
        nudge_y=0.2,
        box.padding = unit(0.9, "lines"),
        point.padding = unit(.4+input$pointSize*0.1, "lines"),show.legend=F
      )
    
    if (input$polar) {
      p <- p + geom_label(data=df_sectors, aes(x=x+0.5, y=height*1.3, label=Sumo, fill=factor(x)), alpha=0.2, color=line_color)
    } else {
      p <- p + geom_label(data=df_sectors, aes(x=x+0.5, y=height*0.95, label=Sumo, fill=factor(x)), alpha=0.2, color=line_color)
      
    }
    #Remove legend for sectors
    p <- p + guides(fill = "none")
    
    #Remove grid
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())
    
    #Set Aspect ratio to unity, relevant for non-polar plot
    p <- p + theme(aspect.ratio=1)
    
    
    p <- p + labs(x = NULL, y=bquote(' '*-Log[10]*'(p-value)'))
    
    
    # Dim non-identified
    # g <- girafe(ggobj = p, options = list(
    #     opts_hover_inv(css = "opacity:0.5;"),
    #     opts_hover(css = "stroke:black;r:4pt"),
    #     opts_tooltip(css= tooltip_css,delay_mouseover = 0, offy = -90,opacity = .75),
    #     opts_zoom(min = 1, max = 4)
    # ))
    # 
    # return(g)
    
    return(p)
    
    })
    
    
    #Save as HTML widget
    # htmlwidgets::saveWidget(file = "PolarVolcano.html", g)
    
    ##### Set width and height of the plot area
    # width <- reactive ({ input$plot_width })
    # height <- reactive ({ input$plot_width }) 
    # 
    # output$coolplot <- renderPlot(width = width, height = height, {
    #     plot(plot_data())
    # })
    # 
    
    output$iplot <- renderGirafe({
        
        if (input$dark) {
            tooltip_css <- "background-color:white;color:black;padding:10px;border-radius:10px;font-family:sans-serif;font-weight:bold"
        } else {
            tooltip_css <- "background-color:black;color:white;padding:10px;border-radius:10px;font-family:sans-serif;font-weight:bold"
        }
        
        w <- 8
        # h <- w*(input$plot_height/input$plot_width)
        
        x <- girafe(code = plot(plot_data()),
                    
                    width_svg = w, height_svg = w,
                    options = list(
                      # opts_hover_inv(css = "opacity:0.5;"),
                        opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),
                        opts_selection(
                            type = "multiple", css = "fill:#FF3333;stroke:black;"),
                        opts_tooltip(css= tooltip_css,offy = -40, opacity = .9),
                        opts_zoom(min = 1, max = 4)
                    ))
        x
    })
    
    ################ List of user-selected hits #########
    df_user <- reactive({
      
      
      df <- as.data.frame(df_upload())
      
      #select based on text input
      usr_selection <- c(input$user_gene_list,input$iplot_selected)
      if (length(usr_selection >0)) {
      df_selected_by_name <- df %>% filter(GENE_NAME %in% usr_selection) %>% filter(!is.na(Diff))  }
      else return (df <- data.frame())
      
      return(df_selected_by_name)
      
      
    })
    
    
    ####### SAVE ########
    output$downloadHTML <- downloadHandler(
      
      filename <- function() {
        paste("PolaRVolcano_", Sys.time(), ".html", sep = "")
      },
      content <- function(file) {
        if (input$dark) {
          tooltip_css <- "background-color:white;color:black;padding:10px;border-radius:10px;font-family:sans-serif;font-weight:bold"
        } else {
          tooltip_css <- "background-color:black;color:white;padding:10px;border-radius:10px;font-family:sans-serif;font-weight:bold"
        }

        saveWidget(widget = 
                     girafe(ggobj = plot_data(), width_svg = 10, height_svg = 10,
                            options = list(
                              opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),
                              opts_selection(
                                type = "multiple", css = "fill:#FF3333;stroke:black;"),
                              opts_tooltip(css= tooltip_css,offy = -40, opacity = .75),
                              opts_zoom(min = 1, max = 4)
                            ))
                   , file = file)
      })
    
    
    output$downloadPlotPDF <- downloadHandler(
      filename <- function() {
        paste("VolcaNoseR_", Sys.time(), ".pdf", sep = "")
      },
      content <- function(file) {
        pdf(file, width = input$plot_width/72, height = input$plot_height/72)
        plot(plot_data())
        
        dev.off()
      },
      contentType = "application/pdf" # MIME type of the image
    )
    
    selected_state <- reactive({
      input$iplot_selected
    })
    output$console <- renderPrint({
      # input$iplot_hovered
      input$iplot_selected
    })
    
    observeEvent( input$iplot_selected,{
      
      #If selected items are new, add them
      if (length(intersect(input$iplot_selected,input$user_gene_list))==0) {
        genelist.selected <<- c(input$user_gene_list,input$iplot_selected)
        updateSelectizeInput(session, "user_gene_list", selected = genelist.selected)
        #if not, invert selection
      } else {
        common <- intersect(input$user_gene_list,input$iplot_selected)
        combined <- union(input$user_gene_list,input$iplot_selected)
        genelist.selected <<- combined[! combined %in% common]
        updateSelectizeInput(session, "user_gene_list", selected = genelist.selected)
        
      }
      
      
      
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
