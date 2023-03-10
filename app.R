# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(mizer)
require(tidyverse)
sim<-readRDS("sim4.rds")
params<-sim@params
params<-upgradeParams(params)
speciesnames<-params@species_params$species
obsy <- readRDS("~/Dropbox/Nila SouthEastAus/SouthEastAustraliaMizer/yieldObs_timeVariant.rds")
obsy<-reshape2::melt(obsy,id.var=1)
names(obsy)<-c("time","sp","value")

# gear functions
{trawl <- function(w,W25_S,W50_S,...){
  
  SR <- W50_S - W25_S
  S1 <- W50_S*log(3)/SR
  S2 <- S1 / W50_S
  sel <-1 / (1 + exp(S1 - S2*w))
  return(sel)
}
  
  gillnet<-function(w,W50_B,sigma_B,...){
    sel<-exp(-(w-W50_B)^2/sigma_B^2)
    return(sel)
  }}


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Exploring Time Series"),
   
     fluidRow(
       column(4, wellPanel(
      #   sliderInput("Rmax", "log10 Maximum Recruitment:", min = round(log10(min(params@species_params$R_max)),0.1), max = round(log10(max(params@species_params$R_max)),0.1), value = -5,
      #                 step = 0.5),
         sliderInput("erepro", "Reproductive Efficiency:", min = 0.0001, max = 1, value = 0.5,
                     step = 0.01),
          sliderInput("catchability", "Catchability:", min =0.0001, max = 1, value = 0.5,
                      step = 0.1),    
         selectInput("species", "Species:",
                     speciesnames)
       )),
       column(6,
              plotOutput("distPlot", width = 600, height = 600)
       )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
  output$distPlot <- renderPlot({
    
    
    #change erepro of selected species
    species_params(params)$erepro[params@species_params$species==input$species] <- input$erepro
    #change Rmax of selected species
   # species_params(params)$R_max[params@species_params$species==input$species] <- input$Rmax
    #change catchability of selected species
    gear_params(params)$catchability[params@species_params$species==input$species] <-input$catchability
    
    #params <- setReproduction(params)
    
    #re-run model
    simt <- project(params, effort=sim@effort)
    
 
    # output adjusted modelled yields and reshape for plotting
    #y <- getYield(simt)
    #y <- reshape2::melt(y)
    
    
    p<-plotYieldGear(simt,species=input$species)  +
      geom_point(data=filter(obsy,sp==input$species),aes(x = time, y = value, 
                                                   colour = sp),size=0.6) +
      #facet_wrap(~sp) +
      scale_y_continuous(name = "Yield [g/year]")  +
      scale_colour_manual(values = sim@params@linecolour) +
      xlim(1980, 2010)
    p
    

   })
}

# Run the application 
shinyApp(ui = ui, server = server)

