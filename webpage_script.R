library(shiny)
library(leaflet)
library(mapview)
library(data.table) #for data
library(sp) #for coordinates 
library(gstat) #for kriging
library(raster)



#C:/project/new
#C:/project/data



# Set working directory to the folder containing this script:
# (assumes you use R studio)
basedir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(basedir, "data", sep=.Platform$file.sep))


file_list <- gsub("\\.csv$","", list.files(pattern="\\.csv$"))

file_list = list.files(pattern="*.csv")

data_list <- vector("list", "length" = length(file_list))

dir = paste(basedir, "new", sep=.Platform$file.sep)

for (i in seq_along(file_list)) {
  filename = file_list[[i]]
  
  number_of_rows <- 10
  lanlongcity = read.csv(filename, header = FALSE, sep=";", nrows = number_of_rows,
                         skip = 0)
  lanlongcity
  
  lanlongcity = lanlongcity[c(1,6)]
  lanlongcity = as.data.frame(t(lanlongcity))
  tempDF <- lanlongcity
  tempDF[] <- lapply(lanlongcity, as.character)
  colnames(lanlongcity) <- tempDF[1, ]
  lanlongcity = lanlongcity[-1, ] 
  tempDF <- NULL
  
  #lanlongcity
  
  weather_data = read.csv(filename, header = TRUE, sep=";", skip = 10)
  weather_data_full = merge(x = lanlongcity, y = weather_data, by = NULL)
  
  # weather_data_full
  
  #install.packages("data.table") 
  #library(data.table)
  weather_data_full <- data.table(weather_data_full)
  
  
  
  
  #weather_data_full[Year == 2018 & Month == 4 & Day == 15]
  weather_specific = weather_data_full[Year == 2018 & Month == 4 & Day == 14 & Hour == 11]
  
  file <- as.character(paste("file", i, ".csv", sep="")) 
  write.csv(weather_specific, file = file.path(dir, file)) 
  
  
}

setwd(dir)
weather_d <- 
  do.call(rbind,
          lapply(list.files(path = dir), read.csv))

weather <- 
  do.call(rbind,
          lapply(list.files(path = dir), read.csv))

weather_data_full <- weather_data <- NULL
data_list <- NULL 
lanlongcity <- NULL
weather_specific <- NULL
setwd(paste(basedir, "raster", sep=.Platform$file.sep))


# Define UI ----
ui <- fluidPage(
  # Header or Title Panel 
  titlePanel(title = h4("WEB PAGE", align="center")),
  sidebarLayout(
    # Sidebar panel
    sidebarPanel(
      # h3("Select Variable", align="center"),
      selectInput("var1", 
                  h3("Select Variable"), 
                  choices = list("Temperature" = "Temperature...2.m.above.gnd.", 
                                 "Wind" = "Wind.Direction...10.m.above.gnd.", 
                                 "Radiation" = "Wind.Speed...10.m.above.gnd."),
                  selected = "Temperature...2.m.above.gnd."),
      
      h3("Select Date"),
      fluidRow(
        
        column(4,
               
               selectInput("Year",
                           "Year",
                           c("All",
                             unique(as.character(weather_d$Year))))
        ),
        column(4,
               selectInput("Month",
                           "Month",
                           c("All",
                             unique(as.character(weather_d$Month))))
        ),
        column(4,
               selectInput("Day",
                           "Day",
                           c("All",
                             unique(as.character(weather_d$Day))))
        )
      ),
      
      fluidRow(
        align="center",
        submitButton("Submit") 
      )
    ),
    
    # Main Panel
    mainPanel(
      tabsetPanel(type="tab", 
                  tabPanel("Map", 
                           leafletOutput("mapplot"),
                           #  mapview:::plainViewOutput("mapplot"),
                           verbatimTextOutput("summary")),
                  tabPanel("Table", tableOutput("data")),
                  tabPanel("Plot",
                           # fluidRow(...)
                           plotOutput("myhist"),
                           plotOutput("plot2")
                  ),
                  tabPanel("Structure")
                  
      )
      
    )
    
  )
)  


# Define server logic ----
server <- function(input, output) {
  
  station <- data.frame(weather,
                        lat = c("LON"),
                        long = c("LAT"))
  station$LON <- as.numeric(station$LON)
  station$LAT <- as.numeric(station$LAT)
  
  coordinates(station) <- ~LON + LAT
  
  # Set the projection. They were latitude and longitude, so use WGS84 long-lat projection
  proj4string(station) <- CRS("+init=epsg:4326")
  
  station <- spTransform(station, CRSobj = CRS("+init=epsg:3857"))
  
  # View the station location using the mapview function
  
  
  
  ori <- SpatialPoints(cbind(8, 54), proj4string =  CRS("+init=epsg:4326")) 
  # Convert the projection of ori
  # Use EPSG: 3857 (Spherical Mercator)
  ori_t <- spTransform(ori, CRSobj = CRS("+init=epsg:3857"))
  
  coordinates(ori_t)
  
  x_ori <- round(coordinates(ori_t)[1, 1]/100) * 100
  y_ori <- round(coordinates(ori_t)[1, 2]/100) * 100
  
  # Define how many cells for x and y axis
  x_cell <-90
  y_cell <- 110
  
  # Define the resolution to be 1000 meters
  cell_size <- 10000
  
  # Create the extent
  ext <- extent(x_ori, x_ori + (x_cell * cell_size), y_ori, y_ori + (y_cell * cell_size)) 
  
  # Initialize a raster layer
  ras <- raster(ext)
  
  # Set the resolution to be
  res(ras) <- c(cell_size, cell_size)
  ras[] <- 0
  
  # Project the raster
  projection(ras) <- CRS("+init=epsg:3857")
  
  
  # Convert to spatial pixel
  st_grid <- rasterToPoints(ras, spatial = TRUE)
  gridded(st_grid) <- TRUE
  st_grid <- as(st_grid, "SpatialPixels")
  
  
  TheVariogram=variogram(Temperature...2.m.above.gnd.~1, data=station)
  #plot(TheVariogram)
  TheVariogramModel <- vgm(model="Sph")
  #plot(TheVariogram, model=TheVariogramModel) 
  
  
  lzn.fit <- fit.variogram(TheVariogram, model=TheVariogramModel)
  #plot(TheVariogram, model=lzn.fit)
  
  
  lzn.kriged <- krige(Temperature...2.m.above.gnd. ~ 1, station, st_grid, model=lzn.fit)
  
  rasterDF <- raster(lzn.kriged)
  projection(rasterDF) <- CRS("+init=epsg:3857")
  #  rasterDF
  
  
  writeRaster(rasterDF, filename="test.tif", format="GTiff", overwrite=TRUE)
  
  Power_Plant <- matrix(c("Ringhals Nuclear Power Plant", 12.110833,57.259722), ncol=3)
  colnames(Power_Plant) <- c('NAME', 'LON', 'LAT')
  Power_Plant <- data.table(Power_Plant)
  Power_Plant$LON <- as.numeric(as.character(Power_Plant$LON))
  Power_Plant$LAT <- as.numeric(as.character(Power_Plant$LAT))
  
  icon.npp <- makeAwesomeIcon(icon = "industry", lib = "fa", markerColor = 'lightgreen', iconColor = 'black')
  
  r <- raster("test.tif")
  pal <- colorNumeric(c("#ffffff", "#e56979", "#ef001f"), values(r),
                      na.color = "transparent")
  
  output$mapplot <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Stamen.TonerLite,
                       options = providerTileOptions(noWrap = TRUE))%>%
      addTiles()%>%
      setView(lng=10,
              lat=55,
              zoom=5)%>%
      addProviderTiles("Esri.WorldImagery",
                       group="Esri WorldImagery")%>%
      addTiles(group="Open Street Map")%>%
      addMouseCoordinates(style = "basic") %>%
      # Add 2 marker groups
      addAwesomeMarkers(data=Power_Plant,
                        lng=~LON,
                        lat=~LAT,
                        group="Power_Plant",
                        icon = icon.npp
      ) %>%
      
      
      addCircleMarkers(data=weather_d,
                       lng=~LON,
                       lat=~LAT,
                       radius=0.2,
                       color="black",
                       fillColor="black",
                       stroke = TRUE,
                       fillOpacity = 0.8,
                       group="Weather Stations",
                       label = ~ as.character(Temperature...2.m.above.gnd.) 
                       
      ) %>%
      addRasterImage(rasterDF,
                     colors = pal,
                     opacity = 0.8,
                     group="Surface"
      ) %>%
      addLegend(pal = pal,
                position = 'topright',
                values = values(r), 
                title = "Surface"
      ) %>%
      # Add the control widget
      addLayersControl(overlayGroups = c("Weather Stations","Power_Plant", "Surface") ,
                       position = 'bottomright',
                       baseGroups = c("Open Street Map","Esri WorldImagery"),
                       options = layersControlOptions(collapsed = TRUE))
  })
  
  output$myhist <- renderPlot({
    hist(weather_d[,input$var1], 
         col= "#EC5965" , 
         main = paste("Histogram of" , input$var1), 
         xlab=paste(input$var1))
  })
  
  output$plot2 <- renderPlot({plot(TheVariogram, model=lzn.fit)
  })
  
  output$data <- renderTable({
    weather_d[,c('CITY','LON','LAT', input$var1)]}, 
    striped = TRUE, 
    hover = TRUE, 
    bordered = TRUE, 
    rownames = TRUE, 
    colnames = TRUE)
}


# Run the app ----
shinyApp(ui = ui, server = server) 