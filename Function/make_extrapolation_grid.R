#'#***************************************************************************
#'@author : Felix Massiot-Granier
#'Adapted from Clara Peron, clara.peron@mnhn.fr  
#'and Cecilia O'Leary and Cole Monnahan | December 2020
#'@update : 2024-10-22
#'@email : Felix.Massiot-Granier@mnhn.fr
#'#***************************************************************************
#'@description :  Creating an extrapolation grid
#'#***************************************************************************
#'@param depth_range vector of min and max depth range that is to be used to predict in the @Extrapol_grid
#'
#'
#' Some reference
#' An example of how to create user-defined extrapolation regions (extents) for VAST.
#' From https://github.com/James-Thorson-NOAA/VAST/wiki/Creating-an-extrapolation-grid

  #cat("\n")
  #cat("2. Make_Extrapolation_Grid")
  #cat("\n")
  
  #ToKeepOrNotToKeep(toKeep)

#  1. Load bathy data

  # Gebco 2024 :
  bathy <- rast("Data/gebco_2025_n-27.1582_s-75.8145_w25.1719_e116.1563.nc")
  
  # Check the CRS of the bathymetry data
  crs(bathy)

# 2. Define depth range -------------------------------------------------------------

  depth_range= c(0, -2500)
  selected_isobaths <- depth_range

#3. Crop the Bathymetry Data to the Area of Interest --------------------

  # Define the cropping extent (xmin, xmax, ymin, ymax)
  crop_extent <- ext(62.5,75,-52,-45)
  
  # Crop the raster to the defined extent
  bathymetry_cropped <- suppressWarnings(crop(bathy, crop_extent))

#4.  Filter the bathymetry  ---------------------------------------------

  # Raster of the sampled depth 0 to 1200 m deep
  filtered_raster <- ifel(bathymetry_cropped >= depth_range[2] & bathymetry_cropped < depth_range[1], bathymetry_cropped, NA)
  
  plot_filtered_raster <- plot(filtered_raster, main = "Filtered Raster (0 to 2500 m)")
  
  # Nice map
    suppressMessages(
      SampledDepthRaster <- ggplot() +
      geom_spatraster(data = filtered_raster,) +
      scale_fill_whitebox_c(
        palette = "deep",
        labels = scales::label_number(suffix = "m"),
        n.breaks = 12,
        guide = guide_legend(reverse = TRUE)
      ) +
      labs(
        fill = "",
        title = paste0("Filtered Raster (",depth_range[1]," to ",-depth_range[2]," m)")
      ) +
      theme_light())
    
    
    
# 5. Define  projection system UTM and epsg_code ---------------------------------------------------------
    
    # Get UTM zone for conversion to UTM projection
    longitude <- sum(st_bbox(filtered_raster )[c("xmin", "xmax")])/2
    
    # Calculate the UTM zone
    utm_zone <- floor((longitude + 180) / 6) + 1
    
    # Determine the EPSG code
    epsg_code <- paste0("EPSG:", ifelse(sum(st_bbox(filtered_raster )[c("ymin", "ymax")])/2 >= 0, 326, 327), utm_zone)
    
    # Alternative would be to use the dedicated Kerguelen EPSG
    # epsg_code="EPSG:3336"
    
# 6. Remove outside EEZ-------------------------------------------------------------
    
    ZEE_Ker <- EEZ_sf[EEZ_sf$Country=="Kerguelen Islands",]
    
    filtered_raster <- mask(filtered_raster, ZEE_Ker)
    
    plot(filtered_raster, main = "Filtered Raster (0 to 1200 m)")
              
# 7. Remove Island ------------------------------------------------------------
    
    # Transform FR_island_sf wgs84 to FR_island_sf projected
    FR_island_sf_UTM <- st_transform(structure(FR_island_sf, proj4string = "EPSG:4326"), epsg_code)
    
    ggplot() +  geom_sf(data = FR_island_sf_UTM)
    
 #7.1 Project and draw buffer around islands --------------------------------------
    
    ## Draw buffer of 5 km around each location
    FR_island_sf_UTM_buffer <- st_buffer(FR_island_sf_UTM, 8000, nQuadSegs = 10)

    ggplot() + xlab("Longitude") + ylab("Latitude")+
      geom_sf(data = FR_island_sf_UTM_buffer ,
              fill = "darkgrey")
    
#7.2 Remove Island values from raster
    # Unproject
    FR_island_sf_buffer <- st_transform(structure(FR_island_sf_UTM_buffer, proj4string =  epsg_code),"EPSG:4326")

    # Rasterize polygon
    FR_island_buffer_raster <- rasterize(FR_island_sf_buffer, filtered_raster, field = NA, background = 1)
    
    # Remove values on the islands : create NA
    filtered_raster_hole <- filtered_raster*FR_island_buffer_raster

    plot(filtered_raster_hole, main = "Filtered Raster (0 to 1200 m)")
  
  suppressMessages(
    map_filtered_raster <- ggplot() +
      geom_spatraster(data = filtered_raster_hole,) +
      scale_fill_whitebox_c(
        palette = "deep",
        labels = scales::label_number(suffix = "m"),
        n.breaks = 12,
        guide = guide_legend(reverse = TRUE)
      ) +
      labs(
        fill = "",
        title = paste0("Filtered Raster (",depth_range[1]," to ",-depth_range[2]," m)")
      ) +
      theme_light())
  
# Transform to UTM
  filtered_raster_hole_UTM <- terra::project(filtered_raster_hole, epsg_code)

# 8. Create a polygon with the island for barrier SDMtmb ---------------------

  Poly_Ker_polygon <- st_as_sf(as.polygons(filtered_raster_hole_UTM)) # Félix ? Ca ne marche pas. Long et pb sur le plot
  
  #Remove the elevation info that we don't need
  Poly_Ker_polygon <- Poly_Ker_polygon[, setdiff(names(Poly_Ker_polygon), "elevation")]

  saveRDS(Poly_Ker_polygon,'Data/Ker_polygon.rds') # sauver l'objet FR_island_sf_UTM_buffer plutôt
  cat("Saved:", here('Data/Ker_polygon.rds'), "\n")
  saveRDS(FR_island_sf_UTM_buffer,'Data/Ker_polygon.rds') # sauver l'objet FR_island_sf_UTM_buffer plutôt
  cat("Saved:", here('Data/Ker_polygon.rds'), "\n")
  
# Previous approach used to get contour from a raster  
  # # Step 3: Generate Contours
  # # Define contour levels (isobaths)
  # contour_intervals <- selected_isobaths
  # 
  # # Generate contour lines for the cropped raster
  # contours <- terra::as.contour(bathymetry_cropped, levels = contour_intervals)
  # contours <- terra::as.contour(bathymetry_cropped, levels = 0)


# 8. Create an extrapolation grid  ----------------------------------------
  
## 8.1. Generate a grid that matches the area we want to predict on --------
  # Size of grid **in meters** (since working in UTM). Controls the resolution of the grid.
  cell_size <- 10000  # 10km # Careful if UTM is in m or km !

  # Make a grid from the filtered_raster_UTM 

  region_grid <- st_make_grid(filtered_raster_hole_UTM, cellsize = cell_size, what = "centers")
  
  # Trun into a vector for terra functions
  region_grid=vect(region_grid)
  extracted_values <- extract(filtered_raster_hole_UTM, region_grid)
  # Add extracted values as an attribute to the grid points
  region_grid$values <- extracted_values[, 2]  # Second column contains the raster values
  
  # Filter points within the raster's valid region (non-NA)
  region_grid <- region_grid[!is.na(region_grid$values), ]

  plot(region_grid)

# ### 4.2. Removing outer  ------------------------------
# #### 4.2.1 Removing Lakes, remote mounts -----------------------------------------
# # It hase already been done manually
# 
# # Lakes
# # removingLakesIsland<-locator()
# # removingOuter1<-locator()
# # removingOuter2<-locator()
# # saveRDS(removingLakesIsland, 'Data/Spatial/removingLakesIsland.rds')
# # saveRDS(removingOuter1, 'Data/Spatial/removingOuter1.rds')
# # saveRDS(removingOuter2, 'Data/Spatial/removingOuter2.rds')
# 
# removingOuter1<-data.frame(readRDS(here('Data/Spatial/removingOuter1.rds')))
# removingOuter2<-data.frame(readRDS(here('Data/Spatial/removingOuter2.rds')))
# removingLakesIsland<-data.frame(readRDS(here('Data/Spatial/removingLakesIsland.rds')))
# 
# removingLakesIsland <- rbind(removingLakesIsland, removingLakesIsland[1, ]) 
# removingOuter1 <- rbind(removingOuter1, removingOuter1[1, ]) 
# removingOuter2 <- rbind(removingOuter2, removingOuter2[1, ]) 
# 
# 
# # Create a multipolygon with locator polygons
# multipolygon <- st_multipolygon(list(
#   list(as.matrix(removingLakesIsland)),  # First polygon
#   list(as.matrix(removingOuter1)), # Second polygon
#   list(as.matrix(removingOuter2)) # third polygon
# ))
# 
# # Convert to an sf object
# multipolygon_sf <- st_sf(
#   id = 1,                              # Attribute column
#   geometry = st_sfc(multipolygon),     # Geometry
#   crs = epsg_code                          # Set CRS (e.g., UTM Zone 42S)
# )
# 
# # Remove the points within the multypoligon 
# region_grid_final <- erase(region_grid, vect(multipolygon_sf))
# 
# # Format the file to be an SF object and change the column name to depth
# region_grid_final <- st_as_sf(region_grid_final)
# names(region_grid_final)[names(region_grid_final) == "values"] <- "depth"


# Format the file to be an SF object and change the column name to depth
region_grid <- st_as_sf(region_grid)
names(region_grid)[names(region_grid) == "values"] <- "depth"


#### 4.2.2 Selecting only the EEZ -----------------------------------------
# 
# polygon_transformed <- suppressMessages(st_transform(EEZ_sf %>% filter(EEZ=="Kerguelen Islands Exclusive Economic Zone"), crs = epsg_code))
# 
# region_grid_final <- region_grid_final[
#   st_within(region_grid_final, polygon_transformed, sparse = FALSE),
# ]


### 4.3. Create and saving plot ---------------------------------------------

FinalGrid <- ggplot() +
  # Plot the multipolygon
  #geom_spatvector(data = multipolygon_sf, fill = "#a6cee3", color = "#1f78b4", alpha = 0.6) +
  # Plot region_grid
  geom_spatvector(data = region_grid, aes(color = "region_grid"), size = 1.5) +
  # Plot region_grid_final
  #geom_spatvector(data = region_grid_final, aes(color = "region_grid_final"), size = 1.5) +
  # Define custom color palette
  scale_color_manual(values = c("region_grid_final" = "#3B9AB2", "region_grid" = "#E1AF00")) +
  # Add title and legend
  labs(
    title = "Final grid and removed outer points",
    subtitle = "",
    color = "Point Type"
  ) +
  # Minimalist theme for a cleaner look
  theme_minimal() +
  # Adjust legend and title styling
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = ""
  )

if(SavePlots) {
  save_plot(plot = FinalGrid,
            device = Device,
            path = FigOutputPath)
}



### 4.4 Finalizing dataset format and saving ---------------------------------------
coords <- st_coordinates(region_grid)  # Extract coordinates
df_grid <- cbind(as.data.frame(region_grid["depth"]), coords)  # Combine with depth

df_grid <- with(df_grid,
<<<<<<< Updated upstream:Rscript/5.1.Make_extrapolation_grid.R
                  data.frame(X=X,
                             Y=Y,
=======
                  data.frame(Lon=X,
                             Lat=Y, 
>>>>>>> Stashed changes:Rscript/2.Make_Extrapolation_Grid.R
                             Id=seq_len(nrow(df_grid)),
                             Area_m2=((cell_size)^2), # Cell size in meters
                             depth=depth),
                  row=1:nrow(df_grid))
                 


# Saving
saveRDS(df_grid,"Data/Spatial/Extrapol_grid.rds")
cat("Saved:", here('Data/Spatial/Extrapol_grid.rds'), "\n")












# 
# # # 5. Create an extrapolation grid  
# # 
# # # Construct the extrapolation grid for VAST using sf package
# # # Size of grid **in meters** (since working in UTM). Controls the resolution of the grid.
# # cell_size <- 10000  # 10km # carefull if UTM is in m or km !
# # 
# # # Make grid
# # region_grid <- st_make_grid(region_polygon, cellsize = cell_size, what = "centers")
# # 
# # plot(region_grid)
# # # Convert region_grid to Spatial Points to SpatialPointsDataFrame
# # region_grid <- as(region_grid, "Spatial")
# # region_grid_sp <- as(region_grid, "SpatialPointsDataFrame")
# # 
# # # Combine shapefile data (region_polygon) with Spatial Points(region_grid_spatial) & place in SpatialPointsDataFrame data
# # # (this provides you with your strata identifier (here called  Id) in your data frame))
# # region_polygon <- as(region_polygon, "Spatial")
# # region_grid_sp@data <- over(region_grid, region_polygon)
# # 

# # 
# # # Convert the Spatial* objects to sf objects
# # region_grid_sf <- st_as_sf(region_grid)
# # region_polygon_sf <- st_as_sf(region_polygon)
# # 
# # # Ensure both objects have the same CRS
# # if (st_crs(region_grid_sf) != st_crs(region_polygon_sf)) {
# #   region_polygon_sf <- st_transform(region_polygon_sf, st_crs(region_grid_sf))
# # }
# # 
# # # Perform a spatial join
# # region_grid_sf <- st_intersection(region_grid_sf, region_polygon_sf)
# # 
# # # Convert back to sp if needed
# # region_grid_sp <- as(region_grid_sf, "Spatial")
# # 
# # ggplot() +
# #   geom_sf(data = region_polygon_sf, fill = "blue", alpha = 0.3, color = "darkblue", size = 0.5) +
# #   geom_sf(data = region_grid_sf, color = "red", size = 0.3, alpha = 0.7) +
# #   theme_minimal() +
# #   labs(
# #     title = "Region Polygon and Grid",
# #     subtitle = "Visualizing Spatial Relationships",
# #     x = "Longitude",
# #     y = "Latitude"
# #   )
# # 
# # 
# # 
# # # # Convert back to lon/lat coordinates (geographic coordinates) as that is what VAST uses
# # region_grid_FIN <- as.data.frame(spTransform(region_grid_sp, crs_UTM))
# # #   
# # # Reformat data  
# # region_df <- with(region_grid_FIN,
# #                   data.frame(Lon=coords.x1,
# #                              Lat=coords.x2, 
# #                              Id=seq_len(nrow(region_df)),
# #                              Area_km2=((cell_size)^2)),
# #                   row=1:nrow(region_grid_FIN))
# # 
# # # Filter out the grid that does not overlap (outside extent)
# # region <- subset(region_df, !is.na(Id))
# # 
# # # This is the final file needed.
# # str(region)
# # plot(region)
# # saveRDS(region, file = "Data/Poker_region.Rds")
# # # saveRDS(region, file = "Data/Poker_hole_region.Rds")
# # 
# # #  plot(region)
# # 
# # 
# # 
# # # 6. Get coordinates of the donut hole for the barrier function
# # # get all area without hole and then apply barrier 
# # # stop at region_polygon <- st_as_sf(region_polygon)
# # 
# # 
# # sp_hole_geo <- st_as_sf(sp_hole)
# # Hole <- as_Spatial(sp_hole_geo)
# # 
# # 
# # # Get UTM zone for conversion to UTM projection
# # Hole@proj4string <- crs_POKER
# # 
# # # Project in projected coordinates (UTM)
# # Hole_polygon <- spTransform(Hole, crs_UTM)
# # Hole_polygon <- st_as_sf(Hole_polygon)
# # 
# # saveRDS(Hole_polygon,'Data/Hole_polygon.rds')
# # 
# 
# 
# # 
# # #7.  Extrapolation of depth from Gebco 
# # library(terra)
# # # Warning from Clara
# # 
# # # Step 1: Load the GEBCO bathymetry data
# # gebco_file <- "Data/GEBCO_10_Dec_2024/gebco_2024_n-44.0991_s-54.8657_w58.7549_e77.6953.nc"  
# # gebco <- rast(gebco_file)
# # 
# # # Step 2: Create a data frame of coordinates (longitude and latitude)
# # coords <- region
# # 
# # # Step 3: Extract depth values
# # # Convert the coordinates into a SpatVector
# # points <- vect(region, geom = c("Lon", "Lat"), crs = crs(gebco))
# # 
# # # Extract depth values
# # depth_values <- extract(gebco, points)
# # 
# # # Combine depth values with the coordinates
# # result <- cbind(region, depth = depth_values[, 2])
# # 
# # saveRDS(result,"Data/extrapol_grid.rds")
# 
# 
# 
# # 

# # # 8. Getting area based data from Gebco 
# # 
# # # Step 1: Define Isobaths and Load Bathymetry Data
# # #selected_isobaths<-c(0,-100,-200,-300,-400,-500,-1000,-2000)
# # selected_isobaths<-c(0,-1200)
# # 
# # # Provide the path to your GEBCO file
# # gebco_file <- "Data/GEBCO_10_Dec_2024/gebco_2024_n-44.0991_s-54.8657_w58.7549_e77.6953.nc"  
# # 
# # # Load the bathymetry data as a terra raster object
# # bathymetry <- rast(gebco_file)
# # 
# # # Check the CRS of the bathymetry data
# # crs(bathymetry)
# # 
# # # Step 2: Crop the Bathymetry Data to the Area of Interest
# # # Define the cropping extent (xmin, xmax, ymin, ymax)
# # crop_extent <- ext(62.5,77,-53,-45)
# # 
# # # Crop the raster to the defined extent
# # bathymetry_cropped <- crop(bathymetry, crop_extent)
# # 
# # # Step 3: Generate Contours
# # # Define contour levels (isobaths)
# # contour_intervals <- selected_isobaths
# # 
# # # Generate contour lines for the cropped raster
# # contours <- terra::as.contour(bathymetry_cropped, levels = contour_intervals)
# # contours <- terra::as.contour(bathymetry_cropped, levels = 0)
# # # Convert contours to an `sf` object
# # isobaths <- st_as_sf(contours) %>%
# #   dplyr::mutate(level = -as.numeric(as.character(level)))
# # 
# # 
# # 
# # # Step 4: Generate Contours for Specific Isobath (e.g., 0m)
# # # Create contours for 0m depth
# # contours <- terra::as.contour(bathymetry_cropped, levels = 0)
# # 
# # # Convert contours to an sf object
# # raster_KerIsland <- st_as_sf(contours) %>%
# #   dplyr::mutate(level = -as.numeric(as.character(level)))
# # 
# # # Plot
# # p <- ggplot()+
# #   ylim(-53,-45)+
# #   xlim(62.5,77)+
# #   
# #   geom_sf(data = isobaths, aes(col = as.factor(level)))+
# #   scale_color_manual("Profondeur",
# #                      values = rep(c("black",wes_palette(5, name = "Zissou1")),
# #                                   length.out = length(selected_isobaths)))+
# #   theme(panel.background = element_rect(fill = "white", color = "white"),
# #         panel.grid = element_blank())+
# #   ggspatial::annotation_scale(location = "br",
# #                               bar_cols = c("grey40","white"),
# #                               line_col = "grey40",
# #                               text_col = "grey40")
# # 
# # 
# # p
# # 
# # 
# # region <- as_Spatial(isobaths)
# # 
# # 
# 
# 
# 
# 
# 




