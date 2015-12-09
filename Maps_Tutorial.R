# SPATIAL DATA & MAPS IN R  - WORKSHOP
# March 31, 2015 
# Matthew Bayly

# INDEX.
#============================================================================================#
# 1. Get digital elevation model of Region, download as tiles and paste together 
	#1A list all DEM files to create & download from R
	#1B (raster) read in raster files & assign a unique name
	#1C (Mosaic) paste individual tiles together into a single raster 
	#1D (writeRaster) Save rasters
#============================================================================================#
# 2. Get shapefiles of political boarders & other features to add to Map
	#2A (readOGR) - load shapefiles into R	
	#2B Check to see that everything is OK with shapefiles & rasters
#============================================================================================#
# 3. Add compass rose
#============================================================================================#
# 4. Simple Plot: Convert DEM to fancy hillshade map & plot
	#4A (raster) load in dem
	#4B (extent) crop down & trim raster for plotting 
	#4C (crop) trim down rasters to specified extent with 'crop'
	#4D (terrain & hillShade) make hillshade layer (to view topographical relief)
	#4E crop & manipulate shapefile for plotting	
	#4F plot at continuous colour intervals for elevation values	
	#4G add on additional features 		
	#4H Another option - Plot elevation in specific band colours with topographical relief
#============================================================================================#
# 5. Adding points to a map (sites, species ect)
	#5A make spatial: 1. coordinates 2. projection
	#5B (extract) maybe you want to extract a values for points from one of the underlying raster layers too 	 
#============================================================================================#
# 6. Raster manipulations (raster calculator)
	#6A basic operations: division, multiplication, addition ect. 
	#6B (reclassify) - set values / reclassify a raster 
	#6C (make a mask) working with a mask to crop another raster
	#6D (aggregate) decrease the resolution of the raster (or summarize) 
	#6E (disaggregate) increase the resolution of the raster (probably to match another layer)
	#6F (rasterize) working between vector and raster data (many summary options)
	#6G (clip raster) - by a polygon		
#============================================================================================#
# 7. Polygon manipulation 
	#7A  make some polygons to work with through this section  
	#7B extract polygons by field attribute (repeated from above)
	#7C (crop by square extent) Create the clipping frame
	#7D (crop by another irregular polygon)
	#7E (join polygons)
	#7F (dissolve) for west_prov dissolve provinces & other features to one
	#7G (overlap) select polyons only when points overlap
#============================================================================================#
# 8. Working with projection 		
	#8A define/assign a projection 
	#8B reproject  and transform vector data
	#8C reproject and transform raster
#============================================================================================#
# 9. Add map insert (little overview map)
	
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

# (Directories) Create & set directories for project
	# Create a new folder called 'R_map_workshop'    
	print("Where do you want the workshop folder to be located?")
	setwd("C:/Users/DW/Desktop") 							# must set for your personal computer
	###############################################
	dir.create(paste(getwd(), "/R_map_workshop", sep=""))	# Create main workshop folder
	setwd(paste(getwd(),"/R_map_workshop", sep="")) 		# must set for your personal computer

###########################################
main <- getwd()												# directory object
dir.create(paste(main, "/layers", sep=""))					# create a folder for gis layers
path.layers <- paste(main, "/layers", sep="")				# directory object
dir.create(paste(main, "/figures", sep=""))					# create a final figures folder
path.fig <- paste(main, "/figures", sep="")
directories <- c("path.layers", "path.fig")
setwd(path.layers)

# (Load libraries)
	install.packages(c("dismo", "rgeos", "XML", "raster", "rgdal", "maptools", "sp")) # uncheck, remove #
	library(XML)
	library(raster)
	library(dismo)
	library(rgdal)
	library(rgeos)
	library(maptools)


	
#============================================================================================#
# 1. Get digital elevation model of Region, download as tiles and paste together 
#============================================================================================#

# Download DEM tiles & paste together into single layer
	#DOWNLOAD DEMS & UNZIP
		# DEM values from http://hydrosheds.cr.usgs.gov/datadownload.php?reqdata=3demg
		require(XML) 

			#1A list all DEM files to create & download from R
					my_tiles <- data.frame() # will become list of items to download.
					tiles <- c("n50w115_dem_grid.zip",
					"n50w120_dem_grid.zip") # would include a big list of tiles here
					for(i in 1:length(tiles)){
					# hyperlink to download data from web page through R
					temp <- paste("http://earlywarning.usgs.gov/hydrodata/sa_dem_3s_grid/NA/", tiles[i], sep="")
					temp <- data.frame(temp)
					my_tiles <- rbind(my_tiles, temp); rm(temp)
					}
					for (i in 1:length(tiles)) {
					# download files 
					download.file(as.character(my_tiles[i, 1]), tiles[i]) 
					}
					files <- list.files(pattern = "dem_grid") # make list of files
			# upzip zip files of DEM
					for (i in 1:length(files)){
					# unzip all files in dir and delete them afterwards
					sapply(list.files(pattern = "*.zip"), unzip) 
					unlink(list.files(pattern = "*.zip"))
					}

	#1B (raster) read in each raster files & assign a unique name
		library(raster)
		# roots to tiles in subfolders
		files_p = gsub(pattern = "_grid.zip", replacement = "", files)
		for (i in 1:length(files)){
			# load in each of the files in the subfolder
			rasty <- raster(paste("./", files_p[i], "/", files_p[i], "/w001001.adf", sep="")) 
			# assign a new name to each tile
			assign(paste("rast", i, sep=""), rasty)
			rm(rasty) # remove old rasty
		}
	
	#1C (Mosaic) paste individual tiles together into a single raster (should take about ~ 5 mins)
		rast1 # raster object attributes
		rast2
		#plot(rast2); plot(rast1, add=T) # basic plot
		#my_big_grid <- mosaic(rast1, rast2, fun=max) # takes too long for tutorial
		 
	#1D (writeRaster) Save raster DEM file to avoid repeating the same steps over again
		setwd(path.layers)
		writeRaster(rast2, filename="big_dem.tif", datatype='INT4S', overwrite=TRUE)
		
	

#============================================================================================#
# 2. Get shapefiles of political boarders & other features to add to Map
#============================================================================================#

# ORIGIONAL SOURCE: http://www.diva-gis.org/gdata
# Country: Canada (Roads & Administrative areas) 
# ***Warning script does automatic download of files to folder****
		# Download Canadian administrative areas & roads 
					urls <- c("http://biogeo.ucdavis.edu/data/diva/adm/CAN_adm.zip",
						"http://biogeo.ucdavis.edu/data/diva/rds/CAN_rds.zip") # the files online you want to download
							names <- c("CAN_adm.zip", "CAN_rds.zip")  
							for (i in 1:length(urls)) {
							download.file(urls[i], names[i]) 
							}
		 
		# Unzip files 
			CAN <- list.files(pattern = "CAN")
							# upzip shapefiles
							for (i in 1:length(CAN)){
							# unzip all files in dir and delete them afterwards
							sapply(list.files(pattern = "*.zip"), unzip)
							unlink(list.files(pattern = "*.zip"))
							}
		
	#2A (readOGR) - load shapefiles into R
		my_root <- path.layers
		admin0 <- readOGR(my_root,'CAN_adm0')
		admin1 <- readOGR(my_root,'CAN_adm1')
		#roads <- readOGR(my_root,'CAN_roads') 						# roads file from above
		#admin2 <- readOGR(my_root,'CAN_adm2') 						# districts 
		#admin3 <- readOGR(my_root,'CAN_adm3')						# municipalities

##################################################################	
	#2B Check to see that everything is OK with shapefiles & rasters
		# (plot) shapefiles 
			#plot(admin0)
			plot(admin1, add=T, border="red")
			#plot(admin2, add=T, border="green")
			#plot(admin3, add=T, border="blue")
			head(admin1, 15) # attributes of each polygon feature 
			class(admin1)	
		# (plot) rasters
			rast2
			projection(rast2); class(rast2)
			plot(rast2, add=T)
			# image() is just a fast fxn for quick visualization of big files 
			# image(rast2) 
			
			
		
#============================================================================================#
# 3. Add compass rose
#============================================================================================#

# quick function just for fun 
#'Compass Rose' from 'sp' library & modified by Jim Lemon @ http://r-sig-geo.2731867.n2.nabble.com/How-to-diplasy-a-compass-rose-on-a-map-td4509034.html
compassRose<-function(x,y,rot=0,cex=1) { 
		 oldcex<-par(cex=cex) 
		 mheight<-strheight("M") 
		 xylim<-par("usr") 
		 plotdim<-par("pin") 
		 xmult<-(xylim[2]-xylim[1])/(xylim[4]-xylim[3])*plotdim[2]/plotdim[1] 
		 point.angles<-seq(0,7*pi/4,by=pi/4)+pi*rot/180 
		 crspans<-rep(c(mheight*3,mheight/2),4) 
		 xpoints<-cos(point.angles)*crspans*xmult+x 
		 ypoints<-sin(point.angles)*crspans+y 
		 polygon(xpoints,ypoints) 
		 txtxpoints<-cos(point.angles[c(1,3,5,7)])*1.33*crspans[1]*xmult+x 
		 txtypoints<-sin(point.angles[c(1,3,5,7)])*1.33*crspans[1]+y 
		 text(txtxpoints,txtypoints,c("E","N","W","S")) 
		 par(oldcex) 
} 
	
	
	
#============================================================================================#
# 4. Simple Plot: Convert DEM to fancy hillshade map & plot
#============================================================================================#

directories; setwd(path.layers)
#4A (raster) load in dem
	dem <- raster('big_dem.tif')

#4B (extent) crop down & trim raster for plotting 
	# make an extent object
		my_box <- extent(-119, -115, 50.1, 54)					# ymin, ymax, xmin, xmax

#4C (crop) trim down rasters to specified extent with 'crop'
		dem <- crop(dem, my_box)								# crop down objects

#4D (terrain & hillShade) make hillshade layer (to view topographical relief)
	slope <- terrain(dem, opt='slope') 							# make slope raster
	aspect <- terrain(dem, opt='aspect')						# make aspect raster
	hill <- hillShade(slope, aspect, 40, 270) 					# make hillshade raster
	#plot(hill, col=grey(10:250/250, alpha=0.5), main='', legend=FALSE)
	
#4E crop & manipulate shapefile for plotting
	# crop by an irregular polygonshape - clip polygon by another polygon 
		head(admin1)											# e.g. subset province = Alberta or BC
		west_prov <- admin1[ admin1$NAME_1 %in% c("Alberta","British Columbia"), ] # basic R subsetting (AB or BC)
		plot(west_prov)
		
	# OR, crop by extent object
		my_box2 = as(my_box, "SpatialPolygons")					# convert extent box to shapefile (rectangle)
		proj4string(my_box2) = projection(admin1)				# assign spatial projection to extent object
		plot(my_box2, add=T, col="red", border="green")
		# temporarily convert from spatial polygon df to spatial polygon
		admin1b <- SpatialPolygons(admin1@polygons,proj4string=admin1@proj4string) 
		clipy <- gIntersection(west_prov, my_box2, byid=TRUE, drop_not_poly=TRUE)
		plot(clipy, add=T, border="blue")
		
#4F plot at continuous colour intervals for elevation values
	plot(hill, col=grey(0:250/250, alpha=0.7), 
			axes=FALSE, box=FALSE, main='', legend=FALSE)		# nice to show topography
	plot(dem, add=T, col=terrain.colors(25, alpha=0.5),			
			main='', legend=FALSE, axes=FALSE, box=FALSE)		# ideally this would be some other layer of interest
			
#4G add on additional features 
		plot(clipy, add=T, border="blue", lwd=1.5)
		# add in compass rose
		#compassRose(longitude, latitude) 
			compassRose(-115.7, 53.6, cex=0.40)  
		# plot on roads too?
			#setwd(path.layers); roads <- readOGR(my_root,'CAN_roads') 						
			#roads <- gIntersection(roads, my_box2, byid=TRUE, drop_not_poly=TRUE)
			#plot(roads, add=T, col="darkgrey", lwd=0.5);compassRose(-115.7, 53.6, cex=0.40)   
		# (legend) make fancy custom legend
			r.range <- c(minValue(dem), maxValue(dem))
			# use legend only - play with settings
			plot(dem, legend.only=TRUE, col=terrain.colors(200, alpha=0.5), legend.width=1, legend.shrink=0.75, axis.args=list(at=seq(r.range[1], r.range[2], 200),labels=seq(r.range[1], r.range[2], 200), cex.axis=0.6),legend.args=list(text='Elevation (m)', side=4, font=2, line=2.5, cex=0.8))	
				
				
###################################################################################
#(4H) Another option - Plot elevation in specific band colours with topographical relief
#	  specific break points according to key elevation ranges 
		breakpoints <- c(0,500,1000,1500,2000,2500,5000)
		colors <- c("red","white","blue","green","pink","lightgrey")
		plot(dem,breaks=breakpoints,col=colors, main='')

		
# (optional) great way to save final plot/ figure from R as pdf!	
# DIRECTORY > PDF > PAR STUFF > PLOT, PLOT PLOT > DEV OFF	
	#setwd(path.fig)
	#pdf(file="Sample_plot1.pdf", width=11, height=8.5)
	#par(mfrow=c(1,2), mar=c(0.5,1,1.5,2))
	#plot(hill, col=grey(0:250/250, alpha=0.7), axes=FALSE, box=FALSE, main='', legend=FALSE); plot(dem, add=T, col=terrain.colors(25, alpha=0.5), main='', legend=FALSE, axes=FALSE, box=FALSE); plot(clipy, add=T, border="blue", lwd=1.5); compassRose(-115.7, 53.6, cex=0.70)  ; r.range <- c(minValue(dem), maxValue(dem)); plot(dem, legend.only=TRUE, col=terrain.colors(200, alpha=0.5), legend.width=1, legend.shrink=0.75, axis.args=list(at=seq(r.range[1], r.range[2], 200),labels=seq(r.range[1], r.range[2], 200), cex.axis=0.6),legend.args=list(text='Elevation (m)', side=4, font=2, line=2.5, cex=0.8))
	#breakpoints <- c(0,500,1000,1500,2000,2500,5000); colors <- c("red","white","blue","green","pink","lightgrey"); plot(dem,breaks=breakpoints,col=colors, main='', axes=FALSE, box=FALSE)
	#dev.off()

	
	
#============================================================================================#
# 5. Adding points to a map (sites, species ect)
#============================================================================================#

	# Usually you would import your own csv here with point locations of sites
	# with latitude and longitude as columns in the data frame, but here will just 
	# make up a few points
	
		lat=c(51.754,53.598,52.152,52.482, 51.1); long=c(-117.01, -118.184, -117.7, -115.297, -118.3); name=c("bananas", "Baggers Canyon", "Milk Run", "Big Eagle", "whistle")
		sites = data.frame(name, lat, long)
		# sites = read.csv("sites.csv") # ideally would just do this
		sites
	
	#5A make spatial: 1. coordinates 2. projection
		coordinates(sites) = ~long + lat
		proj4string(sites) = CRS(projection(dem))

		plot(dem); plot(sites, pch=19, add=T)
		text(sites$long,(sites$lat+0.15),labels=sites$name, cex=0.5)

	#5B (extract) maybe you want to extract a values for points from one of the underlying raster layers too 
		els <- extract(dem, sites) 
		# add to points layer 
			# need to temporarily convert back to df from points
			class(sites); sites <- data.frame(sites); class(sites)
			sites$els <- els; sites; 
			# make spatial again
			coordinates(sites) = ~long + lat; proj4string(sites) = CRS(projection(dem))
		# plot on values
			plot(dem); plot(sites, pch=19, add=T)
			text((sites$long),(sites$lat+0.15),labels=paste(sites$name, ", ", sites$els, " m", sep=""), cex=0.5)
	
	 
	 
#============================================================================================#
# 6. Raster manipulations (raster calculator)
#============================================================================================#
	# basic tutorial of key overview functions
	# modified from: Dana Tomlin
	# source: http://rstudio-pubs-static.s3.amazonaws.com/1057_1f7e9ac569644689b7e4de78c1fece90.html
	# trim down original raster to speed up tutorial 
			tutorial_e <- extent(-116.73, -116.71, 51.845, 51.855)
			# crop with extent object
			rasty <- crop(dem, tutorial_e) 
			plot(rasty, axes=FALSE, box=FALSE) # small file/zoomed in
			rasty; summary(values(rasty)) # adding values helps with summary
			par(mfrow=c(1,1))
		
	#6A basic operations: division, multiplication, addition ect. 
		rasty_a <- rasty/100; rasty_b <- log(rasty + 0.5)
		rasty_c <- rasty/rasty # create a mask make values = 1 
		# can try using cellStats for summary statistics
		rasty_d <- (rasty-cellStats(rasty,mean))/cellStats(rasty,sd) # could even normalize ect. 		
		par(mfrow=c(2,3))
			plot(rasty, axes=FALSE, box=FALSE, main="Original"); plot(rasty_a, axes=FALSE, box=FALSE, main="a. Divide"); plot(rasty_b, axes=FALSE, box=FALSE, main="b. Transform"); plot(rasty_c, axes=FALSE, box=FALSE, main="c. Mask"); plot(rasty_d, axes=FALSE, box=FALSE, main="d. Normalized")
		
	#6B (reclassify) - set values / reclassify a raster 
		# e.g. all values <2200 become 5; all values over 2200 but less than 2600 become "NA"... and all values greater than 2600 become 13. 
		m <- c(0, 2200, 5,  2200, 2600, "NA",  2600, 5000, 13)# from 0, to 2200, become=5, .... 2200, to 2600, become=12...
		rclmat <- matrix(m, ncol=3, byrow=TRUE)
		reclassed <- reclassify(rasty, rclmat) # reclassify!
		
	#6C (make a mask) working with a mask to crop another raster
		mask <- reclassed/reclassed
		summary(values(mask))
		# all values in mask are '1' and others are NA so 
		# e.g.  2465 * NA = NA     and    2639 * 1 = 2639m
		applied_mask <- mask*rasty 
		par(mfrow=c(2,2))
			plot(rasty, axes=FALSE, box=FALSE, main="1. Original"); plot(reclassed, axes=FALSE, box=FALSE, main="2. reclassified"); plot(mask, axes=FALSE, box=FALSE, main="3. the mask"); plot(applied_mask, axes=FALSE, box=FALSE, main="4. applied mask")

	#6D (aggregate) decrease the resolution of the raster (or summarize) 
		r_agg_me = aggregate(rasty, fact=4, fun=mean, expand=TRUE) # basic mean 
		r_agg_mn = aggregate(rasty, fact=5, fun=min, expand=TRUE)   # based on minimum
		r_agg_h = aggregate(rasty, fact=c(1,10), fun=max, expand=TRUE)  # veticle only (or disproportionately)  
		par(mfrow=c(2,2)); plot(rasty, axes=FALSE, box=FALSE, main="Original"); plot(r_agg_me, axes=FALSE, box=FALSE, main="Agg mean"); plot(r_agg_mn, axes=FALSE, box=FALSE, main="Agg min"); plot(r_agg_h, axes=FALSE, box=FALSE, main="Vert only")

	#6E (disaggregate) increase the resolution of the raster (probably to match another layer)
		r_dis <- disaggregate(rasty, fact=c(10, 10))
		# just have to add some values to see dif 
		test <- rasty; values(test)= seq(1:ncell(test))
		values(r_dis)= seq(1:ncell(r_dis))
		par(mfrow=c(1,2)); plot(test, axes=FALSE, box=FALSE, main="Original"); plot(r_dis, axes=FALSE, box=FALSE, main="disaggregate")
			# not really what disaggregate actually does (no interpolation) but just best way to visualize it. 
	 
	#6F (rasterize) working between vector and raster data (many summary options)
		# make random points for this exercise 
			x <- runif(150, -116.73, -116.71); y <- runif(150, 51.845, 51.855); name <- seq(1, 150, by=1); rand_pts = data.frame(name, y, x)
			head(rand_pts)
				# make spatial 
				# make spatial: 1. coordinates 2. projection
					coordinates(rand_pts) = ~x + y
					proj4string(rand_pts) = CRS(projection(dem))
					par(mfrow=c(1,1)); plot(rasty); plot(rand_pts, add=T)
			rast_pts <- rasterize(rand_pts, rasty, name, fun=sum) 			# count of points in cells 
			rast_ptsb <- rast_pts/rast_pts 									# presence or absence of point in cell
			rast_ptsc <- rast_ptsb*rasty									# use origional values (also possible with fun arguments
			# see the 'fun' arguments list for 'rasterize' & be creative
			par(mfrow=c(1,3)); plot(rast_pts, axes=FALSE, box=FALSE, main="sum of pts"); plot(rast_ptsb, axes=FALSE, box=FALSE, main="a mask"); plot(rast_ptsc, axes=FALSE, box=FALSE, main="orig el vals"); 
			points(rand_pts) # add points back in 
	
	#6G (clip raster) - by a polygon
				# create a loppy star to cut raster with 
				x=c(-116.725, -116.723,	-116.726, -116.7226, -116.7209,	-116.7194, -116.71588, -116.71864,-116.71680,-116.7213, -116.725)
				y=c(51.846, 51.848, 51.850, 51.850, 51.853, 51.850, 51.850, 51.847, 51.8453, 51.8470, 51.846)
				star <- cbind(x, y); p = Polygon(star)
				ps = Polygons(list(p),1)
				loppy_star = SpatialPolygons(list(ps)); plot(rasty); plot(loppy_star, add=T); class(loppy_star)
			# should be a mask
			star_clip <- rasterize(loppy_star, rasty, fun='sum'); plot(star_clip)
			# then get values of raster
			star_clip <- star_clip*rasty; plot(star_clip)		
	
			
		
#============================================================================================#
# 7. Polygon manipulation 
#============================================================================================#
# best to look at packages 'maptools' and 'rgeos'
	par(mfrow=c(1,1))
	
	#7A  make some polygons to work with through this section  
		par(mfrow=c(1,1)); plot(west_prov) # see above
		points(sites)
		# XL buffer around points will act as next polygon
		s_buff <- circles(sites, d=280000, lonlat=TRUE); plot(s_buff, add=T, border="red")
		library(rgeos)
		s_buff <- gUnaryUnion(s_buff@polygons) # make s_buff polygon
		s_buff <- SpatialPolygonsDataFrame(s_buff,data=as.data.frame("buffer")) # and then spatial polygon dataframe
		projection(s_buff) <- projection(west_prov)
		# basic layout

	#7B extract polygons by field attribute (repeated from above)
		#head(admin1)											# e.g. subset province = Alberta or BC
		#west_prov <- admin1[ admin1$NAME_1 %in% c("Alberta","British Columbia"), ] # basic R subsetting (AB or BC)
		#plot(west_prov)
		
	#7C (crop by square extent) Create the clipping frame
		CP <- as(extent(-123.8, -112.1, 46.5, 54), "SpatialPolygons")
		proj4string(CP) <- CRS(proj4string(west_prov))
		plot(west_prov); plot(CP, add=T, border="green")
		clip_A <- gIntersection(west_prov, CP, byid=TRUE)
		plot(clip_A)
		clip_B <- gIntersection(s_buff, CP, byid=TRUE)
		plot(clip_B, add=T, border="red")
	#7D (crop by another irregular polygon)
		clip_C <- gIntersection(west_prov, s_buff, byid=TRUE)

	#7E (join polygons)
		join_obj <- gUnion(clip_A, clip_B, byid=TRUE)
		plot(join_obj)
		# getting messy

	#7F (dissolve) for west_prov dissolve provinces & other features to one
		diss_att <- gUnaryUnion(join_obj) # see the id option to select an attribute field
		plot(clip_A, border="red", lty=3, lwd=.7)	
		plot(clip_B, border="blue", add=T, lwd=2)	
		plot(diss_att, lwd=2) 							
			# consider here as to whether or not you want to keep objects as a
			# 'spatialpolygondataframe' or just if 'spatialpolygon' will still work?
			# tends to be more difficult to find appropriate functions for spdf
		
	#7G (overlap) select polyons only when points overlap
		# not included in this tutorial, but great examples here: https://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon
		# would use the 'over' command 
		# over('spatialpolygon', 'df of points', fn = sum) 						
		# then would select by polygon attribute in df ( > 0) ect. 
		
###### review plots
	par(mfrow=c(2,3))
		plot(admin1, main="orig"); plot(CP, add=T, border="blue", lwd=0.5)
		plot(west_prov, main="by attribute"); plot(s_buff, add=T, border="red", lwd=0.5); plot(CP, add=T, border="blue", lwd=0.5)
		plot(clip_A, main="clip by frame"); plot(clip_B, add=T, border="red", lwd=0.5)
		plot(join_obj, main="joined")
		plot(diss_att, main="dissolved")

	
	
#============================================================================================#
# 8. Working with projections
#============================================================================================#		
	# go to the layers folder, loop up the '.PRJ' files, open in 
	# a text editor program (check projection). 
	
	# check files projection
	projection(dem)
	projection(admin1)
		# notice both WGS84 (with some other minor differences)
		# projections are a text string, subtle differences above don't really matter here 
		# don't confuse 'projections' with 'datums'
		# building & defining your own projections (see)
		# https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf
		# try using ERSI's projection engine or some other tool if really having trouble here. 
		
	#Moat common / most useful projections 
		# Unprojected/geographic
			prj.wgs = "+proj=longlat +ellps=WGS84"
		# Albers equal area 
			prj.aea = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0"
		# Lambert conformal conic 
			prj.lcc = "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
		# Azimuthal Equidistant 
			prj.not_sure <- "+proj=aeqd  +R=6371000 +lat_0=51 +lon_0=7" #????? not sure 
	
##################################
	#8A define/assign a projection 
		# r needs to know what projection something is in initially before we can reproject & transform. 
		projection(west_prov)
		projection(sites)
		# already assigned
			# if we just assign the projection as lcc the points are not transformed 
			# they are just reassigned ~ THIS IS INCORRECT!
			projection(sites) <- prj.lcc # this would be incorrect (only reassigned)
			projection(sites)
			projection(sites) <- prj.wgs # fixed, assigned to wgs
			projection(sites)
		
	#8B reproject  and transform vector data
		# plot origional 
		par(mfrow=c(2,2))
		plot(west_prov); plot(sites, add=T)
		
		# convert to lcc
		west_prov.lcc <- spTransform(west_prov, CRS(prj.lcc))
		sites.lcc <- spTransform(sites, CRS(prj.lcc))
		plot(west_prov.lcc); plot(sites.lcc, add=T)

		# convert to aea
		west_prov.aea <- spTransform(west_prov, CRS(prj.aea))
		sites.aea <- spTransform(sites, CRS(prj.aea))
		plot(west_prov.aea); plot(sites.aea, add=T)

		# convert to dist(?)
		west_prov.not_sure <- spTransform(west_prov, CRS(prj.not_sure))
		sites.not_sure <- spTransform(sites, CRS(prj.not_sure))
		plot(west_prov.not_sure); plot(sites.not_sure, add=T)
	
	#8C reproject and transform raster
		# should avoid when possible
		# involves resampling across cells
		?projectRaster # but see
	
	
#============================================================================================#
# 9. Add map insert (overview map)
#============================================================================================#

# Might want to add in an insert/overview map to show where your study took place. 
# We want to plot this map insert at a larger scale over top of our original map 
# and have it show the extent of our study area with a box. 
# We will accomplish this by tweeking the par() parameters, plot margins, to create overlapping 
# maps. Then we will plot in our overview map and add a small box of the study area 
# using two different extent objects.  


# Adding a map insert is a bit strange 
# We will tweek the par() parameters to create overlaping maps at different scales. 


par(mai = c(1,1,1,1))
plot(clip_A, bty = "l")

par(mai = c(3.5,3.5,1,1), bg="yellow")
par(new = TRUE) # keep this in 
plot(admin1)
plot(clip_A, border="red", add=T)

#
##
###
####
#####
######
#####
####
###
##
#
# the end





	
	
#============================================================================================#
# XX. extra useful functions 
#============================================================================================#
