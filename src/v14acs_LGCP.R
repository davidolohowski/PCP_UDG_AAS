library(tidyverse)
library(spatial)
library(raster)
library(zoo)
library(spatstat)
library(excursions)
library(INLA)
library(inlabru)
library(spatialEco)
library(RColorBrewer)
library(rgeos)
library(R.utils)
library(doParallel)
library(foreach)
library(viridis)
library(tikzDevice)

v_acs <- read_csv('data/v14acs_GC.csv')

X <- c(1, 4097, 4215, 217)
Y <- c(1, 100, 4242, 4042)
region <- Polygon(cbind(X,Y))
region <- SpatialPolygons(list(Polygons(list(region),'region')))
region <- SpatialPolygonsDataFrame(region, data.frame(id = region@polygons[[1]]@ID, row.names = region@polygons[[1]]@ID))
plot(region)
points(v_acs[,3:4])

grid <- makegrid(region, n = 20000)
grid <- SpatialPoints(grid, proj4string = CRS(proj4string(region)))
grid <- crop(grid, region)

gal1 <- as.data.frame(grid)

names(gal1) <- c('x', 'y')

gal1$z <- ((gal1$x - 3368)*cos(10*pi/180) - (gal1$y-2716)*sin(10*pi/180))^2 + ((gal1$x - 3368)*sin(10*pi/180) + (gal1$y-2716)*cos(10*pi/180))^2/1.25^2

ggplot(gal1, aes(x,y)) + geom_raster(aes(fill = exp(-z/570^2))) + coord_fixed() + geom_point(data = v_acs, aes(x,y)) + geom_contour(aes(z =exp(-z/570^2)))

coordinates(gal1) <- ~x+y

gridded(gal1) <- T

f.gal1 <- function(x, y) {
  spp <- SpatialPoints(data.frame(x = x, y = y), proj4string = fm_sp_get_crs(gal1))
  proj4string(spp) <- fm_sp_get_crs(gal1)
  
  v <- over(spp, gal1)
  if (any(is.na(v$z))) {
    v$z <- inlabru:::bru_fill_missing(gal1, spp, v$z)
  }
  return(v$z)
}

spv <- SpatialPoints(v_acs[,3:4])

v_mesh <- inla.mesh.2d(loc.domain = region@polygons[[1]]@Polygons[[1]]@coords, max.edge = c(100, 1000))

plot(v_mesh)

v_spde <- inla.spde2.pcmatern(mesh = v_mesh, alpha = 2,
                              prior.range = c(200, 0.5),
                              prior.sigma = c(1.5, 0.5))

alpha_fun_1 <- function(u){
  qexp(pnorm(u), rate = 1)
}

R_fun_1 <- function(u){
  qlnorm(pnorm(u), meanlog = 6.5, sdlog = 0.25)
}

cmp <- coordinates ~ gal1(f.gal1(x,y), model = "offset") +
  beta1(1, model = 'linear') + 
  R_internal_1(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
  alpha_internal_1(1, model="linear", mean.linear=0, prec.linear=1) +
  f(coordinates, model = v_spde) + Intercept(1)

form <- coordinates ~ Intercept + 
  beta1*exp(-(gal1/R_fun_1(R_internal_1)^2)^alpha_fun_1(alpha_internal_1)) + f

spv <- SpatialPointsDataFrame(spv, data = data.frame(id = 1:nrow(as.data.frame(spv))))

system.time(
  fit <- bru(components = cmp, family = 'cp', data = spv,
                                    samplers = region, domain = list(coordinates = v_mesh), 
                                    formula = form,
             options = list(bru_initial=list(gal = 1),
                            control.inla=list(int.strategy = 'ccd'),
                            control.compute=list(config=TRUE),
                            control.predictor = list(compute = TRUE)))
)

summary(fit)

X <- c(1, 4300, 4300, 1)
Y <- c(1, 1, 4300, 4300)
S <- Polygon(cbind(X,Y))
S <- SpatialPolygons(list(Polygons(list(S),'region')))
S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))

set.seed(1246)
lambda <- predict(fit, pixels(v_mesh, nx = 300, ny = 300, mask = S), ~ exp(f), seed = 1246)
lambda <- as.data.frame(lambda)

tikz(file = "v14acs_intensity_Us.tex", standAlone=T, width = 5, height = 4.5*5/6)
ggplot(lambda, aes(x/4300,y/4300)) + 
  geom_raster(aes(fill = mean)) + scale_fill_viridis_c(name = '$\\exp(\\mathcal{U}(s)))$') + 
  geom_point(data = v_acs/4300, aes(x,y), size = 0.1) +
  xlab('x') + ylab('y') + coord_fixed() + theme_minimal()
dev.off()
system('pdflatex v14acs_intensity_Us.tex')
