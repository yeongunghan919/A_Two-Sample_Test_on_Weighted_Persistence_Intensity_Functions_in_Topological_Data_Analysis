library(TDA)
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
# generate data
sig=c(0.01,0.02,0.03,0.04,0.05,0.0) # noise level, the last noise is not noise 
npc=50 # number of point clouds per set
nset=100 # number of sets


set.seed(47) # seed

#function of homo Poisson point process on Torus 
Poisson_homo_torus = function(lamb,R,r){
  total_intensity = 4*(pi**2)*R*r
  n = rpois(1,lamb*total_intensity)
  return(torusUnif(n,a=r,c=R))
}
torus_lamb=5
torus_R=2
torus_r=1
# generate datax
torus_pc=torus_pd=list()
for (ii in 1:length(sig)){
  torus_pc[[ii]]=torus_pd[[ii]]=list()
  for (jj in 1:(npc*nset)) {
    # generate point cloud
    not_noise_data=Poisson_homo_torus(lamb=torus_lamb,R=torus_R,r=torus_r)
    torus_pc[[ii]][[jj]]=not_noise_data+matrix(rnorm(prod(dim(not_noise_data)),0,sig[ii]),dim(not_noise_data)[1],dim(not_noise_data)[2])
    # compute persistence diagram
    torus_diag=ripsDiag(torus_pc[[ii]][[jj]],maxdimension=1,maxscale=3)
    df= as.data.frame(unclass(torus_diag$diagram))
    colnames(df) = NULL 
    torus_pd[[ii]][[jj]]=unlist(df,use.names=FALSE) #we transform each diagram to data.frame for transforming json 
  }
}



saveRDS(torus_pd, "Rdata/npc=50_nset=100torus_pd.Rdata")
library(jsonlite)

# JSON
writeLines(toJSON(torus_pd),"json/npc=50_nset=100torus_pd.json")


######################################################################
library(TDA)
library(jsonlite)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Radius perturbations; 0 corresponds to the null model
radius_shift <- c(0.01,0.012,0.014,0.016,0.0)

npc  <- 100L  # number of point clouds in each simulation
nset <- 500L  # number of Monte Carlo simulations

set.seed(47)

# Homogeneous Poisson point process on a torus
Poisson_homo_torus <- function(lamb, R, r) {
  surface_area <- 4 * pi^2 * R * r
  n_points <- rpois(1L, lambda = lamb * surface_area)
  
  torusUnif(n = n_points, a = r, c = R)
}

torus_lamb <- 5
torus_R    <- 2
torus_r    <- 1

# Structure:
# torus_pd[[radius shift]][[simulation]][[point cloud]]
torus_pd <- vector("list", length(radius_shift))

for (ii in seq_along(radius_shift)) {
  
  current_R <- torus_R + radius_shift[ii]
  torus_pd[[ii]] <- vector("list", nset)
  
  for (bb in seq_len(nset)) {
    
    simulation_pd <- vector("list", npc)
    
    for (jj in seq_len(npc)) {
      
      # No Gaussian noise: change only the major radius
      point_cloud <- Poisson_homo_torus(
        lamb = torus_lamb,
        R = current_R,
        r = torus_r
      )
      
      torus_diag <- ripsDiag(
        X = point_cloud,
        maxdimension = 1,
        maxscale = 3
      )
      
      # Same serialization format as the original code
      simulation_pd[[jj]] <- unname(
        as.numeric(torus_diag$diagram)
      )
    }
    
    torus_pd[[ii]][[bb]] <- simulation_pd
    
    if (bb %% 25L == 0L) {
      message(
        sprintf(
          "R = %.2f: %d/%d simulations completed",
          current_R, bb, nset
        )
      )
    }
  }
}

dir.create("Rdata", recursive = TRUE, showWarnings = FALSE)
dir.create("json", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  torus_pd,
  file = "Rdata/npc=100_nset=500_torus_radius_pd.rds"
)




##########################
######################################################################
library(TDA)
library(jsonlite)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Radius perturbations; 0 corresponds to the null model
radius_shift <-c(0.01,0.016,0.026,0.036,0.0)

npc  <- 100L  # number of point clouds in each simulation
nset <- 500L  # number of Monte Carlo simulations

set.seed(47)

# Homogeneous Poisson point process on a torus
Poisson_homo_torus <- function(lamb, R, r) {
  surface_area <- 4 * pi^2 * R * r
  n_points <- rpois(1L, lambda = lamb * surface_area)
  
  torusUnif(n = n_points, a = r, c = R)
}

torus_lamb <- 5
torus_R    <- 2
torus_r    <- 1

# Structure:
# torus_pd[[radius shift]][[simulation]][[point cloud]]
torus_pd <- vector("list", length(radius_shift))

for (ii in seq_along(radius_shift)) {
  
  current_R <- torus_R + radius_shift[ii]
  torus_pd[[ii]] <- vector("list", nset)
  
  for (bb in seq_len(nset)) {
    
    simulation_pd <- vector("list", npc)
    
    for (jj in seq_len(npc)) {
      
      # No Gaussian noise: change only the major radius
      point_cloud <- Poisson_homo_torus(
        lamb = torus_lamb,
        R = current_R,
        r = torus_r
      )
      
      torus_diag <- ripsDiag(
        X = point_cloud,
        maxdimension = 1,
        maxscale = 3
      )
      
      # Same serialization format as the original code
      simulation_pd[[jj]] <- unname(
        as.numeric(torus_diag$diagram)
      )
    }
    
    torus_pd[[ii]][[bb]] <- simulation_pd
    
    if (bb %% 25L == 0L) {
      message(
        sprintf(
          "R = %.3f: %d/%d simulations completed",
          current_R, bb, nset
        )
      )
    }
  }
}

dir.create("Rdata", recursive = TRUE, showWarnings = FALSE)


saveRDS(
  torus_pd,
  file = "Rdata/npc=100_nset=500_torus_radius_00260036.rds"
)
###############
previous <- readRDS("Rdata/npc=100_nset=500_torus_radius_pd.rds")
torus_pd_new <- readRDS("Rdata/npc=100_nset=500_torus_radius_00260036.rds")

stopifnot(length(previous) == 5L, length(torus_pd_new) == 2L)

# order: 0.01, 0.016, 0.026, 0.036, 0.0
torus_pd_final <- c(
  previous[c(1L, 4L)],
  torus_pd_new,
  previous[5L]
)

stopifnot(all(lengths(torus_pd_final) == 500L))

saveRDS(
  torus_pd_final,
  file = "Rdata/npc=100_nset=500_torus_radius_combined.rds"
)

#################Adding Gaussian Noise ################

library(TDA)
library(jsonlite)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Radius perturbations; 0 corresponds to the null model
radius_shift <- c(0.01, 0.016, 0.026)#c(0.01, 0.016, 0.026, 0.036, 0.0)

npc  <- 50L
nset <- 500L
sigma <- 0.05

set.seed(47)

# Homogeneous Poisson point process on a torus
Poisson_homo_torus <- function(lamb, R, r) {
  surface_area <- 4 * pi^2 * R * r
  n_points <- rpois(1L, lambda = lamb * surface_area)
  
  torusUnif(n = n_points, a = r, c = R)
}

torus_lamb <- 5
torus_R    <- 2
torus_r    <- 1

# Structure:
# torus_pd[[radius shift]][[simulation]][[point cloud]]
torus_pd <- vector("list", length(radius_shift))

for (ii in seq_along(radius_shift)) {
  
  current_R <- torus_R + radius_shift[ii]
  torus_pd[[ii]] <- vector("list", nset)
  
  for (bb in seq_len(nset)) {
    
    simulation_pd <- vector("list", npc)
    
    for (jj in seq_len(npc)) {
      
      point_cloud <- Poisson_homo_torus(
        lamb = torus_lamb,
        R = current_R,
        r = torus_r
      )
      
      # Add independent N(0, sigma^2) noise to every coordinate
      noise <- matrix(
        rnorm(length(point_cloud), mean = 0, sd = sigma),
        nrow = nrow(point_cloud),
        ncol = ncol(point_cloud)
      )
      point_cloud <- point_cloud + noise
      
      torus_diag <- ripsDiag(
        X = point_cloud,
        maxdimension = 1,
        maxscale = 3
      )
      
      simulation_pd[[jj]] <- unname(
        as.numeric(torus_diag$diagram)
      )
    }
    
    torus_pd[[ii]][[bb]] <- simulation_pd
    
    if (bb %% 25L == 0L) {
      message(
        sprintf(
          "R = %.3f: %d/%d simulations completed",
          current_R, bb, nset
        )
      )
    }
  }
}

dir.create("Rdata", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  torus_pd,
  file = "Rdata/npc=100_nset=500_torus_radius_withnoise.rds"
)


#################Adding Gaussian Noise with 0.0 shift ################

library(TDA)
library(jsonlite)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Radius perturbations; 0 corresponds to the null model
radius_shift <- c(0.0)#c(0.01, 0.016, 0.026, 0.036, 0.0)

npc  <- 50L
nset <- 500L
sigma <- 0.05

set.seed(47)

# Homogeneous Poisson point process on a torus
Poisson_homo_torus <- function(lamb, R, r) {
  surface_area <- 4 * pi^2 * R * r
  n_points <- rpois(1L, lambda = lamb * surface_area)
  
  torusUnif(n = n_points, a = r, c = R)
}

torus_lamb <- 5
torus_R    <- 2
torus_r    <- 1

# Structure:
# torus_pd[[radius shift]][[simulation]][[point cloud]]
torus_pd <- vector("list", length(radius_shift))

for (ii in seq_along(radius_shift)) {
  
  current_R <- torus_R + radius_shift[ii]
  torus_pd[[ii]] <- vector("list", nset)
  
  for (bb in seq_len(nset)) {
    
    simulation_pd <- vector("list", npc)
    
    for (jj in seq_len(npc)) {
      
      point_cloud <- Poisson_homo_torus(
        lamb = torus_lamb,
        R = current_R,
        r = torus_r
      )
      
      # Add independent N(0, sigma^2) noise to every coordinate
      noise <- matrix(
        rnorm(length(point_cloud), mean = 0, sd = sigma),
        nrow = nrow(point_cloud),
        ncol = ncol(point_cloud)
      )
      point_cloud <- point_cloud + noise
      
      torus_diag <- ripsDiag(
        X = point_cloud,
        maxdimension = 1,
        maxscale = 3
      )
      
      simulation_pd[[jj]] <- unname(
        as.numeric(torus_diag$diagram)
      )
    }
    
    torus_pd[[ii]][[bb]] <- simulation_pd
    
    if (bb %% 25L == 0L) {
      message(
        sprintf(
          "R = %.3f: %d/%d simulations completed",
          current_R, bb, nset
        )
      )
    }
  }
}

dir.create("Rdata", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  torus_pd,
  file = "Rdata/npc=100_nset=500_torus_radius_withnoise_00shift.rds"
)

prethree=readRDS("Rdata/npc=100_nset=500_torus_radius_withnoise.rds")
fourth=readRDS("Rdata/npc=100_nset=500_torus_radius_withnoise_0036.rds")
last=readRDS("Rdata/npc=100_nset=500_torus_radius_withnoise_00shift.rds")

# order: 0.01, 0.016, 0.026, 0.036, 0.0
torus_pd_noise_final <- c(
  prethree[1:3],
  fourth[1L],
  last[1L]
)

stopifnot(all(lengths(torus_pd_noise_final) == 500L))
length(torus_pd_noise_final)
saveRDS(
  torus_pd_noise_final,
  file = "Rdata/npc=100_nset=500_torus_radius_with_noise.rds"
)
