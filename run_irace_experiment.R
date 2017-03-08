# This is an R-script version of README.Rmd, which is to be used when actually
# running the experiment
#
# This file was generated using knitr::purl("README.Rmd"), followed by 
# uncommenting the relevant code blocks and a suitable name-change.

#===============
### Install manuscript version of the MOEADr package
# Uncomment below to install the exact version used in the manuscript
# require(devtools)
# devtools::install_github("fcampelo/MOEADr/MOEADr@Manuscript-Version")

#===============
### Initial setup
suppressPackageStartupMessages(require(irace))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(smoof))
suppressPackageStartupMessages(require(MOEADr))
if (packageVersion("MOEADr") != "0.1.0.0") {
  stop("Wrong MOEADr version: please install the manuscript version using the code block above")
}

# Build scenario
scenario              <- irace::defaultScenario()
nc                    <- parallel::detectCores() - 1
scenario$parallel     <- nc # Number of cores to be used by irace
scenario$seed         <- 123456 # Seed for the experiment
scenario$targetRunner <- "targetRunner" # Runner function (def. below)
scenario$targetRunnerRetries <- 5 # Retries if targetRunner fails to run
scenario$maxExperiments      <- 728 # Tuning budget

# Read tunable parameter list from file
parameters <- readParameters("./Experiments/Irace tuning/parameters.txt")


#===============
### Training instances
fname   <- paste0("UF_", 1:10)
dims    <- c(20:29,
             31:39,
             41:49,
             51:60)
allfuns <- expand.grid(fname, dims)

scenario$instances <- paste0(allfuns[,1], "_", allfuns[,2])


#===============
### targetRunner function for _irace_
target.runner <- function(experiment, scenario){
  conf <- experiment$configuration
  inst <- experiment$instance

  #=============================================
  # Assemble moead input lists
  ## 1. Problem
  fdef    <- unlist(strsplit(inst, split = "_"))
  myfun   <- smoof::makeUFFunction(dimensions = as.numeric(fdef[3]),
                                   id         = as.numeric(fdef[2]))
  
  fattr   <- attr(myfun, "par.set")
  problem <- list(name       = myfun,
                  xmin       = fattr$pars$x$lower,
                  xmax       = fattr$pars$x$upper,
                  m          = attr(myfun, "n.objectives"))
  
  ##===============
  ## 2. Decomp
  decomp <- list(name = conf$decomp.name)
  if (problem$m == 2){ # <-- 2 objectives
    if(decomp$name == "sld") decomp$H <- 99 # <-- yields N = 100
    if(decomp$name == "Uniform") decomp$N <- 100
  } else { # <-- 3 objectives
    if(decomp$name == "sld") decomp$H <- 16 # <-- yields N = 153
    if(decomp$name == "Uniform") decomp$N <- 150
  }
  
  ##===============
  ## 3. Neighbors
  neighbors <- list(name    = conf$neighbor.name,
                    T       = conf$T,
                    delta.p = conf$delta.p)
  
  ##===============
  ## 4. Aggfun
  aggfun <- list(name = conf$aggfun.name)
  if (aggfun$name == "PBI") aggfun$theta <- conf$aggfun.theta
  
  ##===============
  ## 5. Update
  update <- list(name       = conf$update.name,
                 UseArchive = conf$UseArchive)
  if (update$name != "standard") update$nr <- conf$nr
  if (update$name == "best") update$Tr <- conf$Tr
  
  ##===============
  ## 6. Scaling
  scaling <- list(name = conf$scaling.name)
  
  ##===============
  ## 7. Constraint
  constraint<- list(name = "none")
  
  ##===============
  ## 8. Stop criterion
  stopcrit  <- list(list(name    = "maxeval",
                         maxeval = 300000))
  
  ##===============
  ## 9. Echoing
  showpars  <- list(show.iters = "dots",
                    showevery  = 25)
  
  ##===============
  ## 10. Variation stack
  variation <- list(list(name = conf$varop1),
                    list(name = conf$varop2),
                    list(name = conf$varop3),
                    list(name = conf$varop4))
  
  for (i in seq_along(variation)){
    if (variation[[i]]$name == "binrec") {
      variation[[i]]$rho <- get(paste0("binrec.rho", i), conf)
    }
    if (variation[[i]]$name == "diffmut") {
      variation[[i]]$basis <- get(paste0("diffmut.basis", i), conf)
      variation[[i]]$Phi   <- NULL
    }
    if (variation[[i]]$name == "polymut") {
      variation[[i]]$etam <- get(paste0("polymut.eta", i), conf)
      variation[[i]]$pm   <- get(paste0("polymut.pm", i), conf)
    }
    if (variation[[i]]$name == "sbx") {
      variation[[i]]$etax <- get(paste0("sbx.eta", i), conf)
      variation[[i]]$pc   <- get(paste0("sbx.pc", i), conf)
    }
    if (variation[[i]]$name == "localsearch") {
      variation[[i]]$type     <- conf$ls.type
      variation[[i]]$gamma.ls <- conf$gamma.ls
      variation[[i]]$trunc.x  <- TRUE
    }
  }
  
  ##===============
  ## 11. Seed
  seed <- conf$seed
  
  #=============================================
  # Run MOEA/D
  out <- moead(problem, decomp,  aggfun, neighbors, variation, update,
               constraint, scaling, stopcrit, showpars, seed)
  
  #=============================================
  # return IGD
  Yref <- as.matrix(read.table(paste0("./Experiments/Irace tuning/pf_data/",
                                      fdef[1], fdef[2], ".dat")))
  return(calcIGD(Y = out$Y, Yref = Yref))
}

## Running the experiment
irace.output <- irace::irace(scenario, parameters)
