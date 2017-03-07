# Learning stuff...

scenario <- irace::defaultScenario()
scenario$targetRunner <- "runMOEADr"
scenario$trainInstancesDir <- ""
scenario$seed <- 20040408
scenario$parallel <- 3

parameters <- readParameters("parameters.txt")

targetRunner <- function(experiment, scenario){
  conf <- experiment$configuration

  # Assemble moead input lists
  ## 1. Problem
  # problem   <- list(name       = experiment$instance,
  #                   xmin       = rep(0, ???),
  #                   xmax       = rep(1, ???),
  #                   m          = ???)

  ## 2. Decomp
  decomp <- list(name = conf$decomp.name)
  if(decomp$name == "sld") decomp$H <- conf$H
  if(decomp$name == "Uniform") decomp$N <- conf$N

  ## 3. Neighbors
  neighbors <- list(name    = conf$neighbor.name,
                    T       = conf$T,
                    delta.p = conf$delta.p)

  ## 4. Aggfun
  aggfun <- list(name = conf$aggfun.name)
  if (aggfun$name == "PBI") aggfun$theta <- conf$aggfun.theta

  ## 5. Update
  update <- list(name       = conf$update.name,
                 UseArchive = conf$UseArchive)
  if (update$name != "standard") update$nr <- conf$nr
  if (update$name == "best") update$Tr <- conf$Tr

  ## 6. Scaling
  scaling <- list(name = conf$scaling.name)

  ## 7. Constraint
  constraint<- list(name = "vbr",
                    type = "ts")

  ## 8. Stop criterion
  stopcrit  <- list(list(name    = "maxeval",
                         maxeval = 10000))

  ## 9. Echoing
  showpars  <- list(show.iters = "dots",
                    showevery  = 25)

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
    }
  }

  ## 11. Seed
  seed <- conf$seed

  out <- moead(problem, decomp,  aggfun, neighbors, variation, update,
               constraint, scaling, stopcrit, showpars, seed)
}
