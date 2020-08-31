#! /usr/bin/env Rscript

require("parallel", quietly = TRUE)
require("snow", quietly = TRUE)
## require("Rmpi", quietly = TRUE) # mpi.quit("no")
# nTasks <- detectCores()
nTasks = 1
cl <- makeCluster(nTasks, type = "MPI")


out = clusterMap(cl = cl,
                 fun = function(x) paste(system2('hostname', stdout = TRUE),".PID",
                                        as.character(Sys.getpid()), sep = ""),
                 as.list(1:nTasks))

print(out)

stopCluster(cl)
## stopCluster.MPIcluster()
## mpi.close.Rslaves()

##clusterEvalQ(cl = cl, Rmpi::mpi.finalize()

## mpi.exit()
## Rmpi::mpi.finalize()
## mpi.close.Rslaves(dellog  =  FALSE)
# mpi.quit("no")
