suppressPackageStartupMessages({
  library(batchtools)
  library(dplyr)
  library(glue)
})

setwd('R/TL/sieve')
source('codes/simulation.R')

# paths -------------------------------------------------------------------


projdir  <- normalizePath(".", mustWork = TRUE)
code_dir <- file.path(projdir, "codes")
rdsdir <- file.path(projdir, "rds")
dir.create(rdsdir, showWarnings = FALSE, recursive = TRUE)


# job table ---------------------------------------------------------------


jobs <- bind_rows(
  tibble(tag = glue("p1={p1v}"), 
         K = 1,  np = np0, nt = nt0, p1 = p1v, p2 = p20),
  tibble(tag = glue("p2={p2v}"), 
         K = 1,  np = np0, nt = nt0, p1 = p10,  p2 = p2v),
  tibble(tag = glue("K={Kv}"), 
         K = Kv, np = np0, nt = nt0, p1 = p10,  p2 = p20),
  tibble(tag = glue("nt={ntv}"), 
         K = 1,  np = np0, nt = ntv, p1 = p10,  p2 = p20),
  tibble(tag = glue("np={npv}"), 
         K = 1,  np = npv, nt = nt0, p1 = p10,  p2 = p20)
)

# workers -----------------------------------------------------------------


run_one <- function(K, np, nt, p1, p2, tag, code_dir, rdsdir, 
                    ntest, Nrep, regime0) {
  source(file.path(code_dir, "sv-cv.R"))
  source(file.path(code_dir, "functions.R"))
  source(file.path(code_dir, "simulation.R"))
  
  outfile <- file.path(rdsdir, paste0(tag, ".rds"))
  if (file.exists(outfile)) return(outfile)  # idempotent
  
  gg_arr <- main(Nrep, K, np, nt, p1, p2, ntest, regime0)
  
  tmp <- paste0(outfile, ".part")
  saveRDS(gg_arr, tmp)
  file.rename(tmp, outfile)
  outfile
}

# Registry (uses .batchtools.conf.R in project root) ----------------------


regdir <- file.path(projdir, "bt-reg")
reg <- makeRegistry(
  file.dir = regdir,
  conf.file = file.path(projdir, ".batchtools.conf.R"),
  packages = c("dplyr", "glue")
)

# Map jobs -> workers (pass only scalars + minimal shared constant --------

ids <- batchMap(
  fun = run_one,
  K = jobs$K, np = jobs$np, nt = jobs$nt, p1 = jobs$p1, p2 = jobs$p2,
  tag = jobs$tag,
  more.args = list(
    code_dir = code_dir, rdsdir = rdsdir, 
    ntest = ntest, Nrep = Nrep, regime0 = regime0
  ),
  reg = reg
)

# Submit (tweak resources here or via defaults in .batchtools.conf --------


submitJobs(
  ids,
  reg = reg,
  resources = list(
    walltime = 60 * 60 * 6, # x hours
    ntasks = 16, ncpus = 15, # Number of required cpus per task
    memory = 8192 * 4 # Memory in megabytes for each cpu
  )
)

waitForJobs(reg = reg)
print(getStatus(reg = reg))

cat("\nSaved files:\n")
print(list.files(rdsdir, pattern = "\\.rds$", full.names = TRUE))