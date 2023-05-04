# Setup session====
library(ggplot2)
library(here)
#library(tidyverse)
source(here("R","read-admb.R"))
source(here("R","plotfns.R"))


# Check on data files--------
mkdir("data")
mkdir("output")
# Run models Goal to run 4 models with increaseing data...
mn       <- c("catch_only", "add_index", "add_catch_age", "q_vary","M_est")

# Show demo of diff function...

# Compile code-----
setwd("model")
system("admb -f bf")

# Run all models ---------
# set names for IO
dn       <- paste0(mn,".dat")
fn       <- paste0(mn,".rep")
n <- length(fn);n
# Lame loop because 10 second effort to use apply failed...
for (i in 1:n)  {
  runmod <- (paste0("./bf -nox -ind ../data/",dn[i]))
  print(runmod)
  system(runmod)
  #copy result
  cp("R_report.rep", paste0("../output/",fn[i]))
}

# Create list of all model's outputs
M        <- lapply(paste0("../output/",fn), read_rep)
names(M) <- mn
names(M[[2]])

# Make a likelihood table
ft <- tab_fit(M); ft
# use flextable to put on 
flextable::flextable(ft |> format(digits=2,trim=FALSE),align=R) 
flextable::flextable(tab_fit(M)) |> flextable::autofit()
# Make a generic plot----
plot_rec(M,runs=c(1,2))
plot_fmort(M,runs=c(1:5))
plot_ind(M,runs=c(2,3,4))
names(M[[4]])
plot_rec(M,runs=c(4,5))
plot_fmort(M,runs=c(4,5))
plot_ssb(M,runs=c(1,5))
# test plotting-------
plot_fmort(M)

plot_ssb(M)
plot_ssb(M,wrap=TRUE)
plot_ssb(M,runs=c(1,2))
plot_ssb(M,runs=c(1,2,3))
plot_ssb(M,runs=c(4,5))
plot_ssb(M,runs=c(2,3))
plot_ssb(M,runs=c(1,3))
plot_ssb(M,runs=c(3,4))
plot_ssb(M,runs=c(1,4))
