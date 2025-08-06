# Package test

# install.packages("pak")
pak::pak("reidsteele2/phenometrics")

# Read in data
bfly = read.csv('bfly.csv')

# Library
library(phenometrics)

# test stuff
comm = community(bfly, quants = c(0.25, 0.75))

cbs = class_by_species(comm)

cby = class_by_year(comm)

comm_mm = comm_mismatch(bfly, quants = c(0.25, 0.75), base_y = 10)

cbs_mm = comm_mm_by_spp(comm_mm)

cby_mm = comm_mm_by_year(comm_mm)
