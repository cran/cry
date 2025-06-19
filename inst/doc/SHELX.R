## ----setup--------------------------------------------------------------------
library(cry)
library(ggplot2)

## ----echo=TRUE----------------------------------------------------------------
datadir <- system.file("extdata",package="cry")
all_files <- list.files(datadir)
print(all_files)

## ----echo=TRUE----------------------------------------------------------------
filename <- file.path(datadir,"shelxc.log")
obj_shelxc <- read_SHELX_log(filename)
class(obj_shelxc)
names(obj_shelxc)

## -----------------------------------------------------------------------------
plot_SHELX(obj_shelxc, var = obj_shelxc$Chi_sq, type = "shelxc",
           title_chart = "Chis ^2") +
  theme_cry()
    

## ----echo=TRUE----------------------------------------------------------------
filename <- file.path(datadir,"shelxd.log")
obj_shelxd <- read_SHELX_log(filename)
class(obj_shelxd)
names(obj_shelxd)

## ----fig_width = 16, fig_height= 14-------------------------------------------
plot_SHELX(filename = obj_shelxd, type = "shelxd") +
  theme_cry()

## ----echo=TRUE----------------------------------------------------------------
## read the two hands log files separately
filename_i <- file.path(datadir,"shelxe_i.log")
obj_shelxe_i <- read_SHELX_log(filename_i)
class(obj_shelxe_i)
names(obj_shelxe_i)
cycle_i <- obj_shelxe_i$CYCLE
class(cycle_i)
names(cycle_i)
FOM_mapCC_i <- obj_shelxe_i$FOM_mapCC
class(FOM_mapCC_i)
names(FOM_mapCC_i)
Site1_i <- obj_shelxe_i$Site1
class(Site1_i)
names(Site1_i)
Site2_i <- obj_shelxe_i$Site2
class(Site2_i)
names(Site2_i)

filename_o <- file.path(datadir,"shelxe_o.log")
obj_shelxe_o <- read_SHELX_log(filename_o)
class(obj_shelxe_o)
names(obj_shelxe_o)
cycle_o <- obj_shelxe_o$CYCLE
class(cycle_i)
names(cycle_i)
FOM_mapCC_0 <- obj_shelxe_o$FOM_mapCC
class(FOM_mapCC_i)
names(FOM_mapCC_i)
Site1_o <- obj_shelxe_o$Site1
class(Site1_i)
names(Site1_i)
Site2_o <- obj_shelxe_o$Site2
class(Site2_i)
names(Site2_i)

## -----------------------------------------------------------------------------
plot_SHELX(filename = obj_shelxe_i, filename_e = obj_shelxe_o,
           type = "shelxe") +
  theme_cry()

