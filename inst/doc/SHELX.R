## ----setup--------------------------------------------------------------------
library(cry)

## ----echo=TRUE----------------------------------------------------------------
datadir <- system.file("extdata",package="cry")
all_files <- list.files(datadir)
print(all_files)

## ----echo=TRUE----------------------------------------------------------------
filename <- file.path(datadir,"shelxc.log")
obj_shelxc <- readSHELXlog(filename)
class(obj_shelxc)
names(obj_shelxc)

## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(obj_shelxc, aes(1/(Res)^2, d_sig)) +
 geom_point() + geom_line() + theme_bw() +
 xlab(expression(h^2 * (ring(A)^-2))) + 
  ylab(expression(Delta*F/sig(Delta*F)))


## ----echo=TRUE----------------------------------------------------------------
filename <- file.path(datadir,"shelxd.log")
obj_shelxd <- readSHELXlog(filename)
class(obj_shelxd)
names(obj_shelxd)

## ---- fig_width = 16, fig_height= 14------------------------------------------
library(ggplot2)

ggplot(obj_shelxd, aes(CCall, CCweak)) +
  geom_point() + theme_bw() +
  xlab('CCall') + ylab('CCweak') 

## ----echo=TRUE----------------------------------------------------------------
## read the two hands log files separately
filename_i <- file.path(datadir,"shelxe_i.log")
obj_shelxe_i <- readSHELXlog(filename_i)
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
obj_shelxe_o <- readSHELXlog(filename_o)
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

