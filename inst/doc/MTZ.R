## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- ech=TRUE----------------------------------------------------------------
library(cry)
library(ggplot2)

## ----echo=TRUE----------------------------------------------------------------
datadir <- system.file("extdata",package="cry")
all_files <- list.files(datadir)
print(all_files)

## ----echo=TRUE----------------------------------------------------------------
filename <- file.path(datadir,"1dei_phases.mtz")
objMTZ <- readMTZ(filename)
class(objMTZ)
names(objMTZ)

## ----echo=TRUE----------------------------------------------------------------
hdr <- objMTZ$header
class(hdr)
names(hdr)

## ----echo=TRUE----------------------------------------------------------------
print(objMTZ$header$PROJECT)
print(objMTZ$header$CRYSTAL)
print(objMTZ$header$DATASET)

## ----echo=TRUE----------------------------------------------------------------
objMTZ$header$COLUMN[,1]
objMTZ$reflections[1:5,]

## ----echo=TRUE----------------------------------------------------------------
# List the different values of H
unique(objMTZ$reflections$H)

# Select all reflections with H=1
idx <- which(objMTZ$reflections$H == 1)
length(idx)  # 373 reflections have H=0

# Save these reflections in a different object
refs <- objMTZ$reflections[idx,]

# Show the first 10 selected reflections
refs[1:10,]

# Find out the range of FP/sigFP for the selected reflections
range(refs$FP/refs$SIGFP)

## ----echo=TRUE,fig.height=5,fig.width=8---------------------------------------
idx <- which(objMTZ$reflections$FP/objMTZ$reflections$SIGFP >= 1)
length(idx) # 8377 reflections have FP/sigFP >= 1

# The reflections can be extracted
refs <- objMTZ$reflections[idx,]

# Histogram of I/sigI with ggplot2
ggplot(data = refs, aes(FP/SIGFP)) +
  geom_histogram(fill = "#5494AB", colour = "#125CD2", binwidth = 3) +
  xlab("FP/sigFP") + ylab("Number of Reflections") +
  labs(title = "Signal-to-noise ratio FP/sigFP greater than 1") +
  theme_cry()

## ----echo=TRUE----------------------------------------------------------------
# Extract cell parameters from header
cpars <- objMTZ$header$CELL
print(cpars)
a <- cpars[1]; b <- cpars[2]; c <- cpars[3]
aa <- cpars[4]; bb <- cpars[5]; cc <- cpars[6]

# Resolution of reflection (1,0,0) (0,1,0) (0,0,1)
hkl_to_reso(1,0,0,a,b,c,aa,bb,cc)
hkl_to_reso(0,1,0,a,b,c,aa,bb,cc)
hkl_to_reso(0,0,1,a,b,c,aa,bb,cc)

# Reflections with higher Miller indices have higher resolutions
hkl_to_reso(10,0,0,a,b,c,aa,bb,cc)
hkl_to_reso(10,0,-20,a,b,c,aa,bb,cc)

## ----echo=TRUE----------------------------------------------------------------
resos <- hkl_to_reso(objMTZ$reflections$H,
                     objMTZ$reflections$K,
                     objMTZ$reflections$L,
                     a,b,c,aa,bb,cc)

# Resolution range for all data
range(resos)

# Select data with resolution between 5 and 9 angstroms
idx <- which(resos >= 5 & resos <= 9)
length(idx)  # Only 314 reflections

## ----echo=TRUE----------------------------------------------------------------
# Copy original data to 2 separate lists
refs <- objMTZ$reflections
hdr <- objMTZ$header
length(refs[,1])  # Number of reflections before cut
print(hdr$NCOL[2])

# Cut data to 5 angstroms resolution
# (see previous chunk of code for resos)
idx <- which(resos >= 5)
refs <- refs[idx,]
length(refs[,1])  # Now there are only 333 observations

# Modify specific parts of header
hdr$NCOL[2] <- length(refs[,1])            # Number of reflections
hdr$RESO <- sort(1/range(resos[idx])^2)    # Resolution range
hdr$DATASET[2,2] <- "Data cut to 5A"  # Dataset
obsmin <- apply(refs,2,min,na.rm=TRUE)
obsmax <- apply(refs,2,max,na.rm=TRUE)
hdr$COLUMN[,3] <- obsmin                   # COLUMN
hdr$COLUMN[,4] <- obsmax

## ----echo=TRUE----------------------------------------------------------------
# Original COLSRC
print(hdr$COLSRC)

# Change date and time
hdr <- change_COLSRC(hdr)

# New COLSRC
print(hdr$COLSRC)

## ----echo=TRUE----------------------------------------------------------------
# Temporary directory for output
tdir <- tempdir()
fname <- file.path(tdir,"new.mtz")

# Write changed data to the new MTZ file
writeMTZ(refs,hdr,fname,title="New truncated dataset")

## ----echo=TRUE----------------------------------------------------------------
# Read data from "new.mtz"
newMTZ <- readMTZ(fname,TRUE)

