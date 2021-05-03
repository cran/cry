## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- ech=TRUE----------------------------------------------------------------
library(cry)

## ----echo=TRUE----------------------------------------------------------------
datadir <- system.file("extdata",package="cry")
all_files <- list.files(datadir)
print(all_files)

## ----echo=TRUE----------------------------------------------------------------
filename <- file.path(datadir,"xds00_ascii.hkl")
objXDS1 <- readXDS_ASCII(filename)
class(objXDS1)
names(objXDS1)

## ----echo=TRUE----------------------------------------------------------------
# 'reflections' is a data frame
print(class(objXDS1$reflections))

# First 5 observations
print(objXDS1$reflections[1:5,])

## ----echo=TRUE----------------------------------------------------------------
# Load the other XDS ascii file
filename <- file.path(datadir,"6vww_xds_ascii_merged.hkl")
objXDS2 <- readXDS_ASCII(filename,message=TRUE)

# This file has been produced by XSCALE
print(objXDS2$processing_info$MERGE)

# First 5 observations
print(objXDS2$reflections[1:5,])

## ----echo=TRUE----------------------------------------------------------------
# Copy reflections data for the first and
# second dataset, to save on later typing
refs1 <- objXDS1$reflections
refs2 <- objXDS2$reflections

# Check whether the set of reflections is unique

# Dataset 1
nrefs1 <- length(refs1[,1])
print(nrefs1)
n_unique1 <- length(unique(refs1[,1:3])[,1])
print(n_unique1)

# Dataset 2
nrefs2 <- length(refs2[,1])
print(nrefs2)
n_unique2 <- length(unique(refs2[,1:3])[,1])
print(n_unique2)

## ----echo=TRUE----------------------------------------------------------------
# Space group of first dataset
print(objXDS1$header$SPACE_GROUP_NUMBER)

# Space group of second dataset
print(objXDS2$header$SPACE_GROUP_NUMBER)

# Extended Hermann-Mauguin symbol for symmetry n 163
SG <- translate_SG(163)$msg
print(SG)

## ----echo=TRUE----------------------------------------------------------------
# Names of data columns
names(refs2)

# Selection based on I/sigI > 2
idx <- which(refs2$IOBS/refs2$SIGMA.IOBS > 2)
print(length(idx))  #  21579 reflections have I/sigI > 2

## ----echo=TRUE,fig.height=3,fig.width=5---------------------------------------
# Create new I/sigI array
isi <- refs2$IOBS/refs2$SIGMA.IOBS

# Histogram 
hh <- hist(isi,breaks=30,plot=FALSE)
cuts <- cut(hh$breaks,c(-Inf,2,Inf))
plot(hh,main="I/sigI",xlab="I/sigI",col=cuts)

## ----echo=TRUE----------------------------------------------------------------
# First reflections in original data ordering
print(refs2[1:10,])
# 

## ----echo=TRUE----------------------------------------------------------------
# First create the ordering index
idx <- order(refs2$IOBS,decreasing=TRUE)

# Now display the first 10 reflections. They correspond to the
# 10 reflections with highest intensity
print(refs2[idx[1:10],])

## ----echo=TRUE----------------------------------------------------------------
# Ordering index
idx <- order(refs2$IOBS/refs2$SIGMA.IOBS,decreasing=TRUE)

# Display first 10 highest I/sigI reflections
print(refs2[idx[1:10],])

## ----echo=TRUE----------------------------------------------------------------
# Add an I/sigI column
refs2 <- cbind(refs2,data.frame(IsigI=refs2$IOBS/refs2$SIGMA.IOBS))

# Now display with the order previously obtained
print(refs2[idx[1:10],])

## ----echo=TRUE----------------------------------------------------------------
summary(refs2)

## ----echo=TRUE,fig.height=3,fig.width=5---------------------------------------
# Initial histogram
mtxt <- "I/sigI and top quartile for I"
hh <- hist(refs2$IOBS/refs2$SIGMA.IOBS,breaks=30,main=mtxt,
           xlab="I/sigI")

# Find top quartile of I
top_qt <- quantile(refs2$IOBS)[4]

# Data with absolute intensity in the top quartile
idx <- which(refs2$IOBS >= top_qt)

# Colour portion of histogram in red
hist(refs2$IOBS[idx]/refs2$SIGMA.IOBS[idx],breaks=30,col=2,add=TRUE)

## ----echo=TRUE,fig.height=3,fig.width=5---------------------------------------
# Cell parameters are needed to calculate resolutions
cpars <- objXDS2$header$UNIT_CELL_CONSTANTS

# Resolutions in angstroms. (resos has been computed previously)
resos <- hkl_to_reso(refs2$H,refs2$K,refs2$L,cpars[1],cpars[2],cpars[3],
                     cpars[4],cpars[5],cpars[6])
summary(resos)

# Select high resolution data ( <= 2 angstroms)
idx <- which(resos <= 2)

# Now draw the histogram
mtxt <- "I/sigI and top quartile for I"
hh <- hist(refs2$IOBS/refs2$SIGMA.IOBS,breaks=30,main=mtxt,
           xlab="I/sigI")
hist(refs2$IOBS[idx]/refs2$SIGMA.IOBS[idx],breaks=30,col=2,add=TRUE)

## ----echo=TRUE----------------------------------------------------------------
# Date
ddd <- as.character(Sys.Date())
stmp <- strsplit(ddd,"-")[[1]]
allMonths <- c("Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec")
objXDS2$header$DATE <- paste0(stmp[3],"-",
            allMonths[as.integer(stmp[2])],"-",stmp[1])

# Generated by
objXDS2$header$GENERATED_BY <- "CRY"

# Data to 2 angstroms resolution
idx <- which(resos >= 2)
objXDS2$reflections <- objXDS2$reflections[idx,]

# Temporary directory for output
tdir <- tempdir()
fname <- file.path(tdir,"new.hkl")

# Write changed data to the new XDS ascii file
writeXDS_ASCII(objXDS2$processing_info,objXDS2$header,
               objXDS2$reflections,fname)

## ----echo=TRUE----------------------------------------------------------------
# Read data from "new.hkl"
newXDS <- readXDS_ASCII(fname,message=TRUE)

