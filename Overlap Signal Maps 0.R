# Overlap Signal Maps: Script 0
# P R Barber, Jan 2023

# read raw data from txt files, make tiff and panel.csv

#BiocManager::install("imcRtools")
library(imcRtools)
library(EBImage)

# wehre are the txt files
folder <- "raw"
raw <- readImagefromTXT(folder)

# ONLY USING THE FIRST FILE FOUND - NEEDS A FOR LOOP
image_name <- names(raw)[1]

img <- raw@listData[[1]]
channel <- sub("Di", "", dimnames(img)[[3]])

# There is Target information in the first line of the text file
con <- file(paste0(folder, "/", image_name, ".txt"), "r")
header <- readLines(con, n = 1)
close(con)
cols <- strsplit(header, "\t")

# which col has the first channel? Should only be one!
first <- grep(channel[1], cols[[1]])
last <- first + length(channel) - 1
# pull out relevant col names
cols <- cols[[1]][first:last]

# Some will have a '-' and these are labelled
named <- grep("-", cols)
names <- sub("^[0-9]{2,3}[A-Za-z]{1,2}-([a-zA-Z0-9]+)\\([a-zA-Z0-9]+\\)$", "\\1", cols)

panel <- data.frame(channel)
panel$name <- channel
panel$name[named] <- names[named]  # replace channel with the names we know
panel$keep <- 0
panel$keep[named] <- 1   # mark those with names to keep

# Save the panel file
filename <- "panel.csv"
write.csv(panel, filename, row.names = F)

# Save the channels to keep in a tiff in the img folder
dir.create("img", showWarnings = F)
filename <- paste0("img/", image_name, ".tiff")
img_keep <- img[,,named]
writeImage(img_keep, filename)

# NB these are 16bit images 0-65535, Steinbock output are 32bit (unknown range)