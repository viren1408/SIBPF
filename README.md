# SIBPF
The code gets steep spectral sources in a given radio image. Various options in the code will be added here soon 

## Crossmatching Subroutine ##

#Step 1 : Extracting sources from the input fits file > Get region from header > Get the NVSS and TGSS counterparts for the image region

#Step 2 : Executing the logic for cross-matching and creating files > 1. Image sourcesxTGSSxNVSS, 2.Image sourcesxTGSS, 3.Image sourcesxNVSS , 4.Unmatched sources(after crossmatching) , unmatched sources(TGSS), unmatched sources(NVSS). Create Ra dec plots of matched sources to show a degree of matching qualitatively plot can be saved for reference.

## Spectral index subroutine ##

#Step 1: Spectral Index Plots and Files > 1. Master spectral index df

#Step 2(optional): include Vizer spectral index at that location or for a data-frame to check for frequencies less than 10 GHz 
