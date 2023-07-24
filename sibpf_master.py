# Spectral Index Based Pulsar Finder #
def Crossmatching(dir,region,file_path_image,file_path_TGSS,show_matches=False,save_xmatch_csv=False):
  import pandas as pd
  import warnings
  warnings.simplefilter(action='ignore', category=FutureWarning)
  from astroquery.vizier import Vizier
  import astropy.units as u
  import astropy.coordinates as coord
  import numpy as np
  import os
  import matplotlib.pyplot as plt
  #Input fits file (Image processed using AIPS): 
  from astropy.io import fits

  # Open the FITS file and get the header

  Image_file = fits.open(file_path_image)
  header = Image_file[0].header
  file_name, file_ext = os.path.splitext(file_path_image)

  # Open the FITS image file and get the header
  with fits.Open (file_path_image) as hdul:
    header = hdul[0].header

  # pixel scale from the header
  cdelt1 = header['CDELT1']
  cdelt2 = header['CDELT2']
  pixel_scale = np.sqrt(cdelt1**2 + cdelt2**2)

  # size of the image in pixels
  nx = header['NAXIS1']
  ny = header['NAXIS2']

  # Calculate the radius of the image in degrees
  radius_deg = (np.sqrt((nx/2)**2 + (ny/2)**2) * pixel_scale)/2 # to be understood why we are getting double value

  # Get the reference RA and Dec from the header
  ref_ra = header['CRVAL1']
  ref_dec = header['CRVAL2']

  #Creating Source Files 

#GMRT sources in the processed image (sources extracted using BDSF)

  import bdsf

  # Specify the file path and name
  file_path = str(dir)+'/'+str(region)+'.sav'

  # Open the file in write mode
  open(file_path, mode='w', newline='')
  # Load the input FITS file
  input_image = (file_path_image)
  save_file = str(dir)+'/'+str(region)+'.sav'
  img = bdsf.process_image(save_file, filename=input_image, quiet=True)
  img.write_catalog(format='csv', catalog_type='srl',clobber = True)
  img.export_image(img_format = 'fits',img_type = 'gaus_model',clobber = True)
  data = pd.read_csv(str(region)+'.pybdsm.srl',index_col = None,skiprows=5)
  observed_sources = pd.DataFrame(data)
  observed_sources.set_index("# Source_id", inplace = True)

#Creating required dataframae for PyBDSF sources with only required columns 
  observed_sources = observed_sources.iloc[:,[0,1,3,5,6]]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#extract the header information from the Pybdsf fits file to corstraint the crossmatch the sources. 
  Image_file = fits.open(str(region)+'.pybdsm_gaus_model.fits')
  header = Image_file[0].header
  bmaj = header['BMAJ']
  bmin = header['BMIN']
  theta_obs = 2*(header['BMAJ']-header['BMIN']) 
#NVSS sources in the region of the target image: 

  Vizier.ROW_LIMIT = -1
  result = Vizier.query_region(coord.SkyCoord(ra=ref_ra, dec=ref_dec,unit=(u.deg, u.deg),frame='icrs'),radius= radius_deg*u.degree,catalog=["NVSS"])
  table = result[0]
  NVSS_main = table.to_pandas()
  NVSS = NVSS_main.iloc[:,[0,1,2,5,6]] #printing only the relevent  columns 

# Convert RA and Dec from sexagesimal units to decimal degrees
  from astropy.coordinates import SkyCoord
  coords = SkyCoord(NVSS['RAJ2000'], NVSS['DEJ2000'], unit=(u.hourangle, u.deg))
  NVSS['RAJ2000'] = coords.ra.degree
  NVSS['DEJ2000'] = coords.dec.degree

  NVSS.rename(columns = {'S1.4':'Total_flux_NVSS'}, inplace = True)
  NVSS.Total_flux_NVSS = NVSS.Total_flux_NVSS*1e-3 #converting to Jy: PyBDsf gives values in Jy
  NVSS.rename(columns = {'e_S1.4':'E_flux_NVSS'}, inplace = True)
  NVSS.E_flux_NVSS = NVSS.E_flux_NVSS*1e-3 #converting to Jy

  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#TGSS sources in the region of the target image:  
  main_TGSS = pd.read_csv(file_path_TGSS,sep='\t')

  TGSS = main_TGSS.iloc[:,[0,1,3,5,6]]
  TGSS.Total_flux = TGSS.Total_flux*1e-3  #converting to Jy
  TGSS.E_Total_flux = TGSS.E_Total_flux*1e-3 #converting to Jy
# Create a boolean mask for sources within the circular region
  distances = np.sqrt((TGSS['RA'] - ref_ra)**2 + (TGSS['DEC'] - (ref_dec))**2)
  mask = distances <= radius_deg

# Extract the sources within the circular region
  TGSS = TGSS[mask]

# Reset the index of the sources DataFrame
  TGSS = TGSS.reset_index(drop=True)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# Defining the logic for the crossmatching 

#code block for cross_verification of sources for chan_1
  data_1 = []
  for k in range (len(observed_sources)):
    for i in range (len(TGSS)):
      ra_1 = observed_sources[' RA'][k]
      dec_1 = observed_sources[' DEC'][k]
      ra_2 = TGSS['RA'][i]
      dec_2 = TGSS['DEC'][i]
      theta = np.degrees(np.arccos(np.sin(np.radians(dec_1))*np.sin(np.radians(dec_2))+np.cos(np.radians(dec_1))*np.cos(np.radians(dec_2))*np.cos(np.radians(ra_1-ra_2)))) #angular_seperation
      if abs(theta) <= theta_obs:  #constraint : 2(B_maj-b_min) , taken from the header of the fits file
        Ra = TGSS['RA'][i]
        Dec = TGSS['DEC'][i]
        TGSS_id = TGSS['Source_name'][i]
        source_id_1 = observed_sources[' Isl_id'][k]
        Ra_1 = observed_sources[' RA'][k]
        Dec_1 = observed_sources[' DEC'][k]
        flux_tgss = TGSS['Total_flux'][i]
        e_flux_tgss = TGSS['E_Total_flux'][i]
        observed_E_flux = observed_sources[' E_Total_flux'][k]
        observed_Total_flux = observed_sources[' Total_flux'][k]
        data_1.append((Ra,Dec,theta,TGSS_id,source_id_1,Ra_1,Dec_1,flux_tgss,e_flux_tgss,observed_Total_flux,observed_E_flux))

#column = ('Ra','Dec','theta(deg)','NVSS_id','source_id_1','Ra_1','Dec_1','flux_nvss','e_flux_nvss','Total_flux_chan_1','e_flux_1')
  column = ('Ra','Dec','Sep(deg)','TGSS_id','source_id_Observed','Ra_observed','Dec_observed','flux_TGSS','e_flux_TGSS','observed_Total_flux','observed_E_flux')
  TGSS_match = pd.DataFrame(data_1,columns=column)

# Define relevent arrays 
  NVSS_flux_values = []
  e_NVSS_flux_values = []
  j_values =[]

####initiate loop for cross-matching####
  for i in range(len(NVSS)): 
    for j in range(len(TGSS_match)):
      ra_1 = NVSS['RAJ2000'][i]
      dec_1 = NVSS['DEJ2000'][i]
      ra_2 = TGSS_match['Ra'][j]
      dec_2 = TGSS_match['Dec'][j]
      theta = np.degrees(np.arccos(np.sin(np.radians(dec_1))*np.sin(np.radians(dec_2))+np.cos(np.radians(dec_1))*np.cos(np.radians(dec_2))*np.cos(np.radians(ra_1-ra_2))))
      if abs(theta) <= 0.0014 :  #constraint : 2(B_maj-b_min)
        NVSS_flux =NVSS['Total_flux_NVSS'][i]
        e_flux_nvss = NVSS['E_flux_NVSS'][i]
        NVSS_flux_values.append((NVSS_flux))
        e_NVSS_flux_values.append((e_flux_nvss))
        j_values.append((j))

### j-index determines the values of the cross-matched sources ####

  Matched_TGSSxNVSS = pd.DataFrame(TGSS_match.loc[j_values])
  Matched_TGSSxNVSS = Matched_TGSSxNVSS.reset_index(drop = True)
  Matched_TGSSxNVSS.insert(11, "flux_NVSS",NVSS_flux_values, True) # insert the flux_tgss
  Matched_TGSSxNVSS.insert(12, "e_flux_NVSS",e_NVSS_flux_values, True) 

  if save_xmatch_csv == True:
    filename = str(dir)+'/'+str(region)+'/Matched_sources_'+str(region)+'.csv'
    Matched_TGSSxNVSS.to_csv(filename, index=False)

  if show_matches == True :
    # Plot the matched sources on a ra dec which would show the matching of different catalogs by different color circles
    fig, ax = plt.subplots(figsize=(8,6))
    ax.scatter(NVSS['RAJ2000'], NVSS['DEJ2000'], color="none", edgecolor="red", label='NVSS')
    ax.scatter(TGSS['RA'], TGSS['DEC'], color="none", edgecolor="blue", label='TGSS')
    ax.scatter(observed_sources[' RA'],observed_sources[' DEC'],color="none", edgecolor="green", label='Observed_sources')
    ax.scatter(Matched_TGSSxNVSS['Ra'],Matched_TGSSxNVSS['Dec'],color="black",s=50, linewidth=0.5, alpha=0.7, label='Matched')
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('Dec (deg)')
    ax.set_title('Matched sources_'+str(region))
    ax.legend()
    ax.grid(True)
    #plt.show()
    file_loc = str(dir)+'/'+str(region)+'/matched_sources.png'
    plt.savefig(file_loc)
    plt.show()

  return Matched_TGSSxNVSS

## Subroutine for Spectral_index value ##

def Spectral_index(matched_sources,dir,region,spectral_hist=True,get_spectral_index =True,get_candidates=True):
  #importing required libraries 
  import numpy as np 
  import pandas as pd 
  #from astropy.table import Table
  from scipy.stats import linregress
  import matplotlib.pyplot as plt
  #Defining the arrays
  frequencies = np.array([1.50*1e08,3.2270*1e08,1.4*1e09]) #in Hz
  logoffreq = []
  flux = []
  e_flux = []
  spectral_in = []

  #Best_fit_line
  def f(x):
    return slope*x + intercept

  x = np.linspace(np.log10(1.0*1e08),np.log10(frequencies[-1]*2), 100)

  #####Initialization of Loop#####

  for i in range(len(matched_sources)): 
    title = matched_sources['TGSS_id'][i]
    Ra = matched_sources['Ra'][i]
    Dec = matched_sources['Dec'][i]
    log_of_flux = np.array([np.log10(matched_sources['flux_TGSS'][i]),np.log10(matched_sources['observed_Total_flux'][i]),np.log10(matched_sources['flux_NVSS'][i])])
    flux_1 = np.array([matched_sources['flux_TGSS'][i],matched_sources['observed_Total_flux'][i],matched_sources['flux_NVSS'][i]])
    e_flux = np.array([matched_sources['e_flux_TGSS'][i],matched_sources['observed_E_flux'][i],matched_sources['e_flux_NVSS'][i]])

    slope, intercept, r_value, p_value, std_err = linregress(np.log10(frequencies),log_of_flux)
    error = std_err / abs(slope)

    spectral_in.append((i,Ra,Dec,title,slope,error))
    

    #print("The spectral index",slope)

    fig, ax = plt.subplots()
    plt.grid()
    plt.scatter(np.log10(frequencies),log_of_flux,marker='.')
    plt.errorbar(np.log10(frequencies),log_of_flux,yerr = 0.434*(e_flux/flux_1),alpha=0.7,ecolor='black',capsize=2,ls = 'none')
    plt.title(title)
    plt.xlabel('Log of Frequency')
    plt.ylabel('Log of Flux')
    plt.tight_layout
    

    ax.plot(x,f(x),'--g',label = "best_fit_\u03B1 = "+str(slope.round(3))+"+-"+str(error.round(3)))
    plt.legend(loc='best')
    
    filename = str(dir)+'/'+str(region)+'/spectral_index_plot_{}.png'.format(i)
   
    plt.savefig(filename)
    #plt.show()
    
    #Close plot
    plt.close()
    

  #making the dataframe of all spectral indices 
  column = ('sr.no','Ra','Dec','TGSS_id','spectral_index','e_spectral_index')    
  spectral_index = pd.DataFrame(spectral_in,columns = column)

  if get_spectral_index == True:
    filename = str(dir)+'/'+str(region)+'/spectral_index_'+str(region)+'.csv'
    spectral_index.to_csv(filename, index=False)


  if get_candidates == True:
    pulsar_candidates = spectral_index[spectral_index['spectral_index'] < -0.9]
    filename = str(dir)+'/'+str(region)+'/Pulsar_candidtaes_'+str(region)+'.csv'
    pulsar_candidates.to_csv(filename, index=False)
    
    #Plot them 
    fig, ax = plt.subplots(figsize=(8,6))
    plt.scatter(spectral_index['sr.no'],spectral_index['spectral_index'],marker = '.')
    plt.scatter(pulsar_candidates['sr.no'],pulsar_candidates['spectral_index'],s=50, linewidth=0.5, alpha=0.7,color='r',label = 'Selected_candidates')
    plt.errorbar(spectral_index['sr.no'],spectral_index['spectral_index'],yerr =spectral_index['e_spectral_index'],alpha=0.7,ecolor='black',capsize=2,ls = 'none')
    plt.axhline(y=-0.9, color='b', linestyle='--', label='Spectral index = -0.9')
    plt.axhspan(ymin=spectral_index['spectral_index'].min(), ymax=-0.9, alpha=0.3, color='g', label='Spectral index < -0.9')
    plt.legend()
    plt.grid()
    plt.title('Pulsar_candidates')
    plt.ylabel('Spectral Index')
    plt.xlabel('Sr.no')
    loc = str(dir)+'/'+str(region)+'/Pulsar_candidates.png'
    plt.savefig(loc)


  if spectral_hist == True:
     ax_1 = spectral_index.hist(column='spectral_index', bins=25, grid=True, figsize=(12,8), color='#86bf91')
     # add a vertical line at spectral index = -0.9
     plt.axvline(x=-0.9, color='b', linestyle='--', label='Spectral index = -0.9')
     # Plot error bars
     #spectral_index.errorbar(bin_centers, counts, yerr=np.sqrt(counts), xerr=spectral_index['e_spectral_index'], fmt='none', ecolor='black')

     # add a translucent shade to the region with spectral index steeper than -0.9
     plt.axvspan(xmin=spectral_index['spectral_index'].min(), xmax=-0.9, alpha=0.3, color='r', label='Spectral index < -0.9')
     
     # set the plot title and labels
     plt.title('Spectral Index Distribution')
     plt.xlabel('Spectral Index')
     plt.ylabel('Counts')

     # add a legend
     plt.legend()
    
     loc = str(dir)+'/'+str(region)+'/spectral_index_hist.png'
     plt.savefig(loc)
     plt.show
  return spectral_index,pulsar_candidates

#Change the directory before executing the code

matched_sources = Crossmatching(dir='/home/viren/work/DATA/Spectral_Analysis',region='G047.6-11.0',file_path_image='/home/viren/work/DATA/Spectral_Analysis/G047.6-11.0/G047.6-11.0.FITS',file_path_TGSS='/home/viren/work/DATA/Spectral_Analysis/G047.6-11.0/TGSSADR1_7sigma_catalog.tsv',show_matches=True,save_xmatch_csv=True)

Spectral_index(matched_sources,dir='/home/viren/work/DATA/Spectral_Analysis',region='G047.6-11.0')

#Getting Pulsars in the region of image 

import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from psrqpy import QueryATNF

# convert the input values to sexagesimal format
c = SkyCoord(ra=ref_ra*u.deg, dec=ref_dec*u.deg, frame='icrs')
ra_sex = c.ra.to_string(unit=u.hour, sep=':')
dec_sex = c.dec.to_string(unit=u.degree, sep=':')

# set the circular boundary centre (RAJ then DECJ) and radius in sexagesimal format
circular_boundary = [ra_sex, dec_sex, radius_deg]

# query the ATNF database
query = QueryATNF(params=['JNAME', 'RAJ', 'DECJ', 'S400', 'S1400'], circular_boundary=circular_boundary)

# create NumPy arrays for each column in the query
jname = np.array(query['JNAME'])
raj = np.array(query['RAJ'])
decj = np.array(query['DECJ'])
s400 = np.array(query['S400'])
s1400 = np.array(query['S1400'])

# calculate spectral index using NumPy operations
spec_idx = np.where((np.isnan(s400)) | (np.isnan(s1400)), np.nan, (np.log10(s400) - np.log10(s1400)) / (np.log10(400) - np.log10(1400)))

# combine the arrays into a Pandas DataFrame
pulsar_spectral_index = pd.DataFrame({
    'JNAME': jname,
    'RAJ': raj,
    'DECJ': decj,
    'S400' : s400,
    's1400' : s1400,
    'spectral_index': spec_idx
})

print(pulsar_spectral_index)
