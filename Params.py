###### PARAMS ###

working_directory='/home/fabrizio/Desktop/SFRD/COSMOS'

opt_bands=['HSC_g','HSC_r','HSC_i','HSC_z','HSC_y'] #Names of the Optical bands
nir_bands=['Y','J','H','Ks'] #Names of the NIR bands
irac_bands=['IRAC1','IRAC2','IRAC3','IRAC4'] #Names of the IRAC channels

opt_fwhms= [0.79,0.75,0.61,0.68,0.68] #PSF FWHM (optical bands)
nir_fwhms=[0.82, 0.79, 0.76, 0.75] #PSF FWHM (NIR bands)
irac_fwhms=['ch1','ch2','ch3','ch4'] #PSF FWHM (IRAC channels)

opt_zpt= [27,27,27,27,27] #Zero-Point (optical bands)
nir_zpt= [30,30,30,30] #Zero-Point (NIR bands)
irac_zpt= [20.48,20.48,20.48,20.48] #Zero-Point (IRAC channels)

opt_factors=[1.4,1.4,1.5,1.4,1.4] #Corrective factor (optical bands)
nir_factors=[2.75,2.6,2.5,2.4] #Corrective factor (NIR bands)
irac_factors=[1,1,1,1] #Corrective factor (IRAC channels)

irac_flux_factors=[1.1,1.1,1.37,1.57] #Aperture corrections (IRAC channels)

irac_shape=(27,27) #Shape of the IRAC maps

input_catalogue= working_directory + '/Farmer/RS_NIRDark.csv' #Path for the input catalogue (ID, RA, DEC)
output_catalogue = working_directory +'/Out_Catalogue_HSC_pub.csv' #Path for the output catalog
PATH_SAVE= working_directory +'/Try/Pdfs_1/' #Path where to save the output plots
NAME= 'COSMOS' #Name for the sub-directory
