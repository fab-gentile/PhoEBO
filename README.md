# PhoEBO - Photometry Extractor for Blended Objects
<p align="center">
<img src= "https://user-images.githubusercontent.com/83715159/223449606-8ca57af2-97e4-454b-a48c-876070f806e1.png" width=450, height=200>
</p>


# Introduction:
 PhoEBO is a code developed by Fabrizio Gentile, Margherita Talia (University of Bologna), and other collaborators to extract accurate photometry of blended objects in low-resolution astronomical images. The code is optimized for the deblending of the Radio-Selected NIRdark galaxies in the COSMOS field, but it can be employed in other similar studies. A full description of the code is presented in Gentile et al (subm), together with its validation on simulated data. For further details, feel free to contact the author at fabrizio.gentile3@unibo.it.

# Disclaimer:
We are still working to transform PhoEBO from "a code that we are using" to "a code available for everyone". Be patient if there are still some bugs and the documentation is not perfect. If you spot any errors or if you'd like to collaborate to improve the code or the documentation, feel free to get in touch with us!

# Instructions:

The main class of PhoEBO (the "image" class) is included in the _Phoebo.py_ script, therefore it can be simply imported in a Python script as:
 ```python
 from PhoEBO import image
 ```
In this GitHub repository, you can also find a _Param.py_ script containing the main parameters of the code and an example script (_Example.py_) to run PhoEBO.

## Step 1: Choose the parameters:

By editing the _Params.py_ file, you can change the parameters guiding the PhoEBO execution. Here you can find a brief description of each of them:

**working_directory**: Directory in which most of the files will be read/written  

**opt_bands**: Names of the optical bands that will be analyzed by PhoEBO  

**nir_bands**: Names of the NIR bands  
**irac_bands**: Names of the IRAC channels  

**opt_fwhms**: PSF FWHM of the (Gaussian) PSFs of the optical bands  
**nir_fwhms**: PSF FWHM of the (Gaussian) PSFs of the NIR bands  
**irac_fwhms**: PSF FWHM of the PSFs of the IRAC channels.   _(NB: To use the original IRAC PSFs as in Gentile+(subm) keep it to ['ch1', 'ch2', 'ch3', 'ch4'])_   

**opt_zpt**: Photometric zero-point (AB mag) of the optical maps  
**nir_zpt**: Photometric zero-point (AB mag) of the NIR maps  
**irac_zpt**: Photometric zero-point (AB mag) of the IRAC maps  

**opt_factors**: Corrective factor for the uncertainties in the optical bands    
**nir_factors**: Corrective factor for the uncertainties in the NIR bands    
**irac_factors**: Corrective factor for the uncertainties in the IRAC bands     
_(NB: See Gentile+(subm) for a description of these parameters)_    

**irac_flux_factors**: Aperture corrections for the IRAC fluxes    

**irac_shape**: Shape (px,px) of the IRAC maps    

**input_catalogue**: Path of the input catalog (see next step)   
**output_catalogue**: Path of the otput catalog (see next step)  
**PATH_SAVE**: Path of the output plots _(NB: The plotting routines are still work-in-progress! Use them carefully!)_  
**NAME**: A name for the directories containing the data (see next step)

## Step 2 : Preparing the data  

1. Create an input catalog as a csv file with three columns (ID, RA, DEC) for each of the galaxies you'd like to analyze with PhoEBO. The RA-DEC should be the coordinates of the radio-birght source for which you'd like to estimate the flux.
1. Inside a Working Directory (WD, hereafter), create a series of sub-directories where to store the images that PhoEBO will analyze. Each directory should contain only the data of a single source: therefore, it should be named "NAME_ID" (e.g. COSMOS_1, COSMOS_2 ...) where the IDs are those included in the **input catalog**;
2. In each sub-directory, place the maps of each source. Each scientific map should be named as its band (e.g. Y.fits, J-fits, IRAC1.fits...). These names are those included in the **opt_bands**, **nir_bands**, and **IRAC_bands** parameters in the _Params.py_ script. You should also upload an uncertainty map for each band: the file should end with __unc.fits_ (e.g. IRAC1_unc.fits). These maps are crucial to optimize the models and to estimate the uncertainties on the extracted fluxes!

## Step 3: Running the code

Running PhoEBO on a single image consists of generating an _image_ object and executing the various functions listed below. A basic example script to run PhoEBO on multiple images can be found in the repository.

# Code Structure

The PhEBO code is built around a main python class called "Image"


