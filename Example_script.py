from astropy.io import fits
import numpy as np
import os
import pandas as pd
import Params as par
from PhoEBO import image


os.chdir(par.working_directory)
df=pd.read_csv(par.input_catalogue)
f=open(par.output_catalogue,'w')
f.write('ID')
for band in par.opt_bands+par.nir_bands+par.irac_bands:
    f.write(',%s,d%s,FLAG_%s'%(band,band,band))
f.write('\n')

ids=df.ID_FG[:2]     
for i,ID in enumerate(df.ID_FG[:2]):
    # try:
    print(ID,'/',len(ids)) 
    os.chdir('./%s_%i'%(par.NAME,ID))
    os.system('cp *.fits ../Try/')
    os.chdir('../Try/')
    
    for band in ['Ks','J','Y','H']:
        data=fits.getdata('%s_wgt.fits'%band)
        hdr=fits.getheader('%s_wgt.fits'%band)
        data=np.sqrt(1/data)
        fits.writeto('%s_unc.fits'%band,data=data,header=hdr,overwrite=True)
                       
    ex=image(['./3GHz.fits','./3GHz_unc.fits','./HSC_i.fits', './Ks.fits', 'Ks.fits'],[0.61,0.75,0.75],[0.15,0.15,0.15])
    ex.load_irac_psfs(MODE='Upload')
    ex.add_bands(par.opt_bands, par.nir_bands, par.irac_bands, par.opt_fwhms, par.nir_fwhms, par.irac_fwhms)
    ex.detect(4,5)
    ex.create_models()
    ex.optimize_models()
    cat=ex.extract_flux([1,1.5,2], par.opt_factors, par.nir_factors, par.irac_factors, par.irac_flux_factors ,subtract=True)
    ex.generate_plots(ID)
    del ex 
    
    f.write('%i'%ID)
    for j in range(len(cat)):
        f.write(',%f,%f,%r'%(cat[j,0],cat[j,1],bool(cat[j,2])))
    f.write('\n')
    
    os.chdir('../')
    # except:
    #     print('ERROR')
    #     os.chdir('../')

f.close()
