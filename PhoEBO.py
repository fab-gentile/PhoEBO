import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import sep
from astropy.stats import sigma_clipped_stats
from astropy.convolution import convolve
import scipy.optimize as opt
from photutils.psf import create_matching_kernel, CosineBellWindow
from photutils.datasets import load_irac_psf
from photutils.psf.matching import resize_psf
from astropy.modeling.models import Gaussian2D
import photutils.aperture as ape
import warnings
from skimage.transform import resize as resize_ski
from phoebo_utils import resize, crop_center, adu_to_ab, ab_to_ujy, adu_to_ab_err, ab_to_ujy_err
import Params as par

warnings.filterwarnings('ignore')



###### MAIN SCRIPT #####
    

class image():
    def __init__(self,path_det,fwhm_det,pixel_scale):
        self.radio_path, self.radio_unc_path, self.opt_path, self.nir_path, self.irac_path = path_det
        
        self.opt_det=fits.getdata(self.opt_path).byteswap().newbyteorder()
        self.nir_det=fits.getdata(self.nir_path).byteswap().newbyteorder()
        self.irac_det=fits.getdata(self.irac_path).byteswap().newbyteorder()
        self.irac_err=fits.getdata(self.irac_path.replace('.fits','_unc.fits')).byteswap().newbyteorder()
        self.radio=fits.getdata(self.radio_path).byteswap().newbyteorder()
        self.radio_unc=fits.getdata(self.radio_unc_path).byteswap().newbyteorder()
        
        self.opt_det_hdr=fits.getheader(self.opt_path)
        self.nir_det_hdr=fits.getheader(self.nir_path)
        self.irac_det_hdr=fits.getheader(self.irac_path)
        
        self.opt_det_shape=self.opt_det.shape
        self.nir_det_shape=self.nir_det.shape
        self.irac_det_shape=self.irac_det.shape
        

        self.opt_det_pixel_scale, self.nir_det_pixel_scale, self.irac_det_pixel_scale = pixel_scale
        
        self.fwhm_opt_det, self.fwhm_nir_det, self.fwhm_irac_det = fwhm_det
        
    def load_irac_psfs(self, MODE):
        if MODE=='Download':
            print('Loading IRAC PSFs')
            self.IR1=load_irac_psf(channel=1,show_progress=True)
            print('Done 1')
            self.IR2=load_irac_psf(channel=2,show_progress=True)
            print('Done 2')
            self.IR3=load_irac_psf(channel=3,show_progress=True)
            print('Done 3')
            self.IR4=load_irac_psf(channel=4,show_progress=True)
            print('Done 4')
        if MODE=='Upload':
            self.IR1=fits.open('./IRAC_PSF/ch1.fits')[0]
            self.IR2=fits.open('./IRAC_PSF/ch2.fits')[0]
            self.IR3=fits.open('./IRAC_PSF/ch3.fits')[0]
            self.IR4=fits.open('./IRAC_PSF/ch4.fits')[0]
            
        
    def add_bands(self, opt_bands, nir_bands, irac_bands, opt_fwhms, nir_fwhms, irac_fwhms):
        self.opt_fwhms=opt_fwhms
        self.nir_fwhms=nir_fwhms
        self.irac_fwhms=irac_fwhms
        
        self.irac_shape=par.irac_shape
        self.irac_psfs=np.zeros([len(irac_bands),159,159])
        
        self.opt_psfs=np.zeros([len(opt_bands),99,99])
        self.nir_psfs=np.zeros([len(nir_bands),99,99])

        
        self.irac_images=np.zeros([len(irac_bands),self.irac_shape[0],self.irac_shape[1]])
        self.irac_noises=np.zeros([len(irac_bands),self.irac_shape[0],self.irac_shape[1]])
        
        self.opt_images=np.zeros([len(opt_bands),self.opt_det_shape[0],self.opt_det_shape[1]])
        self.opt_noises=np.zeros([len(opt_bands),self.opt_det_shape[0],self.opt_det_shape[1]])
        
        self.nir_images=np.zeros([len(nir_bands),self.nir_det_shape[0],self.nir_det_shape[1]])
        self.nir_noises=np.zeros([len(nir_bands),self.nir_det_shape[0],self.nir_det_shape[1]])
        
        self.opt_bands, self.nir_bands, self.irac_bands = opt_bands, nir_bands, irac_bands
        
        for i, (band,fwhm) in enumerate(zip(self.irac_bands,self.irac_fwhms)):
            self.irac_images[i,:,:]=fits.getdata('./'+band+'.fits')
            self.irac_images[i,:,:][np.where(np.isnan(self.irac_images[i,:,:]))]=np.nanmean(self.irac_images[i,:,:])
            self.irac_hdr=fits.getheader('./'+band+'.fits')
            self.irac_noises[i,:,:]=fits.getdata('./'+band+'_unc.fits')
            
            if fwhm=='ch1': 
                ch1_hdu = self.IR1  
            elif fwhm=='ch2':
                ch1_hdu = self.IR2
            elif fwhm=='ch3':
                ch1_hdu = self.IR3
            elif fwhm=='ch4':
                ch1_hdu = self.IR4
                
            psf = ch1_hdu.data  
            scale=ch1_hdu.header['SECPIX']
            psf1= resize_psf(psf, scale, self.irac_det_pixel_scale)
            
            y, x = np.mgrid[0:psf1.shape[0], 0:psf1.shape[0]]
            gm1 = Gaussian2D(amplitude=100,x_stddev=self.fwhm_irac_det/2.355/self.irac_det_pixel_scale,y_stddev=self.fwhm_irac_det/2.355/self.irac_det_pixel_scale,x_mean=psf1.shape[0]//2,y_mean=psf1.shape[0]//2)
            g1 = gm1(x, y)
            g1 /= g1.sum()
            
            window = CosineBellWindow(alpha=0.35)
            psf=create_matching_kernel(g1, psf1, window=window)
            if not(psf.shape[0]%2): psf=psf[1:,1:]
            psf=crop_center(psf,159,159)
            self.irac_psfs[i,:,:]=psf
            
            
        for i, (band,fwhm) in enumerate(zip(self.nir_bands,self.nir_fwhms)):
            self.nir_images[i,:,:]=resize(fits.getdata('./'+band+'.fits'),fits.getheader('./'+band+'.fits'),self.nir_det_hdr)
            self.nir_images[i,:,:][np.where(np.isnan(self.nir_images[i,:,:]))]=np.nanmean(self.nir_images[i,:,:])
            self.nir_noises[i,:,:]=fits.getdata('./'+band+'_unc.fits')
            if fwhm-self.fwhm_nir_det>0:
                y, x = np.mgrid[0:99*3,0:99*3]
                gm2 = Gaussian2D(amplitude=100,x_stddev=np.sqrt((fwhm/2.355/self.nir_det_pixel_scale/3)**2-(self.fwhm_nir_det/2.355/ self.nir_det_pixel_scale/3)**2),y_stddev=np.sqrt((fwhm/2.355/ self.nir_det_pixel_scale/3)**2-(self.fwhm_nir_det/2.355/ self.nir_det_pixel_scale/3)**2),x_mean=99*1.5,y_mean=99*1.5)
                g2 = gm2(x, y)
                g2 /= g2.sum()
                g2=resize_psf(g2, self.nir_det_pixel_scale/3,self.nir_det_pixel_scale)
                self.nir_psfs[i,:,:]=g2
                
        for i, (band,fwhm) in enumerate(zip(self.opt_bands,self.opt_fwhms)):
            self.opt_images[i,:,:]=resize(fits.getdata('./'+band+'.fits'),fits.getheader('./'+band+'.fits'),self.nir_det_hdr)
            self.opt_images[i,:,:]=fits.getdata('./'+band+'.fits')
            self.opt_images[i,:,:][np.where(np.isnan(self.opt_images[i,:,:]))]=np.nanmean(self.opt_images[i,:,:])
            self.opt_noises[i,:,:]=fits.getdata('./'+band+'_unc.fits')
            if fwhm-self.fwhm_opt_det>0:
                y, x = np.mgrid[0:99*3,0:99*3]
                gm2 = Gaussian2D(amplitude=100,x_stddev=np.sqrt((fwhm/2.355/self.opt_det_pixel_scale/3)**2-(self.fwhm_opt_det/2.355/self.opt_det_pixel_scale/3)**2),y_stddev=np.sqrt((fwhm/2.355/self.opt_det_pixel_scale/3)**2-(self.fwhm_opt_det/2.355/self.opt_det_pixel_scale/3)**2),x_mean=99*1.5,y_mean=99*1.5)
                g2 = gm2(x, y)
                g2 /= g2.sum()
                g2=resize_psf(g2,self.opt_det_pixel_scale/3,self.opt_det_pixel_scale)
                self.opt_psfs[i,:,:]=g2
            
    def detect(self,thresh,npix):
        
        # IRAC part
        
        mean, median, std=sigma_clipped_stats(self.irac_det)
        _, self.irac_seg=sep.extract(self.irac_det, thresh=mean+std*thresh,minarea=npix,segmentation_map=True,deblend_cont=1,deblend_nthresh=1)
        
        rad=resize_ski(self.radio,self.irac_det.shape)
        rad_unc=resize_ski(self.radio_unc,self.irac_det.shape)
        
        
        _, self.rad_seg=sep.extract(rad, err=rad_unc, thresh=5.5 ,minarea=5 ,segmentation_map=True,deblend_cont=1e-10,deblend_nthresh=1000)

        s=self.rad_seg.shape[0]
        p=self.rad_seg[s//2,s//2]
        self.rad_seg[np.where(self.rad_seg!=p)]=0
        self.rad_seg[np.where(self.rad_seg==p)]=1
        self.rad_seg*=(np.amax(self.irac_seg)+1)
        self.irac_seg=np.amax(np.stack([self.rad_seg,self.irac_seg]),axis=0)
        
        self.irac_seg[:20,:]=0
        self.irac_seg[:,:20]=0
        self.irac_seg[80:,:]=0
        self.irac_seg[:,80:]=0
        
        all_id=np.unique(self.irac_seg)
        self.n_contaminants_irac=len(all_id)-2
        
        
        self.source_irac_seg=np.copy(self.irac_seg)
        source_id=self.source_irac_seg[self.source_irac_seg.shape[0]//2,self.source_irac_seg.shape[1]//2]
        self.source_irac_seg[np.where(self.source_irac_seg!=source_id)]=0
        self.source_irac_seg=self.source_irac_seg/np.amax(self.source_irac_seg)
        self.source_irac=self.source_irac_seg*self.irac_det
        
        conts_id=all_id[np.where((all_id!=source_id)&(all_id!=0))]
        self.contaminants_irac_seg=np.zeros([self.n_contaminants_irac,self.irac_det_shape[0],self.irac_det_shape[1]])
        self.contaminants_irac=np.zeros([self.n_contaminants_irac,self.irac_det_shape[0],self.irac_det_shape[1]])
        for n,i in enumerate(conts_id):
            cont=np.copy(self.irac_seg)
            cont[np.where(cont!=i)]=0
            self.contaminants_irac_seg[n,:,:]=cont
            self.contaminants_irac_seg[n,:,:]=self.contaminants_irac_seg[n,:,:]/np.amax(self.contaminants_irac_seg[n,:,:])
            self.contaminants_irac[n,:,:]=self.irac_det*self.contaminants_irac_seg[n,:,:]
            
        # NIR part
        
        mean, median, std=sigma_clipped_stats(self.nir_det)
        _, self.nir_seg=sep.extract(self.nir_det,thresh=mean+thresh*std,minarea=npix,segmentation_map=True,deblend_cont=1,deblend_nthresh=1)
        
        rad=resize_ski(self.radio,self.nir_det.shape)
        rad_unc=resize_ski(self.radio_unc,self.nir_det.shape)
        
        _, self.rad_seg=sep.extract(rad, err=rad_unc, thresh=5.5 ,minarea=5 ,segmentation_map=True,deblend_cont=1e-10,deblend_nthresh=1000)

        s=self.rad_seg.shape[0]
        p=self.rad_seg[s//2,s//2]
        self.rad_seg[np.where(self.rad_seg!=p)]=0
        self.rad_seg[np.where(self.rad_seg==p)]=1
        self.rad_seg*=(np.amax(self.nir_seg)+1)
        self.nir_seg=np.amax(np.stack([self.rad_seg,self.nir_seg]),axis=0)
        
        self.nir_seg[:20,:]=0
        self.nir_seg[:,:20]=0
        self.nir_seg[80:,:]=0
        self.nir_seg[:,80:]=0
        
        all_id=np.unique(self.nir_seg)
        self.n_contaminants_nir=len(all_id)-2
        
        
        self.source_nir_seg=np.copy(self.nir_seg)
        source_id=self.source_nir_seg[self.source_nir_seg.shape[0]//2,self.source_nir_seg.shape[1]//2]
        self.source_nir_seg[np.where(self.source_nir_seg!=source_id)]=0
        self.source_nir_seg=self.source_nir_seg/np.amax(self.source_nir_seg)
        self.source_nir=self.source_nir_seg*self.nir_det
        
        conts_id=all_id[np.where((all_id!=source_id)&(all_id!=0))]
        self.contaminants_nir_seg=np.zeros([self.n_contaminants_nir,self.nir_det_shape[0],self.nir_det_shape[1]])
        self.contaminants_nir=np.zeros([self.n_contaminants_nir,self.nir_det_shape[0],self.nir_det_shape[1]])
        for n,i in enumerate(conts_id):
            cont=np.copy(self.nir_seg)
            cont[np.where(cont!=i)]=0
            self.contaminants_nir_seg[n,:,:]=cont
            self.contaminants_nir_seg[n,:,:]=self.contaminants_nir_seg[n,:,:]/np.amax(self.contaminants_nir_seg[n,:,:])
            self.contaminants_nir[n,:,:]=self.nir_det*self.contaminants_nir_seg[n,:,:]
            
        # OPT part
        
        mean, median, std=sigma_clipped_stats(self.opt_det)
        _, self.opt_seg=sep.extract(self.opt_det,thresh=mean+thresh*std,minarea=npix,segmentation_map=True,deblend_cont=1,deblend_nthresh=1)
        
        rad=resize_ski(self.radio,self.opt_det.shape)
        rad_unc=resize_ski(self.radio_unc,self.opt_det.shape)
        
        _, self.rad_seg=sep.extract(rad, err=rad_unc, thresh=5.5 ,minarea=5 ,segmentation_map=True,deblend_cont=1e-10,deblend_nthresh=1000)
 
        s=self.rad_seg.shape[0]
        p=self.rad_seg[s//2,s//2]
        self.rad_seg[np.where(self.rad_seg!=p)]=0
        self.rad_seg[np.where(self.rad_seg==p)]=1
        self.rad_seg*=(np.amax(self.opt_seg)+1)
        self.opt_seg=np.amax(np.stack([self.rad_seg,self.opt_seg]),axis=0)
        
        self.opt_seg[:20,:]=0
        self.opt_seg[:,:20]=0
        self.opt_seg[80:,:]=0
        self.opt_seg[:,80:]=0
        
        all_id=np.unique(self.opt_seg)
        self.n_contaminants_opt=len(all_id)-2
        
        
        self.source_opt_seg=np.copy(self.opt_seg)
        source_id=self.source_opt_seg[self.source_opt_seg.shape[0]//2,self.source_opt_seg.shape[1]//2]
        self.source_opt_seg[np.where(self.source_opt_seg!=source_id)]=0
        self.source_opt_seg=self.source_opt_seg/np.amax(self.source_opt_seg)
        self.source_opt=self.source_opt_seg*self.opt_det
        
        conts_id=all_id[np.where((all_id!=source_id)&(all_id!=0))]
        self.contaminants_opt_seg=np.zeros([self.n_contaminants_opt,self.opt_det_shape[0],self.opt_det_shape[1]])
        self.contaminants_opt=np.zeros([self.n_contaminants_opt,self.opt_det_shape[0],self.opt_det_shape[1]])
        for n,i in enumerate(conts_id):
            cont=np.copy(self.opt_seg)
            cont[np.where(cont!=i)]=0
            self.contaminants_opt_seg[n,:,:]=cont
            self.contaminants_opt_seg[n,:,:]=self.contaminants_opt_seg[n,:,:]/np.amax(self.contaminants_opt_seg[n,:,:])
            self.contaminants_opt[n,:,:]=self.opt_det*self.contaminants_opt_seg[n,:,:]
        
        
        
    def create_models(self):
        
    #IRAC part
    
        self.source_irac_mod1=np.zeros([len(self.irac_bands),self.irac_det_shape[0],self.irac_det_shape[1]])
        self.source_irac_mod=np.zeros([len(self.irac_bands),par.irac_shape[0],par.irac_shape[1]])
        
        for i,band in enumerate(self.irac_bands):
            self.source_irac_mod1[i,:,:]=convolve(self.source_irac,self.irac_psfs[i])
            self.source_irac_mod[i,:,:]=resize(self.source_irac_mod1[i,:,:],self.irac_det_hdr,self.irac_hdr)
            self.source_irac_mod[i,:,:]/=np.amax(self.source_irac_mod[i,:,:])
            
        self.contaminants_irac_mod1=np.zeros([self.n_contaminants_irac, len(self.irac_bands),self.irac_det_shape[0],self.irac_det_shape[1]])
        self.contaminants_irac_mod=np.zeros([self.n_contaminants_irac,len(self.irac_bands),par.irac_shape[0],par.irac_shape[1]])

        for i in range(self.n_contaminants_irac):
            for i_band in range(len(self.irac_bands)):
                self.contaminants_irac_mod1[i,i_band,:,:]=convolve(self.contaminants_irac[i],self.irac_psfs[i_band])
                self.contaminants_irac_mod[i,i_band,:,:]=resize(self.contaminants_irac_mod1[i,i_band,:,:],self.irac_det_hdr,self.irac_hdr)
                self.contaminants_irac_mod[i,i_band,:,:]/=np.amax(self.contaminants_irac_mod[i,i_band,:,:])
                
    #NIR part
    
        self.source_nir_mod1=np.zeros([len(self.nir_bands),self.nir_det_shape[0],self.nir_det_shape[1]])
        
        for i,band in enumerate(self.nir_bands):
            if self.nir_fwhms[i]-self.fwhm_nir_det==0: 
                self.source_nir_mod1[i,:,:]=self.source_nir
            else:
                self.source_nir_mod1[i,:,:]=convolve(self.source_nir,self.nir_psfs[i],normalize_kernel=False)
            self.source_nir_mod1[i,:,:]/=np.amax(self.source_nir_mod1[i,:,:])
            
        self.contaminants_nir_mod1=np.zeros([self.n_contaminants_nir, len(self.nir_bands),self.nir_det_shape[0],self.nir_det_shape[1]])

        for i in range(self.n_contaminants_nir):
            for i_band in range(len(self.nir_bands)):
                if self.nir_fwhms[i_band]-self.fwhm_nir_det==0:       
                    self.contaminants_nir_mod1[i,i_band,:,:]=self.contaminants_nir[i]
                else:
                    self.contaminants_nir_mod1[i,i_band,:,:]=convolve(self.contaminants_nir[i],self.nir_psfs[i_band],normalize_kernel=False)
                self.contaminants_nir_mod1[i,i_band,:,:]/=np.amax(self.contaminants_nir_mod1[i,i_band,:,:])
                
    #OPT part
    
        self.source_opt_mod1=np.zeros([len(self.opt_bands),self.opt_det_shape[0],self.opt_det_shape[1]])
        
        for i,band in enumerate(self.opt_bands):
            if self.opt_fwhms[i]-self.fwhm_opt_det==0: 
                self.source_opt_mod1[i,:,:]=self.source_opt
            else:
                self.source_opt_mod1[i,:,:]=convolve(self.source_opt,self.opt_psfs[i],normalize_kernel=False)
            self.source_opt_mod1[i,:,:]/=np.amax(self.source_opt_mod1[i,:,:])
            
        self.contaminants_opt_mod1=np.zeros([self.n_contaminants_opt, len(self.opt_bands),self.opt_det_shape[0],self.opt_det_shape[1]])

        for i in range(self.n_contaminants_opt):
            for i_band in range(len(self.opt_bands)):
                if self.opt_fwhms[i_band]-self.fwhm_opt_det==0:       
                    self.contaminants_opt_mod1[i,i_band,:,:]=self.contaminants_opt[i]
                else:
                    self.contaminants_opt_mod1[i,i_band,:,:]=convolve(self.contaminants_opt[i],self.opt_psfs[i_band],normalize_kernel=False)
                self.contaminants_opt_mod1[i,i_band,:,:]/=np.amax(self.contaminants_opt_mod1[i,i_band,:,:])
                
                
    def optimize_models(self):
        
        def residual(s,img,source,contaminants,noise):
            c=np.copy(contaminants)
            sou=np.copy(source)
            mean,median,std=sigma_clipped_stats(img)
            sou=s[0]*sou
            for i in range(len(contaminants)):
                c[i]=contaminants[i]*s[i+1]
            cont=np.nansum(c,axis=0)
            mod=cont+sou
            res=np.nansum(((mod-img)/noise)**2)
            return res
    #IRAC Part
        self.images_irac_sub=np.zeros([len(self.irac_bands),par.irac_shape[0],par.irac_shape[1]])
        self.irac_flags=[]
        for i in range(len(self.irac_bands)):
            x0=np.zeros(self.n_contaminants_irac+1)
            x0[0]=np.amax(self.source_irac)/np.amax(self.source_irac_mod[i,:,:])
            x0[1:]=[sigma_clipped_stats(self.irac_images[i])[2] for _ in range(self.n_contaminants_irac)]
            bound=[(0,None) for _ in range(len(x0))]
            result=opt.minimize(residual,x0=x0,args=(self.irac_images[i],self.source_irac_mod[i],self.contaminants_irac_mod[:,i,:,:],self.irac_noises[i]),options={'maxiter':1000}, bounds=bound)
            for j in range(self.n_contaminants_irac):
                self.contaminants_irac_mod[j,i,:,:]=self.contaminants_irac_mod[j,i,:,:]*result.x[j+1]
            self.source_irac_mod[i]*=result.x[0]
            self.images_irac_sub[i]=self.irac_images[i]-np.nansum(self.contaminants_irac_mod[:,i,:,:],axis=0)
            self.irac_flags.append(result.success)
            
    #NIR Part
        self.images_nir_sub=np.zeros([len(self.nir_bands),self.nir_det_shape[0],self.nir_det_shape[1]])
        self.nir_flags=[]
        for i in range(len(self.nir_bands)):
            x0=np.zeros(self.n_contaminants_nir+1)
            x0[0]=np.amax(self.source_nir)/np.amax(self.source_nir_mod1[i,:,:])
            x0[1:]=[sigma_clipped_stats(self.nir_images[i])[2] for _ in range(self.n_contaminants_nir)]
            bound=[(0,None) for _ in range(len(x0))]
            result=opt.minimize(residual,x0=x0,args=(self.nir_images[i],self.source_nir_mod1[i],self.contaminants_nir_mod1[:,i,:,:],self.nir_noises[i]),options={'maxiter':1000}, bounds=bound)
            for j in range(self.n_contaminants_nir):
                self.contaminants_nir_mod1[j,i,:,:]=self.contaminants_nir_mod1[j,i,:,:]*result.x[j+1]
            self.source_nir_mod1[i]*=result.x[0]
            self.images_nir_sub[i]=self.nir_images[i]-np.nansum(self.contaminants_nir_mod1[:,i,:,:],axis=0)
            self.nir_flags.append(result.success)
            
    #OPT Part
        self.images_opt_sub=np.zeros([len(self.opt_bands),self.opt_det_shape[0],self.opt_det_shape[1]])
        self.opt_flags=[]
        for i in range(len(self.opt_bands)):
            x0=np.zeros(self.n_contaminants_opt+1)
            x0[0]=np.amax(self.source_opt)/np.amax(self.source_opt_mod1[i,:,:])
            x0[1:]=[sigma_clipped_stats(self.opt_images[i])[2] for _ in range(self.n_contaminants_opt)]
            bound=[(0,None) for _ in range(len(x0))]
            result=opt.minimize(residual,x0=x0,args=(self.opt_images[i],self.source_opt_mod1[i],self.contaminants_opt_mod1[:,i,:,:],self.opt_noises[i]),options={'maxiter':1000}, bounds=bound)
            for j in range(self.n_contaminants_opt):
                self.contaminants_opt_mod1[j,i,:,:]=self.contaminants_opt_mod1[j,i,:,:]*result.x[j+1]
            self.source_opt_mod1[i]*=result.x[0]
            self.images_opt_sub[i]=self.opt_images[i]-np.nansum(self.contaminants_opt_mod1[:,i,:,:],axis=0)
            self.opt_flags.append(result.success)
        fits.writeto('sub.fits',data=self.images_opt_sub[0],overwrite=True)

                
                
    def extract_flux(self, radius, opt_factors, nir_factors, irac_factors,irac_flux_factors,subtract):
        self.tab=np.ones([(len(self.opt_bands)+len(self.nir_bands)+len(self.irac_bands)),3])
        
        #opt
        ann=ape.CircularAnnulus([self.opt_images[0].shape[0]//2,self.opt_images[0].shape[0]//2],r_in=radius[1]/self.opt_det_pixel_scale,r_out=radius[2]/self.opt_det_pixel_scale)
        circ=ape.CircularAperture([self.opt_images[0].shape[0]//2,self.opt_images[0].shape[0]//2], r=radius[0]/self.opt_det_pixel_scale)
        for i in range(len(self.opt_bands)):
            aps=[circ,ann]
            phot_table=ape.aperture_photometry(self.images_opt_sub[i], aps, error=self.opt_noises[i])
            bkg_mean=phot_table['aperture_sum_1']/ann.area
            bkg_sum=bkg_mean*circ.area
            if subtract:
                self.tab[i,0]=phot_table['aperture_sum_0'].value[0]-bkg_sum.value[0]
                self.tab[i,1]=phot_table['aperture_sum_err_0'].value[0]
                self.tab[i,2]=self.opt_flags[i]
            else:
                self.tab[i,0]=phot_table['aperture_sum_0'].value[0]
                self.tab[i,1]=phot_table['aperture_sum_err_0'].value[0]
                self.tab[i,2]=self.opt_flags[i]
                
        
        #nir
        ann=ape.CircularAnnulus([self.nir_images[0].shape[0]//2,self.nir_images[0].shape[0]//2],r_in=radius[1]/self.nir_det_pixel_scale,r_out=radius[2]/self.nir_det_pixel_scale)
        circ=ape.CircularAperture([self.nir_images[0].shape[0]//2,self.nir_images[0].shape[0]//2], r=radius[0]/self.nir_det_pixel_scale)
        for i in range(len(self.nir_bands)):
            aps=[circ,ann]
            phot_table=ape.aperture_photometry(self.images_nir_sub[i], aps, error=self.nir_noises[i])
            bkg_mean=phot_table['aperture_sum_1']/ann.area
            bkg_sum=bkg_mean*circ.area
            if subtract:
                self.tab[len(self.opt_bands)+i,0]=phot_table['aperture_sum_0'].value[0]-bkg_sum.value[0]
                self.tab[len(self.opt_bands)+i,1]=phot_table['aperture_sum_err_0'].value[0]
                self.tab[len(self.opt_bands)+i,2]=self.nir_flags[i]
            else:
                self.tab[len(self.opt_bands)+i,0]=phot_table['aperture_sum_0'].value[0]
                self.tab[len(self.opt_bands)+i,1]=phot_table['aperture_sum_err_0'].value[0]
                self.tab[len(self.opt_bands)+i,2]=self.nir_flags[i]
            
        #IRAC
        ann=ape.CircularAnnulus([self.irac_images[0].shape[0]//2,self.irac_images[0].shape[0]//2],r_in=radius[1]/0.6,r_out=radius[2]/0.6)
        circ=ape.CircularAperture([self.irac_images[0].shape[0]//2,self.irac_images[0].shape[0]//2], r=radius[0]/0.6)
        for i in range(len(self.irac_bands)):
            aps=[circ,ann]
            phot_table=ape.aperture_photometry(self.images_irac_sub[i], aps, error=self.irac_noises[i])
            bkg_mean=phot_table['aperture_sum_1']/ann.area
            bkg_sum=bkg_mean*circ.area
            if subtract:
                self.tab[len(self.nir_bands)+len(self.opt_bands)+i,0]=phot_table['aperture_sum_0'].value[0]-bkg_sum.value[0]
                self.tab[len(self.nir_bands)+len(self.opt_bands)+i,1]=phot_table['aperture_sum_err_0'].value[0]
                self.tab[len(self.nir_bands)+len(self.opt_bands)+i,2]=self.irac_flags[i]
            else:
                self.tab[len(self.nir_bands)+len(self.opt_bands)+i,0]=phot_table['aperture_sum_0'].value[0]
                self.tab[len(self.nir_bands)+len(self.opt_bands)+i,1]=phot_table['aperture_sum_err_0'].value[0]
                self.tab[len(self.nir_bands)+len(self.opt_bands)+i,2]=self.irac_flags[i]
        
        self.new=np.ones(self.tab.shape)*-99
        
        #opt
        for i, (factor, zpt) in enumerate(zip(opt_factors,par.opt_zpt)):
            self.new[i,0]=adu_to_ab(self.tab[i,0],zpt)
            self.new[i,0]=ab_to_ujy(self.new[i,0])
            if self.tab[i,0]>0:
                self.new[i,1]=adu_to_ab_err(self.tab[i,1],self.tab[i,0])
                self.new[i,1]=ab_to_ujy_err(self.new[i,1],self.new[i,0])*factor
            else:
                self.new[i,1]=adu_to_ab(self.tab[i,1],zpt)
                self.new[i,1]=ab_to_ujy(self.new[i,1])*factor
            self.new[i,2]=self.tab[i,2]
            
        for i, (factor,zpt) in enumerate(zip(nir_factors,par.nir_zpt)):
            i=len(self.opt_bands)+i
            self.new[i,0]=adu_to_ab(self.tab[i,0],zpt)
            self.new[i,0]=ab_to_ujy(self.new[i,0])
            if self.tab[i,0]>0:
                self.new[i,1]=adu_to_ab_err(self.tab[i,1],self.tab[i,0])
                self.new[i,1]=ab_to_ujy_err(self.new[i,1],self.new[i,0])*factor
            else:
                self.new[i,1]=adu_to_ab(self.tab[i,1],zpt)
                self.new[i,1]=ab_to_ujy(self.new[i,1])*factor
            self.new[i,2]=self.tab[i,2]

                    
        for i, (factor,flux_factor,zpt) in enumerate(zip(irac_factors,irac_flux_factors,par.irac_zpt)):
            i=len(self.opt_bands)+len(self.nir_bands)+i
            self.new[i,0]=adu_to_ab(self.tab[i,0],zpt)
            self.new[i,0]=ab_to_ujy(self.new[i,0])*flux_factor
            if self.tab[i,0]>0:
                self.new[i,1]=adu_to_ab_err(self.tab[i,1],self.tab[i,0])
                self.new[i,1]=ab_to_ujy_err(self.new[i,1],self.new[i,0])*factor
            else:
                self.new[i,1]=adu_to_ab(self.tab[i,1],zpt)
                self.new[i,1]=ab_to_ujy(self.new[i,1])*flux_factor
            self.new[i,2]=self.tab[i,2]
        
        return self.new
                
                
    
    def generate_plots(self,ID):
        fig, axs= plt.subplots(7,7,figsize=(8.27,11.69))
        for i in range(7):
            for j in range(7):
                axs[i,j].axis('off')
        axs[0,2].imshow(self.nir_images[0],cmap='jet')
        axs[0,2].set(title='Chi2')
        axs[0,3].imshow(self.radio,cmap='jet')
        axs[0,3].set(title='COSMOS_%i\n\n 3GHz'%ID)
        axs[0,4].imshow(self.nir_seg,cmap='jet')
        axs[0,4].set(title='Seg')


        for i,band in enumerate(self.opt_bands[:2]):
            m,M=np.amin(self.opt_images[i]),np.amax(self.opt_images[i])
            axs[i+1,0].imshow(self.opt_images[i], vmin=m, vmax=M,cmap='jet')
            axs[i+1,1].imshow(np.nansum([np.nansum(self.contaminants_opt_mod1[:,i,:,:],axis=0),self.source_opt_mod1[i]],axis=0), vmin=m, vmax=M,cmap='jet')
            axs[i+1,2].imshow(self.images_opt_sub[i], vmin=m, vmax=M,cmap='jet')
            if self.new[i,0]>0:
                axs[i+1,1].set(title=r'%s: (%.2f $\pm$ %.2f) $\mu$Jy'%(band,self.new[i,0],self.new[i,1]))
            else:
                axs[i+1,1].set(title=r'%s: Upp.Lim. %.2f $\mu$Jy'%(band,self.new[i,1]))
                
        for i,band in enumerate(self.opt_bands[3:]):
            n=i+3
            m,M=np.amin(self.opt_images[n]),np.amax(self.opt_images[n])
            axs[i+1,4].imshow(self.opt_images[n], vmin=m, vmax=M,cmap='jet')
            axs[i+1,5].imshow(np.nansum([np.nansum(self.contaminants_opt_mod1[:,n,:,:],axis=0),self.source_opt_mod1[n]],axis=0), vmin=m, vmax=M,cmap='jet')
            axs[i+1,6].imshow(self.images_opt_sub[n], vmin=m, vmax=M,cmap='jet')
            if self.new[n,0]>0:
                axs[i+1,5].set(title=r'%s: (%.2f $\pm$ %.2f) $\mu$Jy'%(band,self.new[n,0],self.new[n,1]))
            else:
                axs[i+1,5].set(title=r'%s: Upp.Lim. %.2f $\mu$Jy'%(band,self.new[n,1]))
                
        for i,band in enumerate(self.nir_bands[:2]):
            m,M=np.amin(self.nir_images[i]),np.amax(self.nir_images[i])
            axs[i+3,0].imshow(self.nir_images[i], vmin=m, vmax=M,cmap='jet')
            axs[i+3,1].imshow(np.nansum(self.contaminants_nir_mod1[:,i,:,:],axis=0)+self.source_nir_mod1[i], vmin=m, vmax=M,cmap='jet')
            axs[i+3,2].imshow(self.images_nir_sub[i], vmin=m, vmax=M,cmap='jet')
            if self.new[len(self.opt_bands)+i,0]>0:
                axs[i+3,1].set(title=r'%s: (%.2f $\pm$ %.2f) $\mu$Jy'%(band,self.new[len(self.opt_bands)+i,0],self.new[len(self.opt_bands)+i,1]))
            else:
                axs[i+3,1].set(title=r'%s: Upp.Lim. %.2f $\mu$Jy'%(band,self.new[len(self.opt_bands)+i,1]))
                
        for i,band in enumerate(self.nir_bands[2:]):
            n=i+2
            m,M=np.amin(self.nir_images[n]),np.amax(self.nir_images[n])
            axs[i+3,4].imshow(self.nir_images[n], vmin=m, vmax=M,cmap='jet')
            axs[i+3,5].imshow(np.nansum(self.contaminants_nir_mod1[:,n,:,:],axis=0)+self.source_nir_mod1[n], vmin=m, vmax=M,cmap='jet')
            axs[i+3,6].imshow(self.images_nir_sub[n], vmin=m, vmax=M,cmap='jet')
            if self.new[len(self.opt_bands)+n,0]>0:
                axs[i+3,5].set(title=r'%s: (%.2f $\pm$ %.2f) $\mu$Jy'%(band,self.new[len(self.opt_bands)+n,0],self.new[len(self.opt_bands)+n,1]))
            else:
                axs[i+3,5].set(title=r'%s: Upp.Lim. %.2f $\mu$Jy'%(band,self.new[len(self.opt_bands)+n,1]))
        
        for i,band in enumerate(self.irac_bands[:2]):
            m,M=np.amin(self.irac_images[i]),np.amax(self.irac_images[i])
            axs[i+5,0].imshow(self.irac_images[i], vmin=m, vmax=M,cmap='jet')
            axs[i+5,1].imshow(np.nansum(self.contaminants_irac_mod[:,i,:,:],axis=0)+self.source_irac_mod[i], vmin=m, vmax=M,cmap='jet')
            axs[i+5,2].imshow(self.images_irac_sub[i], vmin=m, vmax=M,cmap='jet')
            if self.new[len(self.opt_bands)+len(self.nir_bands)+i,0]>0:
                axs[i+5,1].set(title=r'%s: (%.2f $\pm$ %.2f) $\mu$Jy'%(band,self.new[len(self.opt_bands)+len(self.nir_bands)+i,0],self.new[len(self.opt_bands)+len(self.nir_bands)+i,1]))
            else:
                axs[i+5,1].set(title=r'%s: Upp.Lim. %.2f $\mu$Jy'%(band,self.new[len(self.opt_bands)+len(self.nir_bands)+i,1]))
                
        for i,band in enumerate(self.irac_bands[2:]):
            n=i+2
            m,M=np.amin(self.irac_images[n]),np.amax(self.irac_images[n])
            axs[i+5,4].imshow(self.irac_images[n], vmin=m, vmax=M,cmap='jet')
            axs[i+5,5].imshow(np.nansum(self.contaminants_irac_mod[:,n,:,:],axis=0)+self.source_irac_mod[n], vmin=m, vmax=M,cmap='jet')
            axs[i+5,6].imshow(self.images_irac_sub[n], vmin=m, vmax=M,cmap='jet')
            if self.new[len(self.opt_bands)+len(self.nir_bands)+n,0]>0:
                axs[i+5,5].set(title=r'%s: (%.2f $\pm$ %.2f) $\mu$Jy'%(band,self.new[len(self.opt_bands)+len(self.nir_bands)+n,0],self.new[len(self.opt_bands)+len(self.nir_bands)+n,1]))
            else:
                axs[i+5,5].set(title=r'%s: Upp.Lim. %.2f $\mu$Jy'%(band,self.new[len(self.opt_bands)+len(self.nir_bands)+n,1]))
        
        fig.savefig(par.PATH_SAVE+'%i.pdf'%ID ,dpi=600, bbox_inches='tight')
        plt.close(fig)
        
        

