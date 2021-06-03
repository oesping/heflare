
import sunpy.map
#from sunpy.instr.aia import aiaprep
from aiapy.calibrate import register, update_pointing, normalize_exposure

import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u

import numpy as np

from astropy.time import Time


# Function to normalize the grid values
def aia_normalize(array):
    """Normalizes numpy arrays into scale 0.0 - 1.0"""
    array_min, array_max = array.min(), array.max()
    return ((array - array_min)/(array_max - array_min))




def aia_submap(filelist,timerange,xran,yran,key='aia',mean=True,normalize=True):
    maplist=sunpy.map.Map(sorted(filelist))
    obstimes=np.array([tmap.date.datetime for tmap in maplist])
    cut_list=np.logical_and(obstimes > timerange.start.datetime, obstimes<timerange.end.datetime)
    cmaplist=np.array(maplist)[cut_list]
    ###
    maps_lv2=[]
    data_lv2=None
    print(len(cmaplist))
    #import pdb
    #pdb.set_trace()
    try:
        for tmap in cmaplist:
            m_norm=tmap
            if key == 'aia':
                m_updated_pointing = update_pointing(tmap)
                m_registered = register(m_updated_pointing)
                m_norm= normalize_exposure(m_registered)
            if key == 'hmi':
                tmap.meta['crder1'] = 0
                tmap.meta['crder2'] = 0
                tmapa=tmap.resample((1024, 1024)*u.pix) ##alig aia
                m_norm=tmapa.rotate(order=3)
            
            top_right = SkyCoord(xran[1]*u.arcsec, yran[1]*u.arcsec, frame=m_norm.coordinate_frame)
            bottom_left = SkyCoord(xran[0] * u.arcsec, yran[0]* u.arcsec, frame=m_norm.coordinate_frame)
            spmap =m_norm.submap(bottom_left,top_right=top_right)
            maps_lv2.append(spmap)
             #rspmap=spmap.rotate(angle=90 * u.deg)
        if  mean and len(maps_lv2) >1:
            ##update map
            prep_map=maps_lv2[0]
            ptime=np.mean(obstimes[cut_list]-timerange.start.datetime)+timerange.start.datetime
            from astropy.time import Time
            header=prep_map.meta
            header['date-obs']=Time(ptime).isot
            datas=[tmp.data for tmp in maps_lv2]
            mean_data=np.array(datas).mean(axis=0)
            if  normalize:
                where_are_NaNs =np.isnan(mean_data)
                where_are_negs=mean_data <0.0
                mean_data[where_are_NaNs]=1.0
                mean_data[where_are_negs]=1.0
                lmean_data=np.log10(mean_data)
                data_lv2=(lmean_data-np.min(lmean_data))/(np.max(lmean_data)-np.min(lmean_data))
            else:
                data_lv2=mean_data
            prep_map= sunpy.map.Map(data_lv2, header)
            return prep_map
        
        else:
            return maps_lv2[0]
    except:
        print("No valid data!")
        return 0


def aia_hmi_align(aiamap,hmimap,sample=[1024,1024]):
    from reproject import reproject_interp
    ##reproject to aia
    map_aia=aiamap.resample((sample[0],sample[1])*u.pix)
    map_hmi=hmimap.resample((sample[0],sample[1])*u.pix)
    output, footprint = reproject_interp(map_hmi, map_aia.wcs, map_aia.data.shape,order="nearest-neighbor")
    out_hmi = sunpy.map.Map(output, map_aia.wcs)
    return out_hmi
    

    
    

def aia_composite(Rmap,Gmap,Bmap,weight=1.0):
    redn = weight*aia_normalize(Rmap.data)
    greenn =weight*aia_normalize(Gmap.data)
    bluen = weight*aia_normalize(Bmap.data)    
    rgbim= np.dstack((redn, greenn, bluen))
    return rgbim
 




