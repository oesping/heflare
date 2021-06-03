from astropy.time import Time
from datetime import datetime,timedelta
from  sunpy.timeseries import TimeSeries


def norh_tcx2csv(odir='.',dtype='tcx',wave='17'):
    '''read tcx data to csv'''
    norhdb='/data1/datahub/raw/NoRH'
    wavetype={'17':'tca',
                 '34':'tcz'}
    norhfile=[]
    print,norhdb
    for year in range(1992,1993,1):
        otime=datetime(year,1,1)    
        while otime.year < year:
            norh_path=norhdb+'/'+dtype+'/'+otime.strftime('/%Y/%m/%d/')
            norh_file=norh_path+wavetype[wave]+otime.strftime('/%y/%m/%d')
            norh =TimeSeries(norh_file[0], source='NoRH') 
            norhdf=norh.to_dataframe()
            norhdf.index.name='time'
            norhdf.rename(columns = {'Correlation Coefficient':'Coefficient'}, inplace = True)
            norhdf.tocsv(norh_file+'csv')
            otime=otime+timedelta(days=1.0)
    print,"Task finished!"


