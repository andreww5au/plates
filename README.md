# plates
Tool for converting Perth Observatory scanned plates into FITS files

The raw data consists of TIFF files - one for the plate itself, another for the
paper envelope around the plate. These TIFF files are large - around 15000x15000 pixels for
the plates (~400 MB each), and 3000x2000 for the envelopes (~2 MB each).

Hand-written metadata about each plate (from the envelope, and from observing log books
written at the time) has been transcribed and collated into an Excel spreadsheet
containing >30 columns of data. Each row describes a single plate scan.

This program reads the Excel spreadsheet, and for each row:
   - Uses the 'plate number' field to find the TIFF files for the plate and envelope.
   - Parses the metadata, working around typos, missing and invalid data, etc.
   - Reads in the plate and envelope TIFF files.
   - Writes a 1000x1000 pixel JPEG version of the plate image, with a 'Perth Observatory'
     watermark, into a 'thumbnail' directory (~75 kB each).
   - Writes a full size JPEG version of the plate image, without a watermark, into a
     'jpeg' directory (~ 3 MB each).
   - Writes a FITS version of the plate, with as much metadata as possible, into a 'fits' 
     directory.
   - Generates an all-sky map image showing all plates scanned so far as regions on the map 
     roughly corresponding to the plate size. At the end of the run, these images are 
     assembled into a video showing the progress of the plates over time.

Note that many of the plates were taken as multiple exposures, with small offsets between
each exposure. This was done to allow accurate astrometry over a wide range of stellar
magnitudes. Unlike CCD images, photographic emulsion has a very poor dynamic range, so 
brighter stars saturate into a featureless blob that is hard to determine an accurate centre
for. Instead, a typical astrometric plate will might have three exposures - for example, 
240 seconds, 120 seconds, and 10 seconds, with a small offset (a few tens of arcseconds)
between them, by opening the shutter multiple times.

The positions of the brightest stars on the plate are measured using the shortest exposure, 
with the least saturation - but this exposure only shows the brighter stars. The positions of
the faintest stars are measured using the longest exposure - but this one has the brightest 
stars hopelessly saturated. The offsets between the exposures are determined statistically
using measurements of all stars visible in more than one of the exposures on that plate.

A typical FITS header looks like:
SIMPLE  = 'T       '           / File does conform to FITS standard             
BITPIX  =                   16 / 16-bit signed integers                         
NAXIS   =                    0 / number of array dimensions                     
EXTEND  =                    T                                                  
CREATOR = 'plates.py v1.0.0 by Andrew.Williams@curtin.edu.au'                   
OBSERVAT= 'Perth Observatory, Astrographic dome.' / Observatory name            
TELESCOP= 'Astro   '           / Telescope name                                 
LATITUDE= '-32:00:27.98'       / Latitude                                       
LONGITUD= '116:08:11.36'       / Longitude East                                 
INSTRUME= 'Astrographic Plate Camera' / Instrument                              
OBSERVER= 'Y       '           / Person guiding the telescope                   
OBJECT  = 'Ast Cat Accepted'   / Target name                                    
GUIDER  = '' / Guider Settings                                                  
SEQNUM  =                 2682 / Plate scanning sequence number                 
PLATE   = '3150    '           / Plate number                                   
PLSIZE  = 'L       '           / Plate Size                                     
EMULSION= 'glass   '           / Plate type                                     
RASTRING= '113500  '                                                            
RA      =    11.58333333333333 / Right Ascension in hours                       
DESTRING= '-3700   '                                                            
DEC     =                -37.0 / Declination in degrees                         
ENV_DATE= '1913,may,23'        / Date written on envelope                       
DATE-OBS= '1913-05-23'         / Date of observation                            
TIME-OBS= '11:38:51.837'       / UTC at start of exposure                       
TIMESYS = 'UTC     '           / Time System is UTC                             
LST     = '11:25:20.005'       / Local Sidereal Time on envelope                
JD      =     2419910.98532219 / Julian Date at start of exposure               
MJD-OBS =    19910.48532218978 / Modified Julian Date at start of exposure      
EXPTIME =                  372 / Total exposure time in seconds                 
LSTART1 = '11:25:20.005'       / LST at start of exposure                       
LEND1   = '11:29:20.005'       / LST at end of exposure                         
LSTART2 = '11:29:30.005'       / LST at start of second exposure                
LEND2   = '11:31:30.005'       / LST at end of second exposure                  
LSTART3 = '11:31:40.005'       / LST at start of third exposure                 
LEND3   = '11:31:53.005'       / LST at end of third exposure                   
EQUINOX =                 1900 / Equinox missing, assuming 1900                 
ALT     =                 84.7 / Altitude in Degrees                            
AZ      =                164.5 / Azimuth in Degrees                             
COMMENT Ast Obs Bk 6Pg046                                                       
COMMENT LST Exposure                                                            
HISTORY Scanned by 'SVK' on 2019-06-07 00:00:00 (week 9)                        
END                                                                   

The FITS files:
    - 
