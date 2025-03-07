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



The FITS files:
    - 
