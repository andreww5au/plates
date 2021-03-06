
import argparse
import csv
import datetime
import glob
import os
import traceback
import warnings

from unidecode import unidecode

import astropy
from astropy.units import deg, hour
from astropy.coordinates import Angle, SkyCoord, EarthLocation
from astropy.io import fits
from astropy.time import Time, TimeDelta

import numpy as np

import matplotlib
from matplotlib import pyplot
import matplotlib.ticker as mticker

from PIL import Image
Image.MAX_IMAGE_PIXELS = 933120000

PR = 1.0   # Plate 'radius' - eg, 1.0 means a 2 x 2 degree square plate.

ASTROGRAPHIC_POS = EarthLocation.from_geodetic(lon="116:08:11.36", lat="-32:00:27.98", height=386.0)

VERSION = '1.0.0'
TIFFDIRS = ['C:/Data/plates/tiff/testing1']   # List of directories to look for TIFF files in
FITSDIR = 'C:/Data/plates/fits'      # Base directory to write FITS files out to
HDRDIR = 'C:/Data/plates/headers'
LOGDIR = 'C:/Data/plates/logs'

MONTHS = {'jan':1, 'january':1,
          'feb':2, 'february':2, 'fen':2,
          'mar':3, 'march':3,
          'apr':4, 'april':4,
          'may':5,
          'jun':6, 'june':6,
          'jul':7, 'july':7, 'ul':7,
          'aug':8, 'august':8,
          'sep':9, 'sept':9, 'september':9, 'seo':9,
          'oct':10, 'october':10,
          'nov':11, 'november':11,
          'dec':12, 'december':12}

XTL = np.arange(0, 24)
XT = (XTL - 12) * np.pi / 12

USAGE = """
Takes a spreadsheet containing plate data, and either converts scanned TIFF
images to FITS files, collects and plots statistics on the plate metadata,
or both.
"""

# Suppress the warnings from astropy about lack of precision for dates earlier than 1950's
warnings.filterwarnings("ignore")


class DateError(ValueError):
    pass


def get_infilename(platenum=''):
    """
    Given a plate number, return a full path/filename to read the TIFF file. If more than
    one file with that plate number is found, use the highest indexed (most recent) scan.

    If no files are found, return None

    :param platenum: plate number
    :return: input file path/name, or None if the file is not found.
    """
    flist = []
    for tiffdir in TIFFDIRS:
        fspec = os.path.join(tiffdir, 'Plate %s???.tif' % platenum)
        flist += glob.glob(fspec)
        fspec = os.path.join(tiffdir, 'Plate%s???.tif' % platenum)
        flist += glob.glob(fspec)
    flist.sort()
    if flist:
        return flist[-1]   # Highest numbered scan of that plate
    else:
        return    # No files found


def get_coverfilename(platenum=''):
    """
    Given a plate number, return a full path/filename to read the plate envelope TIFF file. If more than
    one file with that plate number is found, use the highest indexed (most recent) scan.

    If no files are found, return None

    :param platenum: plate number
    :return: cover envelope input file path/name, or None if the file is not found.
    """
    flist = []
    for tiffdir in TIFFDIRS:
        fspec = os.path.join(tiffdir, 'Plate %sCover???.tif' % platenum)
        flist += glob.glob(fspec)
        fspec = os.path.join(tiffdir, 'Plate%sCover???.tif' % platenum)
        flist += glob.glob(fspec)
    flist.sort()
    if flist:
        return flist[-1]   # Highest numbered scan of that plate cover
    else:
        return    # No files found


def get_outfilename(tiffname='', platenum=''):
    """
    Given a plate number, return a full path/filename to write the output FITS file.

    Use the last directory component from the original file name as a subdirectory to create the
    FITS file in.

    :param tiffname: Full path/filename of original TIFF file
    :param platenum: plate number
    :return: output file path/name
    """
    if not tiffname:
        return os.path.join(FITSDIR, 'ungrouped', 'plate%s.fits' % platenum)
    dirlist = os.path.split(os.path.dirname(tiffname))
    if len(dirlist) > 1:
        subdir = dirlist[-1]    # the directory that the tiff file is in
    else:
        subdir = ''

    return os.path.join(FITSDIR, subdir, 'plate%s.fits' % platenum)


def get_headerfilename(tiffname='', platenum=''):
    """
    Given a plate number, return a full path/filename to write the bare FITS header.

    Use the last directory component from the original file name as a subdirectory to create the
    FITS header file in.

    :param tiffname: Full path/filename of original TIFF file
    :param platenum: plate number
    :return: output file path/name
    """
    if not tiffname:
        return os.path.join(HDRDIR, 'ungrouped', 'plate%s.txt' % platenum)
    dirlist = os.path.split(os.path.dirname(tiffname))
    if len(dirlist) > 1:
        subdir = dirlist[-1]    # the directory that the tiff file is in
    else:
        subdir = ''

    return os.path.join(HDRDIR, subdir, 'plate%s.txt' % platenum)


def isdone(tiffname='', platenum=''):
    """
    Returns true if the given plate number has already been processed

    :param tiffname: Full path/filename of original TIFF file
    :param platenum:
    :return: True if the given plate number has already been processed
    """
    return os.path.exists(get_outfilename(tiffname=tiffname, platenum=platenum))


def get_observatory(telescope=''):
    """
    Given a telescope name, return the name of its observatory
    :param telescope: Telescope name from scanned plate
    :return: Observatory name to use in OBSERVAT card in FITS header
    """
    if telescope.lower().strip() == 'astro':
        return 'Perth Observatory, Astrographic dome.', ASTROGRAPHIC_POS, 'Astrographic Plate Camera'
    else:
        return 'Perth Observatory', ASTROGRAPHIC_POS, 'Astrographic Plate Camera'


def tstring_to_timerec(basis, year, month, day, tstring, obslocation=None):
    """
    Takes a time string in decimal hours, and returns an astropy.time.Time object for the given location.
    :param basis: time basis - one of 'lst', 'utc', or 'awst'
    :param year: Date of observation
    :param month:
    :param day:
    :param tstring: local sidereal time - either in hours, as a string ('12:34:56'), or as a fraction of a day ('0.475925925925926')
    :param obslocation: astropy.EarthLocation object
    :return: an astropy.time.Time object
    """
    # print(tstring)
    tstring_hours = Angle(tstring, unit=hour).hour
    # print(tstring_hours)
    obsdate = datetime.datetime(year, month, day)
    obsdateT = Time(obsdate, location=obslocation)
    # print(obsdateT)
    # print(obsdateT.sidereal_time('apparent').hour)
    if basis == 'lst':
        utc_starthour = (tstring_hours - obsdateT.sidereal_time('apparent').hour) * 23.93447 / 24.0
        if utc_starthour < 0:
            utc_starthour += 23.93447
        # print(utc_starthour)
    elif basis == 'utc':
        utc_starthour = tstring_hours
    elif basis == 'awst':
        utc_starthour = tstring_hours - 8.0   # For Perth, assume time zone is always 8 hours.
        if utc_starthour < 0:
            utc_starthour += 24.0
            obsdateT -= TimeDelta(86400, format='sec')    # If AWST < 8.0, it's the previous UT date.
    else:
        # print("Invalid time basis '%s'" % basis)
        return
    obsdateT += TimeDelta(utc_starthour * 3600, format='sec')
    obsdateT.location = obslocation
    # print(obsdateT)
    # print()
    return obsdateT


def parsehex(hexstring='', ra=True):
    """Parse a string guaranteed to be in XX:YY:ZZ.zzz format, but
       make sure that the parts are value deg/hours, minutes and seconds.

       If all looks OK, return a float with the value.
       If there's a problem, return an error string.
    """
    parts = hexstring.split(':')
    p1, p2, p3 = '0', '0', '0'
    if len(parts) == 2:
        p1, p2 = tuple(parts)
    elif len(parts) == 3:
        p1, p2, p3 = tuple(parts)
    else:
        return "Too many colons in value: %s" % hexstring

    try:
        p1v = int(p1)
    except ValueError:
        return "Invalid hours component in value: %s" % hexstring

    try:
        p2v = float(p2)
    except ValueError:
        return "Invalid minutes component in value: %s" % hexstring

    try:
        p3v = float(p3)
    except ValueError:
        return "Invalid seconds component in value: %s" % hexstring

    if ra:
        if p1v >= 24:
            return "RA hours > 24: %s" % hexstring
    else:
        if p1v > 90:
            return "abs(Dec degrees) > 90: %s" % hexstring

    if p2v >= 60.0:
        return "minutes > 60: %s" % hexstring

    if p3v >= 60.0:
        return "seconds > 60: %s" % hexstring

    val = p1v + (p2v / 60.0) + (p3v / 3600.0)
    if ra and (val >= 24.0):
        return "RA value > 24.0: %s" % hexstring
    elif (not ra) and (val > 90.0):
        return "abs(Dec) value > 90.0: %s" % hexstring
    return val


def parseval(valstring='', ra=True):
    """Parse an ra/dec value in one of a dozen stupid formats:
    RAs:
        12:34
        12:34:56
        12:34:56.7
        12:34.56
        12.3
        12
        1234
        1234.5
        123456
        123456.7
    Decs:
        All above, with or without a sign (+/-) in front.
    """
    s = valstring.strip()

    # Strip off the sign character, if present
    sign = '+'
    if s[0] in ['+', '-']:
        sign = s[0]
        s = s[1:].strip()

    if sign == '+':
        sgn = 1
    else:
        if ra:
            return "Unexpected negative value for RA: %s" % valstring
        sgn = -1

    rv = None

    if ':' in s:   # 12:34, 12:34:56, 12:34:56.7, 12:34.56
        rv = parsehex(hexstring=valstring, ra=ra)

    if rv is None:
        if '.' in s:
            bdp, adp = s.split('.')
            if len(bdp) <= 2:    # 1.2, 12.3
                try:
                    rv = float(s)
                except:
                    return "Invalid float conversion for %s" % valstring
            elif len(bdp) == 4:  # 1234.5
                rv = parsehex("%s:%s.%s" % (bdp[:2], bdp[2:], adp), ra=ra)
            elif len(bdp) == 6:  # 123456.7
                rv = parsehex("%s:%s:%s.%s" % (bdp[:2], bdp[2:4], bdp[4:], adp), ra=ra)
            else:
                rv = "weird number of chars before decimal point: %s" % valstring
        else:
            if len(s) <= 2:  # 1, 12
                try:
                    rv = float(s)
                except ValueError:
                    rv = "Invalid float conversion for %s" % valstring
            elif len(s) == 4:   # 1234
                rv = parsehex("%s:%s" % (s[:2], s[2:4]), ra=ra)
            elif len(s) == 6:   # 123456
                rv = parsehex("%s:%s:%s" % (s[:2], s[2:4], s[4:6]), ra=ra)
            else:
                rv = "weird number of characters without decimal point: %s" % valstring

    if type(rv) is not float:
        return rv
    else:
        if ra and (rv >= 24.0):
            return "RA value > 24.0: %s" % valstring
        elif (not ra) and (rv > 90.0):
            return "abs(Dec) value > 90.0: %s" % valstring
        else:
            return sgn * rv


def parsedate(datestring):
    """Date nominally looks like '1902,jul,22'
    """
    datestring = datestring.replace('.', ',')   # Replace any dots with commas, because some rows use them instead
    datestring = datestring.replace(':', ',')   # Replace any colons with commas as well
    bits = datestring.split(',')
    bits = [x for x in bits if x]    # Drop any empty components - eg 1902,,jun,22
    if len(bits) > 3:
        return None    # Can't fix this, too many components.

    if len(bits) == 3:    # Catch '1902,jun,22', '1902, jun, 22', etc - dots or commas with or without spaces before/after them
        return bits[0].strip(), bits[1].strip().lower(), bits[2].strip()

    datestring = datestring.replace(',', ' ')   # Any commas are now replaced with spaces as a separator
    bits = datestring.split(' ')
    bits = [x for x in bits if x]    # Drop any empty components - eg 1902,,jun 22
    if len(bits) > 3:
        return None    # Can't fix this, too many components.

    if len(bits) == 3:  # Catch '1902 jun 22', '1902 jun,22', etc - spaces, or dots or commas with or without spaces before/after them
        return bits[0].strip(), bits[1].strip().lower(), bits[2].strip()

    if len(bits) == 2:
        b0, b1 = bits[0].strip(), bits[1].strip()
        if (len(b0) >= 5) and (b0[:-3].isdigit()) and (b0[-3:].lower() in MONTHS):   # Catch 1902jun,22
            return b0[:-3], b0[-3:].lower().strip(), b1.strip()
        if (len(b1) >= 4) and (b1[:3].lower() in MONTHS) and (b1[3:].isdigit()):     # Catch 1902,jun22
            return b0.strip(), b1[:3].lower().strip(), b1[3:].strip()

    if len(bits) == 1:   # Try to handle cases with no spaces or commas, eg '1902jun22'
        ybit = ''
        mbit = ''
        dbit = ''
        bad = False
        for x in datestring.strip():
            if x == ' ':
                continue
            elif x.isdigit():
                if mbit:
                    dbit += x    # First group of digits before any letters
                else:
                    ybit += x    # Second group of digits, after the month section
            else:
                if ybit and (not dbit):
                    mbit += x    # First group of letters, after exactly one group of digits
                else:
                    bad = True    # non-digits before any digits, or after the second group of digits
        if not bad:
            return ybit, mbit.lower().strip(), dbit

    return None    # Catch anything else


def dostats(ra=None, dec=None, envdate='', analysis=''):
    global covmap, covcount
    if (ra is None) or (dec is None):
        return
    radeg = ra.deg
    decdeg = dec.deg
    ras = np.arange(int((radeg - PR) * 10), int((radeg + PR) * 10), dtype=np.int32)
    ras[(ras < 0)] += 3600
    ras[(ras >= 3600)] -= 3600

    decs = np.arange(int((decdeg - PR) * 10) + 900, int((decdeg + PR) * 10) + 900, dtype=np.int32)
    decs = decs[decs > 0]
    decs = decs[decs < 1800]

    area = np.meshgrid(ras, decs)
    covmap[tuple(area)] = 255
    covcount[tuple(area)] += 1


def do_plate(row=None, dofits=False, analysis=''):
    """
    Given a single row from the spreadsheet, process the file associated with that row.

    If dofits is True, write out converted FTIS files, if not, just do stats.

    :param row: A tuple of values from one row of the CSV file.
    :param dofits: If True, write out converted fits files.
    :param analysis: String that determines what analysis/plotting will be done
    :return: tuple with (sequence number, gotRADec, readTIFF, wroteFITS)
    """
    (seqnum,  # A: Integer sequence number
     scanner,  # B: Initials of the person doing the plate scan
     scandate,  # C: Date scanned (DD-MM-YYYY)
     platenum,  # D: String plate number (eg A23, 2435, etc)

     scannedyn,  # E: Y if the plate has been scanned, otherwise N
     platesize,  # F: L or S

     telescope,  # G: Telescope name (string)
     envdate,  # H: Date plate taken (YYYY,mmm,DD) where mmm is 'jan', 'feb', etc - eg '1901,may,23'
     equinoxstr,  # I: Coordinate equinox for ra/dec - if blank, use 1900 and save warning msg
     envra,  # J: RA, HHMMSS, unknown equinox (probably either B1900 or B1950)
     envdec,  # K: Dec, [+/-]DDMM[SS], unknown equinox (probably either B1900 or B1950)
     envobject,  # L: Object name string
     envobjtype,  # M: Object type string

     bookexpbasis,  # N: One of LST, UT, ST, WAST or AWST. If blank, use LST and save warning msg
     bookexptimes,  # O: could be float, int, or a list of space seperated numbers,
     # any of which might have parentheses around them. Eg '10 (20) 10'.
     # Note that the space might be missing, eg '10(20) 10'
     books1, booke1,  # P, Q: Exposure start and end times in LST, eg '8:56:32'
     books2, booke2,  # R, S:
     books3, booke3,  # T, U:
     books4, booke4,  # V, W:

     envguider,  # X: Initials of the person guiding the telescope
     envguidesettings,  # Y: Settings for the offset plate guider
     envha,  # Z: Hour Angle of the target (at exposure start?). Usually blank.
     envew,  # AA: Telescope tube is East (E) or West (W) of pier. Usually blank.
     envtemp,  # AB: Temperature in deg C (usually blank). Usually blank.
     envtransp,  # AC: Transp (float? No idea what it means). Usually blank.
     envscint,  # AD: Scint - presumably scintillation, in arcsec? Usually blank.
     envemulsion,  # AE: Emulstion type (string)
     envbox,  # AF: Box. Usually blank. (string)
     envfilter,  # AG: Filter. Usually blank. (string)
     envdevelop,  # AH: Develop. Usually blank. (string)
     envaperture,  # AI: Aperture. Usually blank. (string)
     envgrating,  # AJ: Grating. Usually blank. (string)
     envremarks,  # AK: Remarks written on envelope (string)

     othermarkings,  # AL: Any other markings (string)
     scannercomments,  # AM: Comments by the person scanning the plate (string)
     scanningweek,   # AN: Week number that the plate was scanned
     calc_ra,  # AO: RA in degrees, calculated from other columns (hand tuned?)
     calc_dec,  # AP: Dec in degrees, calculated from other columns (hand tuned?)
     *extras,    # Catch all for an extra columns
     ) = tuple(row)

    # Convert any strings that might have unicode characters into ASCII, making an attempt
    # to choose ASCII characters that are as close to the unicode ones as possible (stripping
    # accents, etc). FITS file headers must be in normal printable ASCII characters.
    envobject = ''.join([x for x in unidecode(envobject) if ' ' <= x <= '~'])
    envobjtype = ''.join([x for x in unidecode(envobjtype) if ' ' <= x <= '~'])
    envdevelop = ''.join([x for x in unidecode(envdevelop) if ' ' <= x <= '~'])
    envgrating = ''.join([x for x in unidecode(envgrating) if ' ' <= x <= '~'])
    envremarks = ''.join([x for x in unidecode(envremarks) if ' ' <= x <= '~'])
    othermarkings = ''.join([x for x in unidecode(othermarkings) if ' ' <= x <= '~'])
    scannercomments = ''.join([x for x in unidecode(scannercomments) if ' ' <= x <= '~'])
    ra = None
    if envra:
        if (len(envra.strip()) > 1) or (envra.strip() == '0'):
            rv = parseval(valstring=envra, ra=True)
            if type(rv) is not float:
                print("Seq# %s error in RA: %s" % (seqnum, rv))
            else:
                ra = Angle(rv, unit=hour)

    dec = None
    if envdec:
        if (len(envdec.strip()) > 1) or (envdec.strip() == '0'):
            rv = parseval(valstring=envdec, ra=False)
            if type(rv) is not float:
                print("Seq# %s error in Dec: %s" % (seqnum, rv))
            else:
                dec = Angle(rv, unit=deg)

    radiff = None
    try:
        calc_ra = float(calc_ra)
        if 0 <= calc_ra <= 360.0:
            if (ra is not None):
                radiff = (ra.hour * 15) - calc_ra
                if radiff > 1.0:
                    print("Seq# %s Warning: RA=%s renders to %1.1f hours, not calc_ra of %1.1f hours" % (seqnum, envra, ra.hour, calc_ra / 15.0))
            else:
                ra = Angle(calc_ra / 15.0, unit=hour)
                print('    Overriding weird RA of %s with value %1.4f hours from calc_ra' % (envra, calc_ra / 15.0))
        else:
            if ra is None:
                print('    ', end='')
            print('Seq# %s error: calc_ra of %1.4f degrees is outside 0 to 360 degree range' % (seqnum, calc_ra))
    except ValueError:
        if (ra is None) and envra.strip():
            print('    Could not find valid RA for %s' % envra)

    decdiff = None
    try:
        calc_dec = float(calc_dec)
        if -90.0 <= calc_dec <= 58.0:
            if (dec is not None):
                decdiff = dec.deg - calc_dec
                if decdiff > 1.0:
                    print("Seq# %s Warning: DEC=%s renders to %1.1f, not calc_dec of %1.1f" % (seqnum, envdec, dec.deg, calc_dec))
            else:
                dec = Angle(calc_dec, unit=deg)
                print('    Overriding weird DEC of %s with value %1.4f from calc_dec' % (envdec, calc_dec))
        else:
            if dec is None:
                print('    ', end='')
            print('Seq# %s error: calc_dec of %1.4f is outside -90 to +58 range' % (seqnum, calc_dec))
    except:
        if (dec is None) and envdec.strip():
            print('    Could not find valid DEC for %s' % envdec)

    dostats(ra=ra, dec=dec, envdate=envdate, analysis=analysis)

    if not dofits:
        # Return tuple is (sequence number, gotRADec, readTIFF, wroteFITS)
        return (seqnum, (ra is not None) and (dec is not None), False, False)

    tiff_filename = get_infilename(platenum=platenum)
    primary_hdu = fits.PrimaryHDU()
    got_tiff = False
    if (tiff_filename is None) or (not os.path.exists(tiff_filename)):
        # print('Seq# %s File %s not found' % (seqnum, tiff_filename))
        pass
    else:
        print('Seq# %s File %s found' % (seqnum, tiff_filename))
        if isdone(tiffname=tiff_filename, platenum=platenum):
            print('Seq# %s File %s already processed' % (seqnum,
                                                         get_outfilename(tiffname=tiff_filename, platenum=platenum)))
        else:
            tiff_img = Image.open(tiff_filename, 'r')
            primary_hdu = fits.PrimaryHDU(np.frombuffer(tiff_img.tobytes(), dtype=np.uint16))
            primary_hdu.data.shape = tiff_img.size
            got_tiff = True

    cover_filename = get_coverfilename(platenum=platenum)
    cover_hdu = None
    if (cover_filename is None) or (not os.path.exists(cover_filename)):
        # print('Seq# %s - Cover envelope scan file %s not found' % (seqnum, cover_filename))
        pass
    elif got_tiff:
        cover_img = Image.open(cover_filename, 'r')
        cover_hdu = fits.ImageHDU(np.frombuffer(cover_img.tobytes(), dtype=np.uint8))
        cover_hdu.data.shape = cover_img.size

    head = primary_hdu.header
    head.set('SIMPLE', 'T', 'File does conform to FITS standard')
    head.set('BITPIX', 16, '16-bit signed integers')
    head.set('CREATOR', 'plates.py v%s by Andrew.Williams@curtin.edu.au' % VERSION)
    obsname, obslocation, instrument = get_observatory(telescope)
    head.set('OBSERVAT', obsname, 'Observatory name')
    head.set('TELESCOP', telescope, 'Telescope name')
    head.set('LATITUDE', obslocation.lat.to_string(sep=':', pad=True), 'Latitude')
    head.set('LONGITUD', obslocation.lon.to_string(sep=':', pad=True), 'Longitude East')
    head.set('INSTRUME', instrument, 'Instrument')
    head.set('OBSERVER', envguider, 'Person guiding the telescope')
    head.set('OBJECT', envobject, 'Target name')
    head.set('GUIDER', envguidesettings, 'Guider Settings')
    head.set('SEQNUM', int(seqnum), 'Plate scanning sequence number')
    head.set('PLATE', platenum, 'Plate number')
    head.set('PLSIZE', platesize, 'Plate Size')
    head.set('EMULSION', envemulsion, 'Plate type')
    head.set('COMMENT', envremarks)
    head.set('COMMENT', othermarkings)
    head.set('HISTORY', "Scanned by '%s' on %s (week %s)" % (scanner, scandate, scanningweek))
    if scannercomments:
        head.set('HISTORY', scannercomments)

    if envra:
        head.set('RASTRING', envra)
    if ra is not None:
        head.set('RA', ra.hour, 'Right Ascension in hours')

    if envdec:
        head.set('DESTRING', envdec)
    if dec is not None:
        head.set('DEC', dec.deg, 'Declination in degrees')

    if envha:
        head.set('ENVHA', envha, 'Hour Angle, from plate envelope')
    if envew:
        head.set('COMMENT', "Telescope is %s of pier" % envew)
    if envtemp:
        head.set('COMMENT', "Temperature %s from envelope" % envtemp)
    if envtransp:
        head.set('COMMENT', "Transparency %s from envelope" % envtransp)
    if envscint:
        head.set('COMMENT', "Scintillation %s from envelope" % envscint)
    if envbox:
        head.set('COMMENT', "Box %s from envelope" % envbox)
    if envfilter:
        head.set('FILTER', envfilter, 'Filter name')
    if envdevelop:
        head.set('COMMENT', "Developed %s from envelope" % envdevelop)
    if envaperture:
        head.set('COMMENT', "Aperture %s from envelope" % envaperture)
    if envgrating:
        head.set('COMMENT', "Grating %s from envelope" % envgrating)

    books1T, books2T, books3T, books4T = None, None, None, None
    booke1T, booke2T, booke3T, booke4T = None, None, None, None

    head.set('ENV_DATE', envdate, 'Date written on envelope')
    try:
        if '&' in envdate:
            # print('Plate taken over multiple days - using first day: %s' % envdate)
            result = parsedate(envdate[:envdate.find('&')])
        elif '-' in envdate:
            # print('Plate taken over multiple days - using first day: %s' % envdate)
            result = parsedate(envdate[:envdate.find('-')])
        elif '/' in envdate:
            # print('Plate taken over multiple days - using first day: %s' % envdate)
            result = parsedate(envdate[:envdate.find('/')])
        else:
            result = parsedate(envdate)

        if result is not None:
            yearstring, monthname, daystring = result
        else:
            raise DateError

        year, month, day = int(yearstring), MONTHS[monthname.lower().strip()], int(daystring)
        tmpstring = bookexpbasis.strip().upper()
        if tmpstring in ['LST', 'ST'] or ('LST' in othermarkings):
            basis = 'lst'
            tmesg_text = 'Local Sidereal Time on envelope'
        elif tmpstring in ['UT', 'UTC', 'GMT']:
            basis = 'utc'
            tmesg_text = 'UTC converted to Local Sidereal Time'
        elif tmpstring in ['AWST', 'WAST']:
            basis = 'awst'
            tmesg_text = 'AWST converted to Local Sidereal Time'
        else:
            basis = 'lst'
            tmesg_text = 'No time zone on envelope, assumed LST'

        if books1:
            if ':' not in books1:
                books1 = Angle(float(books1) * 24, unit=hour).to_string(sep=':', precision=0)
            books1T = tstring_to_timerec(basis=basis,
                                         year=year, month=month, day=day,
                                         tstring=books1,
                                         obslocation=obslocation)
            head.set('DATE-OBS', '%4d-%02d-%02d' % (year, month, day), 'Date of observation')
            head.set('TIME-OBS', books1T.iso.split()[1], 'UTC at start of exposure')
            head.set('TIMESYS', 'UTC', 'Time System is UTC')
            head.set('LST', books1T.sidereal_time('apparent').to_string(unit=hour,
                                                                        sep=':',
                                                                        pad=True),
                     tmesg_text)
            head.set('JD', books1T.jd, 'Julian Date at start of exposure')
            head.set('MJD-OBS', books1T.mjd, 'Modified Julian Date at start of exposure')
            if booke1:
                if ':' not in booke1:
                    booke1 = Angle(float(booke1) * 24, unit=hour).to_string(sep=':', precision=0)
                booke1T = tstring_to_timerec(basis=basis,
                                             year=year, month=month, day=day,
                                             tstring=booke1,
                                             obslocation=obslocation)
                exptime = int(booke1T.gps - books1T.gps + 0.5)

                if books2 and booke2:
                    if ':' not in books2:
                        books2 = Angle(float(books2) * 24, unit=hour).to_string(sep=':', precision=0)
                    if ':' not in booke2:
                        booke2 = Angle(float(booke2) * 24, unit=hour).to_string(sep=':', precision=0)
                    books2T = tstring_to_timerec(basis=basis,
                                                 year=year, month=month, day=day,
                                                 tstring=books2,
                                                 obslocation=obslocation)
                    booke2T = tstring_to_timerec(basis=basis,
                                                 year=year, month=month, day=day,
                                                 tstring=booke2,
                                                 obslocation=obslocation)
                    exptime += int(booke2T.gps - books2T.gps + 0.5)

                if books3 and booke3:
                    if ':' not in books3:
                        books3 = Angle(float(books3) * 24, unit=hour).to_string(sep=':', precision=0)
                    if ':' not in booke3:
                        booke3 = Angle(float(booke3) * 24, unit=hour).to_string(sep=':', precision=0)
                    books3T = tstring_to_timerec(basis=basis,
                                                 year=year, month=month, day=day,
                                                 tstring=books3,
                                                 obslocation=obslocation)
                    booke3T = tstring_to_timerec(basis=basis,
                                                 year=year, month=month, day=day,
                                                 tstring=booke3,
                                                 obslocation=obslocation)
                    exptime += int(booke3T.gps - books3T.gps + 0.5)

                if books4 and booke4:
                    if ':' not in books4:
                        books4 = Angle(float(books4) * 24, unit=hour).to_string(sep=':', precision=0)
                    if ':' not in booke4:
                        booke4 = Angle(float(booke4) * 24, unit=hour).to_string(sep=':', precision=0)
                    books4T = tstring_to_timerec(basis=basis,
                                                 year=year, month=month, day=day,
                                                 tstring=books4,
                                                 obslocation=obslocation)
                    booke4T = tstring_to_timerec(basis=basis,
                                                 year=year, month=month, day=day,
                                                 tstring=booke4,
                                                 obslocation=obslocation)
                    exptime += int(booke4T.gps - books4T.gps + 0.5)

                head.set('EXPTIME', exptime, 'Total exposure time in seconds')
    except DateError:
        print('Seq# %s missing or invalid date: %s' % (seqnum, envdate))
        head.set('DATE-OBS', envdate, 'Date of observation')
    except:
        head.set('DATE-OBS', envdate, 'Date of observation')
        print("Seq# %s error in parsing date/times: %s" % (seqnum, traceback.format_exc()))

    if books1T and booke1T:
        head.set('LSTART1', books1T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at start of exposure')
        head.set('LEND1', booke1T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at end of exposure')

    if books2T and booke2T:
        head.set('LSTART2', books2T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at start of second exposure')
        head.set('LEND2', booke2T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at end of second exposure')

    if books3T and booke3T:
        head.set('LSTART3', books3T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at start of third exposure')
        head.set('LEND3', booke3T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at end of third exposure')

    if books4T and booke4T:
        head.set('LSTART4', books4T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at start of fourth exposure')
        head.set('LEND4', booke4T.sidereal_time('apparent').to_string(unit=hour, sep=':', pad=True), 'LST at end of fourth exposure')

    stime = books1T
    midtime = stime   # Fall back to start time, if we don't have an endtime
    if books1T is not None:
        if booke4T is not None:
            etime = booke4T
        elif booke3T is not None:
            etime = booke3T
        elif booke2T is not None:
            etime = booke2T
        else:
            etime = booke1T

        if etime is not None:
            midtime = stime + (etime - stime) / 2.0   # Mid-point of all shutter open times.

    try:
        equinox = int(equinoxstr)
        if equinox == 1900:
            eq_msg = 'Coordinate equinox 1900'
            frame = astropy.coordinates.FK4(equinox='B1900', obstime=midtime)
        elif equinox == '1950':
            eq_msg = 'Coordinates in B1950'
            frame = astropy.coordinates.FK4(equinox='B1950', obstime=midtime)
        elif equinox == '2000':
            eq_msg = 'Coordinates in J2000'
            frame = astropy.coordinates.FK4(equinox='J2000', obstime=midtime)
        else:
            eq_msg = 'Equinox given as %d' % equinox
            frame = astropy.coordinates.FK4(equinox='B1950', obstime=midtime)
    except:
        equinox = 1900
        eq_msg = 'Equinox missing, assuming 1900'
        frame = astropy.coordinates.FK4(equinox='B1900', obstime=midtime)

    head.set('EQUINOX', equinox, eq_msg)

    if (ra is not None) and (dec is not None) and (midtime is not None):
        target = SkyCoord(ra=ra, dec=dec,
                          frame=frame,
                          location=obslocation)
        target_app = target.transform_to('altaz')
        head.set('ALT', int(target_app.alt.deg * 10) / 10, 'Altitude in Degrees')
        head.set('AZ', int(target_app.az.deg * 10) / 10, 'Azimuth in Degrees')
        if target_app.alt.deg < 0:
            print("Seq# %s has Alt=%1.1f, Az=%1.1f and is below horizon" % (seqnum,
                                                                            target_app.alt.deg,
                                                                            target_app.az.deg))
            head.set('COMMENT', 'WARNING - Altitude of %1.1f deg is below horizon' % target_app.alt.deg)

    if cover_hdu:
        hdulist = fits.HDUList([primary_hdu, cover_hdu])
    else:
        hdulist = fits.HDUList([primary_hdu])

    if got_tiff:
        hdulist.writeto(get_outfilename(tiffname=tiff_filename, platenum=platenum), overwrite=True)
    hdulist[0].header.tofile(get_headerfilename(tiffname=tiff_filename, platenum=platenum),
                             sep='\r\n',
                             padding=False,
                             overwrite=True)
    # Return tuple is (sequence number, gotRADec, readTIFF, wroteFITS)
    return (seqnum, (ra is not None) and (dec is not None), True, True)


def do_all(fname='', dofits=False, analysis=''):
    """
    Given a CSV file name, process all the TIFF files that have not yet been processed.

    :param fname:
    :param dofits: If True, write out converted fits files.
    :param analysis: String that determines what analysis/plotting will be done
    :return:
    """
    results = []
    with open(fname, newline='') as csvfile:
        spamreader = csv.reader(csvfile)
        rownum = 0
        for row in spamreader:
            rownum += 1
            (seqnum, *_) = tuple(row)
            if not seqnum.strip().isdigit():
                continue  # Skip header lines
            results.append(do_plate(row=row, dofits=dofits, analysis=analysis))
    return results


def genplots(count=0, count_radec=0, count_tiff=0, count_fits=0):
    imgmap = Image.frombytes(mode='L', size=(1800, 3600), data=covmap.astype('uint8') * 4).transpose(Image.ROTATE_90)
    imgmap.save('C:/Data/Plates/covmap.png')
    imgcount = Image.frombytes(mode='L', size=(1800, 3600), data=covcount.astype('uint8') * 4).transpose(
        Image.ROTATE_90)
    imgcount.save('C:/Data/Plates/covcount.png')

    # figure out the mapping of pixel coords to lat/lon for the image
    lon = np.linspace(np.pi, -np.pi, covcount.shape[0])
    lat = np.linspace(-np.pi / 2., np.pi / 2., covcount.shape[1])
    Lon, Lat = np.meshgrid(lon, lat)

    c = pyplot.get_cmap(name='jet')
    c.set_bad('white')

    fig = pyplot.figure(figsize=(16, 8), dpi=100)
    ax = fig.add_subplot(111, projection="mollweide")
    ax.set_title('Astrographic plate coverage: %d plates, %d with valid RA/Dec' % (count, count_radec))
    ax.set_xticks(XT)
    ax.set_xticklabels(XTL)
    ax.set_ylabel('Declination', fontsize='large')
    data = covcount.transpose().astype(np.float32)
    data[data == 0] = np.NaN
    im = ax.pcolormesh(Lon, Lat, data, cmap=c, vmax=20)

    # add a color bar
    cb = pyplot.colorbar(im, ax=ax, shrink=0.6)
    formatter = mticker.FuncFormatter(lambda y, _:'{:d}'.format(int(y)) + ('+' if int(y) == 20 else ''))
    cb.ax.yaxis.set_major_formatter(formatter)
    cb.set_label('Number of plates', fontsize=16)

    # Set a grid and make the axes and tick markers larger and more bold
    ax.grid(True, lw=1.5, ls=':')
    pyplot.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    pyplot.savefig('c:/data/plates/covcount-map.png', dpi=100)


# TODO - Write a copy of the header to a file on its own, so that Matt can send them to me
#      - Do proper logging, so I can see more about what's happening in the log file.
#      - Decide what to do about a return value for the do_plate() function


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument('fnames', default='', nargs='+',
                        help='CSV file/s to use as input')
    parser.add_argument('--dofits', default=False, action='store_true',
                        help='If specified, convert any TIFF files found to FITS files')
    parser.add_argument('--analysis', default='',
                        help='If specified, carry out the specified stats analysis')
    parser.add_argument('--fitsdir', default=FITSDIR,
                        help='Directory to write FITS files to - default %s' % FITSDIR)
    options = parser.parse_args()

    covmap = np.zeros(shape=(3600, 1800), dtype=np.int32)
    covcount = np.zeros(shape=(3600, 1800), dtype=np.int32)

    print(options.fnames)
    fnames = []
    for fname in options.fnames:
        fnames += glob.glob(fname)   # Windows shell doesn't do wildcard expansion, so do it here.

    results = []
    for fname in fnames:
        results += do_all(fname, dofits=options.dofits, analysis=options.analysis)

    count = 0
    count_radec = 0
    count_tiff = 0
    count_fits = 0
    for row in results:
        (sequnum, gotRADec, readTIFF, wroteFITS) = row
        count += 1
        if gotRADec:
            count_radec += 1
        if readTIFF:
            count_tiff += 1
        if wroteFITS:
            count_fits += 1

    genplots(count=count, count_radec=count_radec, count_tiff=count_tiff, count_fits=count_fits)
