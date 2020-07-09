
import csv

import numpy

from PIL import Image

import astropy
from astropy import units
from astropy.coordinates import Angle, SkyCoord

import matplotlib
from matplotlib import pyplot

PR = 1.0   # Plate 'radius' - eg, 1.0 means a 2 x 2 degree square plate.

covmap = numpy.zeros(shape=(3600, 1800), dtype=numpy.int32)
covcount = numpy.zeros(shape=(3600, 1800), dtype=numpy.int32)

FNAME = 'C:/Data/Plates/AllCoords-20190910.csv'


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
            elif len(bdp) == 3:   # 123.4
                rv = parsehex("0%s:%s" % (bdp[0], bdp[1:]), ra=ra)
            elif len(bdp) == 4:  # 1234.5
                rv = parsehex("%s:%s.%s" % (bdp[:2], bdp[2:], adp), ra=ra)
            elif len(bdp) == 5:  # 12345.6
                rv = parsehex("0%s:%s:%s" % (s[0], s[1:3], s[3:]), ra=ra)
            elif len(bdp) == 6:  # 123456.7
                rv = parsehex("%s:%s:%s.%s" % (bdp[:2], bdp[2:4], bdp[4:], adp), ra=ra)
            elif (len(bdp) == 7) and (bdp[0] == '0'):  # 0123456.7
                rv = parsehex("%s:%s:%s.%s" % (bdp[1:3], bdp[3:5], bdp[5:], adp), ra=ra)
            else:
                rv = "weird number of chars before decimal point: %s" % valstring
        else:
            if len(s) <= 2:  # 1, 12
                try:
                    rv = float(s)
                except ValueError:
                    rv = "Invalid float conversion for %s" % valstring
            elif len(s) == 3: # 123
                rv = parsehex("0%s:%s" % (s[0], s[1:]), ra=ra)
            elif len(s) == 4:   # 1234
                rv = parsehex("%s:%s" % (s[:2], s[2:4]), ra=ra)
            elif len(s) == 5: # 12345
                rv = parsehex("0%s:%s:%s" % (s[0], s[1:3], s[3:]), ra=ra)
            elif len(s) == 6:   # 123456
                rv = parsehex("%s:%s:%s" % (s[:2], s[2:4], s[4:6]), ra=ra)
            elif (len(s) == 7) and (s[0] == '0'):  # 0123456
                rv = parsehex("%s:%s:%s" % (s[1:3], s[3:5], s[5:]), ra=ra)
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


def do_all(fname=FNAME):
    """
    Given a CSV file name, plot all the plate outlines.

    :param fname:
    :return:
    """
    with open(fname, newline='') as csvfile:
        spamreader = csv.reader(csvfile)
        rownum = 0
        for row in spamreader:
            rownum += 1
            if rownum < 2:
                continue  # Skip header line
            rastring, decstring = tuple(row)
            rastring = rastring.strip()
            decstring = decstring.strip()

            ra = None
            radeg = None
            if (len(rastring) > 1) or (rastring == '0'):
                rv = parseval(valstring=rastring, ra=True)
                if type(rv) is not float:
                    print("Line %d error in RA: %s" % (rownum, rv))
                else:
                    ra = Angle(rv, unit=units.hour)
            if ra is not None:
                radeg = ra.deg

            dec = None
            decdeg = None
            if (len(decstring) > 1) or (decstring == '0'):
                rv = parseval(valstring=decstring, ra=False)
                if type(rv) is not float:
                    print("Line %d error in Dec: %s" % (rownum,rv))
                else:
                    dec = Angle(rv, unit=units.deg)
            if dec is not None:
                decdeg = dec.deg

            if radeg and decdeg:
                ras = numpy.arange(int((radeg - PR) * 10), int((radeg + PR) * 10), dtype=numpy.int32)
                ras[(ras < 0)] += 3600
                ras[(ras >= 3600)] -= 3600

                decs = numpy.arange(int((decdeg - PR) * 10) + 900, int((decdeg + PR) * 10) + 900, dtype=numpy.int32)
                decs = decs[decs > 0]
                decs = decs[decs < 1800]

                area = numpy.meshgrid(ras, decs)
                covmap[tuple(area)] = 255
                covcount[tuple(area)] += 1


if __name__ == '__main__':
    do_all()
    imgmap = Image.frombytes(mode='L', size=(1800, 3600), data=covmap.astype('uint8') * 4).transpose(Image.ROTATE_90)
    imgmap.save('C:/Data/Plates/covmap.png')
    imgcount = Image.frombytes(mode='L', size=(1800, 3600), data=covcount.astype('uint8') * 4).transpose(Image.ROTATE_90)
    imgcount.save('C:/Data/Plates/covcount.png')

    # figure out the mapping of pixel coords to lat/lon for the image
    lon = numpy.linspace(numpy.pi, -numpy.pi, covcount.shape[0])
    lat = numpy.linspace(-numpy.pi / 2., numpy.pi / 2., covcount.shape[1])
    Lon, Lat = numpy.meshgrid(lon, lat)

    c = pyplot.get_cmap(name='jet')
    c.set_bad('white')

    fig = pyplot.figure(figsize=(16, 8), dpi=100)
    ax = fig.add_subplot(111, projection="mollweide")
    data = covcount.transpose().astype(numpy.float32)
    data[data == 0] = numpy.NaN
    im = ax.pcolormesh(Lon, Lat, data, cmap=c, vmax=20)

    # add a color bar
    cb = pyplot.colorbar(im, ax=ax, shrink=0.6)
    cb.set_label('Number of plates', fontsize=16)

    # Set a grid and make the axes and tick markers larger and more bold
    ax.grid(True, lw=1.5, ls=':')
    pyplot.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    pyplot.savefig('c:\data\plates\covcount-map.png', dpi=100)