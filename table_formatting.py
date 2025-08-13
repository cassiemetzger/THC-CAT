from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits, ascii
from astropy.table import Table, QTable
import numpy as np
import math
from decimal import Decimal, ROUND_HALF_UP, getcontext

def get_THC_name(allwise, matches): 
    ra = allwise['ra'][matches['WISE idx']]
    dec = allwise['dec'][matches['WISE idx']]
    c = SkyCoord(ra = ra*u.degree, dec = dec*u.degree, frame='fk5')
    strings = c.to_string('hmsdms')
    for i in range(len(strings)): 
        time_part, declination_part = strings[i].split()
        hours = time_part[0:2]
        minutes = time_part[3:5]
        degrees = declination_part[1:3]
        sign = declination_part[0]
        mins = declination_part[4:6]
        convert = float(mins)/60 
        mins = str(convert) 
        mins = mins[2:3]
        strings[i] = f'J{hours}{minutes}{sign}{degrees}{mins}'
    return strings


def format_sigfigs(num, sigfig):
    if num is None or (isinstance(num, float) and math.isnan(num)):
        return "nan"

    # Ensure context precision is high enough
    getcontext().prec = 10

    # Convert to Decimal for consistent rounding
    d = Decimal(str(num))
    rounded = d.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)

    # Format with exactly 2 digits after the decimal point
    return f"{rounded:.{sigfig}f}"
def create_fits_table(x): 
    tbl = QTable(x)
    tbl['ra'] = tbl['ra']*u.degree
    tbl['dec'] = tbl['dec']*u.degree
    tbl['w1mpro'] = tbl['w1mpro']*u.mag 
    tbl['w1sigmpro'] = tbl['w1sigmpro']*u.mag
    tbl['w2mpro'] = tbl['w2mpro']*u.mag
    tbl['w2sigmpro'] = tbl['w2sigmpro']*u.mag
    tbl['w3mpro'] = tbl['w3mpro']*u.mag
    tbl['w3sigmpro'] = tbl['w3sigmpro']*u.mag
    tbl['w1flux'] = tbl['w1flux']*u.mJy
    tbl['w1sigflux'] = tbl['w1sigflux']*u.mJy
    tbl['ML_FLUX_1'] = tbl['ML_FLUX_1']*u.erg/(u.cm**2*u.s)
    tbl['ML_FLUX_ERR_1'] = tbl['ML_FLUX_ERR_1']*u.erg/(u.cm**2*u.s)
    tbl['ML_FLUX_P4'] = tbl['ML_FLUX_P4']*u.erg/(u.cm**2*u.s)
    tbl['ML_FLUX_ERR_P4'] = tbl['ML_FLUX_ERR_P4']*u.erg/(u.cm**2*u.s)
    tbl['Radio flux'] = tbl['Radio flux']*u.mJy
    tbl['Radio freq'] = tbl['Radio freq']*u.GHz
    tbl['Alpha IRX'] = tbl['Alpha IRX']
    tbl['Alpha RIR'] = tbl['Alpha RIR']
    tbl['r_mag'] = tbl['r_mag']*u.mag
    table_hdu = fits.BinTableHDU(tbl) 
    header = table_hdu.header
    header['EXTNAME'] = 'SOURCE_CATALOG'
    header['AUTHOR'] = 'Cassie Metzger' 
    header['DATE'] = '2025-04-03' 
    header['WAVELENG'] = 'Multi'
    header['DATASRC'] = 'eROSITA + WISE + Assorted radio catalogs'
    header['DOI'] = '10.48550/arXiv.2501.12520' 
    header['WAVELENG'] = 'Multi'
    header['SELPRO'] = 'Sources have WISE magnitudes less than 14.3, 13.8, and 12.2 in the 3.4um, 4.6um, and 12um bands, respectively, fall within 5'' of an eRASS1 source, fall within 0.21 of the best fit line y = 0.84x -1.05, have alpha_IRX >= -1, are not extended in the X-ray wavelength, are not extragalactic, have a radio counterpart with 5'', and have alpha_RIR >= -0.43'
    tbl.write('./CATALOGS/THC_catalog.fits', format = 'fits', overwrite = True)
    print('Catalog saved as THC_catalog.fits')