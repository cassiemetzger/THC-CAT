import numpy as np
import astropy 
from astropy.io import fits
from astropy.table import Table, Column, vstack, MaskedColumn, QTable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.io import ascii
from astroquery.vizier import Vizier
import os
import subprocess
import tarfile 
import astroquery 
from astroquery.heasarc import Heasarc
from astroquery.gaia import Gaia
from upload_catalogs import load_allwise, load_erass1, load_radio_catalogs, load_tevcat, load_analysis_catalogs
from table_formatting import get_THC_name, format_sigfigs, create_fits_table

def distance (m,b, x,y): 
   return np.abs(m * x + b-y)/np.sqrt(m**2 + 1**2)

def main(): 
    #constants 

    sol = 2.997925E14
    w1=sol/3.368
    fc_w1=0.9921
    F_w1 = 306.682
    e_energy = 1.6022e-16 #J 
    h = 6.626e-34 #Js
    e_v1 = e_energy/h
    w_v = 3E8/3.4E-6
    cs = 3E8

    #data download

    #wise_download()
    #erass1_download()
    #data upload
    allwise = load_allwise()
    erass1 = load_erass1()

    #matching sources 

    c = SkyCoord(ra = erass1['RA']*u.degree, dec = erass1['DEC']*u.degree)
    catalog = SkyCoord(ra = allwise['ra']*u.degree, dec = allwise['dec']*u.degree)
    idx, d2d, d3d = match_coordinates_sky(c, catalog)
    m = d2d.arcsec <= 5 #we choose a 5 arcsec matching radius
    matches = erass1[m]
    matches["WISE idx"] = idx[m]
    matches["WISE sep"] = d2d.arcsec[m]
    numbers = matches["WISE idx"]
    unique, counts = np.unique(numbers, return_counts=True)
    duplicates = unique[counts > 1] 
    if(len(duplicates) > 0):
        mask = np.array([False if i in duplicates else True for i in matches[f'WISE idx']])
        for val in duplicates: 
            idx = np.where(matches[f'WISE idx'] == val)[0]
            dis1 = matches[f'WISE sep'][idx[0]]
            dis2 = matches[f'WISE sep'][idx[1]]
            if(dis1 < dis2): 
                mask[idx[0]] = True
                mask[idx[1]] = False
            if(dis2 < dis1): 
                mask[idx[0]] = False
                mask[idx[1]] = True
            if(dis1 == dis2): 
                mask[idx[0]] = False
                mask[idx[1]] = True
        matches = matches[mask]
    print(f"{len(matches)} IR-X-ray source matches found.")

    #distance to best fit line

    x = allwise['w2mpro'][matches['WISE idx']]- allwise['w3mpro'][matches['WISE idx']]
    y = allwise['w1mpro'][matches['WISE idx']] - allwise['w2mpro'][matches['WISE idx']]
    m = 0.83809562126163
    b = -1.045969175938493
    matches['distance'] = distance(m, b, x, y)
    m = matches['distance'] < 0.21 
    matches = matches[m]
    print(f"{len(matches)} sources following best fit line cut.")

    #xir cut 

    matches['w1flux'] = F_w1/fc_w1*10**(-allwise['w1mpro'][matches['WISE idx']]/2.5)
    w1errh = allwise['w1mpro'][matches['WISE idx']] - allwise['w1sigmpro'][matches['WISE idx']]
    w1errl = allwise['w1mpro'][matches['WISE idx']] + allwise['w1sigmpro'][matches['WISE idx']]
    err_w1_h = F_w1/fc_w1*10**(-w1errh/2.5)-matches['w1flux']
    err_w1_l = F_w1/fc_w1*10**(-w1errl/2.5)-matches['w1flux']
    matches['w1sigflux'] = 0.5*(np.abs(err_w1_h)+ np.abs(err_w1_l))
    matches['ML_FLUX_1_JY'] = 1E+23 * matches['ML_FLUX_1'] * 1/(e_v1)
    matches['Alpha IRX'] = -np.log10(matches['w1flux']/matches['ML_FLUX_1_JY'])/np.log10(w_v/e_v1)
    m = matches['Alpha IRX'] <= 1
    matches = matches[m]

    #removing extended sources & galactic sources

    print(f"{len(matches)} sources following IRX cut.")
    m = matches["EXT_LIKE"] != 0 
    matches['IAUNAME'] = [name[7:] for name in matches['IAUNAME']]
    extended= Table([allwise['name'][matches['WISE idx'][m]], matches['IAUNAME'][m]], names = ('WISEA', 'eRASS1'))
    extended.write("./CATALOGS/TABLEA2.fits", format='fits', overwrite=True)
    m = matches['EXT_LIKE'] == 0 
    matches = matches[m]
    ra = allwise['ra'][matches['WISE idx']]
    dec = allwise['dec'][matches['WISE idx']]
    c = SkyCoord(ra*u.degree, dec*u.degree, frame = 'icrs')
    long = c.galactic.b.degree
    m = np.abs(long) > 10 
    matches = matches[m]
    print(f"{len(matches)} sources following galactic cut & extended source cut.")

    # radio matching 

    #load in SUMSS, FIRST, and NVSS catalogs 
    sumss, first, nvss, RADIO_instrument_inf = load_radio_catalogs()
    c = SkyCoord(ra = allwise['ra'][matches['WISE idx']]*u.degree,dec = allwise['dec'][matches['WISE idx']]*u.degree, frame = 'icrs')
    #initializing radio columns 
    matches['Radio flux'] =[np.nan]*len(matches)
    matches['Radio freq'] = [np.nan]*len(matches)
    matches['Radio catalog'] = Column(['None']*len(matches), dtype = 'U15')
    #starting with SUMSS
    catalog = SkyCoord(ra = sumss['RA']*u.degree, dec = sumss['DEC']*u.degree, frame = 'icrs')
    idx, d2d, d3d = match_coordinates_sky(c, catalog)
    m = d2d.arcsec <=5 #we choose a 5 arcsec matching radius
    matches['Radio flux'][m] = sumss['INT_FLUX_36_CM'][idx[m]]
    matches['Radio freq'][m] = 843E6 
    matches['Radio catalog'][m] = 'SUMSS'
    #now FIRST 
    catalog = SkyCoord(ra = first['RA']*u.degree, dec = first['DEC']*u.degree, frame = 'icrs')
    idx, d2d, d3d = match_coordinates_sky(c, catalog)
    m = d2d.arcsec <=5
    matches['Radio flux'][m] = first['FLUX_20_CM'][idx[m]]
    matches['Radio freq'][m] = 1.4E9
    matches['Radio catalog'][m] = 'FIRST'
    #finally, NVSS 
    catalog = SkyCoord(ra = nvss['RA']*u.degree, dec = nvss['DEC']*u.degree, frame = 'icrs')
    idx, d2d, d3d = match_coordinates_sky(c, catalog)
    m = d2d.arcsec <=5
    matches['Radio flux'][m] = nvss['FLUX_20_CM'][idx[m]]
    matches['Radio freq'][m] = 1.4E9
    matches['Radio catalog'][m] = 'NVSS'
    #search RADIO master catalog for remaining sources that dont have a flux in these catalogs 
    m = np.isnan(matches['Radio flux'])
    idx = np.array(range(len(matches)))
    idx = idx[m]
    radius = 5 * u.arcsec
    for i in idx: 
        coord = SkyCoord(ra= allwise['ra'][matches['WISE idx'][i]]*u.degree, dec = allwise['dec'][matches['WISE idx'][i]]*u.degree, frame = 'icrs')
        try: 
            table = Heasarc.query_region(coord, mission = "RADIO", radius = radius)
            table['SEARCH_OFFSET_'] = np.array([float(x[:5]) for x in table['SEARCH_OFFSET_']])
            if(np.min(table['SEARCH_OFFSET_']) <=0.0833333): 
                target = np.argmin(table['SEARCH_OFFSET_'])
                if(table['FLUX_20_CM'][target] != 0): 
                    matches['Radio flux'][i] = table['FLUX_20_CM'][target]
                    matches['Radio freq'][i] = cs/(20/100)
                    matches['Radio catalog'][i] = table['DATABASE_TABLE'][target]
                if(table['FLUX_6_CM'][target] != 0 and table['FLUX_20_CM'][target] == 0):
                    matches['Radio flux'][i] = table['FLUX_6_CM'][target]
                    matches['Radio freq'][i] = cs/(6/100)
                    matches['Radio catalog'][i] = table['DATABASE_TABLE'][target]
                if(table['FLUX_20_CM'][target] == 0 and table['FLUX_6_CM'][target] == 0):
                    matches['Radio flux'][i] = table['FLUX_OTHER'][target]
                    s = table['DATABASE_TABLE'][target]
                    matches['Radio catalog'][i] = s.strip()
                    m = RADIO_instrument_inf['Database'] == matches['Radio catalog'][i]
                    matches['Radio freq'][i] = (float(RADIO_instrument_inf['Selected Frequency (MHz)'][m]) *1E6)
        except Exception as e:
            continue
    m = np.isnan(matches['Radio flux'])
    missing_radio = Table([allwise['name'][matches['WISE idx'][m]], matches['IAUNAME'][m]], names = ('WISEA', 'eRASS1'))
    missing_radio.write("./CATALOGS/TABLEA3.fits", format='fits', overwrite=True)
    m = (np.isnan(matches['Radio flux']) == False)
    print(f"{len(matches[m])} sources following radio flux cut.")
    matches = matches[m]
    matches['Radio flux'] = [flux*1E-3 for flux in matches['Radio flux']] #mJy to Jy
    matches['Alpha RIR'] = -np.log10(matches['Radio flux']/matches['w1flux'])/np.log10(matches['Radio freq']/w_v)
    m = matches['Alpha RIR'] <= 0.43
    print(f"{len(matches[m])} sources following alpha_rir cut.")
    matches = matches[m]
    #exclude TeVCAT sources 
    tevcat = load_tevcat()
    c = SkyCoord(ra = allwise['ra'][matches['WISE idx']]*u.degree, dec = allwise['dec'][matches['WISE idx']]*u.degree, frame = 'icrs')
    catalog =SkyCoord(ra= tevcat['Simbad RA deg']*u.degree, dec = tevcat['Simbad DEC deg']*u.degree, frame = 'icrs') #we recommend using Simbad to retrieve the optical counterparts for the TeVCAT sources for better localization
    idx, d2d, d3d = match_coordinates_sky(c, catalog)
    m = d2d.arcsec > 5 #we want to exclude TeVCAT matches, we consider a match to be within 5 arcsec 
    matches = matches[m]
    print(f"{len(matches)} candidate sources found.")

    #constructing final table 
    thc_name = get_THC_name(allwise, matches)
    thc = Table([thc_name, allwise['name'][matches['WISE idx']], matches['IAUNAME'], allwise['ra'][matches['WISE idx']], allwise['dec'][matches['WISE idx']], allwise['w1mpro'][matches['WISE idx']], allwise['w1sigmpro'][matches['WISE idx']], allwise['w2mpro'][matches['WISE idx']], allwise['w2sigmpro'][matches['WISE idx']], allwise['w3mpro'][matches['WISE idx']], allwise['w3sigmpro'][matches['WISE idx']], matches['w1flux'], matches['w1sigflux'], matches['ML_FLUX_1'], matches['ML_FLUX_ERR_1'], matches['ML_FLUX_P4'], matches['ML_FLUX_ERR_P4'], matches['Radio flux'], matches['Radio freq'], matches['Radio catalog'], matches['Alpha IRX'], matches['Alpha RIR']], names = ("THC", "WISEA", "1eRASS", "ra", "dec", "w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro", "w3mpro", "w3sigmpro", "w1flux", "w1sigflux", "ML_FLUX_1", "ML_FLUX_ERR_1", "ML_FLUX_P4", "ML_FLUX_ERR_P4","Radio flux", "Radio freq", "Catalog", "Alpha IRX", "Alpha RIR"))
    thc['r_mag'] = [np.nan]*len(thc)
    for i in range(len(thc)): 
        c = SkyCoord(ra = thc['ra'][i]*u.degree, dec = thc['dec'][i]*u.degree, frame = 'icrs')
        try: 
            j = Gaia.cone_search_async(c, u.Quantity(5, u.arcsec))
            r = j.get_results()
            sel = np.argmin(r['dist'])
            thc['r_mag'][i] = r['phot_rp_mean_mag'][sel]
        except Exception as e:
            print(f"Error retrieving Gaia data for source {thc['ra'][i]}, {thc['dec'][i]}: {e}")
            continue
    #multiwavelength analysis
    fhl2, fhl3, fgl4, lac4, cgh1, bzcat, hsp3, whsp1, whsp2, c20, dt3, dt5, l25, mar25, mas13 = load_analysis_catalogs()
    c = SkyCoord(ra = thc['ra']*u.degree, dec = thc['dec']*u.degree, frame = 'icrs')
    #2fhl 
    catalog = SkyCoord(ra = fhl2['RAJ2000']*u.degree, dec = fhl2['DEJ2000']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.deg <= fhl2['Pos_err_95'][idx]
    thc['2FHL Assoc'] = Column(['']*len(thc), dtype = 'U17')
    thc['2FHL Assoc'][mask] = fhl2['Source_Name'][idx[mask]]
    #3fhl
    catalog = SkyCoord(ra = fhl3['RAJ2000']*u.degree, dec = fhl3['DEJ2000']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.deg <= fhl3['Conf_95_SemiMajor'][idx]
    thc['3FHL Assoc'] = Column(['']*len(thc), dtype = 'U18')
    thc['3FHL Assoc'][mask] = fhl3['Source_Name'][idx[mask]]
    #4fgl
    catalog = SkyCoord(ra = fgl4['RAJ2000']*u.degree, dec = fgl4['DEJ2000']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.deg <= fgl4['Conf_95_SemiMajor'][idx]
    thc['4FGL Assoc'] = Column(['']*len(thc), dtype = 'U18')
    thc['4FGL Assoc'][mask] = fgl4['Source_Name'][idx[mask]]
    thc['4FGL Assoc'] = [name.replace(' ', '') for name in thc['4FGL Assoc']]
    thc['4FGL Assoc'] = [name.replace('LJ', 'L J') for name in thc['4FGL Assoc']]
    #4lac
    catalog = SkyCoord(ra = lac4['RACdeg']*u.degree, dec = lac4['DECdeg']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.arcsec <= 5.0
    thc['4LAC Assoc'] = Column(['']*len(thc), dtype = 'U28')
    thc['4LAC Assoc'][mask] = lac4['Assoc1'][idx[mask]]
    #1cgh
    catalog = SkyCoord(ra = cgh1['RA_1CGH']*u.degree, dec = cgh1['DEC_1CGH']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.arcsec <= 5.0
    thc['1CGH Assoc'] = Column(['']*len(thc), dtype = 'U28')
    thc['1CGH Assoc'][mask] = cgh1['Counterpart_name'][idx[mask]]
    #bzcat
    catalog = SkyCoord(ra = bzcat['RA']*u.degree, dec = bzcat['DEC']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.arcsec <= 5
    thc['BZCAT Assoc'] = Column(['']*len(thc), dtype = 'U15')
    thc['BZCAT OTYPE'] = Column(['']*len(thc), dtype = 'U24')
    thc['BZCAT Assoc'][mask] = [name[10:] for name in bzcat['NAME'][idx[mask]]]
    thc['BZCAT OTYPE'][mask] = bzcat['OBJECT_TYPE'][idx[mask]]
    #3hsp
    catalog = SkyCoord(ra = hsp3['RAJ2000']*u.degree, dec = hsp3['DEJ2000']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.arcsec <= 5
    thc['3HSP Assoc'] = Column(['']*len(thc), dtype = 'U20')
    thc['3HSP Assoc'][mask] = hsp3['Name'][idx[mask]]
    #whsp
    catalog = SkyCoord(ra = whsp1['_RA']*u.degree, dec = whsp1['_DE']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.arcsec <= 5
    thc['1WHSP Assoc'] = Column(['']*len(thc), dtype = 'U15')
    thc['1WHSP Assoc'][mask] = whsp1['_1WHSP'][idx[mask]]
    catalog = SkyCoord(ra = whsp2['_RA']*u.degree, dec = whsp2['_DE']*u.degree, frame = 'icrs')
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    mask = d2d.arcsec <= 5
    thc['2WHSP Assoc'] = Column(['']*len(thc), dtype = 'U15')
    thc['2WHSP Assoc'][mask] = whsp2['_2WHSPJ'][idx[mask]]
    #costamante2020
    thc['C20'] = Column(['N']*len(thc), dtype = 'U1')
    for i in range(len(thc)): 
        target = thc['BZCAT Assoc'][i] 
        target = target.replace(" ", "")
        if(target != "None"): 
            target1 = target[0:9]
            for j in range(len(c20)): 
                c_t = c20[j]
                c_t1 = c_t[0:9]
                if(target1 == c_t1): 
                    target2 = target[10:14]
                    c_t2 = c_t[10:14]
                    if(target2 == c_t2): 
                        thc['C20'][i] = 'Y'
    #dabrusco2023
    thc['D19'] = Column(['N']*len(thc), dtype = 'U1')
    for i in range(len(thc)): 
        target = thc['WISEA'][i]
        target1 = target[0:10]
    for j in range(len(dt3)): 
        catalog = dt3['WISE'][j]
        catalog1 = catalog[0:10]
        if(target1 == catalog1): 
           target2 = target[11:21]
           catalog2 = catalog[11:21]
           if(target2 == catalog2): 
               thc['D19'][i] = 'Y'
    for k in range(len(dt5)): 
        catalog = dt5['WISE'][k]
        catalog1 = catalog[0:10]
        if(target1 == catalog1): 
            target2 = target[11:21]
            catalog2 = catalog[11:21]
            if(target2 == catalog2): 
                thc['D19'][i] = 'Y'
    #láinez2025
    thc['L25'] = Column(['N']*len(thc), dtype = 'U1')
    for i in range(len(thc)): 
        target = thc['4FGL Assoc'][i] 
        target1 = target[5:12]
        for j in range(len(l25)): 
            l_t = l25[j]
            l_t1 = l_t[0:7]
            if(target1 == l_t1): 
                target2 =target[13:17]
                l_t2 = l_t[8:12]
                if(target2 == l_t2): 
                    thc['L25'][i] = 'Y'
    #marchesi2025
    thc['MAR25'] = Column(['N']*len(thc), dtype = 'U1')
    for i in range(len(thc)): 
        target = thc['BZCAT Assoc'][i]
        if(target != "None"):
            target = target.replace(" ", "")
            target1 = target[0:9]
            for j in range(len(mar25)): 
                catalog = mar25['ID_5BZCAT'][j]
                catalog1 = catalog[0:9]
                if(target1 == catalog1): 
                    target2 = target[10:14]
                    catalog2 = catalog[10:14]
                    if(target2 == catalog2): 
                        thc['MAR25'][i] = 'Y'
    #massaro2013
    thc['MAS13'] = Column(['N']*len(thc), dtype = 'U1')
    for i in range(len(thc)): 
        name = thc['THC'][i]
        if(name in mas13['WISE']): 
            thc['MAS13'][i] = 'Y'
            idx = np.where(mas13['WISE'] == name)[0][0]

    #write final table
    thc['ra'] = [float(format_sigfigs(num, 5)) for num in thc['ra']]
    thc['dec'] = [float(format_sigfigs(num, 5)) for num in thc['dec']]
    thc['w1mpro'] = [float(format_sigfigs(num, 3)) for num in thc['w1mpro']]
    thc['w1sigmpro'] = [float(format_sigfigs(num, 3)) for num in thc['w1sigmpro']]
    thc['w2mpro'] = [float(format_sigfigs(num, 3)) for num in thc['w2mpro']]
    thc['w2sigmpro'] = [float(format_sigfigs(num, 3)) for num in thc['w2sigmpro']]
    thc['w3mpro'] = [float(format_sigfigs(num, 3)) for num in thc['w3mpro']]
    thc['w3sigmpro'] = [float(format_sigfigs(num, 3)) for num in thc['w3sigmpro']]
    thc['w1flux'] = thc['w1flux'] *1E3 # convert to mJy
    thc['w1flux'] = [float(format_sigfigs(num, 3)) for num in thc['w1flux']]
    thc['w1sigflux'] = thc['w1sigflux'] *1E3 # convert to mJy
    thc['w1sigflux'] = [float(format_sigfigs(num, 3)) for num in thc['w1sigflux']]
    thc['ML_FLUX_1'] = [float(format_sigfigs(num, 3)) for num in thc['ML_FLUX_1']]
    thc['ML_FLUX_ERR_1'] = [float(format_sigfigs(num, 3)) for num in thc['ML_FLUX_ERR_1']]
    thc['ML_FLUX_P4'] = [float(format_sigfigs(num, 3)) for num in thc['ML_FLUX_P4']]
    thc['ML_FLUX_ERR_P4'] = [float(format_sigfigs(num, 3)) for num in thc['ML_FLUX_ERR_P4']]
    thc['Radio flux'] = thc['Radio flux'] *1E3 # convert to mJy
    thc['Radio flux'] = [float(format_sigfigs(num, 2)) for num in thc['Radio flux']]
    thc['Alpha RIR'] = [float(format_sigfigs(num, 3)) for num in thc['Alpha RIR']]
    thc['Alpha IRX'] = [float(format_sigfigs(num, 3)) for num in thc['Alpha IRX']]
    thc['r_mag'] = [float(format_sigfigs(num, 3)) for num in thc['r_mag']]
    thc['Radio freq'] = thc['Radio freq'] * 1E-9 # convert to GHz 
    thc['Radio freq'] = [float(format_sigfigs(num, 3)) for num in thc['Radio freq']]
    thc.sort('ra')
    create_fits_table(thc)

main()