from astropy.io import fits, ascii
from astropy.table import Table, Column, vstack, MaskedColumn, QTable
import subprocess
import os 

#WISE 
def download_allwise(): 
    print(f"Downloading ALLWISE data from IRSA...")
    try:
        os.system("python wise_download.py") 
    except Exception as e:
        print(f"Error running wise_download.py: {e}")
def load_allwise(): 
    print(f"Loading ALLWISE data...")
    try: 
       hdu_list = fits.open("./CATALOGS/ALLWISE/ALLWISE.fits")
       allwise = Table(hdu_list[1].data)
       hdu_list.close()
       print(f"ALLWISE data loaded successfully.")
       return allwise
    except Exception as e: 
        print(f"Error loading ALLWISE data: {e}")
        return None
#eRASS1 
def download_erass1(): 
    print(f"Downloading eROSITA DR1 data from MPE...")
    try: 
        url = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Main.tar.gz"
        output_file = "./CATALOGS/eRASS1.tar.gz"
        subprocess.run(["curl", "-L", "-o", output_file, url], check=True)
    except Exception as e:
        print(f"Error downloading eRASS1 tar.gz file: {e}")
    try: 
        subprocess.run(['gunzip', output_file], check = True)
    except Exception as e:
        print(f"Error unzipping eRASS1 tar.gz file: {e}")
    try:
        tar_path = "./CATALOGS/eRASS1.tar"
        with tarfile.open(tar_path, 'r') as tar:
            tar.extractall(path="./CATALOGS")
        os.remove("./CATALOGS/eRASS1.tar")
        print(f"eRASS1 data downloaded and extracted.")
    except Exception as e:
        print(f"Error extracting eRASS1 tar file: {e}")
def load_erass1():
    print(f"Loading eRASS1 data...")
    try: 
        hdu_list = fits.open("./CATALOGS/eRASS1_Main.v1.1.fits")
        erass1 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"eRASS1 data loaded successfully.")
        return erass1
    except Exception as e: 
        print(f"Error loading eRASS1 data: {e}")
        return None
#RADIO catalogs 
def load_sumss():
    print(f"Loading SUMSS data...")
    try: 
        hdu_list = fits.open("./CATALOGS/radio/SUMSS.fits")
        sumss = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"SUMSS data loaded successfully.")
        return sumss
    except Exception as e: 
        print(f"Error loading SUMSS data: {e}")
        return None 
def load_first(): 
    print(f"Loading FIRST data...")
    try: 
        hdu_list = fits.open("./CATALOGS/radio/FIRST.fits")
        first = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"FIRST data loaded successfully.")
        return first
    except Exception as e: 
        print(f"Error loading FIRST data: {e}")
        return None   
def load_nvss(): 
    print(f"Loading NVSS data...")
    try: 
        hdu_list = fits.open("./CATALOGS/radio/NVSS.fits")
        nvss = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"NVSS data loaded successfully.")
        return nvss
    except Exception as e: 
        print(f"Error loading NVSS data: {e}")
        return None
def load_radio_catalogs():
    sumss = load_sumss()
    first = load_first()
    nvss = load_nvss()
    RADIO_instrument_inf = ascii.read("./CATALOGS/radio/radiocatalog_frequencies.txt", delimiter = ',')
    return sumss, first, nvss, RADIO_instrument_inf
def load_tevcat(): 
    print(f"Loading TeVCat data...")
    try: 
        tevcat = ascii.read("./CATALOGS/gammaray/TeVCAT_HBLs.csv")
        return tevcat
    except Exception as e: 
        print(f"Error loading TeVCat data: {e}")
        return None
def load_2fhl(): 
    print(f"Loading 2FHL data...")
    try: 
        hdu_list = fits.open("./CATALOGS/gammaray/2FHL.fits")
        fhl2 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"2FHL data loaded successfully.")
        return fhl2
    except Exception as e: 
        print(f"Error loading 2FHL data: {e}")
        return None
def load_3fh(): 
    print(f"Loading 3FHL data...")
    try: 
        hdu_list = fits.open('./catalogs/gammaray/3FHL.fits')
        fhl3 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"3FHL data loaded successfully.")
        return fhl3
    except Exception as e: 
        print(f"Error loading 3FHL data: {e}")
        return None
def load_4FGL():
    print(f"Loading 4FGL data...")
    try: 
        hdu_list = fits.open('./catalogs/gammaray/4FGL_DR4.fits')
        fgl4 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"4FGL data loaded successfully.")
        return fgl4
    except Exception as e: 
        print(f"Error loading 4FGL data: {e}")
        return None
def load_4lac(): 
    print(f"Loading 4LAC data...")
    try: 
        hdu_list = fits.open('./catalogs/gammaray/4LAC-DR3.fits')
        lac4 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"4LAC data loaded successfully.")
        return lac4
    except Exception as e: 
        print(f"Error loading 4LAC data: {e}")
        return None
def load_1cgh():
    print(f"Loading 1CGH data...")
    try: 
        hdu_list = fits.open('./catalogs/gammaray/1CGH_Preliminary_V2.0.fits')
        cgh1 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"1CGH data loaded successfully.")
        return cgh1
    except Exception as e: 
        print(f"Error loading 1CGH data: {e}")
        return None
def load_bzcat(): 
    print(f"Loading BZCAT data...")
    try: 
        hdu_list = fits.open('./catalogs/multiwavelength/bzcat.fits')
        bzcat = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"BZCAT data loaded successfully.")
        return bzcat
    except Exception as e: 
        print(f"Error loading BZCAT data: {e}")
        return None
def load_3hsp(): 
    print(f"Loading 3HSP data...")
    try: 
        hdu_list = fits.open('./catalogs/multiwavelength/3HSP.fits')
        hsp3 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"3HSP data loaded successfully.")
        return hsp3
    except Exception as e: 
        print(f"Error loading 3HSP data: {e}")
        return None
def load_whsp(): 
    print(f"Loading WHSP 1 and 2 data...")
    try: 
        hdu_list = fits.open('./catalogs/multiwavelength/1whsp.fits')
        whsp1 = Table(hdu_list[1].data)
        hdu_list.close()
        hdu_list = fits.open('./catalogs/multiwavelength/2whsp.fits')
        whsp2 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"WHSP 1 and 2 data loaded successfully.")
        return whsp1, whsp2
    except Exception as e: 
        print(f"Error loading WHSP data: {e}")
        return None, None
def load_costamante2020(): 
    print(f"Loading C20 data...")
    try: 
        C20 = ascii.read('./catalogs/multiwavelength/C20.txt')
        C20 = C20[0][:]
        return C20 
    except Exception as e:
        print(f"Error loading C20 data: {e}")
        return None
def load_dabrusco2019(): 
    print(f"Loading D19 Table 3 data...")
    try: 
        hdu_list = fits.open('./catalogs/multiwavelength/dabruscot3.fits')
        dabruscot3 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"D19 Table 3 data loaded successfully.")
        hdu_list = fits.open('./catalogs/multiwavelength/dabruscot5.fits')
        dabruscot5 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"D19 Table 5 data loaded successfully.")
        return dabruscot3, dabruscot5
    except Exception as e: 
        print(f"Error loading D19 data: {e}")
        return None, None
def load_lainez2025(): 
    print(f"Loading L25 data...")
    try: 
        L25 = ascii.read('./catalogs/multiwavelength/L25.txt')
        L25 = L25[0][:]
        return L25
    except Exception as e:
        print(f"Error loading L25 data: {e}")
        return None
def load_marchesi2025(): 
    print(f"Loading MAR25 data...")
    try:
        hdu_list = fits.open('./catalogs/multiwavelength/Marchesi_Xray_Bl_AA_25.fits')
        MAR25 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"MAR25 data loaded successfully.")
        return MAR25
    except Exception as e: 
        print(f"Error loading MAR25 data: {e}")
        return None
def load_massaro2013(): 
    print(f"Loading MAS13 data...")
    try: 
        hdu_list = fits.open('./catalogs/multiwavelength/massaro_selection.fits')
        MAS13 = Table(hdu_list[1].data)
        hdu_list.close()
        print(f"MAS13 data loaded successfully.")
        return MAS13
    except Exception as e: 
        print(f"Error loading MAS13 data: {e}")
        return None
def load_analysis_catalogs(): 
    fhl2 = load_2fhl()
    fhl3 = load_3fh()
    fgl4 = load_4FGL()
    lac4=load_4lac()
    cgh1 = load_1cgh()
    bzcat = load_bzcat()
    hsp3= load_3hsp()
    whsp1, whsp2 = load_whsp()
    c20 = load_costamante2020()
    dt3, dt5 = load_dabrusco2019()
    l25= load_lainez2025()
    mar25 = load_marchesi2025()
    mas13 = load_massaro2013()
    return fhl2, fhl3, fgl4, lac4, cgh1, bzcat, hsp3, whsp1, whsp2, c20, dt3, dt5, l25, mar25, mas13