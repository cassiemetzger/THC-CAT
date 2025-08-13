import subprocess
import glob
import os
import numpy as np
import sys 
import astropy
from astropy.io import fits
from astropy.table import Table, Column, vstack, MaskedColumn, QTable
d_url = "http://irsa.ipac.caltech.edu/data/download/wise-allwise/"
file_names = ["README.txt", "wise-allwise-cat.bz2.sizes.txt"]
answer = input("Which ALLWISE files would you like to download? Either enter 'all' to download all 48 files (you will need 349 GB of storage to store the compressed files and 1141 GB to store the uncompressed files), or enter the number of the files you wish to download in the following format: 01,02,03,...48 \n")
download_full_files = False #turned off by default to avoid breaking computer :( set to True when you're confident you have the space & time to download the full catalog
columns_to_keep = [0, 1, 2, 16, 17, 20, 21, 24, 25] #edit to include the columns you want to keep 
header = ["designation", "ra", "dec", "w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro", "w3mpro", "w3sigmpro"]
lines = 0
if answer == 'all':
    numbers = [f"{i:02d}" for i in range(1, 49)]
    for i in range(len(numbers)): 
        file_names.append(f'wise-allwise-cat-part{numbers[i]}.bz2')
    for file_name in file_names:
        print(f"Downloading file {file_name}...")
        url = d_url + file_name
        output_file = "./CATALOGS/ALLWISE/" + file_name
        if download_full_files:
            try:
                subprocess.run([
                "curl", "-L", "-o", output_file, url
                ], check=True)
                print(f"Full download completed: {output_file}")
                subprocess.run(["bunzip2", output_file], check=True)
                name = file_name
                num = name[21:23]
                os.rename(f"./CATALOGS/ALLWISE/{file_name[:24]}", f"./CATALOGS/ALLWISE/part{num}.txt")
                os.remove(output_file)  # Remove the original file after decompression
                with open(f"./CATALOGS/ALLWISE/part{num}.txt", 'r') as infile, open(f"./CATALOGS/ALLWISE/part{num}_filtered.txt", 'w') as outfile:
                    outfile.write('|'.join(header)+ '\n') 
                    for line in infile: 
                        lines +=1
                        line = line.rstrip('\n')
                        fields = line.split('|')
                        if(len(fields) == 299):
                            for i in columns_to_keep[1:]:
                                if(fields[i] == ''): 
                                    fields[i] = np.nan
                                    fields[i] = float(fields[i])
                                else: 
                                    fields[i] = float(fields[i])
                            if(fields[16]< 14.3 and fields[20]< 13.8 and fields[24]< 12.2):
                                     if((fields[20] - fields[24]) > 1.25 and (fields[20]-fields[24])< 2.4 and (fields[16]-fields[20])> 0.05 and (fields[16]-fields[20] < 0.9)):
                                         selected = [str(fields[i]) for i in columns_to_keep if i < len(fields)]
                                         outfile.write('|'.join(selected) + '\n')
                        else: 
                            print(f"File download incomplete")
                            sys.exit()
                    os.remove(output_file)
            except subprocess.CalledProcessError as e:
                print(f"Error running curl: {e}")
        else: 
            range_header = "0-10485759" #downloads only the first 10 MB of each file
            try: 
                subprocess.run([
                "curl", "-L", "-r", range_header, "-o", output_file, url
                ], check=True)
                print(f"Partial download completed: {output_file}")
                cmd = ["bzip2recover", output_file]
                try: 
                    subprocess.run(cmd, check=True)
                    print("File recovered successfully")
                    recovered_files = sorted(glob.glob('./CATALOGS/ALLWISE/rec*.bz2'))
                    os.remove(output_file)  # Remove the original file after recovery
                    for rf in recovered_files:
                         subprocess.run(["bunzip2", "-k", rf], check=True)
                         decompressed_file = rf[:-4]
                         new_name = decompressed_file + ".txt"
                         os.rename(decompressed_file, new_name)  
                         os.remove(rf)
                    folder = './CATALOGS/ALLWISE'
                    name = file_name
                    num = name[21:23]
                    output_file = os.path.join(folder, f'part{num}.txt')
                    rec_files = sorted(glob.glob(os.path.join(folder, 'rec*')))
                    with open(output_file, 'wb') as outfile:
                        for fname in rec_files:
                            print(f"Adding and removing {fname} ...")
                            with open(fname, 'rb') as infile:
                                outfile.write(infile.read())
                            os.remove(fname)
                    outfile.write('|'.join(header)+ '\n')  # Add header to the output file
                    with open(output_file, 'r') as infile, open(os.path.join(folder, f'/part{num}_filtered.txt'), 'w') as outfile:
                        outfile.write('|'.join(header)+ '\n') 
                        for line in infile: 
                            lines +=1
                            line = line.rstrip('\n')
                            fields = line.split('|')
                            if(len(fields) == 299):
                                for i in columns_to_keep[1:]:
                                    if(fields[i] == ''): 
                                        fields[i] = np.nan
                                        fields[i] = float(fields[i])
                                    else: 
                                        fields[i] = float(fields[i])
                                if(fields[16]< 14.3 and fields[20]< 13.8 and fields[24]< 12.2):
                                     if((fields[20] - fields[24]) > 1.25 and (fields[20]-fields[24])< 2.4 and (fields[16]-fields[20])> 0.05 and (fields[16]-fields[20] < 0.9)):
                                         selected = [str(fields[i]) for i in columns_to_keep if i < len(fields)]
                                         outfile.write('|'.join(selected) + '\n')
                    os.remove(output_file)
                except subprocess.CalledProcessError as e:
                    print(f"Error running bzip2recover: {e}")
            except subprocess.CalledProcessError as e:
                print(f"Error running curl: {e}")
else: 
    numbers = answer.split(',')
    print(numbers)
    for i in range(len(numbers)): 
        file_names.append(f'wise-allwise-cat-part{numbers[i]}.bz2')
    for file_name in file_names:
        print(f"Downloading file {file_name}...")
        url = d_url + file_name
        output_file = "./CATALOGS/ALLWISE/" + file_name
        if download_full_files:
            try:
                subprocess.run([
                "curl", "-L", "-o", output_file, url
                ], check=True)
                print(f"Full download completed: {output_file}")
                subprocess.run(["bunzip2", output_file], check=True)
                name = file_name
                num = name[21:23]
                os.rename(f"./CATALOGS/ALLWISE/{file_name[:24]}", f"./CATALOGS/ALLWISE/part{num}.txt")
                os.remove(output_file)  # Remove the original file after decompression
                with open(f"./CATALOGS/ALLWISE/part{num}.txt", 'r') as infile, open(f"./CATALOGS/ALLWISE/part{num}_filtered.txt", 'w') as outfile:
                    outfile.write('|'.join(header)+ '\n') 
                    for line in infile:
                        lines +=1 
                        line = line.rstrip('\n')
                        fields = line.split('|')
                        if(len(fields) == 299):
                            for i in columns_to_keep[1:]:
                                if(fields[i] == ''): 
                                    fields[i] = np.nan
                                    fields[i] = float(fields[i])
                                else: 
                                    fields[i] = float(fields[i])
                            if(fields[16]< 14.3 and fields[20]< 13.8 and fields[24]< 12.2):
                                     if((fields[20] - fields[24]) > 1.25 and (fields[20]-fields[24])< 2.4 and (fields[16]-fields[20])> 0.05 and (fields[16]-fields[20] < 0.9)):
                                         selected = [str(fields[i]) for i in columns_to_keep if i < len(fields)]
                                         outfile.write('|'.join(selected) + '\n')
                        else: 
                            print(f"File download incomplete")
                            sys.exit()
                    os.remove(output_file)
            except subprocess.CalledProcessError as e:
                print(f"Error running curl: {e}")
        else: 
            range_header = "0-10485759" #downloads only the first 10 MB of each file
            try: 
                subprocess.run([
                "curl", "-L", "-r", range_header, "-o", output_file, url
                ], check=True)
                print(f"Partial download completed: {output_file}")
                cmd = ["bzip2recover", output_file]
                try: 
                    subprocess.run(cmd, check=True)
                    print("file recovered successfully")
                    recovered_files = sorted(glob.glob('./CATALOGS/ALLWISE/rec*.bz2'))
                    for rf in recovered_files:
                         subprocess.run(["bunzip2", "-k", rf], check=True)
                         decompressed_file = rf[:-4]
                         new_name = decompressed_file + ".txt"
                         os.rename(decompressed_file, new_name)  
                         os.remove(rf)
                    os.remove(output_file)  # Remove the original file after recovery
                    folder = './CATALOGS/ALLWISE'
                    name = file_name
                    num = name[21:23]
                    output_file = os.path.join(folder, f'part{num}.txt')
                    rec_files = sorted(glob.glob(os.path.join(folder, 'rec*')))
                    with open(output_file, 'wb') as outfile:
                        for fname in rec_files:
                            #print(f"Adding and removing {fname} ...")
                            with open(fname, 'rb') as infile:
                                outfile.write(infile.read())
                            os.remove(fname)
                    with open(output_file, 'r') as infile, open(f'{folder}/part{num}_filtered.txt', 'w') as outfile:
                        outfile.write('|'.join(header)+ '\n')  # Add header to the output file
                        print(f"Filtering {output_file} ...")
                        for line in infile: 
                            lines+=1
                            line = line.rstrip('\n')
                            fields = line.split('|')
                            if(len(fields) == 299):
                                for i in columns_to_keep[1:]:
                                    if(fields[i] == ''): 
                                        fields[i] = np.nan
                                        fields[i] = float(fields[i])
                                    else: 
                                        fields[i] = float(fields[i])
                                if(fields[16]< 14.3 and fields[20]< 13.8 and fields[24]< 12.2):
                                     if((fields[20] - fields[24]) > 1.25 and (fields[20]-fields[24])< 2.4 and (fields[16]-fields[20])> 0.05 and (fields[16]-fields[20] < 0.9)):
                                         selected = [str(fields[i]) for i in columns_to_keep if i < len(fields)]
                                         outfile.write('|'.join(selected) + '\n')
                    os.remove(output_file)
                except subprocess.CalledProcessError as e:
                    print(f"Error running bzip2recover: {e}")
            except subprocess.CalledProcessError as e:
                print(f"Error running curl: {e}")
print(f"{lines} IR sources processed.")
#concatenate all filtered files into one file
print("Concatenating filtered files into one FITS file...")
input_dir = "./CATALOGS/ALLWISE/"
output_fits = os.path.join(input_dir, "ALLWISE.fits")
files = sorted(glob.glob(os.path.join(input_dir, "part*_filtered.txt")))
temp_txt = os.path.join(input_dir, "ALLWISE.txt")
with open(temp_txt, "w") as outfile:
    first_file = True
    for filename in files:
        with open(filename,"r") as infile: 
            for i, line in enumerate(infile): 
                if first_file or i > 0: 
                    outfile.write(line)
        first_file = False
        os.remove(filename)
table = Table.read(temp_txt, format='ascii.basic', delimiter='|')
table.write(output_fits, format='fits', overwrite=True)
print(f"FITS file created: {output_fits}")
os.remove("./CATALOGS/ALLWISE/ALLWISE.txt")  # Remove the temporary text file
