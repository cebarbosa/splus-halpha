# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 18:06:40 2021

@author: 55119
"""

import glob, os, sys, numpy as np, matplotlib.pyplot as plt
from astropy.table import QTable
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import numpy as np
import splusdata
CUPOM: OFFEXTRA OFFEXTRA

path= ["C:/Users/55119/Dropbox/splus-halpha (1)/FCC_halpha/Fornax10201.pix"]

for x in range(len(path)):
    
    os.chdir(path[x])
    SPLUS-s29s31-50.22486407163766--37.134297664400115-256_R.fits (1)
    result = []
    hAlphaFlux = []
    bands = ['_F660_','_I_','_R_']
    i = 0 
    coef_F660 = 125.3
    coef_R = 1419.0
    fnu_R = 0
    for file in glob.glob("*_swp.fits"):
        
        for i in range(len(bands)):
            if bands[i] in file:
                mainOutput = sys.stdout
        
                fileName = "FILENAME  =  " + file
                openFile = fits.open(file)
                header = openFile[1].header
                data = openFile[1].data
                header = str(header)
                fileNameList = []
                fileNameList.append(fileName)
                headerList = []
                headerList.append(header)
                dataList = []
                dataList.append(data)
                result.append(fileNameList + headerList + dataList)
    
    for i in range(len(result)):
        
        
        fileName2 = str(result[i])
        fileName2.strip()
        fileName2 = fileName2.split('FILENAME')[1]
        fileName2 = fileName2.replace(" ", "")
        fileName2 = fileName2.split('=')[1]
        fileName2 = fileName2.split('.fits', 1)[0]
        
        
        magZP = str(result[i])
        magZP.strip()
        magZP = magZP.split('MAGZP')[1]
        magZP = magZP.replace(" ", "")
        magZP = magZP.split('=')[1]
        magZP = magZP.split('/', 1)[0]
        magZP = float(magZP)
        
        data = result[i][2]

        fnu = np.dot(data,(np.power(10, (-0.4 * (magZP + 48.6)))))
        
        
        fnu_F660 = fnu[1]
        fnu_R = fnu[0]
        fnu_I= fnu[2]
        
        
        fluxTwoBands = np.dot(coef_F660,((fnu_F660 - fnu_R)/(1-(coef_F660/coef_R))))
       # fluxThreeBands = ((fnu_R- fnu_I)-((alpha_R-alpha_I)/(alpha_F660 - alpha_I))*(fnu_F660 - fnu_I))/(beta_F66)*(alpha_I - alpha_R )- (        fluxThreeBands = ((fnu_R- fnu_I)-((alpha_R-alpha_I)/(alpha_F660 - alpha_I))*((fnu_F660 -fnu_I))/((beta_F660)*(alpha_I -alpha_R )-(beta_R))        
        hAlphaFlux.append("FILENAME  =  " + fileName2 + ".fits                 " + "MAGZP  =  " + str(magZP) + "                 F_NU  =  " + str(fnu) )
        
        
        vmax = np.nanpercentile(fluxTwoBands, 95)
        vmin = np.nanpercentile(fluxTwoBands, 80)
        plt.imshow(fluxTwoBands, origin="lower", vmax=vmax, vmin=vmin)
        plt.colorbar()
        plt.show()
        
   #      #############   For 3 Filters    ####################
    
   #      COEFS=fits.open("file:///c:/Users/55119/Dropbox/splus-halpha (1)/tables/coeffs.fits")
   
   #      t = Table.read(COEFS)
    
   #      alpha_F660=t['alpha_x'][0]
   #      beta_F660=t['beta_x'][0]
   #      delta_F660=t['delta_x'][0]

   #      alpha_I=t['alpha_x'][1]
   #      beta_I=t['beta_x'][1]
   #      delta_I=t['delta_x'][1]

   #      alpha_R=t['alpha_x'][2]
   #      beta_R=t['beta_x'][2]
   #      delta_R=t['delta_x'][2]

   #      fluxThreeBands = (((fnu_R- fnu_I)-((alpha_R-alpha_I)/(alpha_F660 - alpha_I))*(fnu_F660 - fnu_I))/((beta_F660)*(alpha_I -alpha_R )-(beta_R)))
    
   # # fluxThreeBands = (((fnu_R- fnu_I)-((t['alpha_x'][2])-(t['alpha_x'][1]))/((t['alpha_x'][0] )- (t['alpha_x'][2]))*(fnu_F660 - fnu_I))/(((t['beta_x'][0]))*((t['alpha_x'][1]) - (t['alpha_x'][2]) )-((t['beta_x'][2))))
   #      vmax = np.nanpercentile(fluxThreeBands, 95)
   #      vmin = np.nanpercentile(fluxThreeBands, 80)
   #      plt.imshow(fluxThreeBands, origin="lower", vmax=vmax, vmin=vmin)
   #      plt.colorbar()
   #      plt.show()
        
     #   print(fileName2)
     #   print(type(fnu))
     #   print(fnu.shape)
      #  print(fluxTwoBands)
        
    #    sys.stdout = outputFile = open("resultado_geral.txt", "w+") #altera as saídas para um arquivo de texto no diretório atual
     #   print (hAlphaFlux) #escreve os resultados no diretório atual
     #   outputFile.close() #fecha a conexão do arquivo gerado
        
      #  sys.stdout = outputFile = open("fnu_F660.txt", "w+") #altera as saídas para um arquivo de texto no diretório atual
     #   print (fnu_F660) #escreve os resultados no diretório atual
      #  outputFile.close() #fecha a conexão do arquivo gerado
        
      #  sys.stdout = outputFile = open("fnu_R.txt", "w+") #altera as saídas para um arquivo de texto no diretório atual
      #  print (fnu_R) #escreve os resultados no diretório atual
        #  outputFile.close() #fecha a conexão do arquivo gerado
        
     #   sys.stdout = outputFile = open("fluxTwoBands.txt", "w+") #altera as saídas para um arquivo de texto no diretório atual
       # print (fluxTwoBands) #escreve os resultados no diretório atual
     #   outputFile.close() #fecha a conexão do arquivo gerado
        
       # sys.stdout = mainOutput #altera a saída do Python para a saída principal (console)
    
Usuário: yw5icrr