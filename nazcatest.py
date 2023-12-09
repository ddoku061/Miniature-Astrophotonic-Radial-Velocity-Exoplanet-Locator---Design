import nazca as nd
import sys
import os
import numpy as np
import nazca.pdk_template as pdk
from numpy.core.fromnumeric import size
import nrcfab as fab
#import nazca.demofab as demo
"""
Honours Project Interim Report
#@Authors:  Doga Dokuz
#2023(c)

**Confidential to NRC**
"""
sys.path.insert(0, 'C:\\Users\\a9\\Desktop\\honours')


import math as math
import nazca.geometries as geom
import string
from nazca.interconnects import Interconnect
import time
start = time.time()

PI = 3.14159
file_name = os.path.basename(__file__)
file_name = file_name[:-3]

from nazca.interconnects import Interconnect

con = Interconnect(width=0.5,radius=10) #make all default interconnects 0.5 um wide and have a bend radius of 10 um
'''everything is in microns'''

#defaults, these can be modified depending on the parameters we want to change
layer = 1 #this is for Applied Nanotools processing
bragg_x_start=0 #where the bragg filters should start
chip_length=4000 #how big the chip is 
bend_radius=50 #this is in case we wanted to make the filters curved to save space 
waveguide_width =0.43 #the width of the waveguide
bragg_period = 0.25 
bragg_duty_cycle = 1
bragg_amplitude = 0.5 #in the simulation it was used 0.002um however, we change them to 0.5 um to be able to see the changes in the design
bragg_length=500 #the length of the bragg filter
bragg_type='triangular' #the shapes of the bragg gratings
test_waveguide_offset=0 
x_offset = -1 #starting offset to be able to connect the fiber cables


#creating a definition for Ca lines to separate the filters
def Ca(layer=1, x_offset = x_offset, bragg_x_start=bragg_x_start, chip_length=4740,bend_radius=50,waveguide_width=0.43,bragg_period1=0.27,bragg_period2=0.25,bragg_period3=0.3,bragg_duty_cycle=0.5,bragg_amplitude=0.2,bragg_length=500,bragg_type='triangular',test_waveguide_offset=test_waveguide_offset,distance = 0):
    fab.filters.bragg_filter_triple(bragg_period1=bragg_period1,bragg_period2=bragg_period2,bragg_period3=bragg_period3,x_offset = x_offset,bragg_x_start=bragg_x_start,chip_length=chip_length,bend_radius=bend_radius, waveguide_width=waveguide_width, bragg_duty_cycle=bragg_duty_cycle,bragg_amplitude=bragg_amplitude, bragg_length=bragg_length,bragg_type='triangular').put(0,distance) #wavelength for 849 nm
   
#creating a definition for Na lines to separate the filters
def Na(layer=1, x_offset = x_offset,bragg_x_start=bragg_x_start, chip_length=4740,bend_radius=50,waveguide_width=0.5,bragg_period=1,bragg_duty_cycle=1,bragg_amplitude=0.5,bragg_length=500,bragg_type='triangular',test_waveguide_offset=test_waveguide_offset,distance = 0):
    fab.filters.bragg_filter_double(layer=1,x_offset = x_offset,bragg_x_start=bragg_x_start,chip_length=4740,bend_radius=50, waveguide_width=0.5,bragg_period1=0.2, bragg_period2= 0.25,bragg_duty_cycle=1,bragg_amplitude=0.5, bragg_length=500,bragg_type='triangular').put(0,distance)
    

#creating CHIP C3
Ca(x_offset = -1, distance = -500)
Ca(x_offset = -1, distance = -800)
Ca(x_offset = -1, distance = -1100)
Ca(x_offset = -1, distance = -1400)
Ca(x_offset = -1, distance = -1700)
Ca(x_offset = -1, distance = -2000)
Ca(x_offset = -1, distance = -2300)
Ca(x_offset = -1, distance = -2600)
Ca(x_offset = -1, distance = -2900)
Ca(x_offset = -1, distance = -3200)
Ca(x_offset = -1, distance = -3500)
Ca(x_offset = -1, distance = -3800)
Ca(x_offset = -1, distance = -4100)
Ca(x_offset = -1, distance = -4400)

#Ccreating CHIP C4
Ca(x_offset = 5000-1, distance = -500)
Ca(x_offset = 5000-1, distance = -800)
Ca(x_offset = 5000-1, distance = -1100)
Ca(x_offset = 5000-1, distance = -1400)
Ca(x_offset = 5000-1, distance = -1700)
Ca(x_offset = 5000-1, distance = -2000)
Ca(x_offset = 5000-1, distance = -2300)
Ca(x_offset = 5000-1, distance = -2600)
Ca(x_offset = 5000-1, distance = -2900)
Ca(x_offset = 5000-1, distance = -3200)
Ca(x_offset = 5000-1, distance = -3500)
Ca(x_offset = 5000-1, distance = -3800)
Ca(x_offset = 5000-1, distance = -4100)
Ca(x_offset = 5000-1, distance = -4400)

#Ccreating CHIP D3
Na(x_offset = -1, distance = -5500)
Na(x_offset = -1, distance = -5800)
Na(x_offset = -1, distance = -6100)
Na(x_offset = -1, distance = -6400)
Na(x_offset = -1, distance = -6700)
Na(x_offset = -1, distance = -7000)
Na(x_offset = -1, distance = -7300)
Na(x_offset = -1, distance = -7600)
Na(x_offset = -1, distance = -7900)
Na(x_offset = -1, distance = -8200)
Na(x_offset = -1, distance = -8500)
Na(x_offset = -1, distance = -8800)
Na(x_offset = -1, distance = -9100)
Na(x_offset = -1, distance = -9400)


# creating CHIP D4
Na(x_offset = 5000-1,distance = -5500)
Na(x_offset = 5000-1, distance = -5800)
Na(x_offset = 5000-1, distance = -6100)
Na(x_offset = 5000-1, distance = -6400)
Na(x_offset = 5000-1, distance = -6700)
Na(x_offset = 5000-1, distance = -7000)
Na(x_offset = 5000-1, distance = -7300)
Na(x_offset = 5000-1, distance = -7600)
Na(x_offset = 5000-1, distance = -7900)
Na(x_offset = 5000-1, distance = -8200)
Na(x_offset = 5000-1, distance = -8500)
Na(x_offset = 5000-1, distance = -8800)
Na(x_offset = 5000-1, distance = -9100)
Na(x_offset = 5000-1, distance = -9400)

#creating an area for the chip 20000 by 20000
chip_1 = fab.DesignArea(length=20000, height=20000)
chip_1.applied_nanotools_edge_multichip(chip_name="ASTRO 13 SOI",write_text=False,write_logos=False, names="D. Dokuz",tiny_chip=True, name_list='vert',second_name=True, logo_position='top right',number_of_chips_x=4,number_of_chips_y=4,total_chip_height=20000, total_chip_width=20000).put(0,0)

#exporting to the gds file
nd.export_gds(filename='mychip')
#just to see how much time it took to run the code
end = time.time()
#print('=================================================================================')
print("Total time to generate mask: " + str(end - start)+ " seconds")
#print('=================================================================================')  