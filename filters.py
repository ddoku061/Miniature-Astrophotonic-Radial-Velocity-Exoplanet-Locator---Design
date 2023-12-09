"""
#-----------------------------------------------------------------------
# This file is part of NRC PDK for Nazca.
#-----------------------------------------------------------------------
#
#@Authors:  Doga Dokuz & Ross Cheriton 
#2023(c)

**NRC PDK additional geometries for filters**
**Confidential to NRC**
"""
from math import sin, cos, acos, sqrt
from itertools import count
from nazca import interconnects
from nazca.cp import x 
import numpy as np
from numpy.core.function_base import linspace
from scipy.interpolate import interp1d
import nazca as nd
import nazca.geometries as geom
import nazca.pdk_template as pdk
import math
from nazca.interconnects import Interconnect
import os

import nrcfab as fab

PI = 3.14159

gds_path_abs=os.path.join(os.path.abspath(__file__ + "\\..\\..\\"),'gdsBB')


#Define the SWG segment
@nd.bb_util.hashme('swg_segment')
def bragg_triangle_segment(bragg_duty_cycle=1,bragg_period=0.5,bragg_amplitude=0.1,layer=1):
    with nd.Cell(hashme=True) as bragg_segment:
        triangle_width=bragg_period*bragg_duty_cycle
        nd.Polygon(layer=layer, points=[(-triangle_width/2,0), (triangle_width/2,0), (0,bragg_amplitude)]).put(0)
    return bragg_segment

@nd.bb_util.hashme('swg_segment_semitriangle')
def bragg_semitriangle_segment(bragg_duty_cycle=1,bragg_period=0.5,bragg_amplitude=0.1,layer=1):
    with nd.Cell(hashme=True) as bragg_segment:
        triangle_width=bragg_period*bragg_duty_cycle
        nd.Polygon(layer=layer, points=[(0,0), (triangle_width/2,0), (0,bragg_amplitude)]).put(0)
    return bragg_segment

@nd.bb_util.hashme('swg_segment_semitriangle2')
def bragg_semitriangle_segment2(bragg_duty_cycle=1,bragg_period=0.5,bragg_amplitude=0.1,layer=1):
    with nd.Cell(hashme=True) as bragg_segment:
        triangle_width=bragg_period*bragg_duty_cycle
        nd.Polygon(layer=layer, points=[(0,0), (triangle_width/2,0), (0,-bragg_amplitude)]).put(0)
    return bragg_segment

@nd.bb_util.hashme('swg_segment_semicircle')
def bragg_semicircle_segment(bragg_duty_cycle=1,bragg_period=0.5,layer=1):
    with nd.Cell(hashme=True) as bragg_segment:
        semicircle_width=bragg_period*bragg_duty_cycle
        nd.Polygon(layer=layer, points=geom.pie(radius=semicircle_width/2, angle=180, N=20)).put(0)
    return bragg_segment

@nd.bb_util.hashme('swg_segment_trapezoid')
def bragg_trapezoid_segment(bragg_duty_cycle=1,bragg_period=0.5,bragg_amplitude=0.1,angle1=45, angle2=45,layer=1):
    with nd.Cell(hashme=True) as bragg_segment:
        trapezoid_width=bragg_period*bragg_duty_cycle
        nd.Polygon(layer=layer, points=geom.trapezoid(length=trapezoid_width, height=bragg_amplitude, angle1=angle1, angle2=angle2)).put(0)
    return bragg_segment

@nd.bb_util.hashme('swg_segment_square')
def bragg_square_segment(bragg_duty_cycle=0.5,bragg_period=0.5,bragg_amplitude=0.1,layer=1):
    with nd.Cell(hashme=True) as bragg_segment:
        trapezoid_width=bragg_period*bragg_duty_cycle
        nd.Polygon(layer=layer, points=geom.box(length=trapezoid_width, width=bragg_amplitude)).put(0)
    return bragg_segment



def bragg_filter(heater_type='amf', heater_pads=True,layer=1, heater_width=10, edge_coupler_length=500, offset_bend_location=550,bragg_x_start=0,chip_length=4000,bend_radius=50, waveguide_width=0.5,bragg_period=0.5, bragg_duty_cycle=1,bragg_amplitude=0.018, angle1=45, angle2=45, bragg_length=100,bragg_type='triangular',test_waveguide_offset=0, arrow=False,xs='Deep',heater_xs='AMF_TiN_heater'):

    """Return a good 1550 nm splitter optimized for 1550 nm
    Args:
        
    Returns:
        Cell: A Cell object with a simple y splitter
    """

    con = Interconnect(xs=xs)
    with nd.Cell(name='bragg_filter') as bragg_filter_cell:
        taper_length=50
        pad_length=150 
        bragg_width = waveguide_width-bragg_amplitude
        number_of_bragg=int(round(bragg_length/bragg_period))

        #Bragg section

        if bragg_type=='triangular':
            #top half triangle and the array
            bragg_semitriangle_segment(layer=layer,bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude/2).put(0,waveguide_width/2+test_waveguide_offset-(bragg_amplitude)/2)
            bragg_triangle_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude).put(0+(bragg_period*bragg_duty_cycle),(waveguide_width-bragg_amplitude)/2+test_waveguide_offset,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])
            #top half triangle at the end
            bragg_semitriangle_segment2(layer=layer,bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude/2).put(bragg_length,waveguide_width/2+test_waveguide_offset-(bragg_amplitude)/2,180)
            

            #bottom half triangle and the array
            bragg_semitriangle_segment2(layer=layer,bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude/2).put(0,-waveguide_width/2+test_waveguide_offset+bragg_amplitude/2)
            bragg_triangle_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude).put(0+(bragg_period*bragg_duty_cycle),-(waveguide_width-bragg_amplitude)/2+test_waveguide_offset,180,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])
            #bottom half triangle at the end
            bragg_semitriangle_segment(layer=layer,bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude/2).put(bragg_length,-waveguide_width/2+test_waveguide_offset+bragg_amplitude/2,180)
            
            
        
        if bragg_type=='semi_circle':
            bragg_semicircle_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period).put(offset_bend_location+2*bend_radius+bragg_x_start,waveguide_width/2+test_waveguide_offset,90,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])
            bragg_semicircle_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period).put(offset_bend_location+2*bend_radius+bragg_x_start,-waveguide_width/2+test_waveguide_offset,270,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])

        if bragg_type=='trapezoid':
            bragg_trapezoid_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude,angle1=angle1, angle2=angle2).put(offset_bend_location+2*bend_radius+bragg_x_start,waveguide_width/2+test_waveguide_offset,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])
            bragg_trapezoid_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude,angle1=angle1, angle2=angle2).put(offset_bend_location+2*bend_radius+bragg_x_start+bragg_period,-waveguide_width/2+test_waveguide_offset,180,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])

        if bragg_type=='square':
            bragg_square_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude).put(offset_bend_location+2*bend_radius+bragg_x_start,waveguide_width/2+bragg_amplitude/2+test_waveguide_offset,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])
            bragg_square_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude).put(offset_bend_location+2*bend_radius+bragg_x_start,-waveguide_width/2-bragg_amplitude/2+test_waveguide_offset,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])

        if bragg_type=='alternating_square':
            bragg_square_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude).put(offset_bend_location+2*bend_radius+bragg_x_start,waveguide_width/2+bragg_amplitude/2+test_waveguide_offset,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])
            bragg_square_segment(layer=layer, bragg_duty_cycle=bragg_duty_cycle,bragg_period=bragg_period,bragg_amplitude=bragg_amplitude).put(offset_bend_location+2*bend_radius+bragg_x_start,-waveguide_width/2-bragg_amplitude/2+test_waveguide_offset,180,array=[int(number_of_bragg), [bragg_period, 0], 1, [0, 1]])

        if heater_type=='amf':
            straight_heater=con.strt(width=heater_width,length=bragg_length,xs=heater_xs,arrow=arrow).put(0,0)
            if heater_pads==True:
                fab.heaters.amf_pads2(pad_width=80,pad_length=150,x_offset=30).put(straight_heater.pin['a0'])
                fab.heaters.amf_pads1(pad_width=80,pad_length=150,x_offset=30).put(straight_heater.pin['b0'])


        #Straight waveguide for the bragg filter
        con.strt(length=bragg_length,width=bragg_width,xs=xs,arrow=arrow).put(0,0)
    
    return bragg_filter_cell

#Bragg filter for Ca II lines in series
def bragg_filter_triple(separation_between_bragg_filter=800,heater_type='amf', heater_pads=True, heater_width=10,  x_offset = -1, offset_bend_location=550,bragg_x_start=0,chip_length=4740,bend_radius=50, waveguide_width=0.43,bragg_period1=0.25, bragg_period2=0.5, bragg_period3=0.5,bragg_duty_cycle=1,bragg_amplitude=0.018, angle1=45, angle2=45, bragg_length=500,bragg_type='triangular',arrow=False,xs='Deep',heater_xs='AMF_TiN_heater'):
    con = Interconnect(xs=xs)
    with nd.Cell(name='bragg_filter_multline') as bragg_filter_cell_multiline_cell:
        
        taper_length=100
        pad_length=150
        con.taper(length=taper_length,width1 = 0.12,width2=waveguide_width,xs=xs,arrow=arrow).put(x_offset,0)    

        #Straight waveguide before the bragg filter (no sidewall patterns here)
        input_and_output_lengths=(chip_length-3*bragg_length-2*taper_length-2*separation_between_bragg_filter)/2
        print(input_and_output_lengths)
        con.strt(length=input_and_output_lengths,width=waveguide_width,xs=xs,arrow=arrow).put()
        #1st filter
        bragg_filter(bragg_length=bragg_length,bragg_period=bragg_period1, bragg_duty_cycle=1,bragg_amplitude=0.018,waveguide_width=waveguide_width,bragg_type='triangular',xs=xs).put(input_and_output_lengths+taper_length+x_offset,0)
        #Straight waveguide between bragg filters
        con.strt(length=separation_between_bragg_filter-1,width=waveguide_width,xs=xs,arrow=arrow).put(input_and_output_lengths+taper_length+bragg_length+x_offset,0)

        #2nd filter
        bragg_filter(bragg_length=bragg_length,bragg_period=bragg_period2, bragg_duty_cycle=1,bragg_amplitude=0.018).put(input_and_output_lengths+bragg_length+taper_length+separation_between_bragg_filter+x_offset-1,0)

        #Straight waveguide between bragg filters
        con.strt(length=separation_between_bragg_filter-1,width=waveguide_width,xs=xs,arrow=arrow).put(input_and_output_lengths+taper_length+2*bragg_length+separation_between_bragg_filter+x_offset-1,0)
        #3rd filter
        bragg_filter(bragg_length=bragg_length,bragg_period=bragg_period3, bragg_duty_cycle=1,bragg_amplitude=0.018).put(input_and_output_lengths+taper_length+2*bragg_length+2*separation_between_bragg_filter+x_offset-2,0)

        con.strt(length=input_and_output_lengths+4,width=waveguide_width,xs=xs,arrow=arrow).put(input_and_output_lengths+taper_length+3*bragg_length+2*separation_between_bragg_filter+x_offset-2,0)
        con.taper(length=taper_length,width2 = 0.12,width1=waveguide_width,xs=xs,arrow=arrow).put() 
        #below is for debugging purposes
        print("*****************************************************************")
        print("Length of input and output waveguide:" + str(input_and_output_lengths))
        print("*****************************************************************")

    return bragg_filter_cell_multiline_cell


#Bragg filters for Na I lines in series
def bragg_filter_double(separation_between_bragg_filter=1320, layer= 1, heater_type ='amf', heater_pads=True, heater_width = 10, x_offset = -1, offset_bend_location=550,bragg_x_start=0,chip_length=4740,bend_radius=50, waveguide_width=0.5,bragg_period1=0.5, bragg_period2=0.5, bragg_duty_cycle=1,bragg_amplitude=0.018, angle1=45, angle2=45, bragg_length=500,bragg_type='triangular',arrow=False,xs='Deep',heater_xs='AMF_TiN_heater'):
    
    con = Interconnect(xs=xs)
    with nd.Cell(name='bragg_filter_multline') as bragg_filter_cell_multiline_cell:
        taper_length=100
        pad_length=150
        con.taper(length=taper_length,width1 = 0.12,width2=waveguide_width,xs=xs,arrow=arrow).put(x_offset,0)

        #Straight waveguide before the bragg filter (no sidewall patterns here)
        input_and_output_lengths=(chip_length-2*bragg_length-separation_between_bragg_filter-2*taper_length)/2
        con.strt(length=input_and_output_lengths,width=waveguide_width,xs=xs,arrow=arrow).put()

        #1st filter
        bragg_filter(bragg_length=500,bragg_period=bragg_period1, bragg_duty_cycle=1,bragg_amplitude=0.018,waveguide_width=waveguide_width,bragg_type='triangular',xs=xs).put(input_and_output_lengths+taper_length+x_offset,0)        #Straight waveguide between bragg filters
        #Straight waveguide between bragg filters
        con.strt(length=separation_between_bragg_filter-1,width=waveguide_width,xs=xs,arrow=arrow).put(input_and_output_lengths+taper_length+bragg_length+x_offset,0)
        
        #2nd filter
        bragg_filter(bragg_length=500,bragg_period=bragg_period2, bragg_duty_cycle=1,bragg_amplitude=0.018).put(input_and_output_lengths+bragg_length+taper_length+separation_between_bragg_filter+x_offset-1,0)
        #Straight waveguide between bragg filters
        con.strt(length=input_and_output_lengths+3,width=waveguide_width,xs=xs,arrow=arrow).put(input_and_output_lengths+taper_length+2*bragg_length+separation_between_bragg_filter+x_offset-1,0)
        
        con.taper(length=taper_length,width2 = 0.12,width1=waveguide_width,xs=xs,arrow=arrow).put() 

    return bragg_filter_cell_multiline_cell