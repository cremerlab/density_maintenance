import numpy as np 

def volume(length, width):
    """
    Computes the total volume of a spherocylinder given a width and total length. 
    """
    return np.pi * (width**2 / 12) * (3 * length - width)

def surface_area(length, width):
    """
    Computes the total surface area of a spherocylinder given a width and length.
    """
    return np.pi * length * width

def envelope_volume(length, width, delta):
    """
    Computes periplasmic volume of a spherocylinder given a total length, width,
    and thickness (delta) of the periplasm. 
    """
    return volume(length, width) - volume(length - 2 * delta, width - 2 * delta)

def surface_area_volume(length, width, delta=0):
    """
    Computes the surface area to volume ratio for a spherocylinder given a 
    length and a width. If a thickness (delta) is supplied, the surface area to
    envelope volume ratio is returned.
    """
    sa = surface_area(length, width)
    if delta == 0:
        vol = volume(length, width)
    else:
        vol = envelope_volume(length, width, delta)
    return sa / vol
    
