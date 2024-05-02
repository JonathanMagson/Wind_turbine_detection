#This version corresponds to the code used in the IPOL article published in july/august 2021, for S2 images.
#Its role is to get sat and sun angles, download image, etc.

from os.path import join, basename
import utm
import datetime
from pysolar.solar import *
import numpy as np
import tifffile
import imagecodecs

################################## Extract sun's angles #######################################
################################## Extract sun's angles #######################################
################################## Extract sun's angles #######################################

#==============================================================================
def get_sun(meta, rep = ""):
# This function gets the sun's angles : azimuth and altitude.
# Its entries are meta, the filename of a Sentinel 2 metadata file ; rep, the repository where you can find meta.
# Its output is a np.array sun_angles=[azisun,altisun] (degrees)

    if rep is None:
        ffile = meta
    else:
        ffile = join(rep, basename(meta))      
        
    with open(meta) as texte:
        lines = texte.readlines()
    rows = lines

    # search for UTM position
    pos = [l for l in rows if "Center" in l]
    pos = pos[0]
    posx, posy = pos.split(')')[0].split('(')[1].split(',')
    posx, posy = eval(posx), eval(posy)

    # conversion into latitude longitude
    info = [l for l in rows if "PROJCRS" in l]
    zone = info[0].split('/')[1].split(' ')[3].strip('",')
    zone = zone[:-3]
    zone_letter = zone[-1]
    zone_number = eval(zone[0:len(zone)-1])
    latitude, longitude = utm.to_latlon(posx, posy, zone_number, zone_letter)

    # search for acquisition time
    granule=[l for l in rows if "granule_date" in l]
    granule=granule[0].split('=')[1]
    aa, mois, jj = granule.split(' ')[0].split('-')
    hh, mm, ss = granule.split(' ')[1].split(':')
    ss = ss[:-1]

    madate = datetime.datetime(int(aa), int(mois), int(jj),
                               hour=int(hh), minute=int(mm), second=int(ss),
                               tzinfo=datetime.timezone.utc)

    # get altitude and azimuth
    altitude = get_altitude(latitude, longitude, madate)
    azimuth = get_azimuth(latitude, longitude, madate)

    return [azimuth,altitude]
#==============================================================================


################################## Extract sat's angles #######################################
################################## Extract sat's angles #######################################
################################## Extract sat's angles #######################################

#==============================================================================
def get_sat(meta, rep = ""):
# This function gets the sat's angles : azimuth and zenith.
# Its entries are meta, the filename of a Sentinel 2 metadata file ; rep, the repository where you can find meta.
# Its output is a np.array sun_angles=[azisat,zensat] (degrees)

    if rep is None:
        ffile = meta
    else:
        ffile = join(rep, basename(meta))      

    with open(meta) as texte:
        lines = texte.readlines()
    rows = lines

    # search for azimuth sat
    pos = [l for l in rows if "satellite_azimuth" in l]
    azimuth = eval(pos[0][53:61])

    # search for zenith sat
    pos = [l for l in rows if "satellite_zenith" in l]
    zenith = eval(pos[0][52:60])

    return [azimuth,zenith]
#==============================================================================


################################## Load sat image #######################################
################################## Load sat image #######################################
################################## Load sat image #######################################

#==============================================================================
def get_image(f, rep = ""):
# This function loads the satellite image from the .tif.
# Its entries are f, the filename of a Sentinel 2 data file ; rep, the repository where you can find f.
# Its output is the np.array of the satellite image

    if rep is None:
        ffile = f
    else:
        ffile = join(rep, basename(f))
    im = tifffile.imread(ffile)
    if len(im.shape) > 2:
        im = im[:,:,0]
    return(im)
#==============================================================================

