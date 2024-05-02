#!/usr/bin/env python3

"""
This function is the runnable one, linking everything, producing the final files for the demo IPOL.
"""

import numpy as np
import argparse
import tifffile
import imagecodecs
import imageio

import utils_eol_compute as uec
import utils_eol_meta_nogdal as uem


def save_NFAs(f, meta, t_NFA=0, t_shadow=25, t_hub=50):
#f : filename (.tif)
#meta : metadata (.txt)
#t_NFA : -logNFA threshold
#t_shadow,t_hub : shadow and hub thresholds
    
    # load the .tif
    im = uem.get_image(f)

    # get sat and sun angles
    sat_angles,sun_angles = uem.get_sat(meta),uem.get_sun(meta)
    azisat,zensat = sat_angles[0],sat_angles[1]
    azisun,altisun = sun_angles[0],sun_angles[1]
    
    # compute -logNFA map
    nfa = uec.main_NFA('S2', im, azisun, altisun, azisat, zensat, diff_s=t_shadow, diff_h=t_hub)

    # save
    tifffile.imsave('nfa.tif',nfa)

    # convert into png & save
    nfapng = (nfa-np.min(nfa))/(np.max(nfa)-np.min(nfa))*256
    imageio.imsave('nfa.png',nfapng)

    # threshold & save
    nfa_detec = 256*(nfa > t_NFA)
    imageio.imsave('nfa_detec.png',nfa_detec)

    # convert input into png & save
    inputpng = (im-np.min(im))/(np.max(im)-np.min(im))*256
    imageio.imsave('input_0.png',inputpng)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('Wind turbine monodate detection using a contrario method (OH)'))
    parser.add_argument('--filename', type=str,
                        help=('Sentinel-2 image filename (.tif)'))
    parser.add_argument('--metaname', type=str,
                        help=('metadata corresponding to the Sentinel-2 file (.txt)'))
    parser.add_argument('--t_NFA', type=int, default=0,
                        help=('-logNFA threshold'))
    parser.add_argument('--t_shadow', type=int, default=25,
                        help=('shadow threshold'))
    parser.add_argument('--t_hub', type=int, default=50,
                        help='hub threshold')
    args = parser.parse_args()

    save_NFAs(f=args.filename, meta=args.metaname, t_NFA=args.t_NFA, t_shadow=args.t_shadow, t_hub=args.t_hub)
   
