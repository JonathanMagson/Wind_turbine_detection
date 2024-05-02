#This version corresponds to the executive code used in the IPOL article published in july/august 2021, for S2 images, monodate.

import sys
import numpy as np
import scipy
import scipy.special

################################## Shadow calculation #######################################
################################## Shadow calculation #######################################
################################## Shadow calculation #######################################

#==============================================================================
def compute_shadow(azimuth, altitude):
# This function calculates the direction and unitary length of the shadow.
# Its entries are two angles :
# - azimuth : sun's azimuth (degrees)
# - altitude : sun's altitude (degrees)

    # convert degrees into radians
    altitude, azimuth = altitude / 180.0 * np.pi, azimuth / 180.0 * np.pi
           
    x = - np.cos(altitude) * np.sin(azimuth)
    y = - np.cos(altitude) * np.cos(azimuth)
    z = - np.sin(altitude)

    return [x, y, z]
#==============================================================================
#==============================================================================
def main_shadow_S2(azisun, altisun,\
               h = 80,resolution = 10,sample_rate = 10, dist_shadow = 15):
# This function calculates the coordinate shifts between central pixel and shadow samples and their neighbors.
# Its entries are
# - azisun : sun's azimuth (degrees)
# - altisun : sun's altitude (degrees)
# - h : tested heights
# - resolution : satellite's resolution (meters)
# - sample_rate : sampling rate (meters)
# - dist_shadow : sample-neighbor distance (meters)
# Its output is an array structured in 6 columns containing the coordinate information as follows : 
#	[row_neighbor1,row_sample,row_neighbor2,column_neighbor1,column_sample,column_neighbor2].
#From each row of this array, we can derive the coordinates of a shadow sample and its two associated neighbors.


    #compute the shadow coordinates    
    [x, y, z] = compute_shadow(azisun, altisun)

    # compute the coordinates in meters in the earth coordinate system such that x grows from West to East and y grows from South to North. These coordinates will correspond to the end of the shadow segment.
    # note that x varies in the same way than the columns of an array (West/left -> East/right), and y varies in the opposite way than the on of the rows (South->North versus top->bottom)
    x, y = h/np.abs(z) * x, h/np.abs(z) * y
        
    # convert in pixels
    x,y = x/resolution, y/resolution
        
    # base of the segment
    cx, cy = 0, 0   

    # end of the shadow segment
    sx, sy = cx + x, cy - y 
    # the minus in -y comes from y growing from South to North in the earth coordinate system, but the rows in an array grow from top to bottom

    # compute unitary orientation vector
    norme_xy = np.sqrt(x*x + y*y)
    norm_x, norm_y = x/norme_xy, -y/norme_xy

    # compute minimum gap between two samples
    samplrate_x, samplrate_y = norm_x*sample_rate/resolution, norm_y*sample_rate/resolution

    # sample the shadow coordinates
    eps = sys.float_info.epsilon    
            #columns : x
    if cx != sx: #tackle the case where the shadow is strictly horizontal or vertical
        cc = [i for i in np.arange(cx,sx,samplrate_x)]
        cc = np.array(cc)
            #rows : y
    if cy != sy:
        rr = [i for i in np.arange(cy,sy,samplrate_y)]
        rr = np.array(rr)
        #tackle the case where the shadow is strictly horizontal or vertical
    if cx == sx:
        lencc = rr.shape[0]
        cc = np.ones((lencc,1))*cx
    if cy == sy:
        lenrr = cc.shape[0]
        rr = np.ones((lenrr,1))*cy
        
    # SHADOW COORDINATES COMPUTED #
    
    # compute the neighbors : 
    fact_xy = norme_xy
    fact_loin = dist_shadow/resolution

    # compute the gap between sample and neighbor. Orthogonal vector of [x,-y] is [y,x].
    petit_decal_x = fact_loin*y/fact_xy
    petit_decal_y = fact_loin*x/fact_xy

    rrN1,ccN1 = rr+petit_decal_y,cc+petit_decal_x
    rrN2,ccN2 = rr-petit_decal_y,cc-petit_decal_x
    # NEIGHBORS COORDINATES COMPUTED #

    # regroup everything in an array of the right shape
    npix = len(rr)

    rrr = np.concatenate((np.transpose(rrN1), np.transpose(rr), np.transpose(rrN2) ))
    rrr = np.transpose(np.reshape(rrr,(3,npix)))
    ccc = np.concatenate((ccN1, cc, ccN2))
    ccc = np.transpose(np.reshape(ccc,(3,npix)))

    tab = np.concatenate((rrr,ccc),axis=1)
    # FINAL ARRAY COMPUTED #

    return(tab)
#==============================================================================


######################################### Hub calculation ######################################
######################################### Hub calculation ######################################
######################################### Hub calculation ######################################

#==============================================================================
def compute_hub(azimuth, zenith):
# This function calculates the direction and unitary length of the shift between bottom and top of tower.
# Its entries are :
# - azimuth : satellite's azimuth (degrees)
# - zenith : satellite's zenith (degrees)

    # conversion degrees into radians    
    zenith, azimuth = zenith / 180.0 * np.pi, azimuth / 180.0 * np.pi    
    
    x = - np.sin(azimuth)
    y = - np.cos(azimuth)
    z = - np.tan(zenith)

    return [x, y, z]
#==============================================================================
#==============================================================================
def get_mesh_triangle(centrex,centrey,nangles = 6,resolution = 10,sample_rate = 10,r_oi = 3):
# This function calculates a triangular mesh of limited size around center coordinates. This mesh will be used to sample the hub pixels.
# Its entries are
# - centrex,centrey : center coordinates
# - nangles : number of sides of the polygon sampling shape
# - resolution : satellite's resolution (m)
# - sample_rate : sampling rate
# - r_oi : size of the region of interest (m)
# Its output is an array containing 2 coordinate columns. Each line gives a sampling set of coordinates.

    # create the stocking lists
    coordx = [centrex]
    coordy = [centrey]

    # compute the farthest coordinates from the center but still in the region of interest
    max_echant = resolution * np.ceil((r_oi+resolution/2)/sample_rate)
    vec_echant = np.arange(0,max_echant*1.0001,sample_rate)/resolution
    
    # compute the coordinates of the mesh, rotating nangles times
    for echanti in vec_echant[1:]:    
        for i in range(nangles):    
                theta = 2*np.pi*i/nangles
                coordpix = [echanti*np.cos(theta)+centrex,echanti*np.sin(theta)+centrey]
                coordx.append((coordpix[0]))
                coordy.append((coordpix[1]))
    
    # wrap up everything            
    coord = np.array((coordx,coordy))
    coord = np.transpose(coord)
    tab = coord
    
    return(tab)
#==============================================================================
#==============================================================================
def main_hub_S2(azisat, zensat,\
               h = 80,resolution = 10,sample_rate = 10,nangles = 6):
# This function calculates the coordinate shifts between central pixel and hub samples.
# Its entries are
# - azisat : satellite's azimuth (degrees)
# - zensat : satellite's zenith (degrees)
# - h : tested heights (meters)
# - resolution : satellite's resolution (meters)
# - sample_rate : sampling rate (meters)
# - nangles : sampling pattern, here = 6
# Its output is the saving of an array containing 2 coordinate columns. The structure is [coordx,coordy].
#The conversion into rows and columns will be [col,row], see main_shadow_S2.

    sat = 'S2'

    # compute the hub coordinates    
    [x, y, z] = compute_hub(azisat, zensat)

    # compute the coordinates in meters in the earth coordinate system such that x grows from West to East and y grows from South to North.
    #These coordinates will correspond to the end of the tower, where the hub stands.
    # note that x varies in the same way than the columns of an array (West/left -> East/right),
    #and y varies in the opposite way than the rows (South->North versus top->bottom)
    x, y = h*np.abs(z) * x, h*np.abs(z) * y
    
    # convert in pixels
    x,y = x/resolution, y/resolution
    
    # base of the tower
    cx, cy = 0, 0

    # top of of the tower
    sx, sy = cx + x, cy - y
    centrex,centrey = sx,sy

    # get the sampled shift coordinates
    coordhub = get_mesh_triangle(centrex,centrey,nangles,resolution,sample_rate,r_oi = 3) #[cc,rr] array

    return(coordhub)
#==============================================================================


#################### bilinear_interpolation tools ##########################
#################### bilinear_interpolation tools ##########################
#################### bilinear_interpolation tools ##########################
# We design a set of tool functions to be used in different cases : array, single point, whole image, ...

#==============================================================================
def bilinear_interpolation(coordinterp, im, opt_glob = 'full'):
# This function interpolates an image in a single set of coordinates.
# Its entries are
# - coordinterp : metadata repository
# - im : the global image
# - opt_glob : option = 'full'. If 'full', interpolates the whole image in the given subpixel shift. 
# Its output is :
# - the whole image shifted by a subpixel shift, if option 'full'.

    if opt_glob == 'full':
        x,y = coordinterp[0],coordinterp[1]
        x1,y1 = int(np.floor(x)),int(np.floor(y))
        x2,y2 = x1+1,y1+1
    
        q11 = im[:-1,:-1]
        q12 = im[:-1,1:]
        q21 = im[1:,:-1]
        q22 = im[1:,1:]

        res1 = (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1))
    
    return res1
#==============================================================================
#==================================================================================   
def interp_plages(im,plagex,plagey,inter):
# This function interpolates an image on (plagex,plagey).
# Its entries are :
# - im : interpolated image
# - plagex,plagey : two growing lists of rate 1
# - inter : interpolation type, can be closest neighbor ('PPV') or bilinear ('bilin')
# Its output is the image interpolated on every pair of coordinates (x,y) with x in plagex, y in plagey.

    if inter == 'PPV':
        xtopleft,ytopleft = round(plagex[0]),round(plagey[0])
        xbotright,tbotright = round(plagex[-1]),round(plagey[-1])
        res = im[xtopleft:xbotright+1,ytopleft:tbotright+1]
    if inter == 'bilin':
        xtopleft,ytopleft = plagex[0],plagey[0]
        xbotright,tbotright = plagex[-1],plagey[-1]
        x1,y1 = int(np.floor(xtopleft)),int(np.floor(ytopleft))
        xend,yend = int(np.floor(xbotright)),int(np.floor(tbotright))
        if abs(xtopleft-round(xtopleft))<10**-9: #tackling calculus approximation
            x1 = int(round(xtopleft))
        if abs(ytopleft-round(ytopleft))<10**-9:
            y1 = int(round(ytopleft))
        if abs(xbotright-round(xbotright))<10**-9:
            xend = int(round(xbotright))
        if abs(tbotright-round(tbotright))<10**-9:
            yend = int(round(tbotright))
        x2,y2 = x1+1,y1+1
        x,y = xtopleft,ytopleft
    
        q11 = im[x1:xend+1,y1:yend+1]
        q12 = im[x1:xend+1,y1+1:yend+2]
        q21 = im[x1+1:xend+2,y1:yend+1]
        q22 = im[x1+1:xend+2,y1+1:yend+2]

        res = (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1))
        
    return(res)
#==================================================================================   


################### Empirical probabilities computation ##################
################### Empirical probabilities computation ##################
################### Empirical probabilities computation ##################

#==============================================================================
def calc_p_emp_shadow(im,xy,xy1,diff_s=0):
# This function estimates the empirical probability for a pixel to pass the shadow test, ie for a random pixel value to be smaller than two neighbors-diff_s, by computing all the differences across the image and taking the proportion passing the test.
# Its entries are :
# - im : image
# - xy : pixel coordinates
# - xy1 : neighbor coordinates
# - diff_s : shadow threshold tightening the test
# Its output is the number of pixels in the image which pass the test divided by the total number of pixels. It approximates a local probability.

    nrow, ncol = im.shape[0],im.shape[1]
    
    # subpixel part of the shift
    decal_x = abs((xy[0]-xy1[0]))%1 
    decal_y = abs((xy[1]-xy1[1]))%1
    
    # pixel part of the shift
    diff_x = int(abs((xy[0]-xy1[0]) - decal_x))
    diff_y = int(abs((xy[1]-xy1[1]) - decal_y))
    
    coordinterp_1 = [decal_x,decal_y]
    coordinterp_2 = [-decal_x,-decal_y]
    opt_glob = 'full'
    
    im_1 = bilinear_interpolation(coordinterp_1, im, opt_glob)
    im_2 = bilinear_interpolation(coordinterp_2, im, opt_glob)
    
    # avoid exiting the image bounds, get the right shape to compare im, im_1, im_2
    if diff_x == 0:
        im_minus = im_1[:,:-(2*diff_y)]
        im_plus = im_2[:,2*diff_y:]
        im_crop = im[:-1,diff_y:-diff_y-1] 
    elif diff_y == 0:
        im_minus = im_1[:-(2*diff_x),:]
        im_plus = im_2[2*diff_x:,:]
        im_crop = im[diff_x:-diff_x-1,:-1]
    else:
        im_minus = im_1[:-(2*diff_x),:-(2*diff_y)]
        im_plus = im_2[2*diff_x:,2*diff_y:]
        im_crop = im[diff_x:-diff_x-1,diff_y:-diff_y-1]

    # get the pixels passing the test
    comp_plus = ((im_crop-im_plus)<-diff_s)*((im_crop-im_minus)<-diff_s)
    
    # count and divide    
    p_emp = sum(sum(comp_plus))/((nrow-2*decal_x)*(ncol-2*decal_y)) 
    
    return p_emp
#==============================================================================
#==================================================================================   
def calc_p_emp_hub(im,resolution,dh,nch = 6,diff_h = 0):
# This function calculates the proportion of pixels in the image which verifies the hub condition, ie pix > pix_i + diff_h
#foreach i in 1..nch (entier), with pix_i located around pix in a regular polygon shape of nch sides, at distance decal.
# Its entries are :
# - im : image
# - resolution : image resolution (meters)
# - dh : distance (meters) between pix and pix_i
# - nch : number of neighbors pix_i
# - diff_h : hub threshold tightening the test 
# Its outuput is the proportion of pixels in the image which pass the hub test. It approximates the corresponding probability.

    # convert dh in pixels
    decal = dh/resolution

    # out-of-bonds safety gap
    d = int(np.ceil(decal))
    
    nx,ny = im.shape[0],im.shape[1]
    
    # set counting variables
    nwin = 0
    ntests = 0
    
    # avoid exiting the image bounds
    imcrop = im[d+1:nx-d-1,d+1:ny-d-1]
    nxc,nyc = imcrop.shape[0],imcrop.shape[1]
    tab_bool = np.ones((nxc,nyc))
    
    # count
    for ii in range(nch):
        theta = np.pi*ii*2/nch
        decalx,decaly = int(decal*np.cos(theta)),int(decal*np.sin(theta))
        imcomp = im[d+1+decalx:nx-d-1+decalx,d+1+decaly:ny-d-1+decaly]
        
        tab_bool = tab_bool*(imcrop > imcomp + diff_h)
        # as soon as a comparison fails, the pixel is deactivated
    
    # wrap up        
    ntests = nxc*nyc
    nwin = np.sum(tab_bool)

    res = nwin/ntests
    
    return(res)
#================================================================================== 


############################## NFA #################################
############################## NFA #################################
############################## NFA #################################

#==============================================================================
def main_NFA(sat,im,azisun, altisun,azisat, zensat,\
             nch = 6,H=[80],r=10,sr=10,ds=15,inter='bilin',dh=30,diff_s=25,diff_h=50):
# This function computes the NFA map of a S2 image.
# Its entries are :
# - sat : satellite (here only 'S2')
# - im : a monocanal sentinel 2 array
# - azisat : satellite's azimuth (degrees)
# - zensat : satellite's zenith (degrees)
# - azisun : sun's azimuth (degrees)
# - altisun : sun's altitude (degrees)
# - nch : number of hub neighbors of a hub sample
# - H : array of tested heights
# - r : satellite resolution
# - sr : sampling rate
# - ds : shadow neighbor distance
# - inter : interpolation type, here only 'bilin' for bilinear
# - dh : hub neighbor distance
# - diff_s : shadow threshold
# - diff_h : hub threshold
# Its output is the NFA map.
    
    # avoid shape compatibility errors
    if len(im.shape) >= 3: 
        im = im[:,:,0]
    (nrow,ncol) = im.shape

    for h in H:

        # get coord csv shadow
        tab_s = main_shadow_S2(azisun, altisun,\
                   h = h,resolution = r,sample_rate = sr, dist_shadow = ds)
        rrN1, rr, rrN2, ccN1, cc, ccN2 = tab_s[:,0],tab_s[:,1],tab_s[:,2],tab_s[:,3],tab_s[:,4],tab_s[:,5]

        # get coord csv hub
        tab_h = main_hub_S2(azisat, zensat,\
                   h = h,resolution = r,sample_rate = sr,nangles = 6)
        hubcc, hubrr = tab_h[:,0],tab_h[:,1]
        
        # get number of sampled pixels
        nk_shadow = len(rr)
        nk_hub = len(hubrr)
        k_max = nk_shadow+nk_hub
    # SHADOW AND HUB SAMPLES COMPUTED #
        
    # compute probabilities        
        
        # compute empirical probabilities
        xy = [rr[0], cc[0]]
        xy1 = [rrN1[0], ccN1[0]]
        mini_ou_maxi = 'mini'
        p_emp_shadow = calc_p_emp_shadow(im,xy,xy1,diff_s)
        p_emp_hub = calc_p_emp_hub(im,r,dh,nch = 6,diff_h = 0)

        # compute weighted average probability                
        pmoy = (nk_shadow*p_emp_shadow + nk_hub*p_emp_hub)/k_max
    # PROBABILITIES COMPUTED #

    # compute K_shadow, the likelihood score for the shadow
        
        nfa = np.ones((im.shape)) # ones will be left too-close-to-the-border pixels

        # move away from the image border, to avoid out of bounds
        distsafety_shadow_x = max(abs(np.concatenate((rr,rrN1,rrN2))))
        distsafety_shadow_y = max(abs(np.concatenate((cc,ccN1,ccN2))))
        distsafety_hub_x = dh/r + max(abs(hubrr))
        distsafety_hub_y = dh/r + max(abs(hubcc))
        ecx = int(np.ceil(max([distsafety_shadow_x,distsafety_hub_x])))
        ecy = int(np.ceil(max([distsafety_shadow_y,distsafety_hub_y])))

        # compute total number of tests
        Ntests = (nrow-2*ecx-1)*(ncol-2*ecy-1)

        # crop the image to keep only the area sufficiently away from the borders
        plagex = np.array(range(ecx+1,nrow-ecx-1))
        plagey = np.array(range(ecy+1,ncol-ecy-1))
        imcrop = im[ecx+1:nrow-ecx-1,ecy+1:ncol-ecy-1]

        # initialize likelihood score array : K_shadow
        K_shadow = np.zeros((plagex.shape[0],plagey.shape[0]))  

        # rather than looping on every image pixels, we loop on every sample shifts.
        #We compare in one go every pixels of the cropped image.
        for i in range(len(rr)): 
            coordsx,coordsy = rr[i],cc[i]
            compx1,compy1, compx2,compy2 = rrN1[i],ccN1[i], rrN2[i],ccN2[i]
            newplagex = plagex + coordsx
            newplagey = plagey + coordsy
            plagecompx1,plagecompy1 = plagex + compx1,plagey + compy1
            plagecompx2,plagecompy2 = plagex + compx2,plagey + compy2

            im_mid = interp_plages(im,newplagex,newplagey,inter)
            im_g = interp_plages(im,plagecompx1,plagecompy1,inter)
            im_d = interp_plages(im,plagecompx2,plagecompy2,inter)

            # increment K_shadow if the pixel passes the test
            add_K = (im_mid < im_d-diff_s)*(im_mid < im_g-diff_s)
            K_shadow = K_shadow + add_K
    # K_shadow COMPUTED #

    # compute K_HUB, the likelihood score for the hub
        # convert dh (m) in pixels
        scaleh = dh/r 
        
        # initialize K_HUB
        K_HUB = np.zeros((plagex.shape[0],plagey.shape[0]))
        
        # same as for the shadow : we loop on every sample shifts, to compare in one go every pixels of the cropped image
        for i in range(len(hubrr)):
            coordsx,coordsy = hubrr[i],hubcc[i]
            newplagex = plagex + coordsx
            newplagey = plagey + coordsy

            # K_INTERM is a stock array, in which pixels are ones at the beginning and become zeros when they fail a comparison test
            K_INTERM = np.zeros((plagex.shape[0],plagey.shape[0]))
            K_INTERM = K_INTERM+1

            # loop on the neighbors
            for j in range(nch):
                theta = j*2*np.pi/nch
                compx1,compy1 = scaleh*np.cos(theta),scaleh*np.sin(theta)
                plagecompx1,plagecompy1 = plagex + compx1,plagey + compy1
                im_mid = interp_plages(im,newplagex,newplagey,inter)
                im_g = interp_plages(im,plagecompx1,plagecompy1,inter)
                COMP = im_mid > (im_g + diff_h)
                
                # avoid shape compatibility problems
                if len(COMP.shape) > 2: #shape de la forme (nx,ny,1)
                    COMP = COMP[:,:,0] 
                                        
                # as soon as a comparison is failed, the pixel is switched off for this sample
                K_INTERM = K_INTERM*COMP

            # increment the likelihood score K_HUB
            K_HUB = K_HUB + K_INTERM
    # K_HUB COMPUTED #

        # get the total likelihood score
        K0 = K_HUB+K_shadow

    # convert likelihood into p-value probability
        
        # compute cumulative probability
        cumpl = 0
        Hval = np.ones(K0.shape)
        for l in range(k_max+1):
            # to avoid float approximation, we must begin by the end
            m = k_max-l
            
            # compute p-value
            pl = scipy.special.binom(k_max, l)*(pmoy**m)*((1-pmoy)**(k_max-m))
            cumpl = cumpl + pl
            
            # stock result
            Hval[K0 == m] = cumpl
    # LIKELIHOOD CONVERTED INTO P-VALUE #

        # NFA computation
        NFA = Ntests*Hval*len(H)

        # transform cropped NFA map into NFA map of the same size of the original image
        valnu = max(np.amax(NFA),Ntests) #value at which set the too-close-to-the-borders pixels
        big_NFA = np.zeros((im.shape[0],im.shape[1])) + valnu
        big_NFA[ecx+1:nrow-ecx-1,ecy+1:ncol-ecy-1] = NFA

        # tackle different height hypothesis
        if h == H[0]:
            stockNFA_H = big_NFA
        else:
            stockNFA_H = np.minimum(stockNFA_H,big_NFA) #keep the best NFA

    # convert NFA into -log(NFA) for visualisation
    NFAfinal = -np.log(stockNFA_H)/np.log(10)
  
    return NFAfinal
#==============================================================================


