import numpy as np
import scipy.ndimage 
import skimage.io 
import skimage.measure
import skimage.morphology
import skimage.segmentation
from .analytical import *
import tqdm
import pandas as pd 

def contour_seg(image, level=0.3, selem='default', perim_bounds=(5, 1E3),
                ip_dist=0.160, ecc_bounds=(0.7, 1), area_bounds=(1, 50),
                return_conts=False, min_int=0.2):
    """
    Identifies contours around dark objects in a phase contrast image.
    Parameters
    ----------
    image: 2d-array
        Phase contrast image of interest.
    level: float
        Level at which to draw contours on black top-hat filtered image.
        Default value is 0.3.
    selem: 2d-array or string
        Structuring element to use for the black top-hat filtering procedure
        Default value is a disk with a diameter of 20 pixels.
    perim_bounds: length 2 tuple
        Lower and upper perimeter bounds of approved objects. This should be
        in units of microns. The default values are 5 and 25 microns for the
        lower and upper bound, respectively.
    ip_dist : float
        Interpixel distance of the image in units of microns per pixel. The
        default value is 0.160 microns per pixel.
    area_bounds : tuple of float
        Upper and lower bounds for selected object areas. These should be
        given in units of square microns.
    ecc_bounds : tuple of float
        Bounds for object eccentricity. Default values are between 0.5 and 1.0.
    return_conts : bool
        If True, the x and y coordinates of the individual contours will be
        returned. Default value is False
    Returns
    -------
    im_lab : 2d-array, int
        Two dimensional image where each individual object is labeled.
    conts : 1d-array
        List of contour coordinates. Each entry of this array comes as
        an x,y pair of arrays. Has the same length as the number of
        contoured objects. This is only returned if `return_conts` is
        True.
    """

    # Apply the white top-hat filter.
    if selem == 'default':
        selem = skimage.morphology.disk(20)

    # Normalize the image.
    image = (image - image.min()) / (image.max() - image.min())

    # Blur and background subtract the image.
    im_blur = skimage.filters.gaussian(image, sigma=5)
    im_sub = image - im_blur

    # Apply the black tophat filter.
    im_filt = skimage.morphology.black_tophat(im_sub, selem)

    # Find the contours and return.
    conts = skimage.measure.find_contours(im_filt, level)

    # Make an empty image for adding the approved objects.
    objs = np.zeros_like(image)

    # Loop through each contour.
    for _, c in enumerate(conts):
        perim = 0
        for j in range(len(c) - 1):
            # Compute the distance between points.
            distance = np.sqrt((c[j + 1, 0] - c[j, 0])**2 +
                               (c[j + 1, 1] - c[j, 1])**2)
            perim += distance * ip_dist

        # Test if the perimeter is allowed by the user defined bounds.
        if (perim > perim_bounds[0]) & (perim < perim_bounds[1]):

            # Round the contours.
            c_int = np.round(c).astype(int)

            # Color the image with the contours and fill.
            objs[c_int[:, 0], c_int[:, 1]] = 1.0

    # Fill and label the objects.
    objs_fill = scipy.ndimage.binary_fill_holes(objs)
    objs_fill = skimage.morphology.remove_small_objects(objs_fill)
    im_lab = skimage.measure.label(objs_fill)

    # Apply filters.
    approved_obj = np.zeros_like(im_lab)
    props = skimage.measure.regionprops(im_lab, image)
    for prop in props:
        area = prop.area * ip_dist**2
        ecc = prop.eccentricity
        if (area < area_bounds[1]) & (area > area_bounds[0]) &\
            (ecc < ecc_bounds[1]) & (ecc > ecc_bounds[0]) &\
                (prop.mean_intensity < min_int):
            approved_obj += (im_lab == prop.label)
    im_lab = skimage.measure.label(approved_obj)

    if return_conts is True:
        return conts, im_lab
    else:
        return [im_lab, im_filt]



def find_zero_crossings(im, selem, thresh):
    """
    This  function computes the gradients in pixel values of an image after
    applying a sobel filter to a given image. This  function is later used in
    the Laplacian of Gaussian cell segmenter (log_segmentation) function. The
    arguments are as follows.
    Parameters
    ----------
    im : 2d-array
        Image to be filtered.
    selem : 2d-array, bool
        Structural element used to compute gradients.
    thresh :  float
        Threshold to define gradients.
    Returns
    -------
    zero_cross : 2d-array
        Image with identified zero-crossings.
    Notes
    -----
    This function as well as `log_segmentation` were written by Justin Bois.
    http://bebi103.caltech.edu/
    """

    # apply a maximum and minimum filter to the image.
    im_max = scipy.ndimage.filters.maximum_filter(im, footprint=selem)
    im_min = scipy.ndimage.filters.minimum_filter(im, footprint=selem)

    # Compute the gradients using a sobel filter.
    im_filt = skimage.filters.sobel(im)

    # Find the zero crossings.
    zero_cross = (((im >= 0) & (im_min < 0)) | ((im <= 0) & (im_max > 0)))\
        & (im_filt >= thresh)

    return zero_cross


# #################
def log_segmentation(im, selem=None, thresh=0.0001, radius=2.0,
                     median_filt=True, clear_border=True, label=False):
    """
    This function computes the Laplacian of a gaussian filtered image and
    detects object edges as regions which cross zero in the derivative.
    Parameters
    ----------
    im :  2d-array
        Image to be processed. Must be a single channel image.
    selem : 2d-array, bool
        Structural element for identifying zero crossings. Default value is
        a 2x2 pixel square.
    radius : float
        Radius for gaussian filter prior to computation of derivatives.
    median_filt : bool
        If True, the input image will be median filtered with a 3x3 structural
        element prior to segmentation.
    selem : 2d-array, bool
        Structural element to be applied for laplacian calculation.
    thresh : float
        Threshold past which
    clear_border : bool
        If True, segmented objects touching the border will be removed.
        Default is True.
    label : bool
        If True, segmented objecs will be labeled. Default is False.
    Returns
    -------
    im_final : 2d-array
        Final segmentation mask. If label==True, the output will be a integer
        labeled image. If label==False, the output will be a bool.
    Notes
    -----
    We thank Justin Bois in his help writing this function.
    https://bebi103.caltech.edu
    """

    # Test that the provided image is only 2-d.
    if len(np.shape(im)) > 2:
        raise ValueError('image must be a single channel!')

    # Determine if the image should be median filtered.
    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = scipy.ndimage.median_filter(im, footprint=selem)
    else:
        im_filt = im
    # Ensure that the provided image is a float.
    if np.max(im) > 1.0:
        im_float = skimage.img_as_float(im_filt)
    else:
        im_float = im_filt

    # Compute the LoG filter of the image.
    im_LoG = scipy.ndimage.filters.gaussian_laplace(im_float, radius)

    # Define the structural element.
    if selem is None:
        selem = skimage.morphology.square(3)

    # Using find_zero_crossings, identify the edges of objects.
    edges = find_zero_crossings(im_LoG, selem, thresh)

    # Skeletonize the edges to a line with a single pixel width.
    skel_im = skimage.morphology.skeletonize(edges)

    # Fill the holes to generate binary image.
    im_fill = scipy.ndimage.morphology.binary_fill_holes(skel_im)

    # Remove small objects and objects touching border.
    im_final = skimage.morphology.remove_small_objects(im_fill)
    if clear_border is True:
        im_final = skimage.segmentation.clear_border(im_final, buffer_size=5)

    # Determine if the objects should be labeled.
    if label is True:
        im_final = skimage.measure.label(im_final)

    # Return the labeled image.
    return im_final


def tophat_filter(image, 
                  large_selem_diam=50, 
                  small_selem_diam=2,
                  threshold='otsu'):
    """
    Performs a series of tophat transforms to amplify phase contrast images for 
    segmentation

    Parameters
    ==========
    image : 2d-array
        The phase contrast image to be filtered. 
    large_selem_diam : int
        The diameter of the large disk structuring element used in generating the 
        black tophat transform. Default is 50 pixels.
    small_selem_diam : int
        The diameter of the small disk structuring element used in generating the 
        white tophat transform. Default is 2 pixels.
    threshold: str or float
        Threshold used for rough segmentation of cells. If of type `float`, 
        provided value is used. If `otsu`, Otsu's method is used. if 'none',
        no thresholding is applied.

    Returns 
    ========
    closed_image : 2d-array
        The filtered and thresholded (if desired) image. 
    """

    # Normalize the image. 
    im_norm = (image - image.min()) / (image.max() - image.min())

    # Perform background subtraction to correct of uneven illumination
    im_blur = skimage.filters.gaussian(image, 30)
    im_sub = im_norm - im_blur

    # Define structuring elements
    lg_selem = skimage.morphology.disk(large_selem_diam)
    sm_selem = skimage.morphology.disk(small_selem_diam)

    # Perform filtering and closing operations
    blk_tophat = skimage.morphology.black_tophat(im_sub, lg_selem)
    wht_tophat = skimage.morphology.white_tophat(blk_tophat, lg_selem)
    closing = scipy.ndimage.grey_closing(wht_tophat, footprint=sm_selem)

    if threshold == 'otsu':
        thresh = skimage.filter.threshold_otsu(closing)
    elif threshold == 'none':
        return closing
    else:
        thresh = skimage.filter.threshold_otsu(closing)
    return closing * (closing > thresh)

def contour_segmentation(image, 
                         filter=True,
                         area_bounds=(10, 100),
                         ecc_bound=0.5,
                         solidity_bound=0.9,
                         ip_dist=0.0167,
                         return_mask=False,
                         return_cells=False,
                         **tophat_kwargs):
    """
    Segments an image using Laplacian of Gaussian then computes smoothed contours 
    of segmented objects. Contour curvature is also computed. 

    Parameters
    ==========
    image : 2d-array
        The image with dark cells to be segmented. 
    filter : bool
        If True, tophat filters will be applied with the supplied kwargs.
    area_bounds : tuple of positive floats 
        The bounds in 2d-projected cell area between which cells are accepted. 
        This must be in the same units as the interpixel distance. Default is
        between 10 and 100 square microns. 
    ecc_bound : positive_float
        Object eccentricity below which objects should be discarded.
        Default is 0.5
    solidity_bound : positive float
        The solidity below chich objects should be discarded. Solidity is 
        defined as the area fraction of the convex hull that is occupied by 
        the object. Default is 0.9.
    ip_dist : float
        The interpixel distance in the same units as `area_bounds`. Default is 
        0.0167 microns per pixel
        TODO: Give more information.
    return_mask : bool
        If True, a mask of segmented objects is returned.
    return_cells : bool
        If True, a dict is compiled containing the cropped segmentation mask 
        and the cropped intensity image. 

    Returns 
    ========
    objects :  pandas DataFrame
        A pandas DataFrame with information regarding  each successfully 
        segmented object.
    mask : 2d-array of ints
        The object segmentation mask after object filtering.  This is only 
        returned if `return_mask == True`.
    cells : dict of dicts
        A dictionary with each segmentation and intensity image for each cell,
        rotated to be aligned. This is only returned if `return_cells==True`.
    """
    if filter:
        _image = tophat_filter(image, **tophat_kwargs)
    else:
        _image = image
    # Perform the laplacian of gaussian segmentation
    log_selem = skimage.morphology.square(2)
    seg = log_segmentation(image, radius=1, selem=log_selem, 
                            thresh=0.001, median_filt=False, 
                            label=False)

    # Clean up the mask and label
    seg = skimage.morphology.binary_erosion(seg, log_selem)
    seg = skimage.morphology.remove_small_holes(seg)
    seg = skimage.morphology.remove_small_objects(seg)
    labeled = skimage.measure.label(seg)

    # Compute properties and create storage objects
    mask = np.zeros_like(image)
    props = skimage.measure.regionprops(labeled)
    cell_images = {}
    objects = pd.DataFrame([])

    # Iterate through each segmented object
    idx = 0
    for p in tqdm.tqdm(props):
        if ((p.area >= area_bounds[0]/ip_dist) & (p.area < area_bounds[1]/ip_dist)) &\
            (p.solidity > solidity_bound) & (p.eccentricity > ecc_bound):
            # Update the mask
            mask += labeled==p.label
            # Crop the original object and rotate.
            padded, _ = pad_bbox(p.bbox, np.shape(labeled), pad=10)
            rot = scipy.ndimage.rotate(labeled[padded]==p.label, -np.rad2deg(p.orientation), order=0) > 0
            rot = skimage.morphology.remove_small_holes(rot)
            rot = skimage.morphology.remove_small_objects(rot)
            relab = skimage.measure.label(rot.astype(int))
            if return_cells:
                rot_int = scipy.ndimage.rotate(_image, -np.rad2deg(p.orientation), order=0)
                rot_props = skimage.measure.regionprops(relab, intensity_image=rot_int)
            else:
                rot_props = skimage.measure.regionprops(relab)
            bbox = rot_props[0].bbox
            rot_pad, _ = pad_bbox(bbox, np.shape(rot), pad=10)
            
            # If an intensity image is desired, also rotate
            if return_cells:
                rot_intensity = scipy.ndimage.rotate(image[padded], -np.rad2deg(p.orientation),  mode='nearest')
                cell_images[idx] = {'intensity_image':rot_intensity[rot_pad],
                                    'segmentation_mask':rot[rot_pad]}                
            # Find contours and perform a uniform filtering of indices
            cont = skimage.measure.find_contours(rot[rot_pad], 0)[0]
            cx = scipy.ndimage.uniform_filter(cont[:, 1], 10, mode='wrap')
            cy = scipy.ndimage.uniform_filter(cont[:, 0], 10, mode='wrap')  

            # Compute the spline.
            tck, u = scipy.interpolate.splprep([cx, cy], per=1, k=5, s=100)
            unew = np.arange(0, 1.0001, 0.0001)
            out = scipy.interpolate.splev(unew, tck)

            # Compute the curvature and assemble the dataframe
            k = compute_curvature(out)
            _df = pd.DataFrame([]) 
            _df['x_coords'] = out[0][:-2]
            _df['y_coords'] = out[1][:-2]
            _df['curvature'] = k
            _df['cell_id'] = idx
            _df.dropna(inplace=True)
            objects = pd.concat([objects, _df], sort=False)
            # Update the cell counter
            idx += 1
            
    # Determine what needs to be returned
    out = [objects]
    if return_mask:
        out.append(mask)
    if return_cells:
        out.append(cell_images)
    if len(out) == 1:
        out = out[0]
    return out


def pad_bbox(bbox, dims, pad=10):
    """
    Compute a slice object that pads a bounding box. 

    Parameters
    ==========
    bbox : tuple
        The coordinates of the bounding box `(min_row, min_col, max_row, max_col)`.
    dims : tuple
        The dimensions of the original image to avoid overshoot of new bounding 
        box
    pad : int
        The number of pixels to add to each side of the bounding box. Default 
        is 10 pixels.
    
    Returns
    =======
    bbox_new :  (numpy slice object, padded_bbox)
        New bbox as a tuple of a slice index and new bbox coordinates)
    """
    x_lower = [bbox[0]-pad if (bbox[0] - pad) >= 0 else bbox[0]][0] 
    x_upper = [bbox[2]+pad if (bbox[2] + pad) <= dims[0]  else bbox[2]][0]
    y_lower = [bbox[1]-pad if (bbox[1] - pad) >= 0 else bbox[1]][0]
    y_upper = [bbox[3]+pad if (bbox[3] + pad) <= dims[1]  else bbox[3]][0]
    padded_bbox = (x_lower, y_lower, x_upper, y_upper)
    return [np.s_[x_lower:x_upper, y_lower:y_upper], padded_bbox] 


def compute_curvature(arr):
    # Pad the entry of each array to allow for boundary conditions
    dx = np.diff(arr[0])
    ddx = np.diff(dx)
    dy = np.diff(arr[1])
    ddy = np.diff(dy)
    k = (dx[:-1] * ddy  - dy[:-1] * ddx)/((dx[:-1]**2  + dy[:-1]**2)**(3/2))     
    return k


def assign_anatomy(data, 
              max_curve=0.5,
              columns = {'groupby' :'cell_id',
                         'curve'   : 'curvature',
                         'x'       : 'spl_x',
                         'y'       : 'spl_y'},
              ip_dist=0.0167 ):

    # Convert supplied curvature threshold to pixel value.
    curve_thresh = ip_dist/max_curve

    df = pd.DataFrame([])  
    for g, d in data.groupby(columns['groupby']):
        # Find the caps
        d = d.copy()
        caps = d[(d[columns['curve']] >= curve_thresh)].copy()
 
        # Determine contour points the cell planes
        bottom_cap_bound = caps[(caps[columns['y']] - caps[columns['y']].mean() ) < 0 ][columns['y']].max()
        upper_cap_bound = caps[(caps[columns['y']] - caps[columns['y']].mean() ) > 0 ][columns['y']].min()

        # Label caps as top and bottom
        caps['component'] = 'bottom'
        caps.loc[d[columns['y'] > bottom_cap_bound], 'component'] = 'top'

        # Find and label edges as left and right
        sides = d[(d[columns['y']] > bottom_cap_bound) & (d[columns['y']] < upper_cap_bound)].copy()
        sides['component'] = 'right'
        sides.loc[(sides[columns['x']] - sides[columns['x']].mean()) < 0, 
                'component'] = 'left'
        df = pd.concat([df, caps, sides], sort=False)
    return df

def measure_biometrics(data,                       
                       peri_width=0.025,
                       ip_dist = 0.0167,
                       columns = {'groupby' :'cell_id',
                                  'curve'   : 'curvature',
                                  'x'       : 'spl_x',
                                  'y'       : 'spl_y',
                                  'component': 'component'}):
    """
    Computes properties of cells with assigned anatomy.  
    """
    #TODO. Better way to calculate this length. maximize distance between sides 
    # connecting poles.
        
    # Compute width properties
    biometrics = pd.DataFrame([])
    for g, d in data.groupby(data[columns['groupby']]):
        # Make measurements
        length = (d[columns['y']].max() - d[columns['y']].min()) * ip_dist
        left = d[columns['component']=='left'][columns['y']].values
        right = d[columns['component']=='right'][columns['x']].values
        width = [np.min(rval - left) * ip_dist for rval in right] 
        width_mean = np.mean(width)
        width_var  = np.var(width)

        # Compute features 
        vol = volume(length, width_mean)
        sa = surface_area(length, width_mean)
        sav = sa / vol
        env_vol = envelope_volume(length, width_mean, peri_width)
        frac_vol = env_vol / vol

        # Assemble dataframe
        _df = pd.DataFrame([])
        _df['length'] = length
        _df['width_mean'] = width_mean
        _df['width_var'] = width_var
        _df['volume'] = vol
        _df['surface_area'] = sa
        _df['surface_to_volume'] = sav
        _df['periplasm_volume'] = env_vol
        _df['periplasm_fractional_volume'] = frac_vol
        _df[columns['groupby']] = g
        biometrics = pd.concat([_df, biometrics])
    return biometrics

