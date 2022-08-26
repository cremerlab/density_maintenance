import numpy as np
import scipy.ndimage
import skimage.io
import skimage.measure
import skimage.morphology
import skimage.segmentation
import skimage.filters
import warnings
from .analytical import *
import tqdm
import pandas as pd
import cv2


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

    blk_tophat = cv2.morphologyEx(im_sub, cv2.MORPH_BLACKHAT, lg_selem)
    #skimage.morphology.black_tophat(im_sub, lg_selem)
    wht_tophat = cv2.morphologyEx(blk_tophat, cv2.MORPH_TOPHAT, lg_selem)
    #skimage.morphology.white_tophat(blk_tophat, lg_selem)
    closing = scipy.ndimage.grey_closing(wht_tophat, footprint=sm_selem)

    if threshold == 'otsu':
        thresh = skimage.filters.threshold_multiotsu(closing)[1]
    elif threshold == 'none':
        return closing
    else:
        thresh = threshold
    return closing * (closing > thresh)


def contour_segmentation(image,
                         filter=True,
                         area_bounds=(1, 1000),
                         ecc_bound=0.5,
                         solidity_bound=0.9,
                         perim_bounds=(2, 15),
                         ip_dist=0.0319,
                         return_mask=False,
                         return_cells=False,
                         intensity_image=None,
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
    perim_bound : tuple of positive floats
        The object perimeter between which cells should be kept. This should be 
        provided in distance units that can be converted to pixel number 
        through the interpixel distance.
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
    cell_df : pandas DataFrame
        A  pandas DataFrame containing cell ID number, segmentation mask, 
        and intensity image.

    """
    if filter:
        _image = tophat_filter(image, **tophat_kwargs)
    else:
        _image = image
    # Perform the laplacian of gaussian segmentation
    log_selem = skimage.morphology.square(2)
    seg = log_segmentation(_image, radius=1, selem=log_selem,
                           thresh=0.001, median_filt=False,
                           label=False)

    # Clean up the mask and label
    seg = skimage.morphology.remove_small_holes(seg)
    seg = skimage.morphology.remove_small_objects(seg)
    labeled = skimage.measure.label(seg)

    # Compute properties and create storage objects
    mask = np.zeros_like(image)
    props = skimage.measure.regionprops(labeled)
    cell_df = pd.DataFrame([])
    objects = pd.DataFrame([])

    # Iterate through each segmented object
    idx = 0
    for p in tqdm.tqdm(props):
        area = p.area_filled * ip_dist**2
        perim = p.perimeter * ip_dist
        if (area >= area_bounds[0]) & (area <= area_bounds[1]) &\
                (p.solidity >= solidity_bound) & (p.eccentricity >= ecc_bound) &\
                ((perim >= perim_bounds[0]) & (perim <= perim_bounds[1])):
            # Update the mask
            mask += labeled == p.label
            # Crop the original object and rotate.
            padded, _ = pad_bbox(p.bbox, np.shape(labeled), pad=10)
            rot = scipy.ndimage.rotate(
                labeled[padded] == p.label, -np.rad2deg(p.orientation), order=0) > 0
            erode_selem = skimage.morphology.square(2)
            rot = skimage.morphology.binary_erosion(rot, erode_selem)
            rot = skimage.morphology.remove_small_holes(rot)
            rot = skimage.morphology.remove_small_objects(rot)
            rot = scipy.ndimage.binary_fill_holes(rot)

            relab = skimage.measure.label(rot.astype(int))
            if return_cells:
                if intensity_image is None:
                    rot_int = scipy.ndimage.rotate(
                        image[padded], -np.rad2deg(p.orientation), order=0, mode='nearest')
                else:
                    rot_int = scipy.ndimage.rotate(
                        intensity_image[padded], -np.rad2deg(p.orientation), order=0, mode='nearest')
                rot_props = skimage.measure.regionprops(
                    relab, intensity_image=rot_int)
            else:
                rot_props = skimage.measure.regionprops(relab)
            if len(rot_props) == 0:
                continue
            bbox = rot_props[0].bbox
            rot_pad, _ = pad_bbox(bbox, np.shape(rot), pad=10)

            # If an intensity image is desired, also rotate
            if return_cells:
                _cell_df = pd.DataFrame([], index=[0])
                _cell_df['cell_id'] = idx
                _cell_df['cell_image'] = [rot_int[rot_pad]]
                _cell_df['mask'] = [rot[rot_pad]]
                cell_df = pd.concat([cell_df, _cell_df], sort=False)

            # Find contours and perform a uniform filtering of indices
            cont = skimage.measure.find_contours(rot[rot_pad], 0)[0]
            cx = scipy.ndimage.uniform_filter(cont[:, 1], 10, mode='wrap')
            cy = scipy.ndimage.uniform_filter(cont[:, 0], 10, mode='wrap')
            if len(cx) < 10:
                continue
            # Compute the spline.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tck, _ = scipy.interpolate.splprep([cx, cy], per=1, k=3, s=80)
                unew = np.arange(0, 1.0001, 0.0001)
                out = scipy.interpolate.splev(unew, tck)

            # Compute the curvature and assemble the dataframe
            k = compute_curvature(out)
            _df = pd.DataFrame([])
            _df['x_coords'] = out[0][:-2]
            _df['y_coords'] = out[1][:-2]
            _df['curvature'] = k / ip_dist
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
        out.append(cell_df)
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
    x_upper = [bbox[2]+pad if (bbox[2] + pad) <= dims[0] else bbox[2]][0]
    y_lower = [bbox[1]-pad if (bbox[1] - pad) >= 0 else bbox[1]][0]
    y_upper = [bbox[3]+pad if (bbox[3] + pad) <= dims[1] else bbox[3]][0]
    padded_bbox = (x_lower, y_lower, x_upper, y_upper)
    return [np.s_[x_lower:x_upper, y_lower:y_upper], padded_bbox]


def compute_curvature(arr):
    # Pad the entry of each array to allow for boundary conditions
    dx = np.diff(arr[0])
    ddx = np.diff(dx)
    dy = np.diff(arr[1])
    ddy = np.diff(dy)
    k = (dx[:-1] * ddy - dy[:-1] * ddx)/((dx[:-1]**2 + dy[:-1]**2)**(3/2))
    return k


def assign_anatomy(data,
                   cap_radius=0.5,
                   columns={'groupby': 'cell_id',
                            'curve': 'curvature',
                            'x': 'x_coords',
                            'y': 'y_coords'}):

    # Convert supplied curvature threshold to pixel value.

    df = pd.DataFrame([])
    for g, d in data.groupby(columns['groupby']):
        # Find the caps
        d = d.copy()
        caps = d[d[columns['curve']] >= 1/cap_radius].copy()

        # Determine contour points the cell planes
        bottom_cap_bound = caps[(
            caps[columns['y']] - caps[columns['y']].mean()) < 0][columns['y']].max()
        upper_cap_bound = caps[(
            caps[columns['y']] - caps[columns['y']].mean()) > 0][columns['y']].min()

        # Label caps as top and bottom
        d['component'] = 'bottom'
        d.loc[d[columns['y']] > bottom_cap_bound, 'component'] = 'top'

        # Find and label edges as left and right
        _caps = d[(d[columns['y']] < bottom_cap_bound) | (
            d[columns['y']] > upper_cap_bound)].copy()
        sides = d[(d[columns['y']] > bottom_cap_bound) & (
            d[columns['y']] < upper_cap_bound)].copy()
        sides['component'] = 'right'
        sides.loc[(sides[columns['x']] - sides[columns['x']].mean()) < 0,
                  'component'] = 'left'
        df = pd.concat([df, _caps, sides], sort=False)
    return df


def measure_biometrics(data,
                       peri_width=0.025,
                       ip_dist=0.0319,
                       columns={'groupby': 'cell_id',
                                'curve': 'curvature',
                                'x': 'x_coords',
                                'y': 'y_coords',
                                'component': 'component'}):
    """
    Computes properties of cells with assigned anatomy.  
    """
    # TODO. Better way to calculate this length. maximize distance between sides
    # connecting poles.

    # Compute width properties
    biometrics = pd.DataFrame([])
    for g, d in data.groupby(data[columns['groupby']]):
        # Make measurements
        length = (d[columns['y']].max() - d[columns['y']].min()) * ip_dist
        left_x = d[d[columns['component']] == 'left'][columns['x']].values
        left_y = d[d[columns['component']] == 'left'][columns['y']].values
        right_x = d[d[columns['component']] == 'right'][columns['x']].values
        right_y = d[d[columns['component']] == 'right'][columns['y']].values
        if (len(left_x) == 0) | (len(left_y) == 0):
            continue

        def hypot(p1, p2):
            return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p2[0])**2)
        width = [np.min(hypot((_x, _y), (left_x, left_y))) *
                 ip_dist for _x, _y in zip(right_x, right_y)]
        width_median = np.median(width)
        width_var = np.var(width)

        # Compute features
        vol = volume(length, width_median)
        sa = surface_area(length, width_median)
        sav = sa / vol
        env_vol = envelope_volume(length, width_median, peri_width)
        frac_vol = env_vol / vol

        # Assemble dataframe
        _df = pd.DataFrame(
            {'length': length,
             'width_median': width_median,
             'width_var': width_var,
             'volume': vol,
             'surface_area': sa,
             'surface_to_volume': sav,
             'periplasm_volume': env_vol,
             'periplasm_fractional_volume': frac_vol,
             columns['groupby']: g}, index=[0])
        biometrics = pd.concat([biometrics, _df])
    return biometrics
