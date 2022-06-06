import numpy as np
import scipy.ndimage 
import skimage.io 
import skimage.measure
import skimage.morphology
import skimage.segmentation

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