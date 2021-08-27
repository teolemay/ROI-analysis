"""
Teophile Lemay; 2021
In this file, I keep a record of the functions I have written for analyzing fluorescence spectroscopy
images
"""

import os
import numpy as np
import tifffile as tf
from scipy import signal
from scipy import ndimage
from skimage import io
import xml.etree.cElementTree as ET
import re
import zipfile
import cv2 as cv
import pandas as pd
import matplotlib.pyplot as plt


#OPENING/ACCESING FILES, FILE INFO
def read_multipage_tiff(tif_path):
    """
    Read multipage tiff file (INO file)
    uses skimage.io.ImageCollection
    flim is FLIM cube
    hyperspectal is hyperspectral cube

    :param tif_path: full path to multipage tiff file
    :return: tuple containing flim, hyperspectral cube
    """
    img = io.ImageCollection(tif_path)
    flim = img[0]
    hyperspectral = img[1]
    return flim, hyperspectral

def get_tiff_metadata(tif_path, as_string=False, show_arrays=False):
    """
    this function will get the metadata from a tiff file (INO file)
    uses the tifffile library to read 

    :param tif_path: full path to tiff file
    :param kwarg as_string: deffault is false, if true function will return metadata as strings instead of dicts
    :paran kwarg show_array: default is False, if true, function will also return image cube arrays
    :return: tupple containging temproal (0) and hyperspectral (1) metadata dictionaries
    if show_arrays == True:
    :return: tupple containging temporal (0) and hyperspectral (1) cube arrays as after metadata
    """
    tif = tf.TiffFile(tif_path)
    flim, spectral = tif.pages[0], tif.pages[1]
    flim_tags = {}
    spec_tags = {}
    for tag in flim.tags.values():
        key, value = tag.name, tag.value
        flim_tags[key] = value
    for tag in spectral.tags.values():
        key, value = tag.name, tag.value
        spec_tags[key] = value
    metadata = (flim_tags, spec_tags)
    if as_string:
        flim_str = ''
        for tag in flim_tags:
            flim_str += f'{tag}: {flim_tags[tag]}\n\n'
        spec_str = ''
        for tag in spec_tags:
            spec_str += f'{tag}: {spec_tags[tag]}\n\n'
        metadata = (flim_str, spec_str)
    if show_arrays:
        arrays = ( flim.asarray(), spectral.asarray() )
        return metadata, arrays
    return metadata

def get_IRF_data(temp_metadata):
    """
    This function will parse through the metadata for INO temporal image cube, and returns an array of IRF values

    :param temp_metadata: metadata for temporal image cube in dict form
    :return IRF_array: array of IRF values in order from ImageDescription
    """
    ParseError = ET.ParseError
    xmlstring = temp_metadata['ImageDescription']
    #remove bad characters
    xmlstring = (lambda x: re.sub(u'[\x00-\x08\x0b\x0c\x0e-\x1F\uD800-\uDFFF\uFFFE\uFFFF]', '', x))(xmlstring)
    try:
        root = ET.fromstring(xmlstring)
    except ParseError:
        fix_here = xmlstring.find('xmlns')+5
        xmlstring = xmlstring[:fix_here] + ':xsi' + xmlstring[fix_here:]
        root = ET.fromstring(xmlstring)
    try:
        hist = root.find('All').find('IRF').find('Histogram')
    except:
        print('Metadata does not have expected format')
        return None
    IRF_array = np.zeros(len(hist)-1)
    for i in range(len(hist)-1):
        IRF_array[i] = hist[i].text
    return IRF_array

def get_ISS_cubes(ISS_path, series_count=1, channels=2, cube_dims=(256, 256, 256)):
    """
    this function will return the ISS image cube for specified ISS-TDFLIM file.

    :param ISS_path: path to ISS file or file name if file in working directory
    :param series_count (optional): count for position/time series. (default is 1)
    :param channels (optional): number of channels (number of different cubes to make). default=2
    :param cube_dims (optional): dimensions of FLIM cube. default= 256x256x256
    
    :return: ISS image cubes in single ndarray. first index of the array designates separate cubes
    """
    slice_size = 2 * channels * cube_dims[0] * cube_dims[1] * cube_dims[2]
    stack = np.empty((series_count, channels, cube_dims[0], cube_dims[1], cube_dims[2]), dtype='uint16') #the dtype could be changed to uint8 if values are always low = less ram used
    with zipfile.ZipFile(ISS_path) as iss_file:
        with iss_file.open('data/PrimaryDecayData.bin', 'r') as bin_file:
            for i in range(series_count):
                slice = bin_file.read(slice_size)
                slice = np.frombuffer(slice, dtype='uint16')
                cubes = slice.reshape(tuple([channels] + [i for i in cube_dims]))
                stack[i] = cubes
    return stack

def get_ISS_metadata(ISS_path):
    """
    get core metadata for ISS file as string

    :param ISS_path: full path to ISS file (or file name if file is in working directory)

    :return: ISS file core metadata
    """
    with zipfile.ZipFile(ISS_path, 'r') as iss:
        corexml = iss.read('dataProps/Core.xml').decode()
    return corexml

def open_ISS(file_path):
    """
    this function will open and get all relevant information from an ISS file.
    cubes array may be up to 5-dimensional: [slice, acceptor/donor, widht, height, depth]

    :param file_path: path to ISS file of interest
    :return: format string, ndarray containing all the cubes, number of cubes, flim time resolution, physical length of a pixel, unit of pixel length, metadata
    """
    metadata = get_ISS_metadata(file_path)
    try:
        #could add values for time series counts, frame counts. easy to add later if needed
        root = ET.fromstring(metadata)
        dimensions = root.find('Dimensions')
        num_channels = int(dimensions.find('ChannelCount').text)
        cube_height = int(dimensions.find('FrameHeight').text)
        cube_width = int(dimensions.find('FrameWidth').text)
        cube_depth = int(root.find('PhotonCountingSettings').find('AdcResolution').text)
        real_width = float(root.find('Boundary').find('FrameWidth').text)
        real_height = float(root.find('Boundary').find('FrameHeight').text)
        laser_freq = float(root.find('PhotonCountingSettings').find('MacroTimeClockFrequency').text)
        dimensions_unit = root.find('CoordUnitType').text
        tac_time_range = float(root.find('PhotonCountingSettings').find('TacTimeRange').text) #nanoseconds
        pixel_dwell_time = float(root.find('PixelDwellTime').text)
        series_count = int(root.find('Dimensions').find('TimeSeriesCount').text) #always 1 unless multiple time series
        position_series_info = root.find('ExtraDimensions').find('PositionSeriesInfo')
        if position_series_info != None:
            series_count = int(position_series_info.find('Count').text)*series_count
    except:
        print('Metadata does not have expected format')
        return None 

    num_cubes = f'{series_count} position/time series; {num_channels} channels; {series_count * num_channels} total cubes'
    flim_time_resolution = (tac_time_range / cube_depth)*1e-9 #picoseconds
    pixel_len = real_width / cube_width #pixels are square so this is all we need
    cubes = get_ISS_cubes(file_path, series_count=series_count, channels=num_channels, cube_dims=(cube_height, cube_width, cube_depth))
    format = f'ISS Opened File Array Format: \n0. Format string \n1. image cubes array \n2. number of cubes \n3. flim time resolution (seconds) \n4. laser frequency \n5. pixel size ({dimensions_unit}) \n6. pixel dwell time (microseconds) \n7. metadata'
    return format, cubes, num_cubes, flim_time_resolution, laser_freq, pixel_len, pixel_dwell_time, metadata

def open_INO(file_path):
    """
    this function will open and get all relevant information from an INO file.

    :param file_path: path to INO file of interest
    :return: format string, ndarray containing both cubes, number of cubes, flim time resolution, IRF array
    """
    #metadata and arrays
    metadata, cubes = get_tiff_metadata(file_path, show_arrays=True)
    num_cubes = 2
    flim_data = metadata[0]
    flim_cube = cubes[0]
    flim_dims = flim_cube.shape
    flim_time_resolution = flim_data['XResolution'][1] / flim_dims[2] * 1e-12  #picoseconds to seconds
    #need to figure out how to definee spacial resolution units
    spectral_data = metadata[1]
    spectral_cube = cubes[1]
    spectral_dims = spectral_cube.shape
    #IRF
    IRF = get_IRF_data(flim_data)
    #pixel information
    ParseError = ET.ParseError
    xmlstring = flim_data['ImageDescription']
    #remove bad characters
    xmlstring = (lambda x: re.sub(u'[\x00-\x08\x0b\x0c\x0e-\x1F\uD800-\uDFFF\uFFFE\uFFFF]', '', x))(xmlstring)
    try:
        root = ET.fromstring(xmlstring)
    except ParseError:
        fix_here = xmlstring.find('xmlns')+5
        xmlstring = xmlstring[:fix_here] + ':xsi' + xmlstring[fix_here:]
        root = ET.fromstring(xmlstring)
    try:
        pixel_len = float(root.find('All').find('OpticalCore').find('PixelSize').text) #micrometer
        pixel_dwell_time = float(root.find('All').find('OpticalCore').find('PixelTimeUs').text)
        laser_freq = float(root.find('All').find('OpticalCore').find('LaserFrequency').text)*1e6 #value is given in MHz, converted here to Hz
    except:
        print('Metadata does not have expected format')
        return None
    format = f'INO Opened File Array Format: \n0. Format string \n1. image cubes array \n2. number of cubes \n3. flim time resolution (seconds) \n4. laser frequency (Hz) \n5. pixel size (micrometers) \n6. pixel dwell time (microseconds) \n7. IRF \n8. metadata'
    return format, cubes, num_cubes, flim_time_resolution, laser_freq, pixel_len, pixel_dwell_time, IRF, metadata

def open_file(file_path):
    """
    this function will open ISS and INO files and extract the relevant data

    :param file_path: path to desired file (.tif or .iss-tdflim)
    :dtype: string
    """
    if file_path[-4:] == '.tif':
        return open_INO(file_path)
    elif file_path[-11:] == '.iss-tdflim':
        return open_ISS(file_path)
    else:
        print('ERROR unexpected file type')

def list_INO_in_folder(INOfolder_path):
    """
    this function makes a list of all properly named INO files in a folder.

    :param INOfolder_path: path to folder

    :return: list of INO file names in folder.
    """
    files = os.listdir(INOfolder_path)
    INO_files = []
    for file in files:
        if file[:6] == 'WellID' and file[-4:] == '.tif':
            INO_files.append(file)
    return INO_files

def separate_wellIDs(INOfolder_path):
    """
    this function takes the path to a folder containing INO files and separates them according to well ID. the function returns a dictionary where the keys are IDs and the values are lists of INO file names with those IDs

    :param INOfolder_path: path to INO folder

    :return: dictionary of well ids and INO files 
    """
    files = os.listdir(INOfolder_path)
    wells = {}
    for file in files:
        if file[:6]  == 'WellID':
            fullID = file[file.find('_')+1:file.find('-')] #underscore _  is before ID; dash - marks the end of the ID
            if fullID in wells:
                wells[fullID].append(file)
            else:
                wells[fullID] = [file]
    return wells

def parse_INO_metadata(file_path):
    """
    This function collects relevant metadata from an INO file.

    :param file_path: path to INO file of interest
    :return: laser frequency (Hz), flim time resolution (s), IRF (array), spectral map (array), pixel_dwell time (us), pixel length (um)
    """
    #metadata and arrays
    metadata = get_tiff_metadata(file_path, show_arrays=False)
    flim_data = metadata[0]
    flim_time_resolution = (flim_data['XResolution'][1] / flim_data['SamplesPerPixel'])*1e-12 #value is given in ps, converted here to seconds
    spectral_data = metadata[1]
    #pixel information
    ParseError = ET.ParseError
    flimxml = flim_data['ImageDescription']
    spectralxml = spectral_data['ImageDescription']
    #remove bad characters
    flimxml = (lambda x: re.sub(u'[\x00-\x08\x0b\x0c\x0e-\x1F\uD800-\uDFFF\uFFFE\uFFFF]', '', x))(flimxml)
    spectralxml = (lambda x: re.sub(u'[\x00-\x08\x0b\x0c\x0e-\x1F\uD800-\uDFFF\uFFFE\uFFFF]', '', x))(spectralxml)
    try:
        flimroot = ET.fromstring(flimxml)
    except ParseError:
        fix_here = flimxml.find('xmlns')+5
        flimxml = flimxml[:fix_here] + ':xsi' + flimxml[fix_here:]
        flimroot = ET.fromstring(flimxml)
    try:
        spectralroot = ET.fromstring(spectralxml)
    except ParseError:
        fix_here = spectralxml.find('xmlns')+5
        spectralxml = spectralxml[:fix_here] + ':xsi' + spectralxml[fix_here:]
        spectralroot = ET.fromstring(spectralxml)
    IRFhist = flimroot.find('All').find('IRF').find('Histogram')
    IRF_array = np.zeros(len(IRFhist)-1)
    for i in range(len(IRF_array)):
        IRF_array[i] = float(IRFhist[i].text)
    DetectWavelens = spectralroot.find('All').find('HyperSpectralCalibrationData').find('Info').find('Spectral').find('DetectorWavelengths')
    spectral_map = np.zeros(len(DetectWavelens)-1)
    for i in range(len(spectral_map)):
        spectral_map[i] = float(DetectWavelens[i].text)
    pixel_len = float(flimroot.find('All').find('OpticalCore').find('PixelSize').text) #micrometer
    pixel_dwell_time = float(flimroot.find('All').find('OpticalCore').find('PixelTimeUs').text)
    laser_freq = float(flimroot.find('All').find('OpticalCore').find('LaserFrequency').text)*1e6 #value is given in MHz, converted here to Hz
    return laser_freq, flim_time_resolution, IRF_array, spectral_map, pixel_dwell_time, pixel_len

def parse_ISS_metadata(metadata):
    """
    this function parses the ISS metadata and collects relevant information for calculating non structured data

    :param metadata: metadata string from ISS file

    :return: laser frequency (Hz), FLIM time resolution (s), pixel_dwell_time (us), pixel_len (um)
    """
    root = ET.fromstring(metadata)
    laser_freq = float(root.find('PhotonCountingSettings').find('MacroTimeClockFrequency').text)
    tac_time_range = float(root.find('PhotonCountingSettings').find('TacTimeRange').text)
    cube_depth = float(root.find('PhotonCountingSettings').find('AdcResolution').text)
    flim_time_resolution = (tac_time_range / cube_depth)*1e-9 #calculated in 
    pixel_dwell_time = float(root.find('PixelDwellTime').text)
    cube_width = int(root.find('Dimensions').find('FrameWidth').text)
    real_width = float(root.find('Boundary').find('FrameWidth').text)
    pixel_len = real_width / cube_width #pixels are square so this is all we need

    return laser_freq, flim_time_resolution, pixel_dwell_time, pixel_len

def get_spectral_map(file_path):
    """
    This function returns the list of spectral wavelengths for 

    :param file_path: path to INO file of interest
    :return: array containing wavelengths for each bin
    """
    #metadata and arrays
    metadata = get_tiff_metadata(file_path, show_arrays=False)
    spectral_data = metadata[1]
    #pixel information
    ParseError = ET.ParseError
    spectralxml = spectral_data['ImageDescription']
    #remove bad characters
    spectralxml = (lambda x: re.sub(u'[\x00-\x08\x0b\x0c\x0e-\x1F\uD800-\uDFFF\uFFFE\uFFFF]', '', x))(spectralxml)
    try:
        spectralroot = ET.fromstring(spectralxml)
    except ParseError:
        fix_here = spectralxml.find('xmlns')+5
        spectralxml = spectralxml[:fix_here] + ':xsi' + spectralxml[fix_here:]
        spectralroot = ET.fromstring(spectralxml)
    DetectWavelens = spectralroot.find('All').find('HyperSpectralCalibrationData').find('Info').find('Spectral').find('DetectorWavelengths')
    spectral_map = np.zeros(len(DetectWavelens)-1)
    for i in range(len(spectral_map)):
        spectral_map[i] = float(DetectWavelens[i].text)
    return spectral_map


#SEGMENTATION
def LoG_kernel(kernel_size=50, sigma=0.25):
    """
    Generage a Laplacian of a Gaussian Kernel
   
    :param kernel_size: [size of kernel], defaults to 51
    :type kernel_size: int, optional
    :param sigma: [Gaussian standard deviation], defaults to 3.25
    :type sigma: float, optional
    :return: [2D numpy array Laplacian of Gaussian ]
    :rtype: [Numpy float32 array]
    """
    log_kernel = np.zeros((kernel_size,kernel_size))
    s2 = sigma*sigma
    s6 = s2 * s2 * s2
    half_size = np.floor(kernel_size/2.0)
    x_arr = np.linspace(-half_size,half_size,kernel_size)
    y_arr = x_arr.copy()
    for i,x in enumerate(x_arr):
        for j,y in enumerate(y_arr):
            h_g = np.exp(-1.0*(x*x + y*y)/(2*s2))
            h = (x*x + y*y - 2.0*s2) * h_g / (2* np.pi*s6)
            log_kernel[i,j]=h
    log_kernel = log_kernel / np.sum(log_kernel)
    return log_kernel

def get_binary_map(flim_img, kernel_size=50, sigma=0.25, bin_map_factor=1, noise_level=20, small_roi_radius=5):
    """
    1- Using a LoG filter, identify signal regions
    2- Using  mean filter create a binary map from convolved image
    3- Remove background noise pixels
    4- remove regions that are smaller than given radius

    :param flim_img: input image
    :type flim_img: numpy 2D array 
    :param noise_level: background noise level in image, defaults to 20
    :type noise_level: int, optional
    :param small_roi_radius: radius of small regions removed, defaults to 5
    :type small_roi_radius: int, optional
    :return: binary segmentation map
    :rtype: uint8 binary [0,255] 2D numpy array
    """
    log = LoG_kernel(kernel_size=kernel_size, sigma=sigma) 
    conv_img = cv.filter2D(np.float32(flim_img),cv.CV_32F,log,None)
    conv_img = cv.normalize(conv_img,None,0,1,cv.NORM_MINMAX,cv.CV_32F)
    conv_img_mean = np.mean(conv_img) * bin_map_factor
    binary_mask = np.ones_like(flim_img)
    binary_mask[conv_img<conv_img_mean]=0
    binary_mask[flim_img < noise_level] = 0 
    ROI_size = small_roi_radius
    disk = cv.getStructuringElement(cv.MORPH_ELLIPSE,(ROI_size,ROI_size))
    opened_binary_mask = cv.morphologyEx(np.float32(binary_mask), cv.MORPH_OPEN, disk)
    opened_binary_mask = np.uint8(opened_binary_mask*255)
    return opened_binary_mask

def get_local_maxima(img,roi_size=10):
    """
    Get local maxima pixels in an image
    using a dilation operation

    :param img: input image
    :type img: 2d numpy array
    :param roi_size: size of structuring element, defaults to 10
    :type roi_size: int, optional
    :return: binary map with single pixel local minima locatoin
    :rtype: 2D numpy array
    """
    float_img = np.float32(img.copy())
    disk = cv.getStructuringElement(cv.MORPH_ELLIPSE, (roi_size,roi_size))
    dilated = cv.dilate(float_img,disk)
    local_max = cv.compare(float_img,dilated,cv.CMP_GE)
    return local_max

def multiseed_watershed_segmentation(binary_map, seed_map):
    """
    Segment binary map using multiseed binary map
    
    :param binary_map : binary map showing cells as foreground
    :param seed_map: binary map containing binary points for seed
    :type binary_map: 2d numpy array (uint8)
    :type binary_roi: 2d numpy array (uint8)
    :return: segmentation map (pixels belonging to each roi have an integer assigned to them)
    :rtype: 2D numpy array
    """
    ret, markers = cv.connectedComponents(seed_map)
    binary_map_3ch = cv.merge((binary_map,binary_map,binary_map))
    segmap = cv.watershed(binary_map_3ch,markers)
    segmap += 1
    segmap = np.uint16(segmap)
    segmap[binary_map==0]=0
    segmap = ndimage.filters.maximum_filter(segmap,size=2)
    segmap[binary_map==0]=0
    return segmap

def watershed_segment_image(flim_img, kernel_size=50, sigma=0.25, bin_map_factor=1, noise_level=20, small_roi_radius=5, roi_size=10, log_binary=True):
    """
    this function groups the watershed segmentation functions together. The function also renames the segments from 1 to number of ROIs
    Parameters are mostly the same as for the individual functions other than log_binary

    :param log_binary: bool. if True, binary map will be made from log of flim image

    :return: segmented ROI map, array of non-zero ROI labels.
    """
    if log_binary:
        noise_level=np.log(noise_level)
        log_img = np.log(flim_img, where=(flim_img != 0))
        binary_map = get_binary_map(log_img, kernel_size=kernel_size, sigma=sigma, bin_map_factor=bin_map_factor, noise_level=noise_level, small_roi_radius=small_roi_radius)
    else:
        binary_map = get_binary_map(flim_img, kernel_size=kernel_size, sigma=sigma, bin_map_factor=bin_map_factor, noise_level=noise_level, small_roi_radius=small_roi_radius)
    local_maxima = get_local_maxima(flim_img, roi_size=roi_size)
    segmented = multiseed_watershed_segmentation(binary_map, local_maxima)
    labels = np.unique(segmented)
    for i, label in enumerate(labels[1:]):
        segmented[segmented == label] = i + 1
    labels = np.unique(segmented)[1:]
    return segmented, labels

def threshold_ROIs(segmented_img, intensity_img, mean_scalar=1):
    """
    this function applies mean thresholding to pixels within each roi of a segmented image

    :param segmented_img: segmented image with different value for each segment
    :param intensity_img: intensity image
    :param mean_scalar
    """
    segmented = segmented_img.copy()
    segments = np.unique(segmented)
    for i, seg in enumerate(segments[1:]): #skip the 0 value background
        row, col = np.where(segmented==seg)
        vals = np.zeros(row.shape)
        for j in range(len(row)):
            vals[j] = intensity_img[row[j], col[j]]
        roi_mean_val = np.mean(vals)
        for j in range(len(row)):
            if intensity_img[row[j], col[j]] < roi_mean_val*mean_scalar:
                segmented[row[j], col[j]] = 0
    return segmented

def otsu_segment_image(image, kernel_size=50, sigma=0.25, small_roi_radius=4, log_binary=False):
    """
    this function quickly segments an image using otsu binarization and connected components

    :param image: image to segment
    :param min_ROI_size: int, optional

    :return: segmented image, ROI labels
    """
    log = LoG_kernel(kernel_size=kernel_size, sigma=sigma) 
    image = cv.filter2D(np.float32(image),cv.CV_32F,log,None)  #I added the LoG kernel here, maybe it will work nicely
    if log_binary:
        image = np.log(image, where=(image!=0))
        maxi = np.max(image)
        image = image / maxi * 255
        rescaled = image.astype(np.uint8)
    else:
        maxi = np.max(image)
        image = image / maxi * 255
        rescaled = image.astype(np.uint8)
    ret, binary = cv.threshold(rescaled, 0, 255, cv.THRESH_BINARY+cv.THRESH_OTSU)
    disk = cv.getStructuringElement(cv.MORPH_ELLIPSE,(small_roi_radius, small_roi_radius))
    opened_binary =  cv.morphologyEx(binary, cv.MORPH_OPEN, disk)
    ret, segmented_img = cv.connectedComponents(opened_binary)
    labels = np.unique(segmented_img)[1:]
    return segmented_img, labels

def segmentation_choice(image, method, kernel_size=50, sigma=0.25, bin_map_factor=1, small_roi_radius=5, noise_level=20, roi_size=10, log_binary=False, threshold=False, mean_scalar=1):
    """
    this function allows the user to choose whether they want to use watershed or otsu for binarization and segmentation.
    care must be taken to provide the proper keyword arguments for the desired segmentation method

    :param image: 2d array. image for segmentation
    :param method: desired segmentation method
    :param small_roi_radius: radius of disk for cv.morph open
    :param noise_level: noise threshold for watershed binary map function
    :param roi_size: roi size for watershed local maxima function
    :param log_binary: bool. if True, a log operation will be applied on the image for the binary map in watershed segmentation

    :return: segmented image and ROI_labels
    """
    if method == 'watershed':
        segmented, labels = watershed_segment_image(image, kernel_size=kernel_size, sigma=sigma, bin_map_factor=bin_map_factor, noise_level=noise_level, small_roi_radius=small_roi_radius, roi_size=roi_size, log_binary=log_binary)
    elif method == 'otsu':
        segmented, labels = otsu_segment_image(image, small_roi_radius=small_roi_radius, log_binary=log_binary)
    else:
        raise ValueError('Unknown segmentation "method" selected. options are "watershed" OR "otsu"')
    if threshold:
        segmented = threshold_ROIs(segmented, image, mean_scalar=mean_scalar)
    return segmented, labels


#LIFETIME ANALYSIS
def com_lifetime(FLIM_cube, IRF, roi_map, roi_labels, time_res=64e-12):
    """
    SIMPLE BUT INEFFICIENT
    this function calculates a lifetime map using center of mass analysis, given an IRF
    :param FLIM_cube: FLIM cube for which an IRF is known
    :param IRF: irf signal
    :param time_res: time resolution of FLIM cube (seconds). default is 64e-12 seconds (= 64 picoseconds)
    
    :return: lifetime map over ROIs
    """
    time_axis = np.arange(0, 512, 1) * time_res
    IRF_COM = np.average(time_axis, weights=IRF)
    time_vals = np.arange(0, FLIM_cube.shape[2]*time_res, time_res)
    time_cube = np.zeros(FLIM_cube.shape) + time_vals
    time_stretched_cube = np.multiply(FLIM_cube, time_cube)
    time_integral = np.sum(time_stretched_cube, axis=2)
    intensity_integral = np.sum(FLIM_cube, axis=2)
    lifetime_map = np.zeros(roi_map.shape)
    for label in roi_labels:
        roi_COMs = np.divide(time_integral, intensity_integral, where=(roi_map==label))
        lifetime_map[roi_map == label] = np.mean(roi_COMs - IRF_COM)
    return lifetime_map

def phasor_coordinates(decay, freq=20e6, delta_t=19.5e-11, harmonic=1):
    """
    This function calculates the phasor coordinates for a given decay and it's parameters by applying the fourier transform

    :param decay: decay signal for which to get phasor coordinates
    :param freq: frequency of laser in Hz (default is 20MHz) 
    :param delta_t: width of time bins in seconds (default is 19.5e-11 seconds = 195ps)
    :param harmonic: default is 1, no higher harmonics

    :return: phasor coordinates g, s
    """
    w = 2 * np.pi * (freq*harmonic)
    time_bins = (np.linspace(1, len(decay), len(decay)) - 0.5) * delta_t
    scale = np.sum(decay)
    g = np.dot(decay, np.cos(w * time_bins)) / scale
    s = np.dot(decay, np.sin(w * time_bins)) / scale
    return g, s

def phasor_tau(g, s, freq=20e6):
    """
    this function calculates lifetime (tau) from the phasor coordinates

    :param g: phasor coordinate g
    :param s: phasor coordinate s
    :param freq: frequency of laser in Hz (default is 20MHz)

    :return: lifetime tau in seconds
    """
    w = 2 * np.pi *freq
    return (1/w) * (s/g)

def phasor_lifetime(decay, freq=20e6, delta_t=19.5e-11, harmonic=1):
    """
    this function uses phasor analysis to calculate lifetime for a decay
    """
    g, s = phasor_coordinates(decay, freq=freq, delta_t=delta_t, harmonic=harmonic)
    tau = phasor_tau(g, s, freq=freq)
    return tau

def oldROI_phasor_lifetime(segmented_img, truncated_img_cube, freq=20e6, delta_t=19.5e-11, harmonic=1):
    """
    this function creates a lifetime map using mean lifetime of each ROI, calculated from phasor coordinates

    :param segmented_img: segmented image
    :param truncated_img_cube: image cube with only wanted slices kept
    :param freq: frequency of laser in Hz. (default is 20MHz)
    :param delta_t: time bin size in seconds (default is 19.5e-11 second = 195ps)
    """
    segmented = segmented_img.copy()
    segments = np.unique(segmented)
    #re-label segments nicely in increments of 1
    for i, seg in enumerate(segments):
        segmented[np.where(segmented==seg)] = i
    segments = np.unique(segmented)
    tau_map = np.zeros(segmented.shape)
    all_taus = np.zeros(len(segments))
    for i, seg in enumerate(segments[1:]):
        row, col = np.where(segmented == seg)
        taus = np.zeros(row.shape)
        for j in range(len(row)): 
            g, s = phasor_coordinates(truncated_img_cube[row[j], col[j]], freq=freq, delta_t=delta_t, harmonic=harmonic)
            taus[j] = phasor_tau(g, s, freq=freq)
        roi_mean_tau = np.mean(taus)
        all_taus[i] = roi_mean_tau
    for k, lifetime in enumerate(all_taus):
        tau_map[segmented == k+1] = lifetime
    return tau_map

def ROI_phasor_lifetime(segmented, labels, truncated_img_cube, freq=20e6, delta_t=19.5e-11, harmonic=1, list_taus=False):
    """
    this function creates a lifetime map using mean lifetime of each ROI, calculated from phasor coordinates

    :param segmented_img: segmented image
    :param truncated_img_cube: image cube with only wanted slices kept
    :param freq: frequency of laser in Hz. (default is 20MHz)
    :param delta_t: time bin size in seconds (default is 19.5e-11 second = 195ps)
    """
    tau_map = np.zeros(segmented.shape)
    all_taus = np.zeros(len(labels))
    for i, seg in enumerate(labels):
        row, col = np.where(segmented == seg)
        taus = np.zeros(row.shape)
        for j in range(len(row)): 
            g, s = phasor_coordinates(truncated_img_cube[row[j], col[j]], freq=freq, delta_t=delta_t, harmonic=harmonic)
            taus[j] = phasor_tau(g, s, freq=freq)
        all_taus[i] = np.mean(taus)
    for k, tau in enumerate(all_taus):
        tau_map[segmented == k+1] = tau
    if list_taus:
        return tau_map, all_taus
    else: 
        return tau_map

def ROI_phasor_data_maps(segmented_img, labels, truncated_img_cube, freq=20e6, delta_t=19.5e-11, harmonic=1):
    """
    this function creates a lifetime map and a phasor coordinate g map, and a phasor coordinate s map, using mean value of each ROI, calculated from phasor coordinates

    :param segmented_img: segmented image
    :param truncated_img_cube: image cube with only wanted slices kept
    :param freq: frequency of laser in Hz. (default is 20MHz)
    :param delta_t: time bin size in seconds (default is 19.5e-11 second = 195ps)

    :return tau_map: 2d lifetime map (seconds)
    :return g_map: 2d map of mean roi g values
    :return s_map: 2d map of mean roi s values
    """
    tau_map = np.zeros(segmented_img.shape)
    g_map = np.zeros(segmented_img.shape)
    s_map = np.zeros(segmented_img.shape)
    all_taus = np.zeros(len(labels))
    all_Gs = np.zeros(len(labels))
    all_Ss = np.zeros(len(labels))
    for i, seg in enumerate(labels):
        row, col = np.where(segmented_img == seg)
        taus = np.zeros(row.shape)
        Gs = np.zeros(row.shape)
        Ss = np.zeros(row.shape)
        for j in range(len(row)): 
            g, s = phasor_coordinates(truncated_img_cube[row[j], col[j]], freq=freq, delta_t=delta_t, harmonic=harmonic)
            Gs[j] = g 
            Ss[j] = s
            taus[j] = phasor_tau(g, s, freq=freq)
        all_taus[i] = np.mean(taus)
        all_Gs[i] = np.mean(Gs)
        all_Ss[i] = np.mean(Ss)
    for k, tau in enumerate(all_taus):
        tau_map[segmented_img == k+1] = tau
        g_map[segmented_img == k+1] = all_Gs[k]
        s_map[segmented_img == k+1] = all_Ss[k]
    return tau_map, g_map, s_map


#NON STRUCTURED DATA
def INO_stack_nonstructured_data_withFRET(INOfolder_path, csv_name, spectralRange0, spectralRange1, spectralRange2, spectralRange3, algorithm, kernel_size, sigma, noise_level, small_roi_radius, roi_size, log_binary, mean_scalar, donorT0Lims, bin_map_factor=1, autoT0=False, ROI_thresholding=True, no_acceptor_lifetime=3.8e-9):
    """
    this function calculates non-structured for an INO image stack

    :param INOfolder_path:
    :param lifetime_bin_lims (optional): limit image cubes to exponential part of decay for lifetime calculation. (default is 20:110)
    :param ROI_thresholding (optional): Bool. if True, mean thresholding will be applied within ROIs. (default is False)
    :param no_acceptor_lifetime (optional): assumed fluorescence lifetime without presence of acceptor in seconds. (default is 3.8e-9 seconds = 3.8 ns)

    :return: dataframe containing: slice id, ROI labels, ROI sizes, donor image intensity, acceptor image intensity, acceptor/donor ratio, lifetimes, FRET efficiency
    """ #dont forget to update function description
    well_IDs_dict = separate_wellIDs(INOfolder_path)

    well_IDs = np.array([])
    slice_ID = np.array([])
    ROI_labels = np.array([])
    ROI_sizes = np.array([])
    T0_Intensity = np.array([])
    donor_ROI_intensity = np.array([])
    spectralRange0_ROI_intensity = np.array([])
    spectralRange1_ROI_intensity = np.array([])
    spectralRange2_ROI_intensity = np.array([])
    spectralRange3_ROI_intensity = np.array([])
    ROI_lifetime = np.array([])
    ROI_G = np.array([])
    ROI_S = np.array([])
    FRET_efficiency = np.array([])

    for key in well_IDs_dict:
        print(f'\nworking on: {key} stack\n') #THIS IS JUST SO I DONT WORRY TOO MUCH IF THIS TAKES A LONG TIME
        INOfiles = well_IDs_dict[key]
        for i, filename in enumerate(INOfiles):
            print(f'working on slice {i} in stack {key}') #same here
            #load image cubes
            cubes_path = INOfolder_path + '\\' + filename
            try: #this try is here in case of broken/correupted INO files just so that we can skip them and keep working.
                donor_cube, spectral_cube = read_multipage_tiff(cubes_path)
                donor_img = np.sum(donor_cube, axis=2)
                donor_decay = np.sum(np.sum(donor_cube, axis=0), axis=0)
                acceptor_img = np.sum(spectral_cube[:, :, spectralRange0[0]:spectralRange0[1]], axis=2) #spectral range 0
                if spectralRange1 != None:
                    spectralRange1Image = np.sum(spectral_cube[:, :, spectralRange1[0]:spectralRange1[1]], axis=2)
                    spectralRange2Image = np.sum(spectral_cube[:, :, spectralRange2[0]:spectralRange2[1]], axis=2)
                    spectralRange3Image = np.sum(spectral_cube[:, :, spectralRange3[0]:spectralRange3[1]], axis=2)
                if i == 0:
                    laser_freq, FLIM_time_res, IRF, spectral_map, pixel_dwell_time, pixel_len = parse_INO_metadata(cubes_path)
                #SEGMENTATION HERE (of donor image)
                segmented_img, labels = segmentation_choice(donor_img, method=algorithm, kernel_size=kernel_size, sigma=sigma, bin_map_factor=bin_map_factor, noise_level=noise_level, small_roi_radius=small_roi_radius, roi_size=roi_size, log_binary=log_binary, threshold=ROI_thresholding, mean_scalar=mean_scalar)                                                                    
                #well ID array
                well_id = np.full(labels.shape, key)
                #SLICE id ARRAY
                ID = np.zeros(labels.shape) + i
                #ROI SIZE, donor/acceptor intensities (mean intensities)
                sizes = np.zeros(labels.shape) 
                donor_intensity = np.zeros(labels.shape)
                spectralRange0Intensity = np.zeros(labels.shape)
                spectralRange1Intensity = np.zeros(labels.shape)
                spectralRange2Intensity = np.zeros(labels.shape)
                spectralRange3Intensity = np.zeros(labels.shape)
                for j, label in enumerate(labels):
                    sizes[j] = np.count_nonzero(segmented_img == label)
                    donor_intensity[j] = np.mean(donor_img[segmented_img == label])
                    spectralRange0Intensity[j] = np.mean(acceptor_img[segmented_img == label])
                    if spectralRange1 != None:
                        spectralRange1Intensity[j] = np.mean(spectralRange1Image[segmented_img == label])
                        spectralRange2Intensity[j] = np.mean(spectralRange2Image[segmented_img == label])
                        spectralRange3Intensity[j] = np.mean(spectralRange3Image[segmented_img == label])
                #ROI LIFETIMES
                lifetimes = np.zeros(labels.shape) 
                Gs = np.zeros(labels.shape)
                Ss = np.zeros(labels.shape)
                decay_start = np.argmax(donor_decay) 
                decay_end = np.argmin(donor_decay[decay_start:]) + decay_start
                lifetime_map, g_map, s_map = ROI_phasor_data_maps(segmented_img, labels, donor_cube[:, :, decay_start+10:decay_end-10], freq=laser_freq, delta_t=FLIM_time_res) 
                for j, label in enumerate(labels):
                    lifetimes[j] = lifetime_map[segmented_img == label][0] #each ROI is already filled with mean lifetimes, this would otherwise make a list of length = the size of the ROI
                    Gs[j] = g_map[segmented_img == label][0]
                    Ss[j] = s_map[segmented_img == label][0]
                #T0 STUFF
                if autoT0:
                    T0_val = np.mean(donor_decay[decay_start - donorT0Lims[0] : decay_start + donorT0Lims[0]+1])
                else:
                    T0_val = np.mean(donor_decay[donorT0Lims[0]:donorT0Lims[1]+1])
                T0_array = np.full(labels.shape, T0_val)
                #FRET EFFICIENCY
                efficiency = (1 - lifetimes/no_acceptor_lifetime) * 100                 
                #ADD ALL TO THE TOTAL RUNNING ARRAYS
                well_IDs = np.append(well_IDs, well_id)
                slice_ID = np.append(slice_ID, ID)
                ROI_labels = np.append(ROI_labels, labels)
                ROI_sizes = np.append(ROI_sizes, sizes)
                T0_Intensity = np.append(T0_Intensity, T0_array)
                donor_ROI_intensity = np.append(donor_ROI_intensity, donor_intensity)
                spectralRange0_ROI_intensity = np.append(spectralRange0_ROI_intensity, spectralRange0Intensity)
                spectralRange1_ROI_intensity = np.append(spectralRange1_ROI_intensity, spectralRange1Intensity)
                spectralRange2_ROI_intensity = np.append(spectralRange2_ROI_intensity, spectralRange2Intensity)
                spectralRange3_ROI_intensity = np.append(spectralRange3_ROI_intensity, spectralRange3Intensity)
                ROI_lifetime = np.append(ROI_lifetime, lifetimes)
                ROI_G = np.append(ROI_G, Gs)
                ROI_S = np.append(ROI_S, Ss)
                FRET_efficiency = np.append(FRET_efficiency, efficiency)
            except Exception as e:
                print(f'FILE error: \n{e}')
                pass
    df = pd.DataFrame({'Well_ID':well_IDs,
        'Slice_ID':slice_ID,
        'ROI_ID':ROI_labels,
        'ROI_size':ROI_sizes,
        'T0_intensity':T0_Intensity,
        'Donor_intensity':donor_ROI_intensity,
        'Acceptor_intensity':spectralRange0_ROI_intensity,
        'Spectral_Range1_intensity':spectralRange1_ROI_intensity,
        'Spectral_Range2_intensity':spectralRange2_ROI_intensity,
        'Spectral_Range3_intensity':spectralRange3_ROI_intensity,
        'Phasor_G':ROI_G,
        'Phasor_S':ROI_S,
        'Donor_lifetime':ROI_lifetime,
        'FRET_efficiency':FRET_efficiency})
    df.dropna() #get rid of all rows containing NaN due to low brightness
    df.to_csv(csv_name)
    print('Process Finished.')
    return df

def ISS_stack_nonstructured_data_withFRET(stack_path, image_stack, csv_name, algorithm, kernel_size, sigma, noise_level, small_roi_radius, roi_size, log_binary, mean_scalar, donorT0Lims, acceptorT0Lims, bin_map_factor=1, autoT0=False, ROI_thresholding=False, no_acceptor_lifetime=3.8e-9):
    """
    this function calculates non-structured for an ISS image stack

    :param stack_path: path to image stack, image stack name if in working directory
    :param lifetime_bin_lims (optional): limit image cubes to exponential part of decay for lifetime calculation. (default is 20:110)
    :param ROI_thresholding (optional): Bool. if True, mean thresholding will be applied within ROIs. (default is False)
    :param no_acceptor_lifetime (optional): assumed fluorescence lifetime without presence of acceptor in seconds. (default is 3.8e-9 seconds = 3.8 ns)

    :return: dataframe containing: slice id, ROI labels, ROI sizes, donor image intensity, acceptor image intensity, acceptor/donor ratio, lifetimes, FRET efficiency
    """ #dont forget to update function description
    metadata = get_ISS_metadata(stack_path)
    laser_freq, FLIM_time_res, pixel_dwell_time, pixel_len = parse_ISS_metadata(metadata)
    #image is directly passed to function in an effort to conserve memory

    slice_ID = np.array([])
    ROI_labels = np.array([])
    ROI_sizes = np.array([])
    donor_T0 = np.array([])
    acceptor_T0 = np.array([])
    donor_ROI_intensity = np.array([])
    acceptor_ROI_intensity = np.array([])
    ROI_lifetime = np.array([])
    ROI_G = np.array([])
    ROI_S = np.array([])
    FRET_efficiency = np.array([])

    for i, slice in enumerate(image_stack):
        print(f'Working on slice {i}...')
        acceptor_cube = slice[0, :, :, :]
        acceptor_img = np.sum(acceptor_cube, axis=2)
        donor_cube = slice[1, :, :, :]
        donor_decay = np.sum(np.sum(donor_cube, axis=0), axis=0)
        acceptor_decay = np.sum(np.sum(acceptor_cube, axis=0), axis=0)
        donor_img = np.sum(donor_cube, axis=2)
        #SEGMENTATION HERE (of donor image)
        segmented_img, labels = segmentation_choice(donor_img, method=algorithm, kernel_size=kernel_size, sigma=sigma, bin_map_factor=bin_map_factor, noise_level=noise_level, small_roi_radius=small_roi_radius, roi_size=roi_size, log_binary=log_binary, threshold=ROI_thresholding, mean_scalar=mean_scalar)
        #SLICE id ARRAY
        ID = np.zeros(labels.shape) + i
        #ROI SIZE, donor/acceptor intensities
        sizes = np.zeros(labels.shape) 
        donor_intensity = np.zeros(labels.shape)
        acceptor_intensity = np.zeros(labels.shape)
        for j, label in enumerate(labels):
            sizes[j] = np.count_nonzero(segmented_img == label)
            donor_intensity[j] = np.mean(donor_img[segmented_img == label])
            acceptor_intensity[j] = np.mean(acceptor_img[segmented_img == label])
        #ROI LIFETIMES
        lifetimes = np.zeros(labels.shape) 
        Gs = np.zeros(labels.shape)
        Ss = np.zeros(labels.shape)
        decay_start = np.argmax(donor_decay) 
        decay_end = np.argmin(donor_decay[decay_start:]) + decay_start
        lifetime_map, g_map, s_map = ROI_phasor_data_maps(segmented_img, labels, donor_cube[:, :, decay_start+10:decay_end-10], freq=laser_freq, delta_t=FLIM_time_res) 
        for j, label in enumerate(labels):
            lifetimes[j] = lifetime_map[segmented_img == label][0] #each ROI is already filled with mean lifetimes, this would otherwise make a list of length = the size of the ROI
            Gs[j] = g_map[segmented_img == label][0]
            Ss[j] = s_map[segmented_img == label][0]
        #T0 STUFF
        acceptor_peak = np.argmax(acceptor_decay)
        if autoT0:
            donorT0_val = np.mean(donor_decay[decay_start - donorT0Lims[0] : decay_start + donorT0Lims[0]+1])
            acceptorT0_val = np.mean(acceptor_decay[acceptor_peak - donorT0Lims[0] : acceptor_peak + donorT0Lims[0]+1])
        else:
            donorT0_val = np.mean(donor_decay[donorT0Lims[0]:donorT0Lims[1]+1])
            acceptorT0_val = np.mean(acceptor_decay[acceptorT0Lims[0]:acceptorT0Lims[1]+1])
        donorT0_array = np.full(labels.shape, donorT0_val)
        acceptorT0_array = np.full(labels.shape, acceptorT0_val)
        #FRET EFFICIENCY
        efficiency = (1 - lifetimes/no_acceptor_lifetime) * 100
        #ADD ALL TO THE TOTAL RUNNING ARRAYS
        slice_ID = np.append(slice_ID, ID)
        ROI_labels = np.append(ROI_labels, labels)
        ROI_sizes = np.append(ROI_sizes, sizes)
        donor_T0 = np.append(donor_T0, donorT0_array)
        acceptor_T0 = np.append(acceptor_T0, acceptorT0_array)
        donor_ROI_intensity = np.append(donor_ROI_intensity, donor_intensity)
        acceptor_ROI_intensity = np.append(acceptor_ROI_intensity, acceptor_intensity)
        ROI_lifetime = np.append(ROI_lifetime, lifetimes)
        ROI_G = np.append(ROI_G, Gs)
        ROI_S = np.append(ROI_S, Ss)
        FRET_efficiency = np.append(FRET_efficiency, efficiency)
    df = pd.DataFrame({'Slice_ID':slice_ID,
    'ROI_ID':ROI_labels,
    'ROI_size':ROI_sizes,
    'Donor_T0_intensity':donor_T0,
    'Acceptor_T0_intensity':acceptor_T0,
    'Donor_intensity':donor_ROI_intensity,
    'Acceptor_intensity':acceptor_ROI_intensity,
    'Donor_lifetime':ROI_lifetime,
    'Phasor_G':ROI_G,
    'Phasor_S':ROI_S,
    'FRET_efficiency':FRET_efficiency})
    df.dropna
    df.to_csv(csv_name)
    print('Process Finished.')
    return df


#OTHER
def fft_convolution(old_array, kernel):
    """
    this function applies nearest neighbor binning to a 2d array 
    this time using the np.signal.fftconvolve 
    this function is preferred when using larger kernel sizes

    :param old_array: the array to resample with nearest neighbour convolution
    :output new_array: numpy array with same dimensions as old_array, but resampled values
    """
    new_array = signal.fftconvolve(old_array, kernel, mode='same')
    return new_array


def make_exp_mask(reference_image_cube, second_deriv_threshold=0.05, min_interval=10):
    """
    this function uses the second derivative of the log of a decay curve to determine where
    it is close to being a nice exponential curve and creates a binary mask for the exponential
    part of the decay curve.

    :param reference_image_cube: image cube for reference sample (ex: fluorescein)
    :param second_diff_threshold: threshold for second derivative (mask will only include the area with secodn derivative near 0)
    :param min_interval: minimum length of interval between the beginning and end of the non-zero portion of the mask
    
    :return: 1d binary mask array
    """
    ref_decay = np.sum(np.sum(reference_image_cube, axis=0), axis=0)
    log_decay = np.log(ref_decay)
    exp_mask = np.zeros_like(ref_decay)
    log_decay = log_decay[np.isfinite(log_decay)] #drop out all the NAN, infs from the array
    start_deriv = np.argmax(log_decay) #assume that good decay curve appears only after greatest decay signal value
    second_deriv = np.gradient(np.gradient(log_decay[start_deriv:])) #apply np.gradient twice (equivalent to derivative for 1darray)
    abs_second_deriv = np.abs(second_deriv)

    #this function makes sure that the interval isn't too short and adjusts the starting point forward as needed to avoid non-smoothness
    #near the beginning of the decay signal.
    def find_interval(abs_deriv, interval_len, threshold, start_point):
        if (len(abs_deriv) - start_point) < interval_len:
            print('\nERROR: No suitable interval found for given interval length and second derivative threshold.\n')
            return None, None
        #start the mask Nat the first element where the second derivative is below the threshold
        start_index = np.where(abs_deriv[start_point:] < threshold)[0][0] + start_point
        #cut off the mask before the first element where the second deriv exceeds threshold value
        end_index = np.argmax(abs_deriv[start_index:] > threshold) + start_index
        if end_index - start_index >= interval_len:
            return start_index, end_index
        else:
            return find_interval(abs_deriv, interval_len, threshold, start_point+1)
    
    start_index, end_index = find_interval(abs_second_deriv, min_interval, second_deriv_threshold, start_point=0)
    exp_mask[start_index+start_deriv:end_index] = 1 #set mask to 1 over the desired interval.
    return exp_mask


def restrict_decay(decay, mask=[], custom_start=None, custom_end=None):
    """
    this function takes a decay signal and applies user defined bounds to it to isolate the desired part of the signal
    for phasor analysis
    :param decay: 1darray of decay signal
    :param mask: binary mask of same shape as decay (should be 1 over desired decay interval.) default=[] -> no mask
    :param custom_start: user defined index at which to start the desired part of the signal if no mask is provided
    :param custom_end: user defined index at which to end the desired part of the signal if no mask is provided.
    
    :return: 1darray containing the decay signal truncated using either the mask or the user defined start, end parameters.
    """
    if len(mask) != 0:
        if mask.shape == decay.shape:
            kept_decay = np.multiply(decay, mask)
            kept_decay = kept_decay[kept_decay != 0]
            return kept_decay
        print('\nERROR: mask shape does not match decay shape\n')
        return None
    kept_decay = decay[custom_start:custom_end]
    return kept_decay


def get_nonstructured_data(acceptor_cube, donor_cube, lifetime_bin_lims=(20,110), AD_noise_level=20, ROI_thresholding=False):
    """
    THIS STILL NEEDS SOME WORK TO MAKE IT AS NICE AS POSSIBLE
    this function produces non-structured data from a acceptor and donor image cubes

    :param acceptor_cube: acceptor image cube
    :param donor_cube: donor image cube
    :param lifetime_bin_lims: limits on image cube, intended to only use the exponential part of the decay curfve for lifetime calculation
    :param AD_noise_level: noise level threshold for A/D ratio image. (default is 20)
    :param ROI_thresholding: Bool. If true, mean thresholding will be applied within ROIs

    :return: ROI labels, ROI sizes, ROI lifetimes, ROI acceptor/donor ratios
    """
    acceptor_img = np.sum(acceptor_cube, axis=2)
    donor_img = np.sum(donor_cube, axis=2)

    #SEGMENTATION HERE (of donor image)
    binary_map = get_binary_map(donor_img)
    local_maxima = get_local_maxima(donor_img)
    segmented_img = multiseed_watershed_segmentation(binary_map, local_maxima)
    #rename ROIS
    segments = np.unique(segmented_img)
    for i, seg in enumerate(segments[1:]): #skip the 0 label for background
        segmented_img[segmented_img == seg] = i+1
    ROI_labels = np.unique(segmented_img)[1:] #skip the 0 label for background

    #THRESHOLDING (optional)
    if ROI_thresholding:
        segmented_img = threshold_ROIs(segmented_img, donor_img)

    #ROI SIZE
    ROI_sizes = np.zeros(ROI_labels.shape) 
    for i, label in enumerate(ROI_labels):
        ROI_sizes[i] = np.count_nonzero(segmented_img == label)
    #A/D RATIOS
    AD_ratio_img = np.divide(np.float32(acceptor_img), np.float32(donor_img), where=(donor_img != 0))
    AD_ratio_img[donor_img < AD_noise_level] = 0
    ROI_ratios = np.zeros(ROI_labels.shape)
    for i, label in enumerate(ROI_labels):
        ROI_ratios[i] = np.mean(AD_ratio_img[segmented_img == label])
    #ROI LIFETIMES
    lifetime_map = ROI_phasor_lifetime(segmented_img, donor_cube[:, :, lifetime_bin_lims[0]:lifetime_bin_lims[1]]) #make this use function parameters instead
    ROI_lifetimes = np.zeros(ROI_labels.shape) 
    for i, label in enumerate(ROI_labels):
        ROI_lifetimes[i] = np.mean(lifetime_map[segmented_img == label])
    
    return ROI_labels, ROI_sizes, ROI_lifetimes, ROI_ratios


def data_binning(xdata, ydata, auto_bins=True, bin_edges=[]):
    """
    this function will bin the fret efficiency values in order to make a nice binding curve plot
    
    :param ratio: array of data for x axis
    :param FRET: array of data for y axis
    :parm auto_bins: bool. if true, bins are set automatically (with more bins earlier in the dataset). if false, bins are defined by bin_edges.
    :param bin_edges: list of bin edges

    :return: bin sizes, mean y per bin, median y per bin, y std deviation per bin, y std error per bin, bin centers, mean x per bin
    """
    x_array = np.array(xdata)
    y_array = np.array(ydata)
    y_bins = []
    x_bins = []
    x_bin_centers = []
    bin_sizes = []
    i = 0         #i should be bin center.
    lim = np.max(x_array)
    if auto_bins:
        while i <= lim:
            if i < lim/4:
                bin_size=lim*0.015
            elif i < lim/2:
                bin_size= lim*0.05
            else:
                bin_size=lim*0.1
            y_vals = y_array[np.logical_and(x_array >= i-bin_size, x_array < i+bin_size)]
            x_vals = x_array[np.logical_and(x_array >= i-bin_size, x_array < i+bin_size)]
            y_bins.append(y_vals)
            x_bins.append(x_vals)
            x_bin_centers.append(i)
            bin_sizes.append(len(y_vals))
            i += 2*bin_size
    elif bin_edges != []:
        for i in range(len(bin_edges) -1):
            center = (bin_edges[i] + bin_edges[i+1])/2
            y_vals = y_array[np.logical_and(x_array >= bin_edges[i], x_array < bin_edges[i+1])]
            x_vals = x_array[np.logical_and(x_array >= bin_edges[i], x_array < bin_edges[i+1])]
            y_bins.append(y_vals)
            x_bins.append(x_vals)
            bin_sizes.append(len(y_vals))
            x_bin_centers.append(center)
    else:
        print('No bins selected')
        return None
    y_means = np.array([])
    y_medians = np.array([])
    y_std_devs = np.array([])
    y_std_errs = np.array([])
    x_means = np.array([])
    for i, bin in enumerate(y_bins):
        y_means = np.append(y_means, np.mean(bin))
        y_medians = np.append(y_medians, np.median(bin))
        y_std_devs = np.append(y_std_devs, np.std(bin))
        y_std_errs = np.append(y_std_errs, ( np.std(bin) / np.sqrt( len(bin) ) ) )
        x_means = np.append(x_means, np.mean(x_bins[i]))
    return bin_sizes, y_means, y_medians, y_std_devs, y_std_errs, x_bin_centers, x_means


def random_cmap(map_len=256, black_background=False):
    """
    this function creates a random color map, useful in segmentation maps

    :param map_len: optional. length of color map. default is 256
    
    :return: random color map.
    """
    from matplotlib import colors
    temp_cmap = np.random.rand(map_len, 3)
    if black_background:
        temp_cmap[0] = 0
    return colors.ListedColormap(temp_cmap)
    

def black_0_cmap(length=256, cmap='plasma'):
    """
    this function takes a given matplotlib colormap and sets the 0 value to black (RGBa)

    :param cmap: str. name of colormap to edit
    
    :return: listed colormap with 0 value set to black
    """
    from matplotlib import cm
    from matplotlib import colors
    oldcmap = cm.get_cmap(cmap)
    newcolors = oldcmap(np.linspace(0, 1, length))
    newcolors[0] = [0, 0, 0, 1]
    return colors.ListedColormap(newcolors)


def overlay_roi_map(donor_img, roi_map, labels):
    """
    this function creates a png of a randomly colored roi_map, overlayed on top of a FLIM donor intensity image.
    
    :param donor_img: intensity image of donor FLIM
    :param roi_map: ROI map for donor image

    :return: creates a png named 'roi-overlay.png'
    """
    rand_cmap = random_cmap(map_len=len(labels))
    segmented_mask = np.ma.masked_where(roi_map == 0, roi_map)
    plt.figure(figsize=(20,20))
    plt.imshow(donor_img, 'gray')
    plt.imshow(segmented_mask, cmap=rand_cmap, interpolation='none')
    plt.xticks([])
    plt.yticks([])
    plt.savefig('roi-overlay.png', bbox_inches='tight')


def closest_index(value, array):
    """
    this function returns the index of the element in an array that is nearest to the specified value.

    :param value: float. value to look for in array
    :param array: numpy array to look for value in

    :return: int. index of element closest to value.
    """
    index = (np.abs(array - value)).argmin()
    return index


def filter_df(df, **filterDict):
    """
    this function will filter a dataframe according to the parameters given by the filterDict **kwarg

    :param df: dataframe to filter
    :param **filterDict: dictionary, filtering parameters to use

    :return: filtered sub-dataframe
    """
    subdf = df.copy()
    cols = list(subdf.columns)
    for key in filterDict:
        value = filterDict[key]
        if value != None:
            if key in cols:
                op = value[:2]
                if key == 'Well_ID':
                    #check well IDs
                    if op == 'ee':
                        limit = value[2:]
                        subdf = subdf[subdf[key] == limit]
                elif key == 'Acceptor/Donor_ratio':
                    limit = float(value[2:])  #maybe want to add something for conversion errors, but error message may be enough
                    if 'Well_ID' in cols:
                        if op == 'gt':
                            subdf = subdf[subdf['Spectral_Range0_intensity']/subdf['Donor_intensity'] > limit]
                        elif op == 'lt':
                            subdf = subdf[subdf['Spectral_Range0_intensity']/subdf['Donor_intensity'] < limit]
                        elif op == 'ge':
                            subdf = subdf[subdf['Spectral_Range0_intensity']/subdf['Donor_intensity'] >= limit]
                        elif op == 'le':
                            subdf = subdf[subdf['Spectral_Range0_intensity']/subdf['Donor_intensity'] <= limit]
                        elif op == 'ee':
                            subdf = subdf[subdf['Spectral_Range0_intensity']/subdf['Donor_intensity'] == limit]
                    else:
                        if op == 'gt':
                            subdf = subdf[subdf['Acceptor_intensity']/subdf['Donor_intensity'] > limit]
                        elif op == 'lt':
                            subdf = subdf[subdf['Acceptor_intensity']/subdf['Donor_intensity'] < limit]
                        elif op == 'ge':
                            subdf = subdf[subdf['Acceptor_intensity']/subdf['Donor_intensity'] >= limit]
                        elif op == 'le':
                            subdf = subdf[subdf['Acceptor_intensity']/subdf['Donor_intensity'] <= limit]
                        elif op == 'ee':
                            subdf = subdf[subdf['Acceptor_intensity']/subdf['Donor_intensity'] == limit]
                else:
                    limit = float(value[2:])  #maybe want to add something for conversion errors, but error message may be enough
                    if op == 'gt':
                        subdf = subdf[subdf[key] > limit]
                    elif op == 'lt':
                        subdf = subdf[subdf[key] < limit]
                    elif op == 'ge':
                        subdf = subdf[subdf[key] >= limit]
                    elif op == 'le':
                        subdf = subdf[subdf[key] <= limit]
                    elif op == 'ee':
                        subdf = subdf[subdf[key] == limit]
    return subdf








