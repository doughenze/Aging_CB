import Mapping
import argparse
import os
import string

import cv2
import mahotas
import geopandas as gpd
import igraph as ig
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import rasterio
import seaborn as sns
import tifffile
import umap
from anndata import AnnData as ad
from matplotlib import patches as mpatches
import scanpy as sc
import anndata
from shapely.affinity import translate
from shapely.geometry import Polygon, MultiPolygon, box, shape
from skimage import img_as_bool, img_as_ubyte
from skimage.measure import find_contours, regionprops, regionprops_table, inertia_tensor, inertia_tensor_eigvals
from skimage.feature import graycomatrix, graycoprops, local_binary_pattern, canny
from skimage.transform import radon
from skimage.filters import gabor
from skimage.morphology import skeletonize, opening, closing, disk
from skimage.segmentation import find_boundaries
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import geojson
import json

import torch
from transformers import AutoFeatureExtractor, AutoImageProcessor, AutoModel

import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure, morphology, filters
from skimage.morphology import skeletonize, convex_hull_image
from scipy.ndimage import label, distance_transform_edt,binary_fill_holes
from scipy.spatial.distance import cdist
from skimage.draw import circle_perimeter, ellipse_perimeter
from skimage.util import img_as_ubyte
from scipy.stats import linregress

import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure, morphology, filters
from skimage.morphology import skeletonize, closing, disk
from scipy.ndimage import distance_transform_edt, label

import pywt
from scipy import fftpack
import gudhi as gd

def find_junctions_and_endpoints(skeleton):
    junctions = np.zeros_like(skeleton)
    endpoints = np.zeros_like(skeleton)
    
    for i in range(1, skeleton.shape[0] - 1):
        for j in range(1, skeleton.shape[1] - 1):
            if skeleton[i, j]:
                neighborhood = skeleton[i-1:i+2, j-1:j+2]
                count = np.sum(neighborhood) - 1
                if count == 1:
                    endpoints[i, j] = 1
                elif count > 2:
                    junctions[i, j] = 1

    return junctions, endpoints

def calculate_total_length(skeleton):
    labeled_skeleton, num_features = label(skeleton)
    regions = measure.regionprops(labeled_skeleton)
    total_length = 0
    branch_lengths = []

    for region in regions:
        coords = region.coords
        length = 0
        for i in range(len(coords) - 1):
            length += np.linalg.norm(coords[i] - coords[i + 1])
        total_length += length
        branch_lengths.append(length)

    mean_branch_length = np.mean(branch_lengths) if branch_lengths else 0

    return total_length, mean_branch_length, len(branch_lengths)

def count_branches_and_points(junctions, endpoints):
    n_branches = np.sum(junctions)
    n_tips = np.sum(endpoints)
    
    return n_branches, n_tips

def calculate_convex_hull_properties(binary_image):
    convex_hull = convex_hull_image(binary_image)
    convex_hull_props = measure.regionprops(measure.label(convex_hull))[0]
    convex_hull_area = convex_hull_props.area
    convex_hull_perimeter = convex_hull_props.perimeter
    convex_hull_bbox = convex_hull_props.bbox
    convex_hull_span_ratio = (convex_hull_bbox[2] - convex_hull_bbox[0]) / (convex_hull_bbox[3] - convex_hull_bbox[1])
    convex_hull_circularity = (4 * np.pi * convex_hull_area) / (convex_hull_perimeter ** 2)
    return convex_hull_area, convex_hull_perimeter, convex_hull_span_ratio, convex_hull_circularity

def calculate_fractal_dimension(binary_image):
    def boxcount(Z, k):
        S = np.add.reduceat(
            np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
            np.arange(0, Z.shape[1], k), axis=1)
        return len(np.where((S > 0) & (S < k*k))[0])

    Z = binary_image > 0
    p = min(Z.shape)
    n = 2**np.floor(np.log2(p))
    n = int(np.log2(n))
    sizes = 2**np.arange(n, 1, -1)

    counts = []
    for size in sizes:
        counts.append(boxcount(Z, size))

    coeffs = np.polyfit(np.log(sizes), np.log(counts), 1)
    return -coeffs[0]

def calculate_euclidean_and_path_distance(skeleton):
    endpoints = np.column_stack(np.where(skeleton))
    if len(endpoints) > 1:
        euclidean_distance = np.linalg.norm(endpoints[0] - endpoints[-1])
        path_distance = np.sum([np.linalg.norm(endpoints[i] - endpoints[i - 1]) for i in range(1, len(endpoints))])
    else:
        euclidean_distance = 0
        path_distance = 0
    return euclidean_distance, path_distance

def calculate_lacunarity(binary_image, box_size=10):
    image = binary_image > 0
    H, W = image.shape
    lacunarities = []

    for i in range(0, H, box_size):
        for j in range(0, W, box_size):
            box = image[i:i+box_size, j:j+box_size]
            if box.size > 0:
                p = np.sum(box) / box.size
                if p > 0 and p < 1:
                    lacunarities.append((p * (1 - p)) / (p**2))
    
    return np.mean(lacunarities)

def calculate_tortuosity(euclidean_distance, path_distance):
    if euclidean_distance == 0:
        return 0
    return path_distance / euclidean_distance

def perform_sholl_analysis(skeleton, soma_coords, max_radius=500, step=20):
    intersections = []
    radii = np.arange(0, max_radius, step)
    for radius in radii:
        circle = np.zeros_like(skeleton, dtype=bool)
        rr, cc = circle_perimeter(soma_coords[0], soma_coords[1], radius, shape=skeleton.shape)
        circle[rr, cc] = True
        intersection_count = np.sum(skeleton & circle)
        intersections.append(intersection_count)
    return radii, intersections

def calculate_sholl_parameters(radii, intersections):
    # Branching Index: Maximum number of intersections
    branching_index = np.max([np.max(intersections[i:] - intersections[i]) for i in range(len(intersections))])
    
    # Critical Radius: Radius at which the number of intersections is maximized
    critical_radius = radii[np.argmax(intersections)]
    
    # Dendritic Maximum: Maximum number of intersections at any radius
    dendritic_maximum = 2*radii[np.argmax(intersections)]
    
    # Sholl Regression Coefficient (k) and Intercept (c)
    log_radii = np.log(radii[1:])  # Skip the first radius to avoid log(0)
    log_intersections = np.log(intersections[1:])  # Corresponding intersections
    regression_result = linregress(log_radii, log_intersections)
    sholl_regression_coefficient = regression_result.slope
    sholl_regression_intercept = regression_result.intercept
    
    # Ramification Index: Ratio of the number of intersections at the critical radius to the number of intersections at the smallest radius
    first_non_zero = next((x for x in intersections if x != 0), np.nan)
    ramification_index = intersections[np.argmax(intersections)] / first_non_zero
    
    # Radius of Influence (ROI): Radius beyond which the number of intersections does not exceed 5% of the maximum number of intersections
    roi_index = np.argmax(np.array(intersections) <= 0.05 * branching_index)
    radius_of_influence = radii[roi_index] if roi_index != 0 else max(radii)
    
    return {
        'Branching Index': branching_index,
        'Critical Radius': critical_radius,
        'Dendritic Maximum': dendritic_maximum,
        'Sholl Regression Coefficient (k)': sholl_regression_coefficient,
        'Sholl Regression Intercept (c)': sholl_regression_intercept,
        'Ramification Index': ramification_index,
        'Radius of Influence': radius_of_influence
    }

def find_soma_by_distance_transform(binary_image):
    # Compute the distance transform of the binary image
    distance_transform = distance_transform_edt(binary_image)
    
    # Skeletonize the binary image
    skeleton = skeletonize(binary_image)
    
    # Find the coordinates of the skeleton points
    skeleton_coords = np.column_stack(np.where(skeleton))
    
    # Find the point on the skeleton with the maximum distance in the distance transform
    max_distance_point = skeleton_coords[np.argmax(distance_transform[skeleton])]
    
    return max_distance_point, distance_transform, skeleton

def segment_soma(binary_image, soma_coords, distance_transform, threshold_ratio=0.5):
    # Threshold the distance transform to segment the soma
    threshold_value = distance_transform[soma_coords[0], soma_coords[1]] * threshold_ratio
    soma_segment = distance_transform >= threshold_value
    
    # Label the segmented regions
    labeled_soma, _ = label(soma_segment)
    
    # Extract the region containing the soma coordinates
    soma_label = labeled_soma[soma_coords[0], soma_coords[1]]
    soma_region = labeled_soma == soma_label
    
    return soma_region

def calculate_soma_properties(soma_region):
    # Measure properties of the soma region
    properties = measure.regionprops(soma_region.astype(int))
    
    if properties:
        prop = properties[0]
        area = prop.area
        perimeter = prop.perimeter
        circularity = (4 * np.pi * area) / (perimeter ** 2) if perimeter > 0 else 0
    else:
        area = 0
        perimeter = 0
        circularity = 0

    return area, perimeter, circularity

def calculate_morphological_moments(binary_image):
    features = {}
    
    # Spatial Moments, Central Moments, and Normalized Moments
    moments = cv2.moments(binary_image.astype(np.uint8))
    
    for key, value in moments.items():
        features[f'spatial_moment_{key}'] = value
        
    for i in range(4):
        for j in range(4):
            features[f'central_moment_mu{i}{j}'] = moments.get(f'mu{i}{j}', 0)
            features[f'normalized_moment_nu{i}{j}'] = moments.get(f'nu{i}{j}', 0)
    
    # Inertia Tensors and Eigenvalues
    inertia_tensor_ = inertia_tensor(binary_image)
    eigenvalues = inertia_tensor_eigvals(binary_image)
    
    for i in range(inertia_tensor_.shape[0]):
        for j in range(inertia_tensor_.shape[1]):
            features[f'inertia_tensor_{i}_{j}'] = inertia_tensor_[i, j]
    
    for i, eigenvalue in enumerate(eigenvalues):
        features[f'inertia_tensor_eigenvalue_{i}'] = eigenvalue
    
    # Hu Moments
    hu_moments = cv2.HuMoments(cv2.moments(binary_image.astype(np.uint8))).flatten()
    
    for i, hu_moment in enumerate(hu_moments):
        features[f'hu_moment_{i+1}'] = hu_moment
    
    # Zernike Moments
    radius = min(binary_image.shape) // 2
    zernike_moments = mahotas.features.zernike_moments(binary_image, radius)
    
    for i, zernike_moment in enumerate(zernike_moments):
        features[f'zernike_moment_{i+1}'] = zernike_moment
    
    return features

def calculate_shape_characteristics(raw_image, binary_image):
    features = {}
    
    # Use binary_image for contour-based features
    contours, _ = cv2.findContours(binary_image.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    contour = contours[0].squeeze()
    
    # Fourier Descriptors
    fourier_result = np.fft.fft(contour)
    for i in range(10):  # First 10 coefficients
        features[f'fourier_{i+1}'] = fourier_result[i].real
    
    # Texture Descriptors - GLCM
    glcm = graycomatrix(raw_image.astype(np.uint8), [1], [0], 256, symmetric=True, normed=True)
    features['glcm_contrast'] = graycoprops(glcm, 'contrast')[0, 0]
    features['glcm_correlation'] = graycoprops(glcm, 'correlation')[0, 0]
    features['glcm_energy'] = graycoprops(glcm, 'energy')[0, 0]
    features['glcm_homogeneity'] = graycoprops(glcm, 'homogeneity')[0, 0]
    
    # Texture Descriptors - LBP
    lbp = local_binary_pattern(raw_image, P=8, R=1, method='uniform')
    for i, value in enumerate(lbp.ravel()):
        features[f'lbp_{i+1}'] = value
        if i == 50:
            break
    
    # Gabor Filters
    gabor_result, _ = gabor(raw_image, frequency=0.6)
    for i, value in enumerate(gabor_result.ravel()):
        features[f'gabor_{i+1}'] = value
        if i == 50:
            break
    
    # Elliptic Fourier Descriptors (EFD)
    def elliptic_fourier_descriptors(contour, order=10):
        #contour = contour[:, 0]
        dxy = np.diff(contour, axis=0)
        dt = np.sqrt((dxy**2).sum(axis=1))
        t = np.concatenate(([0], np.cumsum(dt)))
        T = t[-1]
        phi = 2 * np.pi * np.outer(t, np.arange(1, order+1)) / T
        #print(contour.shape)
        a = 2 / T * np.dot(np.cos(phi).T, contour[:, 0])
        b = 2 / T * np.dot(np.sin(phi).T, contour[:, 0])
        c = 2 / T * np.dot(np.cos(phi).T, contour[:, 1])
        d = 2 / T * np.dot(np.sin(phi).T, contour[:, 1])
        return np.concatenate([a, b, c, d]).ravel()
    
    efd_result = elliptic_fourier_descriptors(contour, order=10)
    for i, value in enumerate(efd_result):
        features[f'efd_{i+1}'] = value
        if i == 50:
            break
    # Curvature Scale Space (approximated)
    def curvature_scale_space(contour):
        t = np.arange(contour.shape[0])
        kappa = np.abs(np.gradient(np.gradient(contour[:, 0], t), t) * np.gradient(contour[:, 1], t) - np.gradient(np.gradient(contour[:, 1], t), t) * np.gradient(contour[:, 0], t)) / ((np.gradient(contour[:, 0], t)**2 + np.gradient(contour[:, 1], t)**2)**1.5)
        return kappa
    
    css_result = curvature_scale_space(contour)
    for i, value in enumerate(css_result):
        features[f'curvature_{i+1}'] = value
        if i == 50:
            break
    # Topological Descriptors - Betti Numbers & Persistent Homology
    skeleton = skeletonize(binary_image)
    features['betti_number'] = np.sum(skeleton)  # Rough approximation
    
    cubical_complex = gd.CubicalComplex(dimensions=skeleton.shape, top_dimensional_cells=skeleton.flatten())
    persistence = cubical_complex.persistence()
    for i, (dim, pers) in enumerate(persistence):
        features[f'persistent_homology_{i+1}_dim'] = dim
        features[f'persistent_homology_{i+1}_birth'] = pers[0]
        features[f'persistent_homology_{i+1}_death'] = pers[1]
        if i == 50:
                break
    
    # Medial Axis Transform
    medial_axis = morphology.medial_axis(binary_image)
    features['medial_axis'] = medial_axis.sum()
    
    # SIFT
    sift = cv2.SIFT_create()
    kp_sift, des_sift = sift.detectAndCompute(raw_image.astype(np.uint8), None)
    if des_sift is not None:
        for i, value in enumerate(des_sift.flatten()):
            features[f'sift_{i+1}'] = value
            if i == 50:
                break
    
    # Wavelet Transform
    coeffs = pywt.wavedec2(raw_image.astype(np.uint8), 'haar', level=2)
    wavelet_coeffs = []
    for coeff in coeffs:
        if isinstance(coeff, tuple):  # Detail coefficients (tuple of arrays)
            for array in coeff:
                wavelet_coeffs.extend(array.flatten())
        else:  # Approximation coefficients (single array)
            wavelet_coeffs.extend(coeff.flatten())

    # Add the flattened coefficients to the features dictionary
    for i, value in enumerate(wavelet_coeffs):
        features[f'wavelet_{i+1}'] = value
        if i == 50:
            break
    # Radon Transform
    theta = np.linspace(0., 180., max(raw_image.shape), endpoint=False)
    sinogram = radon(raw_image, theta=theta)
    for i, value in enumerate(sinogram.flatten()):
        features[f'radon_{i+1}'] = value
        if i == 50:
            break
    return features

def calculate_morphometric_features(cropped_binary_image, regionprop, cropped_iba1_image, dapi_im):
    # Skeletonize the binary image
    binary_skeleton = skeletonize(cropped_binary_image)

    # Identify junction points in the skeleton
    _, endpoints = find_junctions_and_endpoints(binary_skeleton)

    soma_coords, distance_transform, _ = find_soma_by_distance_transform(cropped_binary_image)
    soma_region = segment_soma(cropped_binary_image, soma_coords, distance_transform, threshold_ratio=0.5)
    skeleton_outside_soma = binary_skeleton & ~soma_region
    junctions , _ = find_junctions_and_endpoints(skeleton_outside_soma)
    

    # Calculate morphometric features
    cell_area = regionprop.area
    cell_perimeter = regionprop.perimeter
    convex_hull_area, convex_hull_perimeter, convex_hull_span_ratio, convex_hull_circularity = calculate_convex_hull_properties(cropped_binary_image > 0)
    cell_solidity = cell_area / convex_hull_area
    cell_convexity = convex_hull_perimeter / cell_perimeter
    cell_roughness = cell_perimeter / convex_hull_perimeter
    cell_circularity = (4 * np.pi * cell_area) / (cell_perimeter ** 2)
    eccentricity = regionprop.eccentricity
    
    soma_area, soma_perimeter, soma_circularity = calculate_soma_properties(soma_region)
    
    total_length, mean_branch_length, number_of_branches = calculate_total_length(skeleton_outside_soma)
    n_branches, n_tips = count_branches_and_points(junctions, endpoints)
    n_bifs = n_branches
    
    # Perform Sholl Analysis
    radii, intersections = perform_sholl_analysis(skeleton_outside_soma, soma_coords)
    sholl_parameters = calculate_sholl_parameters(radii, intersections)

    # Calculate fractal dimension
    fractal_dimension = calculate_fractal_dimension(cropped_binary_image)

    # Calculate Euclidean and path distance
    euclidean_distance, path_distance = calculate_euclidean_and_path_distance(skeleton_outside_soma)

    # Calculate lacunarity
    lacunarity = calculate_lacunarity(cropped_binary_image)

    # Calculate tortuosity
    tortuosity = calculate_tortuosity(euclidean_distance, path_distance)

    # Calculate normalized intensity of the IBA1 channel
    intensity_block = regionprops(cropped_binary_image, cropped_iba1_image)
    region_intensity = max(intensity_block, key=lambda x: x.area).intensity_mean
    normalized_intensity = region_intensity / np.mean(cropped_iba1_image)
    
    intensity_block_dapi = regionprops(cropped_binary_image, dapi_im)
    region_intensity_dapi = max(intensity_block_dapi, key=lambda x: x.area).intensity_mean
    normalized_intensity_dapi = region_intensity_dapi / np.mean(dapi_im)
    
    soma_ratio_of_cell = cell_area/soma_area
    
    # Compile results into a dictionary
    results = {
        'Cell Area': cell_area,
        'Cell Perimeter': cell_perimeter,
        'Convex Hull Area': convex_hull_area,
        'Convex Hull Perimeter': convex_hull_perimeter,
        'Cell Solidity': cell_solidity,
        'Cell Convexity': cell_convexity,
        'Cell Roughness': cell_roughness,
        'Cell Circularity': cell_circularity,
        'Convex Hull Span Ratio': convex_hull_span_ratio,
        'Convex Hull Circularity': convex_hull_circularity,
        'Soma Area': soma_area,
        'Soma Perimeter': soma_perimeter,
        'Soma Circularity': soma_circularity,
        'Soma Ratio': soma_ratio_of_cell,
        'Skeleton Length': total_length,
        'Mean Branch Length': mean_branch_length,
        'Number of Branches': number_of_branches,
        'Number of Branching Points': n_bifs,
        'Number of Terminal Points': n_tips,
        'Branching Index': sholl_parameters['Branching Index'],
        'Dendritic Maximum': sholl_parameters['Dendritic Maximum'],
        'Ramification Index': sholl_parameters['Ramification Index'],
        'Radius of Influence': sholl_parameters['Radius of Influence'],
        'Fractal Dimension': fractal_dimension,
        'Euclidean Distance': euclidean_distance,
        'Path Distance': path_distance,
        'Tortuosity': tortuosity,
        'Lacunarity': lacunarity,
        'Eccentricity': eccentricity,
        'Normalized Intensity': normalized_intensity,
        'DAPI Score': normalized_intensity_dapi
    }
    return results

def embed_image(image, extractor, ex_model, width=224, height=224):
    # resize the image to go inside of the transformer
    image = cv2.resize(image, dsize=(width,height), interpolation=cv2.INTER_NEAREST)
    
    # shift it to RGB
    image = cv2.cvtColor(image, cv2.COLOR_GRAY2RGB)
    
    image = np.transpose(image, (2, 0, 1))
    
    image_batch = np.expand_dims(image, axis=0)
    
    image_tensor = torch.tensor(image).float()
    
    image_tensor = image_tensor.to(next(ex_model.parameters()).device)

    # Extract features
    inputs = extractor(images=image_tensor, return_tensors="pt")
    inputs = {k: v.to(next(ex_model.parameters()).device) for k, v in inputs.items()}  # Move to the same device as the model
    outputs = ex_model(**inputs)
    
    return outputs.pooler_output

def regionprop_to_subimage(image_np, regionprop, margin=50):
    # Get the bounding box coordinates
    minr, minc, maxr, maxc = regionprop.bbox

    # Add a margin
    minr = max(minr - margin, 0)
    minc = max(minc - margin, 0)
    maxr = min(maxr + margin, image_np.shape[0])
    maxc = min(maxc + margin, image_np.shape[1])

    # Extract the subimage
    subimage_np = image_np[minr:maxr, minc:maxc]

    return subimage_np

def generate_vectors_for_batch(adata_filtered, base_path, batch_id, proc, mod):
    data_list = []

    transform_file = f'{base_path}{batch_id}/images/micron_to_mosaic_pixel_transform.csv'
    transform_df = pd.read_table(transform_file, sep=' ', header=None)
    transformation_matrix = transform_df.values

    test = adata_filtered[adata_filtered.obs.batchID == batch_id].copy()

    test.obs['x'] = test.obs['x'] * transformation_matrix[0, 0] + transformation_matrix[0, 2]
    test.obs['y'] = test.obs['y'] * transformation_matrix[1, 1] + transformation_matrix[1, 2]

    test_im = Mapping.load_tiff_image(base_path + batch_id + '/binary_image.tif')
    print(f"Raw Image Loaded for {batch_id}")
    dapi_im = Mapping.load_tiff_image(base_path + batch_id + '/images/mosaic_DAPI_z3.tif')

    height, width = test_im.shape[:2]

    for i in tqdm(range(len(test.obs)), desc=f'Processing {batch_id}'):
        x_ax = round(test.obs.iloc[i, :].x)
        y_ax = round(test.obs.iloc[i, :].y)
        
        # Adjusted coordinates to ensure ROI stays within bounds
        x_start = max(0, x_ax - 500)
        x_end = min(width, x_ax + 500)
        y_start = max(0, y_ax - 500)
        y_end = min(height, y_ax + 500)

        # Adjust the point tuple based on the new ROI
        adjusted_x_ax = 500 if (x_ax - 500 >= 0 and x_ax + 500 <= width) else x_ax - x_start
        adjusted_y_ax = 500 if (y_ax - 500 >= 0 and y_ax + 500 <= height) else y_ax - y_start

        point = (adjusted_x_ax, adjusted_y_ax)

        sub_image = test_im[y_start:y_end, x_start:x_end]

        sub_dapi = dapi_im[y_start:y_end, x_start:x_end]


        subtract = cv2.fastNlMeansDenoising(sub_image)
        pre = cv2.adaptiveThreshold((255 - subtract), 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                    cv2.THRESH_BINARY, 201, 2)

        opened = opening(255 - pre, disk(3))
        pre = closing(opened, disk(3))

        filled_image = binary_fill_holes(pre)

        labeled_array, num_features = label(filled_image)

        #point = (500, 500)

        regions = regionprops(labeled_array)
        
        largest_component_label = None
        max_area = 0

        for region in regions:
            if region.area > 150:
            # Find the boundary of the component
                boundary = find_boundaries(labeled_array == region.label, mode='outer')
                boundary_coords = np.column_stack(np.where(boundary))

            # Calculate the distance from the point to the boundary
                distances = np.linalg.norm(boundary_coords - np.array(point), axis=1)
                min_dist = np.min(distances)

                if min_dist < 100:
                    if region.area > max_area:
                        max_area = region.area
                        largest_component_label = region.label

        if largest_component_label is not None:
            isolated_component = labeled_array == largest_component_label

            label_image = np.zeros_like(labeled_array)
            label_image[isolated_component] = 1
            testing_prop = regionprops(label_image)
            raw_trimmed_sub = regionprop_to_subimage(sub_image, testing_prop[0], margin=50)
            norm_label = cv2.normalize(label_image, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
            seg_trimmed_sub = regionprop_to_subimage(norm_label, testing_prop[0], margin=50)
            

            features = calculate_morphometric_features(label_image, testing_prop[0], sub_image, sub_dapi)
            features['embed_raw'] = json.dumps(embed_image(raw_trimmed_sub,proc,mod).tolist()[0])
            #norm_label = cv2.normalize(label_image, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
            features['embed_segment'] = json.dumps(embed_image(seg_trimmed_sub,proc,mod).tolist()[0])
            features['Name'] = test.obs.iloc[i, :].Name
            
            data_list.append(features)

        else:
            features = {key: np.nan for key in [
                'Cell Area', 'Cell Perimeter', 'Convex Hull Area', 'Convex Hull Perimeter',
                'Cell Solidity', 'Cell Convexity', 'Cell Roughness', 'Cell Circularity',
                'Convex Hull Span Ratio', 'Convex Hull Circularity', 'Soma Area', 'Soma Perimeter',
                'Soma Circularity', 'Soma Ratio', 'Skeleton Length', 'Mean Branch Length',
                'Number of Branches', 'Number of Branching Points', 'Number of Terminal Points',
                'Branching Index', 'Critical Radius', 'Dendritic Maximum', 'Fractal Dimension',
                'Euclidean Distance', 'Path Distance', 'Tortuosity', 'Lacunarity', 'Eccentricity','Normalized Intensity', 'DAPI Score', 'embed_raw', 'embed_segment'
            ]}
            features['Name'] = test.obs.iloc[i, :].Name
            data_list.append(features)

    # Convert the list to a DataFrame
    df = pd.DataFrame(data_list)
    return df


def main(adata_path, base_path, batch_id, output_dir):
    ad_parent = sc.read_h5ad(adata_path)
    ad_parent = ad_parent[ad_parent.obs.subclass_label_transfer == 'Microglia NN']
    
    # This line for using just the ABC labels
    ad_parent = ad_parent[ad_parent.obs.subclass_label_confidence > 0.8]
    
    processor = AutoImageProcessor.from_pretrained('facebook/dinov2-base')
    model = AutoModel.from_pretrained('facebook/dinov2-base')
    DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    feature_vectors = generate_vectors_for_batch(ad_parent, base_path, batch_id, processor, model)
    output_file = os.path.join(output_dir, f'feature_vectors_{batch_id}.csv')
    feature_vectors.to_csv(output_file, index=False)
    print(f"Feature vectors for {batch_id} saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some images and generate feature vectors for a specific batch.')
    parser.add_argument('adata_path', type=str, help='Path to the .h5ad file')
    parser.add_argument('base_path', type=str, help='Base directory path for the images')
    parser.add_argument('batch_id', type=str, help='Batch ID to process')
    parser.add_argument('output_dir', type=str, help='Directory to save the output CSV files')
    args = parser.parse_args()
    main(args.adata_path, args.base_path, args.batch_id, args.output_dir)
