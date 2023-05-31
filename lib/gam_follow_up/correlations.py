import numpy as np
from scipy.stats import pearsonr

from gamtools import cosegregation, matrix


chroms = ['chr{}'.format(c) for c in range(1,20)]


def keep_diagonals(in_matrix, diag_from, diag_to=None):
    
    if diag_to is None:
        diag_to = len(in_matrix)
        
    out_matrix = np.zeros_like(in_matrix) * np.nan
    
    for d in range(diag_from, diag_to):
        indices = matrix.kth_diag_indices(out_matrix, d)
        out_matrix[indices] = in_matrix[indices]
        
    return out_matrix


def clean_nas(arr1, arr2):
    
    arr1[arr1 == -1] = np.NaN
    arr2[arr2 == -1] = np.NaN
    
    neither_na = np.logical_not(np.logical_and(np.isnan(arr1), np.isnan(arr2)))
    
    arr1 = np.nan_to_num(arr1[neither_na])
    arr2 = np.nan_to_num(arr2[neither_na])
    
    return arr1, arr2


def clean_matrices(mat_1, mat_2):
    
    flat_1 = mat_1.flatten()
    flat_2 = mat_2.flatten()
    
    flat_1, flat_2 = clean_nas(flat_1, flat_2)
    
    return flat_1, flat_2


def correlate_matrices(matrix_1, matrix_2):
    
    mat_1_vals, mat_2_vals = clean_matrices(matrix_1, matrix_2)
    
    return pearsonr(mat_1_vals, mat_2_vals)[0]


def correlate_from_seg_tables(seg_table1, seg_table2, matrix_regions=[(0, None)]):

    out_correlations = [[] for i in matrix_regions]

    for chrom in chroms:

        chrom_matrix_1 = cosegregation.get_dprime(seg_table1, chrom)

        chrom_matrix_2 = cosegregation.get_dprime(seg_table2, chrom)

        for i, (d_start, d_stop) in enumerate(matrix_regions):

            this_corr= correlate_matrices(keep_diagonals(chrom_matrix_1, d_start, d_stop),
                                          keep_diagonals(chrom_matrix_2, d_start, d_stop))

            out_correlations[i].append(this_corr)
            
    return out_correlations
