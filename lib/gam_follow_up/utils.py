import pandas as pd
import numpy as np

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return [list(grp) for grp in zip(*args)]

def get_merged(original_segmentation, merge=2):
    column_indexes = list(range(len(original_segmentation.columns)))

    # Shuffle column indexes in place
    np.random.shuffle(column_indexes)

    # Group them in threes
    index_sets = list(grouper(column_indexes, merge))

    # Merge random columns of the segmentation in groups of three
    return pd.concat(
        [original_segmentation.iloc[:,indexes].any(axis=1).astype(int) for indexes in index_sets],
        axis=1)

def split_df(input_df, subsample_n=None):
    if subsample_n is None:
        subsample_n = len(input_df.columns) / 2
    assert (subsample_n * 2) <= len(input_df.columns)
    
    selected_cols = np.random.choice(input_df.columns, (subsample_n * 2), replace=False)
    selected_cols_1 = selected_cols[:subsample_n]
    selected_cols_2 = selected_cols[subsample_n:]
    
    assert len(selected_cols_1) == len(selected_cols_2)
    
    return input_df[selected_cols_1], input_df[selected_cols_2]

def chain_and(*conditions):
    conditions = list(conditions)
    all_true = conditions.pop()
    while len(conditions):
        new_cond = conditions.pop()
        all_true = np.logical_and(all_true, new_cond)
    return all_true
