import os

base_folder = os.path.realpath(os.path.join(__file__, '../../..'))

data_raw_folder = os.path.join(base_folder, 'data', 'raw')
data_processed_folder = os.path.join(base_folder, 'data', 'processed')
data_intermediate_folder = os.path.join(base_folder, 'data', 'intermediate')
figures_folder = os.path.join(base_folder, 'out', 'figures')

def in_data_raw(path):
    return os.path.join(data_raw_folder, path)

def in_data_processed(path):
    return os.path.join(data_processed_folder, path)

def in_data_intermediate(path):
    return os.path.join(data_intermediate_folder, path)

def in_figures(path):
    return os.path.join(figures_folder, path)

