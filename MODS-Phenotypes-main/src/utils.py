from IPython.display import display as display
from multiprocessing import Pool
from pathlib import Path
from tqdm.auto import tqdm
import argparse
import numpy as np
import pandas as pd
import pickle
import pyarrow as pa
import pyarrow.parquet as pq
import os
from concurrent.futures import ThreadPoolExecutor
import random
from src.log import setup_logger
logger = setup_logger("utils")

def read_parquet_file(file_path, schema=None, pbar=None):
    """
    Read a single parquet file and optionally apply a given schema.
    
    This function reads the content of a single parquet file from the given path 
    and applies the specified schema, if provided. If a progress bar (tqdm) 
    instance is provided, the function updates the progress bar after reading 
    the file.
    
    Parameters:
    - file_path (str or Path): The path to the parquet file to be read.
    - schema (pyarrow.Schema, optional): A pyarrow schema to be applied when reading the file.
    - pbar (tqdm, optional): A tqdm progress bar instance to be updated after reading the file.
    
    Returns:
    - pyarrow.Table: A table containing the data from the read parquet file.
    """
    try:
        table = pq.read_table(file_path, schema=schema)
        if pbar:
            pbar.update(1)
        return table
    except pa.lib.ArrowInvalid as e:
        if "Failed to parse string: '(null)'" in str(e):
            logger.error(f"Error reading file: {file_path}. Error: {str(e)}")
            # TODO: Log the problematic columns or any other details you need.
        else:
            # Reraise any other unexpected errors.
            raise
        return None

def read_parquet_files_in_parallel(directory, schema=None, max_workers=None):
    """
    Read multiple parquet files from a directory in parallel and return a concatenated table.
    
    This function searches for all parquet files in the given directory, reads 
    them in parallel using multiple threads, and returns a concatenated table 
    of all the read files. A tqdm progress bar is displayed to show the progress 
    of reading. An optional schema can be applied to each file during the reading.
    
    Parameters:
    - directory (str or Path): The directory containing the parquet files to be read.
    - schema (pyarrow.Schema, optional): A pyarrow schema to be applied when reading each file.
    - max_workers (int, optional): The maximum number of threads to be used for parallel reading.
                                   If None, the default number of threads will be used.
    
    Returns:
    - pyarrow.Table: A concatenated table containing the data from all the read parquet files.
    """
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.parquet')]
    
    with tqdm(total=len(files), desc="Reading files", position=0, leave=True) as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(lambda f: read_parquet_file(f, schema, pbar), files))
    
    # Filter out the None values, which indicate files that raised an error.
    valid_results = [res for res in results if res is not None]
    
    return pa.concat_tables(valid_results)

def make_arrow_schema(schema_config):
    """
    Constructs a pyarrow schema based on the provided configuration.
    
    This function processes the given schema configuration, identifies fields 
    that have a data type of "LIST(TIMESTAMP[NS])", and converts them to 
    the appropriate pyarrow type. Fields with other data types are added 
    directly to the schema.
    
    Parameters:
    - schema_config (dict): A dictionary representing the desired schema. 
                            Keys are field names and values are data type descriptors.
    
    Returns:
    - pa.Schema: A pyarrow schema constructed based on the provided configuration.
    """
    datetime_arrays = set()

    for k,v in schema_config.items():
        if v == "LIST(TIMESTAMP[NS])":
            datetime_arrays.add(k)
        else:
            pass

    for entry_value in datetime_arrays:
        schema_config.pop(entry_value, None)

    arrow_schema = pa.schema(schema_config)

    for field_name in datetime_arrays:
        arrow_schema = arrow_schema.append(pa.field(field_name, pa.list_(pa.date64())))
    
    return arrow_schema

def load_pickle(pickle_fp: str) -> dict:
    """
    Load a pickle file.

    Parameters:
    - pickle_fp (str): Path to the pickle file.

    Returns:
    - dict: Unpickled object.
    """
    with open(pickle_fp, "rb") as file:
        pkl = pickle.load(file, encoding="bytes")
    return pkl


def zipsort(columns: list, features: list) -> list:
    """
    Sort and zip two lists.

    Parameters:
    - columns (list): First list.
    - features (list): Second list.

    Returns:
    - list: Zipped result.
    """
    return zip(sorted(columns), sorted(features))


def flatten(xss: list) -> list:
    """
    Flatten a list of lists.

    Parameters:
    - xss (list): Nested list.

    Returns:
    - list: Flattened list.
    """
    return [x for xs in xss for x in xs]


def time_differentials(start_dts: pd.Series, end_dts: pd.Series) -> pd.Series:
    """
    Calculate the time differential in hours for two DateTime pd.Series.

    Parameters:
    - start_dts (pd.Series): Start dates.
    - end_dts (pd.Series): End dates.

    Returns:
    - pd.Series: Time differentials.
    """
    deltas = end_dts - start_dts
    return deltas.dt.days * 24 + deltas.dt.seconds / 3600.0


def find_files(source_path: Path, ext: str = ".pickle") -> list:
    """
    Find files with a specific extension in a directory.

    Parameters:
    - source_path (Path): Directory to be searched.
    - ext (str, optional): Extension of files to find. Default is ".pickle".

    Returns:
    - list: List of file paths.
    """
    fps = [fp for fp in source_path.iterdir() if ext in str(fp)]
    return fps


def find_pickle_paths(file_path, years, sample_rate: float = 1.0) -> list:
    years = set(range(years['start_year'], years['stop_year'])) - set(years['skip_years'])
    pickle_paths = flatten(
        run_imap_multiprocessing(
            find_files, ([Path(file_path)
                         / str(yr) for yr in years]), len(years)
        )
    )
    
    if sample_rate:
        pickle_paths = random.sample(pickle_paths, int(len(pickle_paths)*sample_rate))
    
    return pickle_paths

def denoise_times(times: pd.Series) -> np.ndarray:
    """
    Remove noise from a series of times.

    Parameters:
    - times (pd.Series): Series of times.

    Returns:
    - np.ndarray: Array of unique, sorted times.
    """
    return np.sort(times.array.unique().dropna())


def load_encounter_pickle(pickle_fp: str) -> tuple:
    """
    Load an encounter pickle and return the pickle and its filename.

    Parameters:
    - pickle_fp (str): Path to the pickle file.

    Returns:
    - tuple: Loaded pickle and the filename without extension.
    """
    with open(pickle_fp, "rb") as file:
        pkl = pickle.load(file, encoding="bytes")
    return pkl, pickle_fp.stem


def parse_dts_str(dts_str: str) -> list:
    """
    Parse datetime strings into a list of datetime objects.

    Parameters:
    - dts_str (str): Datetime string.

    Returns:
    - list: List of parsed datetime objects.
    """
    dts_list = dts_str.split(" ")
    dts_parsed = []
    for dt in dts_list:
        dts_parsed.append(
            pd.to_datetime(
                dt.replace("[", "").replace("]", "").replace("'", ""),
                format="%Y-%m-%d %H:%M:%S.%f",
            )
        )
    return dts_parsed


def restricted_float(x: float) -> float:
    """
    Ensure that a float value is within [0.0, 1.0].

    Parameters:
    - x (float): Input float value.

    Returns:
    - float: Validated float value.
    """
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


def flatten_list(nested_list: list) -> list:
    """
    Flatten a potentially nested list.

    Parameters:
    - nested_list (list): Nested list.

    Returns:
    - list: Flattened list.
    """
    flattened_list = []
    for item in nested_list:
        if isinstance(item, list):
            flattened_list.extend(flatten_list(item))
        else:
            flattened_list.append(item)
    return flattened_list


def run_imap_multiprocessing(func, argument_list: list, num_processes: int) -> list:
    """
    Run a function in parallel using multiprocessing and display progress.

    Parameters:
    - func (callable): Function to run.
    - argument_list (list): Arguments to pass to the function.
    - num_processes (int): Number of processes for multiprocessing.

    Returns:
    - list: List of results from the function.
    """
    pool = Pool(processes=num_processes)
    result_list_tqdm = []
    for result in tqdm(
        pool.imap(func=func, iterable=argument_list), total=len(argument_list)
    ):
        result_list_tqdm.append(result)
    return result_list_tqdm
