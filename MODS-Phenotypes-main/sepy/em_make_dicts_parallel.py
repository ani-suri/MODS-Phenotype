import sys
import pandas as pd

# import numpy as np
import sepyDICT as sd
import pickle
import time
from pathlib import Path
import logging
from functools import partial
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from multiprocessing import Pool, cpu_count


## (NOTE FROM GLENN)
## NONE OF THIS DATA SHOULD BE SAVED TO PICKLES
## FOR GODS SAKE PLEASE USE PARQUET INSTEAD !!!
# import pyarrow.parquet as pq
#
# def write_data_to_parquet(data, path):
#     """Utility function to write data to Parquet."""
#     table = pa.Table.from_pandas(data)
#     pq.write_table(table, path)


def process_csn(yearly_instance, pickle_write_path: Path, csn: int):
    """Process CSN and write to pickle."""
    file_name = pickle_write_path / f"{csn}.pickle"

    # Instantiate class for single encounter
    try:
        instance = sd.sepyDICT(yearly_instance, csn)

    except KeyError as ke:
        logger.error(f"Encountered KeyError for CSN {csn}: {ke}")
        return  # ends the processing and moves on

    # Create encounter dictionary
    inst_dict = instance.encounter_dict

    # Create a pickle file for encounter
    with open(file_name, "wb") as picklefile:
        pickle.dump(inst_dict, picklefile)

    logger.debug("MkDct- Instance created")
    summary_functions = [
        sofa_summary,
        sepsis3_summary,
        sirs_summary,
        sepsis2_summary,
        enc_summary,
        comorbidity_summary,
    ]

    for summary_func in summary_functions:
        try:
            summary_func(csn, instance)
        except Exception as e:
            logger.error(f"MkDct- {summary_func.__name__} in CSN {csn}: {e}")


def write_data_to_csv(data, path, index=False, index_label=None):
    """Utility function to write data to CSV."""
    data.to_csv(path, index=index, index_label=index_label)


def add_summary_to_appended_list(
    csn: int,
    data,
    appended_list: list,
    reset_index: bool = False,
    rename_dict: dict = None,
):
    """Utility function to append data to a list after processing."""
    if reset_index:
        data = data.reset_index()
    if rename_dict:
        data = data.rename(columns=rename_dict)
    data["csn"] = csn
    appended_list.append(data)


def sofa_summary(csn: int, instance: sd.sepyDICT, appended_list: list):
    """Add SOFA scores to summary."""
    add_summary_to_appended_list(
        csn,
        instance.encounter_dict["sofa_scores"],
        appended_list,
        reset_index=True,
        rename_dict={"index": "time_stamp"},
    )


def sepsis3_summary(csn: int, instance: sd.sepyDICT, appended_list: list):
    """Add SEPSIS 3 times to summary."""
    add_summary_to_appended_list(
        csn, instance.encounter_dict["sep3_time"], appended_list
    )


def sirs_summary(csn: int, instance: sd.sepyDICT, appended_list: list):
    """Add SIRS scores to summary."""
    add_summary_to_appended_list(
        csn,
        instance.encounter_dict["sirs_scores"],
        appended_list,
        reset_index=True,
        rename_dict={"index": "time_stamp"},
    )


def sepsis2_summary(csn: int, instance: sd.sepyDICT, appended_list: list):
    """Add SEPSIS 2 times to summary."""
    add_summary_to_appended_list(
        csn, instance.encounter_dict["sep2_time"], appended_list
    )


def enc_summary(csn: int, instance: sd.sepyDICT, appended_list: list):
    """Add encounter summaries."""
    enc_summary_dict = {
        **instance.flags,
        **instance.static_features,
        **instance.event_times,
    }
    enc_summary_df = pd.DataFrame(enc_summary_dict, index=[0]).set_index(["csn"])
    appended_list.append(enc_summary_df)


def comorbidity_summary(
    csn: int, instance: sd.sepyDICT, quan_deyo_dict: dict, quan_elix_dict: dict
):
    """Add comorbidity summaries."""
    quan_deyo_dict[csn] = instance.quan_deyo_ICD10_PerCSN.icd_count
    quan_elix_dict[csn] = instance.quan_elix_ICD10_PerCSN.icd_count


def main():
    start = time.perf_counter()

    # Initialize empty vars

    # ICD10
    quan_deyo_ICD10_dict = {}
    quan_elix_ICD10_dict = {}

    # other summaries
    appended_sofa_scores = []
    appended_sep3_time = []
    appended_sirs_scores = []
    appended_sep2_time = []
    appended_enc_summaries = []

    try:
        logger.debug(sys.version)

        csn_list_file_name = Path(sys.argv[1])  # imports file location for CSN list
        logger.info(f"MkDct- The CSN List location: {csn_list_file_name}")

        pickle_path = Path(sys.argv[2])  # sets directoy for yearly pickles
        logger.info(f"MkDct- The pickle directory: {pickle_path}")

        output_path = Path(sys.argv[3])  # sets directory for output
        logger.info(f"MkDct- The output directory: {output_path}")

        num_processes = int(sys.argv[4])  # total number of tasks/processes
        logger.info(f"MkDct- The tot num of processes: {num_processes}")

        bash_year = int(sys.argv[5])  # year num provided in bash
        logger.info(f"MkDct- The import year is: {bash_year}")

    except:
        logger.error("MkDct- There was an error importing one of the arguments.")
        logger.error(
            f"MkDct- You are trying to load the following CSN list {sys.argv[1]}"
        )
        logger.error(f"MkDct- You are trying to use this many processors {sys.argv[2]}")

    #####################################################################
    ########### Creat Encounter List Based on Processor Assignment
    #####################################################################
    # Cohort selector
    ed = 0
    in_pt = 1
    icu = 0
    adult = 1
    vent_row = 0
    vent_start = 0

    # reads the list of csns
    csn_df = pd.read_csv(csn_list_file_name, sep="|", header=0)

    #  only keep csns that meet specified year
    csn_df = csn_df[(csn_df.in_pt == in_pt) & (csn_df.adult == adult)]
    logger.info(f"MkDct- A total of {csn_df.shape[0]} encounters were selected")
    csn_df = csn_df[["csn", "year"]]

    # drop duplicates
    csn_df = csn_df.drop_duplicates()
    csn_df = csn_df[csn_df.year == bash_year]
    total_num_enc = len(csn_df)
    logger.info(f"MkDct- The year {bash_year} has {total_num_enc} encounters.")

    # TODO: GET RID OF CHUNKING CODE!
    # breaks encounter list into chunks, selects correct chunk based on process num
    # chunk_size = int(total_num_enc / num_processes)
    # logger.info(f"MkDct- The ~chunk size is {chunk_size}")

    # split list
    # list_of_chunks = np.array_split(csn_df, num_processes)
    # logger.info(
    #     f"MkDct- The list of chunks has {len(list_of_chunks)} unique dataframes."
    # )

    # uses processor assignment to select correct chunk
    # process_list = list_of_chunks[processor_assignment]["csn"]
    # logger.info(f"MkDct- The process_list head:\n {process_list.head()}")

    # Just needs the CSNs, no chunking
    process_list = csn_df["csn"].tolist()

    # select correct pickle by year
    pickle_name = pickle_path / ("em_y" + str(bash_year) + ".pickle")
    logger.info(f"MkDct- The following pickle is being read: {pickle_name}")

    try:
        # reads the IMPORT class instance (i.e.  1 year of patient data)
        pickle_load_time = time.perf_counter()  # time to load pickle
        with open(pickle_name, "rb") as handle:
            yearly_instance = pickle.load(handle)
            logger.info(
                f"MkDct-Pickle from year {bash_year} was loaded in {time.perf_counter()-pickle_load_time}s."
            )

        logger.info("-----------LOADED YEARLY PICKLE FILE!!!!---------------")

        # if success, make a dir for this year's encounters
        pickle_write_path = output_path / str(bash_year)
        Path.mkdir(pickle_write_path, exist_ok=True)
        logger.info(
            f"MkDct-Directory for year {bash_year} was set to {pickle_write_path}"
        )

        # make empty list to handle csn's with errors
        error_list = []
        start_csn_creation = time.perf_counter()  # times calc's by year

        #################################################
        ############ Make Dicts by CSN
        #################################################

        process_csn_func = partial(process_csn, yearly_instance, pickle_write_path)

        # TODO: Add tqdm progress bar for the pool and a total count, then delete these loggers
        # logger.info(f'MkDct- The current pt csn is: {csn}, which is {count} of {chunk_size} for year {bash_year}')
        # logger.info(f'MkDct- Enounter {count} of {chunk_size} is complete!')
        with Pool(num_processes) as pool:
            # Create tqdm instance
            pbar = tqdm(
                total=len(process_list), desc="Processing CSNs", dynamic_ncols=True
            )

            # Use imap instead of map
            for _ in pool.imap(process_csn_func, process_list):
                pbar.update(1)

            pbar.close()

        #################################################
        ############ Export Sepsis Summary
        #################################################

        # Create directories
        base_dir = output_path / "em_sepsis_summary"
        dirs_to_create = [
            "sofa_summary",
            "sep3_summary",
            "sirs_summary",
            "sep2_summary",
            "encounter_summary",
            "error_summary",
            "quan_deyo_ICD10_summary",
            "quan_elix_ICD10_summary",
        ]

        for dir_name in dirs_to_create:
            Path.mkdir(base_dir / dir_name, exist_ok=True, parents=True)

        # Write data to files
        unique_file_id = f"{processor_assignment}_{bash_year}"
        write_data_to_csv(
            pd.concat(appended_enc_summaries),
            base_dir / "encounter_summary" / f"encounters_summary_{unique_file_id}.csv",
            index=True,
        )
        write_data_to_csv(
            pd.DataFrame(error_list, columns=["csn", "error"]),
            base_dir / "error_summary" / f"error_list_{unique_file_id}.csv",
        )
        write_data_to_csv(
            pd.concat(appended_sofa_scores),
            base_dir / "sofa_summary" / f"sofa_summary_{unique_file_id}.csv",
        )
        write_data_to_csv(
            pd.concat(appended_sep3_time),
            base_dir / "sep3_summary" / f"sepsis3_summary_{unique_file_id}.csv",
        )
        write_data_to_csv(
            pd.concat(appended_sirs_scores),
            base_dir / "sirs_summary" / f"sirs_summary_{unique_file_id}.csv",
        )
        write_data_to_csv(
            pd.concat(appended_sep2_time),
            base_dir / "sep2_summary" / f"sepsis2_summary_{unique_file_id}.csv",
        )
        write_data_to_csv(
            pd.DataFrame.from_dict(quan_deyo_ICD10_dict).T,
            base_dir
            / "quan_deyo_ICD10_summary"
            / f"quan_deyo_ICD10_summary_{unique_file_id}.csv",
            index=True,
            index_label="csn",
        )
        write_data_to_csv(
            pd.DataFrame.from_dict(quan_elix_ICD10_dict).T,
            base_dir
            / "quan_elix_ICD10_summary"
            / f"quan_elix_ICD10_summary_{unique_file_id}.csv",
            index=True,
            index_label="csn",
        )

        #################################################

        logger.info(
            f"MkDct- Time to create write encounter pickles for {bash_year} was {time.perf_counter()-start_csn_creation} (s)"
        )

    except FileNotFoundError:
        logger.error(f"MkDct- Could not find or open the pickle for year {bash_year}.")
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        raise


if __name__ == "__main__":
    main()
