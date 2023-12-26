# -*- coding: utf-8 -*-
"""
Elite Data Hacks
@author: Christopher S. Josef, MD
@email: csjosef@krvmail.com
"""

# import fill_variable as fv
from sepyAGG import labAGG
import time
import pandas as pd
import numpy as np

from functools import reduce

# libs that can probably be removed
# import itertools
# import os


class sepyDICT:
    ########################################################################################################
    ####################                  Class Variables                                ###################
    ########################################################################################################
    v_vital_col_names = [
        "temperature",
        "daily_weight_kg",
        "height_cm",
        "sbp_line",
        "dbp_line",
        "map_line",
        "sbp_cuff",
        "dbp_cuff",
        "map_cuff",
        "pulse",
        "unassisted_resp_rate",
        "spo2",
        "end_tidal_co2",
        "o2_flow_rate",
    ]

    # List of all lab names (some years might not have all listed labs)
    v_numeric_lab_col_names = [
        "anion_gap",
        "base_excess",
        "bicarb_(hco3)",
        "blood_urea_nitrogen_(bun)",
        "calcium",
        "calcium_adjusted",
        "calcium_ionized",
        "chloride",
        "creatinine",
        "gfr",
        "glucose",
        "magnesium",
        "osmolarity",
        "phosphorus",
        "potassium",
        "sodium",
        # CBC
        "haptoglobin",
        "hematocrit",
        "hemoglobin",
        "met_hgb",
        "platelets",
        "white_blood_cell_count",
        "carboxy_hgb",
        # Hepatic
        "alanine_aminotransferase_(alt)",
        "albumin",
        "alkaline_phosphatase",
        "ammonia",
        "aspartate_aminotransferase_(ast)",
        "bilirubin_direct",
        "bilirubin_total",
        "fibrinogen",
        "inr",
        "lactate_dehydrogenase",
        "lactic_acid",
        "partial_prothrombin_time_(ptt)",
        "prealbumin",
        "protein",
        "prothrombin_time_(pt)",
        "thrombin_time",
        "transferrin",
        # Pancreatic
        "amylase",
        "lipase",
        # Cardiac
        "b-type_natriuretic_peptide_(bnp)",
        "troponin",
        # ABG
        "carboxy_hgb",
        "fio2",
        "partial_pressure_of_carbon_dioxide_(paco2)",
        "partial_pressure_of_oxygen_(pao2)",
        "ph",
        "saturation_of_oxygen_(sao2)",
        # Other
        "d_dimer",
        "hemoglobin_a1c",
        "parathyroid_level",
        "thyroid_stimulating_hormone_(tsh)",
        # Inflammation
        "crp_high_sens",
        "procalcitonin",
        "erythrocyte_sedimentation_rate_(esr)",
    ]

    v_string_lab_col_names = [
        # PCR Testing
        "c_diff",
        "covid",
        "mtp",
    ]

    v_all_lab_col_names = v_numeric_lab_col_names + v_string_lab_col_names

    # Glasgow Comma Scale Cols
    v_gcs_col_names = [
        "gcs_total_score",
        "gcs_verbal_score",
        "gcs_eye_score",
        "gcs_motor_score",
    ]

    # Bed Location Cols
    v_bed_info = [
        "bed_location_start",
        "bed_location_end",
        "bed_unit",
        "bed_room",
        "bed_id",
        "bed_label",
        "hospital_service",
    ]

    # Vasopressor cols
    v_vasopressor_names = [
        "norepinephrine",
        "epinephrine",
        "dobutamine",
        "dopamine",
        "phenylephrine",
        "vasopressin",
    ]

    # Vasopressor units
    v_vasopressor_units = [
        "norepinephrine_dose_unit",
        "epinephrine_dose_unit",
        "dobutamine_dose_unit",
        "dopamine_dose_unit",
        "phenylephrine_dose_unit",
        "vasopressin_dose_unit",
    ]

    # Vasopressor Dose by weight
    v_vasopressor_dose = [
        "norepinephrine_dose_weight",
        "epinephrine_dose_weight",
        "dobutamine_dose_weight",
        "dopamine_dose_weight",
        "phenylephrine_dose_weight",
        "vasopressin_dose_weight",
    ]

    # Vassopressors Paired by Name
    v_vasopressor_col_names = [
        "norepinephrine",
        "norepinephrine_dose_unit",
        "norepinephrine_dose_weight",
        "epinephrine",
        "epinephrine_dose_unit",
        "epinephrine_dose_weight",
        "dobutamine",
        "dobutamine_dose_unit",
        "dobutamine_dose_weight",
        "dopamine",
        "dopamine_dose_unit",
        "dopamine_dose_weight",
        "phenylephrine",
        "phenylephrine_dose_unit",
        "phenylephrine_dose_weight",
        "vasopressin",
        "vasopressin_dose_unit",
        "vasopressin_dose_weight",
    ]

    # Vent Col
    v_vent_col_names = ["Status"]

    v_vent_positive_vars = [
        "vent_mode",
        "vent_rate_set",
        "vent_tidal_rate_set",
        "vent_tidal_rate_exhaled",
        "peep",
    ]  #'fio2'

    # Blood Pressure Cols
    v_bp_cols = ["sbp_line", "dbp_line", "map_line", "sbp_cuff", "dbp_cuff", "map_cuff"]

    # SOFA Cols
    v_sofa_max_24h = [
        "SOFA_coag",
        "SOFA_coag_24h_max",
        "SOFA_renal",
        "SOFA_renal_24h_max",
        "SOFA_hep",
        "SOFA_hep_24h_max",
        "SOFA_neuro",
        "SOFA_neuro_24h_max",
        "SOFA_cardio",
        "SOFA_cardio_24h_max",
        "SOFA_resp",
        "SOFA_resp_24h_max",
        "hourly_total",
        "hourly_total_24h_max",
        "delta_24h",
        "delta_24h_24h_max",
    ]

    ########################################################################################################
    ####################                  Instance Variables                             ###################
    ########################################################################################################

    # The following function accepts a list of patient mrns and creates a dictionary
    def __init__(self, imported, csn):
        print(f"SepyDICT- Creating sepyDICT instance for {csn}")
        filter_date_start_time = time.time()

        ##############################
        self.csn = csn
        # set the pat id based on the encounter; take first incase multiple encounters
        try:
            self.pat_id = imported.df_encounters.loc[csn, ["pat_id"]].iloc[0].item()
        except:
            self.pat_id = imported.df_encounters.loc[csn, ["pat_id"]].iloc[0]

        # get filtered dfs for each patient encounter
        self.try_except(imported, self.pat_id, "demographics")

        # print(self.demographics_perCSN)
        self.try_except(imported, csn, "encounters")
        self.try_except(imported, csn, "labs")
        self.try_except(imported, csn, "vasopressor_meds")
        self.try_except(imported, csn, "anti_infective_meds")
        self.try_except(imported, csn, "vitals")
        self.try_except(imported, csn, "vent")
        self.try_except(imported, csn, "gcs")
        self.try_except(imported, csn, "cultures")
        self.try_except(imported, csn, "beds")
        self.try_except(imported, csn, "procedures")
        self.try_except(imported, csn, "diagnosis")
        self.try_except(imported, csn, "quan_deyo_ICD10")
        self.try_except(imported, csn, "quan_elix_ICD10")
        print("SepyDICT- Now making dictionary")
        self.make_dict_elements(imported)
        print("SepyDICT- Now calcuating Sepsis-2")
        self.run_SEP2()
        print("SepyDICT- Now calcuating Sepsis-3")
        self.run_SEP3()
        print("SepyDICT- Now writing dictionary")
        self.write_dict()
        print(
            f"SepyDICT- Selecting data and writing this dict by CSN took {time.time() - filter_date_start_time}(s)."
        )

    def try_except(self, imported, csn, name):
        filt_df_name = name + "_PerCSN"
        df_name = "df_" + name

        try:
            # print(getattr(imported, df_name).loc[[str(csn)],:])
            if name == "demographics":
                setattr(
                    self, filt_df_name, getattr(imported, df_name).loc[[str(csn)], :]
                )
            else:
                setattr(self, filt_df_name, getattr(imported, df_name).loc[[csn], :])
            # print(f'The {name} file was imported')

        except Exception as e:
            # print(e)
            empty_df = getattr(imported, df_name).iloc[0:0]
            empty_df.index.set_names(getattr(imported, df_name).index.names)

            setattr(self, filt_df_name, empty_df)
            print(f"The were no {name} data for csn {csn}")

    def flag_dict(self):
        self.flags = {}

        # ID numbers
        self.flags["csn"] = self.csn
        self.flags["pt_id"] = self.pat_id

        # vent flags
        self.flags["y_vent_rows"] = 0
        self.flags["y_vent_start_time"] = 0
        self.flags["y_vent_end_time"] = 0
        self.flags["vent_start_time"] = pd.NaT

    def static_features_dict(self):
        #######################################
        # static_features: Patient demographic & encounter features that will not change during admisssion
        #######################################
        # from encounters file
        self.static_features = {}
        # some patients have >1 encounter row; have to take 1st row
        self.static_features["admit_reason"] = self.encounters_PerCSN.iloc[0, :][
            "admit_reason"
        ]
        self.static_features["ed_arrival_source"] = self.encounters_PerCSN.iloc[0, :][
            "ed_arrival_source"
        ]
        self.static_features["total_icu_days"] = self.encounters_PerCSN.iloc[0, :][
            "total_icu_days"
        ]
        self.static_features["total_vent_days"] = self.encounters_PerCSN.iloc[0, :][
            "total_vent_days"
        ]
        self.static_features["total_hosp_days"] = self.encounters_PerCSN.iloc[0, :][
            "total_hosp_days"
        ]
        self.static_features["discharge_status"] = self.encounters_PerCSN.iloc[0, :][
            "discharge_status"
        ]
        self.static_features["discharge_to"] = self.encounters_PerCSN.iloc[0, :][
            "discharge_to"
        ]
        self.static_features["encounter_type"] = self.encounters_PerCSN.iloc[0, :][
            "encounter_type"
        ]
        self.static_features["age"] = self.encounters_PerCSN.iloc[0, :]["age"]
        # some patients have >1 demographic row; have to take 1st row
        # print(self.demographics_PerCSN)
        self.static_features["gender"] = self.demographics_PerCSN.iloc[0, :]["gender"]
        self.static_features["gender_code"] = self.demographics_PerCSN.iloc[0, :][
            "gender_code"
        ]
        self.static_features["race"] = self.demographics_PerCSN.iloc[0, :]["race"]
        self.static_features["race_code"] = self.demographics_PerCSN.iloc[0, :][
            "race_code"
        ]
        self.static_features["ethnicity"] = self.demographics_PerCSN.iloc[0, :][
            "ethnicity"
        ]
        self.static_features["ethnicity_code"] = self.demographics_PerCSN.iloc[0, :][
            "ethnicity_code"
        ]
        # self.static_features ['last4_ssn'] = self.demographics_PerCSN.iloc[0,:]['last4_ssn']

    def event_times_dict(self):
        #######################################
        # event_times: Key event times during a patients admission not otherwise specified
        #######################################
        # print(self.encounters_PerCSN.ed_presentation_time)
        self.event_times = {}
        self.event_times["ed_presentation_time"] = self.encounters_PerCSN.iloc[0, :][
            "ed_presentation_time"
        ]
        self.event_times["hospital_admission_date_time"] = self.encounters_PerCSN.iloc[
            0, :
        ]["hospital_admission_date_time"]
        self.event_times["hospital_discharge_date_time"] = self.encounters_PerCSN.iloc[
            0, :
        ]["hospital_discharge_date_time"]
        self.event_times["start_index"] = min(
            self.encounters_PerCSN.iloc[0, :]["hospital_admission_date_time"],
            self.encounters_PerCSN.iloc[0, :]["ed_presentation_time"],
        )
        # Wait time
        self.flags["ed_wait_time"] = (
            self.event_times["hospital_admission_date_time"]
            - self.event_times["ed_presentation_time"]
        ).total_seconds() / 60

        # bed_df = self.beds_PerCSN

    def build_super_table_index(self):
        # this is index is used in the creation of super_table

        start_time = self.event_times["start_index"]
        end_time = self.event_times["hospital_discharge_date_time"]
        self.super_table_time_index = pd.date_range(start_time, end_time, freq="H")

    def cultures_df(self):
        #######################################
        # cultures selects unique for the encounter
        #######################################
        selected_culture_cols = [
            "proc_code",
            "proc_desc",
            "component_id",
            "component",
            "loinc_code",
            "specimen_collect_time",
            "order_time",
            "order_id",
            "lab_result_time",
            "result_status",
            "lab_result",
        ]

        self.cultures_staging = self.cultures_PerCSN[selected_culture_cols]

    def antibiotics_df(self):
        self.abx_staging = self.anti_infective_meds_PerCSN

    def bin_labs(self):
        df = self.labs_PerCSN

        if df.empty:
            # drop the multi index and keep only collection time
            df.index = df.index.get_level_values("collection_time")
            # create new index with super table time index
            self.labs_staging = df.reindex(self.super_table_time_index)

        else:
            df = df.reset_index("collection_time")
            self.labs_staging = (
                df.resample(
                    "60min",
                    on="collection_time",
                    origin=self.event_times["start_index"],
                )
                .agg(labAGG)
                .reindex(self.super_table_time_index)
            )
            self.labs_staging.columns = [x[0] for x in self.labs_staging.columns]

    def bin_vitals(self):
        df = self.vitals_PerCSN

        if df.empty:
            # drop the multi index and keep only collection time
            # df.index = df.index.get_level_values('recorded_time')

            # create new index with super table time index
            self.vitals_staging = df.reindex(self.super_table_time_index)

        else:
            self.vitals_staging = (
                df.resample(
                    "60min", on="recorded_time", origin=self.event_times["start_index"]
                )
                .mean()
                .reindex(self.super_table_time_index)
            )

    def bin_gcs(self):
        df = self.gcs_PerCSN

        if df.empty:
            df = df.drop(columns=["recorded_time"])
            self.gcs_staging = df.reindex(self.super_table_time_index)

        else:
            df = df.resample(
                "60min", on="recorded_time", origin=self.event_times["start_index"]
            ).min()
            df = df.drop(columns=["recorded_time"])
            self.gcs_staging = df.reindex(self.super_table_time_index)

    def bin_vent(self):
        df = self.vent_PerCSN

        if df.empty:
            # No vent times were found so return empty table with
            # all flags remain set at zero
            df = pd.DataFrame(
                columns=["vent_status", "fio2"], index=self.super_table_time_index
            )
            # vent_status and fio2 will get joined to super table later
            self.vent_status = df.vent_status
            self.vent_fio2 = df.fio2

        else:
            # check to see there is a start & stop time
            vent_start = df[df.vent_start_time.notna()].vent_start_time.values
            vent_stop = df[df.vent_stop_time.notna()].vent_stop_time.values

            # If no vent start time then examin vent_plus rows
            if vent_start.size == 0:
                # flags for vent start/stop already set to 0; set the presence of rows below
                self.flags["y_vent_rows"] = 1

                # identify rows that are real vent vals (i.e. no fio2 alone)
                df["vent_status"] = np.where(
                    df[self.v_vent_positive_vars].notnull().any(axis=1), 1, 0
                )

                # check if there are any "real" vent rows; if so
                if df["vent_status"].sum() > 0:
                    # print(df[df['vent_status']>0].vent_status)
                    self.flags["vent_start_time"] = df[
                        df["vent_status"] > 0
                    ].recorded_time.iloc[0]

                df = (
                    df[["recorded_time", "vent_status", "fio2"]]
                    .resample(
                        "60min",
                        on="recorded_time",
                        origin=self.event_times["start_index"],
                    )
                    .first()
                    .reindex(self.super_table_time_index)
                )

                df[["vent_status", "fio2"]] = df[["vent_status", "fio2"]].ffill(limit=6)
                self.vent_status = df.vent_status
                self.vent_fio2 = df.fio2

            # If there is a vent start, but no stop; add 6hrs to start time
            elif vent_stop.size == 0:
                # flag identifies the presence of vent rows, and start time
                self.flags["y_vent_rows"] = 1
                self.flags["y_vent_start_time"] = 1
                self.flags["vent_start_time"] = vent_start[0]
                # flags['y_vent_end_time'] is already set to zero

                df["vent_status"] = np.where(df.notnull().any(axis=1), 1, 0)
                df = (
                    df[["recorded_time", "vent_status", "fio2"]]
                    .resample(
                        "60min",
                        on="recorded_time",
                        origin=self.event_times["start_index"],
                    )
                    .first()
                    .reindex(self.super_table_time_index)
                )
                df[["vent_status", "fio2"]] = df[["vent_status", "fio2"]].ffill(limit=6)
                self.vent_status = df.vent_status
                self.vent_fio2 = df.fio2
                
            else:
                # flag identifies the presence of vent rows, and start time
                self.flags["y_vent_rows"] = 1
                self.flags["y_vent_start_time"] = 1
                self.flags["y_vent_end_time"] = 1
                self.flags["vent_start_time"] = vent_start[0]

                index = pd.Index([])
                vent_tuples = zip(vent_start, vent_stop)

                for pair in set(vent_tuples):
                    index = index.append(pd.date_range(pair[0], pair[1], freq="H"))

                vent_status = pd.DataFrame(
                    data=([1.0] * len(index)), columns=["vent_status"], index=index
                )
                # sets column to 1 if vent was on
                self.vent_status = (
                    vent_status.resample(
                        "60min", origin=self.event_times["start_index"]
                    )
                    .mean()
                    .reindex(self.super_table_time_index)
                )

                self.vent_fio2 = (
                    df[["recorded_time", "fio2"]]
                    .resample(
                        "60min",
                        on="recorded_time",
                        origin=self.event_times["start_index"],
                    )
                    .mean()
                    .reindex(self.super_table_time_index)
                )

    ########################################################################################################

    def bin_vasopressors(self):
        df = self.vasopressor_meds_PerCSN
        vas_cols = (
            self.v_vasopressor_names + self.v_vasopressor_units + ["med_order_time"]
        )
        df = df[vas_cols]

        if df.empty:
            # drop unecessary cols
            df = df.drop(columns=["med_order_time"])

            # if no vasopressers then attach index to empty df
            self.vasopressor_meds_staging = df.reindex(self.super_table_time_index)

        else:
            df = df.resample(
                "60min", on="med_order_time", origin=self.event_times["start_index"]
            ).max()
            # drop unecessary cols
            df = df.drop(columns=["med_order_time"])

            self.vasopressor_meds_staging = df.reindex(self.super_table_time_index)

    ########################################################################################################

    def make_super_table(self):
        dfs = [
            self.vitals_staging,
            self.labs_staging,
            self.gcs_staging,
            self.vent_status,
            self.vasopressor_meds_staging,
            self.bed_status,
        ]

        # merge eveything into super table
        self.super_table = reduce(
            lambda left, right: pd.merge(
                left, right, left_index=True, right_index=True
            ),
            dfs,
        )

        # if there is a vent then update fio2 with vent fio2 vals
        try:
            self.super_table.update(self.vent_fio2, overwrite=False)
        except:
            pass

    def assign_bed_location(self):
        df = self.beds_PerCSN
        # these columns have the flags for bed status
        bed_category_names = ["icu", "imc", "ed", "procedure"]
        # makes an empty dataframe
        bed_status = pd.DataFrame(columns=bed_category_names)

        for i, row in df.iterrows():
            # makes an hourly index from bed strat to bed end
            index = pd.date_range(
                row["bed_location_start"], row["bed_location_end"], freq="H"
            )

            # makes a df for a single bed with the index and bed category values
            single_bed = pd.DataFrame(
                data=np.repeat([row[bed_category_names].values], len(index), axis=0),
                columns=bed_category_names,
                index=index,
            )
            # adds all beds to single df
            bed_status = bed_status.append(single_bed)
        bed_status = bed_status[~bed_status.index.duplicated(keep="first")]

        # this is bed status re_indexed with super_table index; gets merged in later
        self.bed_status = bed_status.reindex(
            self.super_table_time_index, method="nearest"
        )

    def comorbid_dict(self, imported):
        # =============================================================================
        ### ICD10 Calcs
        # =============================================================================
        #         self.ahrq_ICD10_PerCSN = self.ahrq_ICD10_PerCSN.reset_index().groupby(['ICD10']).first().\
        #                                 groupby(['ahrq']).agg(
        #                                 icd_count = pd.NamedAgg(column="csn", aggfunc="count"),
        #                                 date_time = pd.NamedAgg(column="dx_time_date", aggfunc="first"))\
        #                                 .reindex(imported.v_ahrq_labels).rename_axis(None)
        #                                 #.agg({'csn':'count', 'dx_time_date':'first'})\
        #
        #
        #         self.elix_ICD10_PerCSN = self.elix_ICD10_PerCSN.reset_index().groupby(['ICD10']).first().\
        #                                 groupby(['elix']).agg(
        #                                 icd_count = pd.NamedAgg(column="csn", aggfunc="count"),
        #                                 date_time = pd.NamedAgg(column="dx_time_date", aggfunc="first"))\
        #                                 .reindex(imported.v_elix_labels).rename_axis(None)
        #                                 #.agg({'csn':'count', 'dx_time_date':'first'})\
        # =============================================================================

        self.quan_deyo_ICD10_PerCSN = (
            self.quan_deyo_ICD10_PerCSN.reset_index()
            .groupby(["ICD10"])
            .first()
            .groupby(["quan_deyo"])
            .agg(
                icd_count=pd.NamedAgg(column="csn", aggfunc="count"),
                date_time=pd.NamedAgg(column="dx_time_date", aggfunc="first"),
            )
            .reindex(imported.v_quan_deyo_labels)
            .rename_axis(None)
        )
        # .agg({'csn':'count', 'dx_time_date':'first'})\

        self.quan_elix_ICD10_PerCSN = (
            self.quan_elix_ICD10_PerCSN.reset_index()
            .groupby(["ICD10"])
            .first()
            .groupby(["quan_elix"])
            .agg(
                icd_count=pd.NamedAgg(column="csn", aggfunc="count"),
                date_time=pd.NamedAgg(column="dx_time_date", aggfunc="first"),
            )
            .reindex(imported.v_quan_elix_labels)
            .rename_axis(None)
        )
        # .agg({'csn':'count', 'dx_time_date':'first'})\

    # =============================================================================
    #         self.ccs_ICD10_PerCSN = self.ccs_ICD10_PerCSN.reset_index().groupby(['ICD10']).first().\
    #                                 groupby(['ccs_label']).agg(
    #                                 icd_count = pd.NamedAgg(column="csn", aggfunc="count"),
    #                                 date_time = pd.NamedAgg(column="dx_time_date", aggfunc="first"))\
    #                                 .reindex(imported.v_ccs_labels).rename_axis(None)
    #                                 #.agg({'csn':'count', 'dx_time_date':'first'})\
    # =============================================================================
    def calc_icu_stay(self):
        if self.bed_status.icu.sum() > 0:
            # mask all zeros (i.e. make nan) if there is a gap <=12hrs between ICU bed times then if fills it; otherwise it's zero
            gap_filled = (
                self.bed_status.mask(self.bed_status.icu == 0).icu.fillna(
                    method="ffill", limit=12
                )
            ) + (
                self.bed_status.mask(self.bed_status.icu == 0).icu.fillna(
                    method="bfill"
                )
                * 0
            )
            self.gap_filled = gap_filled
            # converts index into a series
            s = gap_filled.dropna().index.to_series()

            # if the delta between index vals is >1hr then mark it a start time
            start_time = s[s.diff(1) != pd.Timedelta("1 hours")].reset_index(drop=True)

            # if the reverse delta between index vals is > -1hr then mark it a end time
            end_time = s[s.diff(-1) != -pd.Timedelta("1 hours")].reset_index(drop=True)

            # makes a df with start, stop tuples
            times = pd.DataFrame(
                {"start_time": start_time, "end_time": end_time},
                columns=["start_time", "end_time"],
            )

            self.event_times["first_icu_start"] = times.iloc[0]["start_time"]

            self.event_times["first_icu_end"] = times.iloc[0]["end_time"]

        # self.event_times ['first_icu'] =  self.beds_PerCSN[self.beds_PerCSN.icu==1].sort_values('bed_location_start').bed_location_start.iloc[0]
        else:
            self.event_times["first_icu_start"] = None
            self.event_times["first_icu_end"] = None

    ########################################################################################################

    def calc_t_susp(self):
        self.abx_order_time = self.abx_staging.med_order_time.unique()

        self.culture_times = self.cultures_staging.order_time.unique()

        hours72 = pd.Timedelta(hours=72)
        hours24 = pd.Timedelta(hours=24)
        hours0 = pd.Timedelta(hours=0)

        # t_susp if t_abx is first

        sus_abx_first = [
            (abx_t, clt_t)
            for abx_t in self.abx_order_time
            for clt_t in self.culture_times
            if (clt_t - abx_t) < hours24 and (clt_t - abx_t) > hours0
        ]

        # t_susp if t_clt is first
        sus_clt_first = [
            (abx_t, clt_t)
            for clt_t in self.culture_times
            for abx_t in self.abx_order_time
            if (abx_t - clt_t) < hours72 and (abx_t - clt_t) > hours0
        ]

        t_susp_list = sus_clt_first + sus_abx_first
        t_suspicion = pd.DataFrame(t_susp_list, columns=["t_abx", "t_clt"])
        t_suspicion["t_suspicion"] = t_suspicion[["t_abx", "t_clt"]].min(axis=1)
        self.t_suspicion = t_suspicion.sort_values("t_suspicion")

    ########################################################################################################
    #### Handles missing weights & heights
    ########################################################################################################

    def fill_height_weight(self, weight_col="daily_weight_kg", height_col="height_cm"):
        """
        Accepts- a patient dictionary and names of weight and height cols.
        Does- 1) First height is back filled to admission
              2) All weights are forward filled until discharge
              3) If no recorded weight during addmisison then patient is assigned an 'average' weight based on gender.
        Returns- An updated version of super_table
        Notes:A height & weight should almost always be recorded in the first 24hrs
        """

        # define path to super_table
        df = self.super_table

        # gender 1=female & 2 = male
        gender = self.static_features["gender_code"]

        # If there is no weight or height substitue in average weight by gender
        if df[weight_col].isnull().all():
            # if pt is a male
            if gender == 2:
                df.iloc[0, df.columns.get_loc(weight_col)] = 89
                df.iloc[0, df.columns.get_loc(height_col)] = 175.3

            # if a pt is a female
            elif gender == 1:
                df.iloc[0, df.columns.get_loc(weight_col)] = 75
                df.iloc[0, df.columns.get_loc(height_col)] = 161.5

            # if a pt gender is undefined then use average of male & female
            else:
                df.iloc[0, df.columns.get_loc(weight_col)] = (89 + 75) / 2
                df.iloc[0, df.columns.get_loc(height_col)] = (175.3 + 161.5) / 2

        # Backfill to admission
        df[weight_col].loc[: df[height_col].first_valid_index()].fillna(
            method="bfill", inplace=True
        )
        df[height_col].loc[: df[height_col].first_valid_index()].fillna(
            method="bfill", inplace=True
        )

        # Fwdfill to discharge
        df[weight_col].fillna(method="ffill", inplace=True)
        df[height_col].fillna(method="ffill", inplace=True)

    ########################################################################################################
    #### Determines best MAP to use
    ########################################################################################################

    def best_map_by_row(self, row):
        """
        Accepts- A patient_dictionary and a row from super_table
        Does- 1)Reviews all blood pressure values per window (i.e. hour)
              2)selects or calculates a the most appropriate value.
        Returns- a map value in a super_table col called 'best_map'
        """

        # function that manually calculates map if it's not done already
        def calc_map(sbp, dbp):
            if pd.isna(sbp) or pd.isna(dbp):
                return float("NaN")
            else:
                return (sbp + (2 * dbp)) / 3

        if ~row[["map_line", "map_cuff"]].isnull().all():
            if abs(row["map_line"] - row["map_cuff"]) < 10:
                best_map = row["map_line"]
            else:
                best_map = row[["map_line", "map_cuff"]].max()

        elif (
            row[["sbp_line", "dbp_line"]].notnull().all()
            or row[["sbp_cuff", "dbp_cuff"]].notnull().all()
        ):
            best_map = max(
                calc_map(row["sbp_line"], row["dbp_line"]),
                calc_map(row["sbp_cuff"], row["dbp_cuff"]),
            )

        else:
            best_map = float("NaN")

        return best_map

    def best_map(
        self,
        v_bp_cols=[
            "sbp_line",
            "dbp_line",
            "map_line",
            "sbp_cuff",
            "dbp_cuff",
            "map_cuff",
        ],
    ):
        """
        Accepts- A patient_dictionary and list of blood pressure columns
        Does- Runs the best_map function for each window (i.e. row) of super_table
        Returns- An updated super_table with best map now included
        """

        # picks or calculates the best map
        self.super_table["best_map"] = self.super_table[v_bp_cols].apply(
            self.best_map_by_row, axis=1
        )

    ########################################################################################################
    #### Calculates P:F & S:F Ratio
    ########################################################################################################

    # Converts FiO2 to decimal if it is not in this form
    def fio2_decimal(self, fio2="fio2"):
        """
        Accepts- Patient dictionary & FiO2 column name
        Does- Checks to see if FiO2 is decimal, if not divides/100
        Returns- FiO2 col that is all decimals (i.e. 0.10 NOT 10%)
        """

        # small function to check if fio2 is decimal by row
        def fio2_row(row, fio2=fio2):
            if row[fio2] <= 1.0:
                return row[fio2]
            else:
                return row[fio2] / 100

        df = self.super_table
        df[fio2] = df.apply(fio2_row, axis=1)

    ########################################################################################################

    # Calculates pf ratio using SpO2 and PaO2 these P:F ratios are saved as new column
    def calc_pf(
        self, spo2="spo2", pao2="partial_pressure_of_oxygen_(pao2)", fio2="fio2"
    ):
        """
        Accepts- Patient dictionary
        Does- 1) Calculates P:F ratio using SpO2 and PaO2
        Returns- two new columns to super_table with P:F ratios
        """
        df = self.super_table

        # if df['fio2'].max() > 1:
        #    df['fio2'] = df['fio2']/100

        df["pf_sp"] = df[spo2] / df[fio2]
        df["pf_pa"] = df[pao2] / df[fio2]
        return

    def single_pressor_by_weight(self, row, single_pressors_name):
        """Accepts a row from an apply function, and a name of a pressor
        Checks the dosing rate and decides if division by weight is needed or not."""

        if single_pressors_name == "vasopressin":
            val = row[single_pressors_name]

        elif row[single_pressors_name + "_dose_unit"] == "mcg/min":
            val = row[single_pressors_name] / row["daily_weight_kg"]

        elif row[single_pressors_name + "_dose_unit"] == "mcg/kg/min":
            val = row[single_pressors_name]

        else:
            val = row[single_pressors_name]
        return val

    def calc_all_pressors(self, v_vasopressor_names=v_vasopressor_names):
        """Accepts- Patient Dictionary, List of Vasopressor names
        Does- Applies the 'single_pressor_by_weight' function to each pressor each pressor
              column, one row at a time .
        Returns- A column for each pressor that is adjusted for weight as needed."""

        df = self.super_table
        for val in v_vasopressor_names:
            df[val + "_dose_weight"] = df.apply(
                self.single_pressor_by_weight, single_pressors_name=val, axis=1
            )

    ########################################################################################################
    #### Cleans up the vasopressors
    ########################################################################################################

    def fill_values(self, labs=None, vitals=None, gcs=None):
        """
        Accepts- Patient Dictionary and list of patient features to fill
        Does- 1. Fwd fills each value for a max of 24hrs
              2. Back fills for a max of 24hrs from admission (i.e. for labs 1hr after admit)
        Returns- Patient Dictionary with filled patient features
        """
        if labs is None:
            v_all_lab_col_names = self.v_all_lab_col_names
        if vitals is None:
            v_vital_col_names = self.v_vital_col_names
        if gcs is None:
            v_gcs_col_names = self.v_gcs_col_names

        numerical_cols = v_all_lab_col_names + v_vital_col_names + v_gcs_col_names

        # Fwdfill to discharge
        self.super_table[numerical_cols] = self.super_table[numerical_cols].ffill(
            limit=24
        )
        # self.super_table[numerical_cols]=self.super_table[numerical_cols].bfill(limit=24)

    def fill_pressor_values(
        self,
        v_vasopressor_names=None,
        v_vasopressor_units=None,
        v_vasopressor_dose=None,
    ):
        """Accepts- 1) Patient Dictionary
                 2) Lists of Initial vasopressor dose, vasopressor units, vasopressor weight based dose
        Does- Forward fills from first non-null value to the last non-null value.
        Returns-
        Notes- The assumption is that the last pressor is the last dose.
        """

        # Uses class variable for function
        if v_vasopressor_names is None:
            v_vasopressor_names = self.v_vasopressor_col_names

        if v_vasopressor_units is None:
            v_vasopressor_units = self.v_vasopressor_units

        if v_vasopressor_dose is None:
            v_vasopressor_dose = self.v_vasopressor_dose

        # create super_table variable
        df = self.super_table

        # fills the value for the initial vasopressor dose
        df[v_vasopressor_names] = df[v_vasopressor_names].apply(
            lambda columns: columns.loc[: columns.last_valid_index()].ffill()
        )

        # fills the vasopressor name
        df[v_vasopressor_units] = df[v_vasopressor_units].apply(
            lambda columns: columns.loc[: columns.last_valid_index()].ffill()
        )

        # fills the weight based vasopressor dose
        df[v_vasopressor_dose] = df[v_vasopressor_dose].apply(
            lambda columns: columns.loc[: columns.last_valid_index()].ffill()
        )

    def calc_comorbidities(self):
        # calculate CCI etc. return a df
        pass

    def calc_worst_pf(self):
        df = self.super_table
        # select worse pf_pa when on vent
        self.flags["worst_pf_pa"] = df[df["vent_status"] > 0]["pf_pa"].min()
        if df[df["vent_status"] > 0]["pf_pa"].size:
            self.flags["worst_pf_pa_time"] = df[df["vent_status"] > 0]["pf_pa"].idxmin(
                axis=1, skipna=True
            )
        else:
            self.flags["worst_pf_pa_time"] = pd.NaT
        # select worse pf_sp when on vent
        self.flags["worst_pf_sp"] = df[df["vent_status"] > 0]["pf_sp"].min()
        if df[df["vent_status"] > 0]["pf_sp"].size:
            self.flags["worst_pf_sp_time"] = df[df["vent_status"] > 0]["pf_sp"].idxmin(
                axis=1, skipna=True
            )
        else:
            self.flags["worst_pf_sp_time"] = pd.NaT

    #######################################################################################################
    def make_dict_elements(self, imported):
        print("INSIDE MAKE_DICT_ELEMENTS")
        # start_dict_time = time.time()
        self.flag_dict()
        # print('flag complete')
        self.static_features_dict()
        # print('static complete')
        self.event_times_dict()
        # print('event complete')
        self.cultures_df()
        # print('cultures complete')
        self.antibiotics_df()
        # print('abx complete')
        self.build_super_table_index()
        # print('super table index complete')
        self.assign_bed_location()
        # print('the beds have been binned')
        self.bin_labs()
        # print('labs complete')
        self.bin_vitals()
        # print('vitals complete')
        self.bin_gcs()
        # print('gcs complete')
        self.bin_vent()
        # print('vent complete')
        self.bin_vasopressors()
        # print('vasopressors complete')
        self.make_super_table()
        # print('super table')
        self.calc_t_susp()
        # print('t susp complete')
        self.fill_height_weight()
        # print('filling height weight complete')
        self.best_map()
        # print('best map selected')
        self.calc_all_pressors()
        # print('vasopressor mg/kg/min calculated')
        # self.fill_values()
        # print('most values filled fwd/bwd')
        # self.fill_pressor_values()
        # print('pressor values filled fwd/bwd')
        self.fio2_decimal()
        # print('fio2 converted to decimal where appropriate')
        self.calc_pf()
        # print('all p:f and s:f ratios calculated')
        self.comorbid_dict(imported)
        # print('the comorbid dictionary is updated')
        self.calc_icu_stay()
        # print('the first icu stay has been calculated')
        self.calc_worst_pf()
        # print('the worst pf has been calc'd and saved')

        # print(f'It took {time.time()-start_dict_time}(s) to process the whole dictionary.')

    def write_dict(self):
        encounter_dict = {}
        encounter_dict["csn"] = self.csn
        encounter_dict["pat_id"] = self.pat_id
        encounter_dict["cultures_PerCSN"] = self.cultures_PerCSN
        encounter_dict["beds_PerCSN"] = self.beds_PerCSN
        encounter_dict["procedures_PerCSN"] = self.procedures_PerCSN
        encounter_dict["diagnosis_PerCSN"] = self.diagnosis_PerCSN
        encounter_dict["flags"] = self.flags
        encounter_dict["static_features"] = self.static_features
        encounter_dict["event_times"] = self.event_times
        encounter_dict["cultures_staging"] = self.cultures_staging
        encounter_dict["abx_staging"] = self.abx_staging
        encounter_dict["labs_staging"] = self.labs_staging
        encounter_dict["vitals_staging"] = self.vitals_staging
        encounter_dict["gcs_staging"] = self.gcs_staging
        encounter_dict["vent_status"] = self.vent_status
        encounter_dict["vent_fio2"] = self.vent_fio2
        encounter_dict["vasopressor_meds_staging"] = self.vasopressor_meds_staging
        encounter_dict["super_table"] = self.super_table

        # Suspicion
        encounter_dict["abx_order_time"] = self.abx_order_time
        encounter_dict["culture_times"] = self.culture_times
        encounter_dict["t_suspicion"] = self.t_suspicion

        # Sepsis 3
        encounter_dict["sofa_scores"] = self.sofa_scores
        encounter_dict["sep3_time"] = self.sep3_time

        # Sepsis 2
        encounter_dict["sirs_scores"] = self.sirs_scores
        encounter_dict["sep2_time"] = self.sep2_time

        # write to the instance
        self.encounter_dict = encounter_dict

    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    ###
    ###             BEGIN SEPSIS CALCS
    ###
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################

    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    ####################        Calc SEPSIS 3                                            ###################
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################

    ########################################################################################################
    ####################                  SOFA FUNCTIONS                                 ###################
    ########################################################################################################

    def SOFA_resp(self, row, pf_pa="pf_pa", pf_sp="pf_sp"):
        """Accepts- class instance, one row from "super_table", "pf" cols
        Does- Calculates Respiratory SOFA score
        Returns- Single value of Respiratory SOFA score"""
        if row[pf_pa] < 100:
            val = 4
        elif row[pf_pa] < 200:
            val = 3
        elif row[pf_pa] < 300:
            val = 2
        elif row[pf_pa] < 400:
            val = 1
        elif row[pf_pa] >= 400:
            val = 0
        else:
            val = float("NaN")
        return val

    def SOFA_resp_sa(self, row, pf_pa="pf_pa", pf_sp="pf_sp"):
        """Accepts- class instance, one row from "super_table", "pf" cols
        Does- Calculates Respiratory SOFA score
        Returns- Single value of Respiratory SOFA score"""
        if row[pf_sp] < 67:
            val = 4
        elif row[pf_sp] < 142:
            val = 3
        elif row[pf_sp] < 221:
            val = 2
        elif row[pf_sp] < 302:
            val = 1
        elif row[pf_sp] >= 302:
            val = 0
        else:
            val = float("NaN")
        return val

    def SOFA_cardio(
        self,
        row,
        dopamine_dose_weight="dopamine_dose_weight",
        epinephrine_dose_weight="epinephrine_dose_weight",
        norepinephrine_dose_weight="norepinephrine_dose_weight",
        dobutamine_dose_weight="dobutamine_dose_weight",
    ):
        """
        Accepts- class instance, one row from "super_table", weight based pressor cols
        Does- Calculates Cardio SOFA score
        Returns- Single value of Cardio SOFA score
        """

        if (
            (row[dopamine_dose_weight] > 15)
            | (row[epinephrine_dose_weight] > 0.1)
            | (row[norepinephrine_dose_weight] > 0.1)
        ):
            val = 4
        elif (
            (row[dopamine_dose_weight] > 5)
            | (
                (row[epinephrine_dose_weight] > 0.0)
                & (row[epinephrine_dose_weight] <= 0.1)
            )
            | (
                (row[norepinephrine_dose_weight] > 0.0)
                & (row[norepinephrine_dose_weight] <= 0.1)
            )
        ):
            val = 3
        elif ((row[dopamine_dose_weight] > 0.0) & (row[dopamine_dose_weight] <= 5)) | (
            row[dobutamine_dose_weight] > 0
        ):
            val = 2
        elif row["best_map"] < 70:
            val = 1

        elif row["best_map"] >= 70:
            val = 0
        else:
            val = float("NaN")
        return val

    def SOFA_coag(self, row):
        if row["platelets"] >= 150:
            val = 0
        elif (row["platelets"] >= 100) & (row["platelets"] < 150):
            val = 1
        elif (row["platelets"] >= 50) & (row["platelets"] < 100):
            val = 2
        elif (row["platelets"] >= 20) & (row["platelets"] < 50):
            val = 3
        elif row["platelets"] < 20:
            val = 4
        else:
            val = float("NaN")
        return val

    def SOFA_neuro(self, row):
        if row["gcs_total_score"] == 15:
            val = 0
        elif (row["gcs_total_score"] >= 13) & (row["gcs_total_score"] <= 14):
            val = 1
        elif (row["gcs_total_score"] >= 10) & (row["gcs_total_score"] <= 12):
            val = 2
        elif (row["gcs_total_score"] >= 6) & (row["gcs_total_score"] <= 9):
            val = 3
        elif row["gcs_total_score"] < 6:
            val = 4
        else:
            val = float("NaN")
        return val

    def SOFA_hep(self, row):
        if row["bilirubin_total"] < 1.2:
            val = 0
        elif (row["bilirubin_total"] >= 1.2) & (row["bilirubin_total"] < 2.0):
            val = 1
        elif (row["bilirubin_total"] >= 2.0) & (row["bilirubin_total"] < 6.0):
            val = 2
        elif (row["bilirubin_total"] >= 6.0) & (row["bilirubin_total"] < 12.0):
            val = 3
        elif row["bilirubin_total"] >= 12.0:
            val = 4
        else:
            val = float("NaN")
        return val

    def SOFA_renal(self, row):
        if row["creatinine"] < 1.2:
            val = 0
        elif (row["creatinine"] >= 1.2) & (row["creatinine"] < 2.0):
            val = 1
        elif (row["creatinine"] >= 2.0) & (row["creatinine"] < 3.5):
            val = 2
        elif (row["creatinine"] >= 3.5) & (row["creatinine"] < 5.0):
            val = 3
        elif row["creatinine"] >= 5.0:
            val = 4
        else:
            val = float("NaN")
        return val

    def SOFA_cardio_mod(
        self,
        row,
        dopamine_dose_weight="dopamine_dose_weight",
        epinephrine_dose_weight="epinephrine_dose_weight",
        norepinephrine_dose_weight="norepinephrine_dose_weight",
        dobutamine_dose_weight="dobutamine_dose_weight",
    ):
        """
        Accepts- class instance, one row from "super_table", weight based pressor cols
        Does- Calculates Cardio SOFA score
        Returns- Single value of Cardio SOFA score
        """

        if (row[epinephrine_dose_weight] > 0.0) & (row[epinephrine_dose_weight] > 0.0):
            val = 4
        elif (row[epinephrine_dose_weight] > 0.0) | (
            row[epinephrine_dose_weight] > 0.0
        ):
            val = 3
        elif (row[dopamine_dose_weight] > 0.0) | (row[dobutamine_dose_weight] > 0):
            val = 2
        elif row["best_map"] < 70:
            val = 1
        elif row["best_map"] >= 70:
            val = 0
        else:
            val = float("NaN")
        return val

    def calc_all_SOFA(self, window=24):
        df = self.super_table
        sofa_df = pd.DataFrame(
            index=self.super_table.index,
            columns=[
                "SOFA_coag",
                "SOFA_renal",
                "SOFA_hep",
                "SOFA_neuro",
                "SOFA_cardio",
                "SOFA_cardio_mod",
                "SOFA_resp",
                "SOFA_resp_sa",
            ],
        )

        sofa_df["SOFA_coag"] = df.apply(self.SOFA_coag, axis=1)
        sofa_df["SOFA_renal"] = df.apply(self.SOFA_renal, axis=1)
        sofa_df["SOFA_hep"] = df.apply(self.SOFA_hep, axis=1)
        sofa_df["SOFA_neuro"] = df.apply(self.SOFA_neuro, axis=1)
        sofa_df["SOFA_cardio"] = df.apply(self.SOFA_cardio, axis=1)
        sofa_df["SOFA_cardio_mod"] = df.apply(self.SOFA_cardio_mod, axis=1)
        sofa_df["SOFA_resp"] = df.apply(self.SOFA_resp, axis=1)
        sofa_df["SOFA_resp_sa"] = df.apply(self.SOFA_resp_sa, axis=1)
        ######## Normal Calcs
        # Calculate NOMRAL hourly totals for each row
        sofa_df["hourly_total"] = sofa_df[
            [
                "SOFA_coag",
                "SOFA_renal",
                "SOFA_hep",
                "SOFA_neuro",
                "SOFA_cardio",
                "SOFA_resp",
            ]
        ].sum(axis=1)

        # Calculate POST 24hr delta in total SOFA Score
        sofa_df["delta_24h"] = (
            sofa_df["hourly_total"]
            .rolling(window=window, min_periods=24)
            .apply(
                lambda x: x.max() - x.min()
                if x.idxmax().value > x.idxmin().value
                else 0
            )
            .tolist()
        )

        # Calculate FIRST 24h delta in total SOFA score
        sofa_df.update(
            sofa_df.loc[sofa_df.index[0:24], ["hourly_total"]]
            .rolling(window=window, min_periods=1)
            .max()
            .rename(columns={"hourly_total": "delta_24h"})
        )

        ######## Modified Calcs
        # Calculate NOMRAL hourly totals for each row
        sofa_df["hourly_total_mod"] = sofa_df[
            [
                "SOFA_coag",
                "SOFA_renal",
                "SOFA_hep",
                "SOFA_neuro",
                "SOFA_cardio_mod",
                "SOFA_resp_sa",
            ]
        ].sum(axis=1)

        # Calculate POST 24hr delta in total SOFA Score
        sofa_df["delta_24h_mod"] = (
            sofa_df["hourly_total_mod"]
            .rolling(window=window, min_periods=24)
            .apply(
                lambda x: x.max() - x.min()
                if x.idxmax().value > x.idxmin().value
                else 0
            )
            .tolist()
        )

        # Calculate FIRST 24h delta in total SOFA score
        sofa_df.update(
            sofa_df.loc[sofa_df.index[0:24], ["hourly_total_mod"]]
            .rolling(window=window, min_periods=1)
            .max()
            .rename(columns={"hourly_total_mod": "delta_24h_mod"})
        )

        # Safe this dataframe into the patient dictionary
        self.sofa_scores = sofa_df

    ####
    #### Hourly Max SOFA IS ON TIME OUTneeds some more attention 5/6/21
    ####
    # =============================================================================
    #     def hourly_max_SOFA(self,
    #                         window = 24):
    #         df = self.sofa_scores
    #
    #         df['SOFA_coag_24h_max'] = df['SOFA_coag'].rolling(window=window, min_periods=2).max().tolist()
    #         df['SOFA_renal_24h_max'] = df['SOFA_renal'].rolling(window=window, min_periods=2).max().tolist()
    #         df['SOFA_hep_24h_max'] = df['SOFA_hep'].rolling(window=window, min_periods=2).max().tolist()
    #         df['SOFA_neuro_24h_max'] = df['SOFA_neuro'].rolling(window=window, min_periods=2).max().tolist()
    #         df['SOFA_cardio_24h_max'] = df['SOFA_cardio'].rolling(window=window, min_periods=2).max().tolist()
    #         df['SOFA_resp_24h_max'] = df['SOFA_resp'].rolling(window=window, min_periods=2).max().tolist()
    #
    #         # hourly sum considering worst in 24hrs
    #         df['hourly_total_24h_max'] = (df[['SOFA_coag_24h_max',
    #                                    'SOFA_renal_24h_max',
    #                                    'SOFA_hep_24h_max',
    #                                    'SOFA_neuro_24h_max',
    #                                    'SOFA_cardio_24h_max',
    #                                    'SOFA_resp_24h_max']].sum(axis=1))
    #         #
    #         df['delta_24h_24h_max'] = df['hourly_total_24h_max'].\
    #         rolling(window=window, min_periods=24).\
    #         apply(lambda x: x.max() - x.min() if x.idxmax().value> x.idxmin().value else 0 ).tolist()
    #
    #         # Set values in first row to zero if NaN
    #         df.iloc[0,:] = df.iloc[0,].fillna(0)
    # =============================================================================

    ########################################################################################################
    ####################        Run all The Sepsis 3 steps                               ###################
    ########################################################################################################

    def run_SEP3(self):
        """Accepts- a SOFAPrep class instance
        Does- Runs all the prep and calc steps for SOFA score calculation
        Returns- A class instance with updated "super_table" and new "sofa_scores" data frame
        """
        # start_sofa_calc = time.time()
        self.calc_all_SOFA()
        # self.hourly_max_SOFA ()
        self.calc_sep3_time()
        self.calc_sep3_time_mod()

        ####Set first sepsis 3 time in the flag dictionary
        # Select the first row that has 3x values
        df = self.sep3_time[self.sep3_time.notna().all(axis=1)].reset_index()
        if df.empty:
            print("No sep3 times to add to flag dict")
            self.flags["first_sep3_susp"] = None
            self.flags["first_sep3_SOFA"] = None
            self.flags["first_sep3_time"] = None
        else:
            print("adding first sep3 times to flag dict")
            self.flags["first_sep3_susp"] = df["t_suspicion"][0]
            self.flags["first_sep3_SOFA"] = df["t_SOFA"][0]
            self.flags["first_sep3_time"] = df["t_sepsis3"][0]

            self.calc_sep3_time_mod()

        # Set first sepsis 3 time in the flag dictionary
        df = self.sep3_time_mod[self.sep3_time_mod.notna().all(axis=1)].reset_index()
        if df.empty:
            print("No sep3_mod times to add to flag dict")
            self.flags["first_sep3_susp_mod"] = None
            self.flags["first_sep3_SOFA_mod"] = None
            self.flags["first_sep3_time_mod"] = None
        else:
            print("adding first sep3_mod times to flag dict")
            self.flags["first_sep3_susp_mod"] = df["t_suspicion"][0]
            self.flags["first_sep3_SOFA_mod"] = df["t_SOFA_mod"][0]
            self.flags["first_sep3_time_mod"] = df["t_sepsis3_mod"][0]

    ########################################################################################################
    ####################        Calc Tsepsis-3                                           ###################
    ########################################################################################################
    def calc_sep3_time(self, look_back=24, look_forward=12):
        # Initialize empty list to hold SOFA times in loops below
        # t_SOFA_list = []

        # Initialize empty df to hold suspicion and sofa times
        sep3_time_df = pd.DataFrame(columns=["t_suspicion", "t_SOFA"])

        # get suspicion times from class
        suspicion_times = (
            self.t_suspicion["t_suspicion"].sort_values().drop_duplicates()
        )

        #### if NO SUSPICION, then get all SOFA >2
        if suspicion_times.empty:
            df = self.sofa_scores
            # get index of times when total change is >= 2
            sofa_times = df[df["hourly_total"] >= 2].index

            if sofa_times.empty:
                pass

            else:
                sofa_times = sofa_times.tolist()[0]

        #### If SUSPICION time is present
        else:
            sofa_times = []
            for suspicion_time in suspicion_times:
                # look back portion of window (i.e. 24hrs before Tsuspicion)
                start_window_time = suspicion_time - pd.Timedelta(hours=look_back)

                # look forward portion of window (i.e. 12hrs after Tsuspicion)
                end_window_time = suspicion_time + pd.Timedelta(hours=look_forward)

                # get all SOFA that had a 2pt change in last 24hrs (this is calculated in SOFA table)
                potential_sofa_times = self.sofa_scores[
                    self.sofa_scores["delta_24h"] >= 2
                ]

                # keep times that are with in a suspicion window
                potential_sofa_times = potential_sofa_times.loc[
                    start_window_time:end_window_time
                ].index.tolist()
                # print("These are potential SOFA Times: {}".format(potential_sofa_times))

                if not potential_sofa_times:
                    sofa_times.append(float("NaN"))
                    # print ("A NaN was appended")
                else:
                    sofa_times.append(potential_sofa_times[0])
                    # print("This SOFA Score was appended: {}".format(potential_sofa_times[0]))

        # this adds Tsofa and Tsusp and picks the min; it's the most basic Tsep calculator
        sep3_time_df["t_suspicion"] = suspicion_times.tolist()
        sep3_time_df["t_SOFA"] = sofa_times
        sep3_time_df["t_sepsis3"] = sep3_time_df.min(axis=1, skipna=False)

        # This adds all the Tsofas that did not become part of a Tsepsis tuple; probably unecessary
        # all_sofa_times = self.sofa_scores[self.sofa_scores['delta_24h'] >= 2].reset_index()
        # sep3_time_df = all_sofa_times['index'].to_frame().merge(sep3_time_df, how='outer', left_on='index',right_on='t_SOFA')
        # sep3_time_df = sep3_time_df.iloc[sep3_time_df['index'].fillna(sep3_time_df['t_suspicion']).argsort()].reset_index(drop=True).drop(columns=['t_SOFA']).rename(columns={'index':'t_SOFA'})

        self.sep3_time = sep3_time_df

    ########################################################################################################
    ####################        Calc Tsepsis-3      MOD                                  ###################
    ########################################################################################################
    def calc_sep3_time_mod(self, look_back=24, look_forward=12):
        # Initialize empty list to hold SOFA times in loops below
        # t_SOFA_list = []

        # Initialize empty df to hold suspicion and sofa times
        sep3_time_df_mod = pd.DataFrame(columns=["t_suspicion", "t_SOFA_mod"])

        # get suspicion times from class
        suspicion_times = (
            self.t_suspicion["t_suspicion"].sort_values().drop_duplicates()
        )

        #### if NO SUSPICION, then get  first SOFA >2
        if suspicion_times.empty:
            df = self.sofa_scores
            # get index of times when total change is >= 2
            sofa_times_mod = df[df["hourly_total_mod"] >= 2].index

            if sofa_times_mod.empty:
                pass

            else:
                sofa_times_mod = sofa_times_mod.tolist()[0]

        #### If SUSPICION time is present
        else:
            sofa_times_mod = []
            for suspicion_time in suspicion_times:
                # look back portion of window (i.e. 24hrs before Tsuspicion)
                start_window_time = suspicion_time - pd.Timedelta(hours=look_back)

                # look forward portion of window (i.e. 12hrs after Tsuspicion)
                end_window_time = suspicion_time + pd.Timedelta(hours=look_forward)

                # =============================================================================
                #                 #hourly SOFA score df windowed to relevant times
                #                 df = self.sofa_scores.loc[start_window_time:end_window_time]
                #
                #                 #Establish SOFA baseline for the windowget first SOFA score
                #                 if start_window_time <= self.event_times['start_index']:
                #                     baseline = 0
                #                 else:
                #                     baseline = df['hourly_total'][0]
                #
                # =============================================================================
                potential_sofa_times_mod = self.sofa_scores[
                    self.sofa_scores["delta_24h_mod"] >= 2
                ].index.tolist()
                # print("These are potential SOFA Times: {}".format(potential_sofa_times))

                if not potential_sofa_times_mod:
                    sofa_times_mod.append(float("NaN"))
                    # print ("A NaN was appended")
                else:
                    sofa_times_mod.append(potential_sofa_times_mod[0])
                    # print("This SOFA Score was appended: {}".format(potential_sofa_times[0]))

        sep3_time_df_mod["t_suspicion"] = suspicion_times.tolist()
        sep3_time_df_mod["t_SOFA_mod"] = sofa_times_mod
        sep3_time_df_mod["t_sepsis3_mod"] = sep3_time_df_mod.min(axis=1, skipna=False)

        all_sofa_times_mod = self.sofa_scores[
            self.sofa_scores["delta_24h_mod"] >= 2
        ].reset_index()
        sep3_time_df_mod = (
            all_sofa_times_mod["index"]
            .to_frame()
            .merge(
                sep3_time_df_mod, how="outer", left_on="index", right_on="t_SOFA_mod"
            )
        )
        sep3_time_df_mod = (
            sep3_time_df_mod.iloc[
                sep3_time_df_mod["index"]
                .fillna(sep3_time_df_mod["t_suspicion"])
                .argsort()
            ]
            .reset_index(drop=True)
            .drop(columns=["t_SOFA_mod"])
            .rename(columns={"index": "t_SOFA_mod"})
        )

        self.sep3_time_mod = sep3_time_df_mod

    ########################################################################################################
    ########################################################################################################
    ########################################################################################################
    ####################        Calc SEPSIS 2                                            ###################
    ########################################################################################################
    ########################################################################################################
    ########################################################################################################

    ########################################################################################################
    ####################        Calc SIRS Score                                          ###################
    ########################################################################################################

    def SIRS_resp(
        self,
        row,
        resp_rate="unassisted_resp_rate",
        paco2="partial_pressure_of_carbon_dioxide_(paco2)",
    ):
        """Accepts- class instance, one row from "super_table", "resp" cols
        Does- Calculates Respiratory SIRS score
        Returns- Single value of Respiratory SIRS score"""
        if row[resp_rate] > 20:
            val = 1
        elif row[paco2] < 32:
            val = 1
        else:
            val = 0
        return val

    def SIRS_cardio(self, row, hr="pulse"):
        """Accepts- class instance, one row from "super_table", "hr" cols
        Does- Calculates Cardiac SIRS score
        Returns- Single value of Cardiac SIRS score"""
        if row[hr] > 90:
            val = 1
        else:
            val = 0
        return val

    def SIRS_temp(self, row, temp="temperature"):
        """Accepts- class instance, one row from "super_table", "temp" cols
        Does- Calculates Temp SIRS score
        Returns- Single value of Temp SIRS score"""
        if row[temp] > 100.4:
            val = 1
        elif row[temp] < 95.8:
            val = 1
        # =============================================================================
        #         ### For Celcius
        #         if row[temp] > 38.0:
        #             val = 1
        #         elif row[temp] < 36.0:
        #             val = 1
        # =============================================================================
        else:
            val = 0
        return val

    def SIRS_wbc(self, row, wbc="white_blood_cell_count"):
        """Accepts- class instance, one row from "super_table", "wbc" cols
        Does- Calculates White Blood Cell Count SIRS score
        Returns- Single value of White Blood Cell Count SIRS score"""
        if row[wbc] > 12.0:
            val = 1
        elif row[wbc] < 4.0:
            val = 1
        # =============================================================================
        #         ## for bands
        #         if row[bands] > 10:
        #             val = 1
        # =============================================================================
        else:
            val = 0
        return val

    def calc_all_SIRS(self, window=24):
        df = self.super_table
        sirs_df = pd.DataFrame(
            index=self.super_table.index,
            columns=["SIRS_resp", "SIRS_cardio", "SIRS_temp", "SIRS_wbc"],
        )

        sirs_df["SIRS_resp"] = df.apply(self.SIRS_resp, axis=1)
        sirs_df["SIRS_cardio"] = df.apply(self.SIRS_cardio, axis=1)
        sirs_df["SIRS_temp"] = df.apply(self.SIRS_temp, axis=1)
        sirs_df["SIRS_wbc"] = df.apply(self.SIRS_wbc, axis=1)

        # Calculate hourly totals for each row
        sirs_df["hourly_total"] = sirs_df.sum(axis=1)

        # Calculate POST 24hr delta in total SIRS Score
        sirs_df["delta_24h"] = (
            sirs_df["hourly_total"]
            .rolling(window=window, min_periods=24)
            .apply(
                lambda x: x.max() - x.min()
                if x.idxmax().value > x.idxmin().value
                else 0
            )
            .tolist()
        )

        # Calculate FIRST 24h delat in total SOFA score
        sirs_df.update(
            sirs_df.loc[sirs_df.index[0:24], ["hourly_total"]]
            .rolling(window=window, min_periods=1)
            .max()
            .rename(columns={"hourly_total": "delta_24h"})
        )

        # Safe this dataframe into the patient dictionary
        self.sirs_scores = sirs_df

    ########################################################################################################
    ####################        Calc Tsepsis-2                                           ###################
    ########################################################################################################

    def calc_sep2_time(self, look_back=24, look_forward=12):
        # Initialize empty df to hold suspicion and SIRS times
        sep2_time_df = pd.DataFrame(columns=["t_suspicion", "t_SIRS"])

        # get suspicion times from class object
        suspicion_times = (
            self.t_suspicion["t_suspicion"].sort_values().drop_duplicates()
        )

        #### if NO SUSPICION, then get all SIRS >2
        if suspicion_times.empty:
            df = self.sirs_scores

            # get index of times when total change is >= 2
            sirs_times = df[df["delta_24h"] >= 2].index

            if sirs_times.empty:
                pass

            else:
                sirs_times = sirs_times.tolist()

        #### If SUSPICION time is present
        else:
            sirs_times = []
            for suspicion_time in suspicion_times:
                # look back portion of window (i.e. 24hrs before Tsuspicion)
                start_window_time = suspicion_time - pd.Timedelta(hours=look_back)

                # look forward portion of window (i.e. 12hrs after Tsuspicion)
                end_window_time = suspicion_time + pd.Timedelta(hours=look_forward)

                potential_sirs_times = self.sirs_scores[
                    self.sirs_scores["delta_24h"] >= 2
                ].index.tolist()

                if not potential_sirs_times:
                    sirs_times.append(float("NaN"))
                    # print ("A NaN was appended")
                else:
                    sirs_times.append(potential_sirs_times[0])
                    # print("This SIRS Score was appended: {}".format(potential_sirs_times[0]))

        sep2_time_df["t_suspicion"] = suspicion_times.tolist()
        sep2_time_df["t_SIRS"] = sirs_times
        sep2_time_df["t_sepsis2"] = sep2_time_df.min(axis=1, skipna=False)

        all_sirs_times = self.sirs_scores[
            self.sirs_scores["delta_24h"] >= 2
        ].reset_index()
        sep2_time_df = (
            all_sirs_times["index"]
            .to_frame()
            .merge(sep2_time_df, how="outer", left_on="index", right_on="t_SIRS")
        )
        sep2_time_df = (
            sep2_time_df.iloc[
                sep2_time_df["t_suspicion"].fillna(sep2_time_df["index"]).argsort()
            ]
            .reset_index(drop=True)
            .drop(columns=["t_SIRS"])
            .rename(columns={"index": "t_SIRS"})
        )

        self.sep2_time = sep2_time_df

    ########################################################################################################
    ####################        Run all The Sepsis 2 steps                               ###################
    ########################################################################################################

    def run_SEP2(self):
        """Accepts- a SOFAPrep class instance
        Does- Runs all the prep and calc steps for SOFA score calculation
        Returns- A class instance with updated "super_table" and new "sofa_scores" data frame
        """
        # start_SEP2_calc = time.time()
        self.calc_all_SIRS()
        self.calc_sep2_time()

        # Set first sepsis 3 time in the flag dictionary
        df = self.sep2_time[self.sep2_time.notna().all(axis=1)].reset_index()
        if df.empty:
            self.flags["first_sep2_susp"] = None
            self.flags["first_sep2_SIRS"] = None
            self.flags["first_sep2_time"] = None
        else:
            self.flags["first_sep2_susp"] = df["t_suspicion"][0]
            self.flags["first_sep2_SIRS"] = df["t_SIRS"][0]
            self.flags["first_sep2_time"] = df["t_sepsis2"][0]

        # print(f'It took {time.time()-start_SEP2_calc}(s) to calculate Sepsis-2 scores.')
