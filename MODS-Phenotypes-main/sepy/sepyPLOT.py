# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 21:00:48 2021

@author: Elite Data Hacks
"""

import sys
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def plot_scores(pt_dictionary, output_path):
    # Set style type
    plt.style.use("seaborn")
    plt.style.context("paper")

    # identify the sofa_score data frame as df
    df_sofa = pt_dictionary["sofa_scores"]
    df_sep3 = pt_dictionary["sep3_time"]
    csn = str(pt_dictionary["csn"])
    flags = {}

    ## Calculate LOS Fraction
    discharge_time = pt_dictionary["event_times"]["hospital_discharge_date_time"]
    admit_time = pt_dictionary["event_times"]["start_index"]
    los_fraction = 0.005 * (discharge_time - admit_time)

    # Set up the figure w/subplots and titles, etc
    num_subplots = 8  # int(len(df_sofa.columns)/2)
    fig, axes = plt.subplots(num_subplots, 1, figsize=(15, 30), sharex=True)

    fig.suptitle(("Hourly SOFA Scores for Patient CSN: " + csn), fontsize=20, y=1.0)
    plt.tight_layout()
    plt.setp(axes, ylim=(0, 10))

    ############### SEP-3 ###############
    # Confirm that there is a Sep 3 Time; if none then set flag dict
    if df_sep3["t_sepsis3"].isnull().all():
        sep_flag = 0
        flags["first_sep3_time"] = None

        ############### SUSPICION ###############
        # confirm that there is a Suspicion time; if none then set flag dict
        if df_sep3["t_suspicion"].isnull().all():
            sus_flag = 0
            flags["first_sep3_susp"] = discharge_time
            # debug text
            print("No Suspicion Times Recorded; Sus now = discharge time")

            ############### SOFA ###############
            # confirm that there is a SOFA time
            if df_sep3["t_SOFA"].isnull().all():
                sofa_flag = 0
                # Since no Susp time; set it to admit_time
                flags["first_sep3_SOFA"] = admit_time
                # debug text
                print("No SOFA Time recorded.")

            else:
                sofa_flag = 1
                flags["first_sep3_SOFA"] = df_sep3["t_SOFA"].dropna().iloc[0]
                # debug txt
                print(f"First SOFA:{format(flags['first_sep3_SOFA'])}")

        else:
            sus_flag = 1
            flags["first_sep3_susp"] = df_sep3["t_suspicion"].dropna().iloc[0]
            # debug text
            print(f"First Sus:{format(flags['first_sep3_susp'])}")

            ############### SOFA ###############
            # confirm that there is a SOFA time
            if df_sep3["t_SOFA"].isnull().all():
                sofa_flag = 0
                # Since no Susp time; set it to admit_time
                flags["first_sep3_SOFA"] = admit_time
                # debug text
                print("No SOFA Time recorded.")

            else:
                sofa_flag = 1
                flags["first_sep3_SOFA"] = df_sep3["t_SOFA"].dropna().iloc[0]
                # debug txt
                print(f"First SOFA:{format(flags['first_sep3_SOFA'])}")

    # If there is a sep3 time set all attributes (Tsofa, Tsepsis, etc.)
    else:
        sep_flag = 1
        df = df_sep3[df_sep3.notna().all(axis=1)].reset_index()
        flags["first_sep3_susp"] = df["t_suspicion"][0]
        flags["first_sep3_SOFA"] = df["t_SOFA"][0]
        flags["first_sep3_time"] = df["t_sepsis3"][0]

    ############### Plot SOFA Lines ###############
    df_sofa.plot(y=["SOFA_coag"], style=".-", ax=axes[0], title="SOFA Coag")
    df_sofa.plot(y=["SOFA_renal"], style=".-", ax=axes[1], title="SOFA Renal")
    df_sofa.plot(y=["SOFA_hep"], style=".-", ax=axes[2], title="SOFA Hep")
    df_sofa.plot(y=["SOFA_neuro"], style=".-", ax=axes[3], title="SOFA Neuro")
    df_sofa.plot(y=["SOFA_resp"], style=".-", ax=axes[4], title="SOFA Resp")
    df_sofa.plot(
        y=["SOFA_cardio", "SOFA_cardio_mod"],
        style=".-",
        ax=axes[5],
        title="SOFA cardio",
    )
    df_sofa.plot(
        y=["hourly_total", "hourly_total_mod"],
        style=".-",
        kind="line",
        ax=axes[6],
        title="Hourly Totals",
    )
    df_sofa.plot(
        y=["delta_24h", "delta_24h_mod"], style=".-", ax=axes[7], title="24h Delta"
    )

    # This is the old code to plot max totals side by side with base elements
    # =============================================================================
    #     df_sofa.plot(y=['SOFA_coag', 'SOFA_coag_24h_max'], kind = 'line', ax=axes[0], title="SOFA Coag")
    #     df_sofa.plot(y=['SOFA_renal', 'SOFA_renal_24h_max'], kind = 'line', ax=axes[1], title="SOFA Renal")
    #     df_sofa.plot(y=['SOFA_hep', 'SOFA_hep_24h_max'], kind = 'line', ax=axes[2], title="SOFA Hep")
    #     df_sofa.plot(y=['SOFA_neuro', 'SOFA_neuro_24h_max'], kind = 'line', ax=axes[3], title="SOFA Neuro")
    #     df_sofa.plot(y=['SOFA_cardio', 'SOFA_cardio_24h_max'], kind = 'line', ax=axes[4], title="SOFA cardio")
    #     df_sofa.plot(y=['SOFA_resp', 'SOFA_resp_24h_max'], kind = 'line', ax=axes[5], title="SOFA Resp")
    #     df_sofa.plot(y=['hourly_total','delta_24h'], kind = 'line', ax=axes[6], title="Hourly Totals")
    #     df_sofa.plot(y=['hourly_total_24h_max', 'delta_24h_24h_max'], kind = 'line', ax=axes[7], title="Hourly Total w/ Worst")
    # =============================================================================

    ############### Create Labels ###############
    #     txtkw_sep = dict(size=14, color = 'r')
    #     txtkw_not_sep = dict(size=14, color = 'y')
    # =============================================================================
    #         flags['first_sep3_susp']
    #         flags['first_sep3_SOFA']
    #         flags['first_sep3_time']
    # =============================================================================

    txt_sus = "Time Suspicion: {}".format(flags["first_sep3_susp"].strftime("%H:%M"))
    txt_sofa = "Time SOFA: {}".format(flags["first_sep3_SOFA"].strftime("%H:%M"))

    first_sus = flags["first_sep3_susp"]
    first_sofa = flags["first_sep3_SOFA"]

    ############### If there is Sepsis Time then do this  ###############

    if sep_flag == 1:
        ############### Tsofa BEFORE Tsuspicion ###############
        if first_sofa < first_sus:
            left = first_sofa
            left_text = "Time SEPSIS: {}".format(first_sofa.strftime("%H:%M"))

            right = first_sus
            right_text = txt_sus

        ###############  Tsuspicion BEFORE Tsofa ###############
        else:
            left = first_sus
            left_text = "Time SEPSIS: {}".format(first_sus.strftime("%H:%M"))

            right = first_sofa
            right_text = txt_sofa

        for i in range(num_subplots):
            # makes the horizontal 2pt SOFA line
            axes[i].axhline(2, color="red", ls="--", alpha=0.3)

            # makes the vertical left line, which is sepsis and is red
            axes[i].axvline(left, color="r", ls=":")

            # makes red, left text 9 hours before vert line
            axes[i].text(
                (left - los_fraction),
                2.5,
                left_text,
                rotation=270,
                horizontalalignment="right",
                size=14,
                color="r",
            )

            # makes the vertical right line, which is T_Sofa or T_Sepsis
            axes[i].axvline(right, color="y", ls=":")

            # makes the left text
            axes[i].text(
                (right + los_fraction),
                2.5,
                right_text,
                rotation=90,
                horizontalalignment="left",
                size=14,
                color="y",
            )

    ############### If there is NO Sepsis Event ###############
    else:
        if sus_flag == 1 and sofa_flag == 1:
            pass  # labels already set, no need to change, keep this condition in case

        elif sus_flag == 1 and sofa_flag == 0:
            txt_sofa = "NO SOFA TIME"

        elif sus_flag == 0 and sofa_flag == 1:
            txt_sus = "NO SUSPICION TIME"
        else:
            txt_sofa = "NO SOFA TIME"
            txt_sus = "NO SUSPICION TIME"

        for i in range(num_subplots):
            # makes the horizontal 2pt SOFA line
            axes[i].axhline(2, color="red", ls="--", alpha=0.3)

            # No Sepsis, so plot SOFA line if there is one
            axes[i].axvline(first_sofa, color="y", ls=":")
            axes[i].text(
                (first_sofa - los_fraction),
                2.5,
                txt_sofa,
                rotation=270,
                horizontalalignment="right",
                size=14,
                color="y",
            )

            axes[i].axvline(first_sus, color="y", ls=":")
            axes[i].text(
                (first_sus + los_fraction),
                2.5,
                txt_sus,
                rotation=90,
                horizontalalignment="left",
                size=14,
                color="y",
            )
            ### Box Indicating NO sepsis
            axes[i].text(
                0.5,
                0.5,
                "NO SEPSIS",
                size=48,
                bbox={"facecolor": "g", "alpha": 0.4, "edgecolor": "none", "pad": 1},
                horizontalalignment="center",
                verticalalignment="center",
                transform=axes[i].transAxes,
            )
    #### Write out figure
    # make the path a path obj
    output_path = Path(output_path)
    # Block out CSN
    # csn ='CENSORED'
    # title of figure
    title = "SOFA_Score_" + csn + ".png"
    fig_path = output_path / title

    plt.savefig(fig_path, bbox_inches="tight", pad_inches=0)
