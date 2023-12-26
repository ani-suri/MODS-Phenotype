import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import date

font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 16}

matplotlib.rc('font', **font)


def index_1D(input_df, service_id, col_name):
    df = input_df.copy()
    df.sort_values([service_id, 'hospital_encounter_current'], ascending=True, inplace=True)
    df[col_name] = (
        df.groupby(service_id)['hospital_encounter_current']
        .rank(method='first')
        .astype(int)
    )
    df.set_index([service_id, col_name], inplace=True)
    df.sort_index(ascending=True, inplace=True)
    return df

def sofa_status_check(df, score_cols, failure_cols):
    df = sofa_failure_check(df, score_cols, failure_cols)
    df['SOFA_no_failure'] = ~df[failure_cols].any(axis=1)
    return df

def sofa_failure_check(df, score_cols, failure_cols):
    for score, failure in zip(score_cols, failure_cols):
        df[failure] = df[score] >= 2
    return df

def normalize_sofa_scores(df, cols):
    # Forward Fill all SOFA columns
    score_cols = [f"{col}_score" for col in cols]
    df[score_cols] = df.loc[:,cols].ffill()

    # Zero score for all remaining SOFA columns after forward fill
    df.loc[:,score_cols] = df.loc[:,score_cols].fillna(0)

    # Cast as INT
    df[score_cols] = df[score_cols].astype(int)
    return df

def plot_sofa_status(df, title):
    today = date.today()
    # TODO plot offset
    plt_data = df.sum()
    labels, counts = plt_data.index, plt_data.values
    sofa_colors = ('b','r','c','m','g','y','k')
    ticks = range(len(counts))
    plt.bar(ticks, counts, align='center', color=sofa_colors)

    ## Configurations
    
    # Presentation
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = [14, 6]
    plt.title(title)
    plt.figtext(0.98, 0.95, today, horizontalalignment='right', size=10, weight='light', in_layout=True)

    
    # Axes
    plt.xticks(ticks, labels,rotation=90)
    plt.xlabel('Type of Organ Failure')
    plt.ylabel('Number of Medical Encounters')
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    width = 0
    spare_width = (1 - width*2)/2
    ax.set_xlim(-spare_width,len(labels)-spare_width)

    # Table
    plt_table = plt.table(cellText=[[str(n) for n in counts]],
                          colLabels=labels,
                          loc='bottom')
    plt_table.scale(1, 1.5)
    plt_table.auto_set_font_size(True)
    plt_table.set_fontsize(11)

    plt.show()
    
def plot_sofa_status_heatmap(data, reference):
    plt_data = data.T.astype(float)
    mask = plt_data.isna()
    fig, ax = plt.subplots(figsize=(16,8))
    # sns.set()
    with sns.axes_style("white"):
        ax = sns.heatmap(plt_data, ax=ax, cmap='vlag', cbar=False, vmin=0, vmax=1, square=True, annot=False, linewidths=.5, mask=mask)
    ax.set_xlim(0, 24)
    # g.set_facecolor('grey')
    # legend = ax.legend()
    # legend.remove()
    ax.set_facecolor('xkcd:light grey')
    plt.xticks(rotation=15)
    plt.tight_layout()
    plt.show()