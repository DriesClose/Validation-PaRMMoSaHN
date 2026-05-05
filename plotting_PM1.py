#! /bin/python

import pandas as pd
import matplotlib.pyplot as plt
import re
from scipy.signal import savgol_filter
import sys

## File inputs

# data_file = sys.argv[1]
# ref_file = sys.argv[2]
# figure_file_prefix = sys.argv[3]
# pos_list_file = sys.argv[4]

data_file = r"D:/growth_experiments_results/results_clostridia/data_excel/260421_biologPM1_Csymbiosum_3.xlsx" # Name of the MultiSkan Excel file (whole pathway)
#ref_file = "250617_essential_N-source_quest.xlsx"
figure_file_prefix = "260421_biologPM1_Csymbiosum_3" # prefix of the plot filename
pos_list_file = figure_file_prefix # filename of the list of positive growth compounds

## Well definitions

wells = {
    'A01': 'Control',
    'A02': 'L-arabinose',
    'A03': 'N-Acetyl-D-Glucosamine',
    'A04': 'D-saccharic acid',
    'A05': 'Succinic acid',
    'A06': 'D-galactose',
    'A07': 'L-aspartic acid',
    'A08': 'L-proline',
    'A09': 'D-alanine',
    'A10': 'D-trehalose',
    'A11': 'D-mannose',
    'A12': 'Dulcitol',
    'B01': 'D-serine',
    'B02': 'D-sorbitol',
    'B03': 'Glycerol',
    'B04': 'L-fucose',
    'B05': 'D-glucuronic acid',
    'B06': 'D-gluconic acid',
    'B07': 'D,L-alpha-glycerol-phosphate',
    'B08': 'D-xylose',
    'B09': 'L-lactic acid',
    'B10': 'Formic acid',
    'B11': 'D-mannitol',
    'B12': 'L-glutamic acid',
    'C01': 'D-glucose-6-phosphate',
    'C02': 'D-galactonic acid-gamma-lactone',
    'C03': 'D,L-malic acid',
    'C04': 'D-ribose',
    'C05': 'Tween 20',
    'C06': 'L-rhamnose',
    'C07': 'D-fructose',
    'C08': 'Acetic acid',
    'C09': 'D-glucose',
    'C10': 'Maltose',
    'C11': 'D-melibiose',
    'C12': 'Thymidine',
    'D01': 'L-asparagine',
    'D02': 'D-aspartic acid',
    'D03': 'D-glucosaminic acid',
    'D04': '1,2-propanediol',
    'D05': 'Tween 40',
    'D06': 'alpha-keto-glutaric acid',
    'D07': 'alpha-keto-butyric acid',
    'D08': 'alpha-methyl-D-galactoside',
    'D09': 'D-lactose',
    'D10': 'Lactulose',
    'D11': 'Sucrose',
    'D12': 'Uridine',
    'E01': 'L-glutamine',
    'E02': 'm-tartaric acid',
    'E03': 'D-glucose-1-phosphate',
    'E04': 'D-fructose-6-phopshate',
    'E05': 'Tween 80',
    'E06': 'alpha-hydroxy-glutaric acid-gamma-lactone',
    'E07': 'alpha-hydroxy-butyric acid',
    'E08': 'beta-methyl-D-glucoside',
    'E09': 'Adonitol',
    'E10': 'Maltotriose',
    'E11': '2-deoxy-adenosine',
    'E12': 'Adenosine',
    'F01': 'Glycyl-L-aspartic acid',
    'F02': 'citric acid',
    'F03': 'myo-inositol',
    'F04': 'D-threonine',
    'F05': 'Fumaric acid',
    'F06': 'Bromo-succinic acid',
    'F07': 'Propionic acid',
    'F08': 'Mucic acid',
    'F09': 'Glycolic acid',
    'F10': 'Glyoxylic acid',
    'F11': 'D-cellobiose',
    'F12': 'Inosine',
    'G01': 'Glycyl-L-glutamic acid',
    'G02': 'Tricarballylic acid',
    'G03': 'L-serine',
    'G04': 'L-threonine',
    'G05': 'L-alanine',
    'G06': 'L-alanyl-glycine',
    'G07': 'Acetoacetic acid',
    'G08': 'N-acetyl-beta-D-mannosamine',
    'G09': 'Mono-methyl-succinate',
    'G10': 'Methyl-pyruvate',
    'G11': 'D-malic acid',
    'G12': 'L-malic acid',
    'H01': 'Glycyl-L-proline',
    'H02': 'p-hydroxy-phenyl-acetic acid',
    'H03': 'm-hydroxy-phenyl-acetic acid',
    'H04': 'Tyramine',
    'H05': 'D-psicose',
    'H06': 'L-Lyxose',
    'H07': 'Glucuronamide',
    'H08': 'Pyruvic acid',
    'H09': 'L-galactonic acid-gamma-lactone',
    'H10': 'D-galacturonic acid',
    'H11': 'Phenylethylamine',
    'H12': '2-amino-ethanol'
    }

### Parameters
## Smoothening
smooth_window = 24
smooth_order = 4

## Metrics
rolling_window = 12
cumsum_cutoff = 0.75
peak_cutoff = 0.01

## Plotting
nrows = 8
ncols = 12
def plot_colour(pos):
    if pos:
        return 'g'
    else:
        return 'r'

## Reading data
def read_data(file):
    def convert_to_hours(seconds):
        return float(seconds)/3600
    
    df = pd.read_excel(file, header = 0, skiprows = 9, engine = "openpyxl", usecols = lambda x: x != 'Reading',
                       converters = {0: convert_to_hours})
    df = df.iloc[:,:-1]
    df = df.rename(columns = {'avg. time [s]': 'time'}).set_index('time')
    return df

df = read_data(data_file)
# df_blanks = read_data(ref_file)

## Processing
# Blanks moeten echt medium-only zijn, geen negative control
# De waardes die nu gebruikt worden, komen uit een eerdere test met P2, dus die zijn niet geldig voor jouw ZMB1.
def subtract_blanks(df, df_blanks = pd.DataFrame()): #HIER ZIT EEN FOUTJE
#    blank_cols = [i for i in df_blanks.columns if i.startswith('B1') or i.startswith('B2')]
#    df_blanks = df_blanks.loc[:, blank_cols]
#    mean_blank_OD = df_blanks.mean(axis = None)
    mean_blank_OD = 0.01 # !Fill your mean medium-only OD value here! (0.01 for ZMB1)
    df_corr = df - mean_blank_OD
    return df_corr

def smoothen(df, smooth_window, smooth_order):
    return df.apply(lambda x: savgol_filter(x, smooth_window, smooth_order, mode="nearest"))

def calculate_fold_change(df, control_label):
    return df.div(df[control_label], axis = 0)

def identify_growth(df_fc):
    df_max_fc = df_fc.max(axis = 0)
    df_pos = df_max_fc >= 2
    df_pos.name = "fold_change"
    return df_pos

def calculate_rolling_std(df, df_ref, col, rolling_window):
    series = df[col]
    series_ref = df_ref[col]
    rolling_difference_std = (series - series_ref).rolling(rolling_window, center = True).std()
    return rolling_difference_std

def detect_droplets(rolling_difference_std, cumsum_cutoff, peak_cutoff):
    rolling_difference_std_cumsum = rolling_difference_std.sum()
    rolling_std_peak = rolling_difference_std.max()
    cumsum_warning = rolling_difference_std_cumsum >= cumsum_cutoff
    peak_warning = rolling_std_peak >= peak_cutoff
    return (cumsum_warning, peak_warning)

df_corr = subtract_blanks(df)
df_corr_smooth = smoothen(df_corr, smooth_window, smooth_order)
df_fc = calculate_fold_change(df_corr_smooth, 'Un0001 (A01)')
df_roll_std = {well: calculate_rolling_std(df_corr, df_corr_smooth, well, rolling_window) for well in df_corr.columns}
df_roll_std = pd.DataFrame.from_dict(df_roll_std)
droplets = {well: detect_droplets(df_roll_std[well], cumsum_cutoff, peak_cutoff) for well in df.columns}
droplets = pd.DataFrame.from_dict(droplets, orient = "index", columns=['cumsum', 'peak'])
df_pos = identify_growth(df_fc)
checks = droplets.merge(df_pos, left_index = True, right_index = True)
checks['compound'] = [wells[re.search(r'\([A-Z][0-9]{2}\)', well).group(0)[1:-1]] for well in checks.index]

## Plotting
time = list(df.index)
fig_od, ax_od = plt.subplots(nrows = nrows, ncols = ncols, sharex = True, sharey = True, figsize = (45,30))
fig_cs, ax_cs = plt.subplots(nrows = nrows, ncols = ncols, sharex = True, sharey = True, figsize = (45,30))
fig_peak, ax_peak = plt.subplots(nrows = nrows, ncols = ncols, sharex = True, sharey = True, figsize = (45,30))
pos_compounds = []

for i, well in enumerate(checks.index):
    compound = checks.loc[well]['compound']
    series = df_corr[well]
    series_smooth = df_corr_smooth[well]
    series_roll_std = df_roll_std[well]
    
    pos = checks.loc[well]['fold_change']
    warning = checks.loc[well][['cumsum', 'peak']].any()
    colour = plot_colour(pos)
    if pos:
        pos_compounds.append(compound)
    
    ax_od_i = ax_od[i//ncols, i%ncols]
    ax_od_i.plot(time, series, colour + warning*"--")
    ax_od_i.plot(time, series_smooth, 'b:')
    ax_od_i.set_title(compound, color = colour)
    
    ax_cs_i = ax_cs[i//ncols, i%ncols]
    pos = checks.loc[well]['cumsum']
    colour = plot_colour(~pos)
    ax_cs_i.plot(time, series_roll_std, colour)
    ax_cs_i.set_title(compound, color = colour)
    
    ax_peak_i = ax_peak[i//ncols, i%ncols]
    pos = checks.loc[well]['peak']
    colour = plot_colour(~pos)
    ax_peak_i.plot(time, series_roll_std, colour)
    ax_peak_i.set_title(compound, color = colour)
    
ax_od[0][0].set_xlabel('Time [h]')
ax_od[0][0].set_ylabel('OD [-]')
ax_cs[0][0].set_xlabel('Time [h]')
ax_cs[0][0].set_ylabel('Roll. diff. std. [-]')
ax_peak[0][0].set_xlabel('Time [h]')
ax_peak[0][0].set_ylabel('Roll. diff. std. [-]')

fig_od.savefig(figure_file_prefix + '_plot.svg')
#fig_cs.savefig(figure_file_prefix + '_roll_std_cumsum.svg')
#fig_peak.savefig(figure_file_prefix + '_roll_std_peak.svg')

pd.Series(pos_compounds).to_csv(pos_list_file, index = False, header = False)


