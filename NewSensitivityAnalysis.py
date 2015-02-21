"""
Sensitivity analysis, using the new search method.
"""
import sys
import itertools

import Sensitivity
import StarData
import matplotlib.pyplot as plt
import numpy as np
import seaborn
import SpectralTypeRelations
import pandas as pd

import Search_slow

MS = SpectralTypeRelations.MainSequence()

def check_sensitivity():
    fileList = []
    for arg in sys.argv[1:]:
        if 1:
            fileList.append(arg)

    badregions = Search_slow.badregions
    interp_regions = Search_slow.interp_regions
    trimsize = Search_slow.trimsize
    prim_vsini = StarData.get_vsini(fileList)

    Sensitivity.Analyze(fileList, prim_vsini,
                        hdf5_file='/media/ExtraSpace/PhoenixGrid/HRS_Grid.hdf5',
                        extensions=True,
                        resolution=None,
                        trimsize=trimsize,
                        badregions=badregions, interp_regions=interp_regions,
                        metal_values=(0.0,),
                        vsini_values=(5,),
                        Tvalues=range(3000, 6000, 100),
                        debug=False,
                        addmode='simple',
                        output_mode='hdf5')


def analyze_sensitivity(hdf5_file='Sensitivity.hdf5'):
    """
    This uses the output of a previous run of check_sensitivity, and makes plots
    :return:
    """
    hdf5_int = Sensitivity.HDF5_Interface(hdf5_file)
    df = hdf5_int.to_df()

    # Group by a bunch of keys that probably don't change, but could
    groups = df.groupby(('star', 'date', '[Fe/H]', 'logg', 'vsini', 'addmode', 'primary SpT'))

    # Plot
    seaborn.set_style('white')
    seaborn.set_style('ticks')
    plot_styles = ['-', '-.', '--', ':']
    ps_cycler = itertools.cycle(plot_styles)
    sig_fig = plt.figure('Significance Summary')
    sig_ax = sig_fig.add_subplot(111)
    rate_fig = plt.figure('Detection Rate')
    rate_ax = rate_fig.add_subplot(111)
    df_list = []
    for key in sorted(groups.groups.keys(), key=lambda l: MS.SpT_To_Number(l[6][:2])):
        g = groups.get_group(key)
        prim_spt = g['primary SpT'].values[0]
        T_groups = g.groupby('temperature')
        T_keys = T_groups.groups.keys()
        Tvals = np.zeros(len(T_keys))
        rate = np.zeros(len(T_keys))
        significance = np.zeros(len(T_keys))
        for i, Tstring in enumerate(sorted(T_keys)):
            numdetected = sum(T_groups.get_group(Tstring)['significance'].notnull())
            med_sig = np.nanmedian(T_groups.get_group(Tstring)['significance'])
            Tvals[i] = float(Tstring)
            rate[i] = 100.0 * numdetected / float(len(T_groups.get_group(Tstring)))
            significance[i] = med_sig

        labelstr = '{} ({})'.format(key[0], prim_spt)
        ls = ps_cycler.next()
        sig_ax.plot(Tvals, significance, ls, lw=2, label=labelstr)
        rate_ax.plot(Tvals, rate, ls, lw=2, label=labelstr)

        # Save the summary in a dataframe
        data = pd.DataFrame(
            data={'star': [key[0]] * rate.size, 'rate': rate, 'significance': significance, 'temperature': Tvals})
        df_list.append(data)

    summary = pd.concat(df_list, ignore_index=True)
    summary.to_csv('Sensitivity_Summary.csv', index=False)

    # Put labels on the plots
    sig_ax.set_xlabel('Temperature (K)', fontsize=15)
    sig_ax.set_ylabel('Median Significance', fontsize=15)
    sig_leg = sig_ax.legend(loc='best', fancybox=True)
    sig_leg.get_frame().set_alpha(0.5)

    rate_ax.set_xlabel('Temperature (K)', fontsize=15)
    rate_ax.set_ylabel('Detection Rate (%)', fontsize=15)
    rate_leg = rate_ax.legend(loc='best', fancybox=True)
    rate_leg.get_frame().set_alpha(0.5)

    # Make a top set of axes that shows spectral type
    sig_top_ax = add_top_axis(sig_ax)
    sig_top_ax.set_xlabel('Spectral Type', fontsize=15)
    rate_top_ax = add_top_axis(rate_ax)
    rate_top_ax.set_xlabel('Spectral Type', fontsize=15)

    plt.show()


def add_top_axis(axis, spt_values=('M5', 'M0', 'K5', 'K0', 'G5', 'G0')):
    # get the full range of the axis.
    xlim = axis.get_xlim()

    # Find the temperatures at each spectral type
    MS = SpectralTypeRelations.MainSequence()
    temp_values = MS.Interpolate('Temperature', spt_values)

    # make the axis
    top = axis.twiny()
    top.set_xticks(temp_values)
    top.set_xlim(xlim)
    top.set_xticklabels(spt_values)
    return top


if __name__ == '__main__':
    if 'analyze' in sys.argv[1]:
        analyze_sensitivity()
    else:
        check_sensitivity()
