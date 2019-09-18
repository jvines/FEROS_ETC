"""ETC.

Exposure time calculator for optic observations in the V band at FEROS

Author: Jose Vines
"""

from __future__ import division, print_function


import astropy.units as u
import scipy as sp

from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource, Range1d, HoverTool

from utils import *

if __name__ == '__main__':
    description = 'Exposure Time Calculator'
    m0 = 12.82
    k = .1  # Extinction
    Z = 3.60994e-09  # V band zero point in erg/cm2/s/A
    ms = 20.0   # bright night sky magnitude.
    # http://www.ls.eso.org/sci/facilities/lasilla/telescopes/d1p5/misc/SkyBrightness.html
    b = 0.0051 * u.nm / u.bin  # spectral width
    b = b.to(u.um / u.bin)
    lam = 555.15 * u.nm  # wavelength ~ 40 order

    f = .14  # Efficiency
    A = 3.8 * u.m**2  # Telescope area
    npx = 3  # number of pixels

    # Airmass
    X = 2.5
    # Calculate energy per photon
    P = energy(lam.to(u.m))

    # Correct magnitude for airmass and extinction
    ms_corrected = correct_mag(ms, X, k)
    # Calculate flux
    Fm = flux(ms_corrected, Z)
    Fm.to(u.joule / u.s / u.m**2 / u.m)
    fm = f
    # Flux per photon
    Ffm = Fm / P
    # Number of photons from background
    B = Ffm * b * fm * A * npx * u.electron

    dark = 0.000305556
    phi = 4.8

    t = sp.linspace(0, 3600, 10000) * u.s
    sn = s_n(S.value, B.value, phi, npx, dark, t.value)

    mags = sp.linspace(6, 20, 10000)
    sn10 = []
    sn15 = []
    sn20 = []
    sn30 = []
    sn45 = []
    sn60 = []
    for m in mags:
        # Correct magnitude for airmass and extinction
        m_corrected = correct_mag(m, X, k)
        # Calculate flux
        F = flux(m_corrected, Z)
        F = F.to(u.joule / u.s / u.m**2 / u.um)
        # Flux per photon
        Ff = F / P
        # Number of photons from object
        S = Ff * b * f * A * u.electron
        times = sp.array([10, 15, 20, 30, 45, 60]) * u.minute.to(u.s)
        sn10.append(s_n(S.value, B.value, phi, npx, dark, times[0]))
        sn15.append(s_n(S.value, B.value, phi, npx, dark, times[1]))
        sn20.append(s_n(S.value, B.value, phi, npx, dark, times[2]))
        sn30.append(s_n(S.value, B.value, phi, npx, dark, times[3]))
        sn45.append(s_n(S.value, B.value, phi, npx, dark, times[4]))
        sn60.append(s_n(S.value, B.value, phi, npx, dark, times[5]))

    # Create Boke plots for a single airmass value
    output_file('S_N_mag_feros_X_{:.2f}.html'.format(X))

    TOOLS = "hover,help"
    TOOLTIPS = [
        ('t', '@texp'),
        ('V', '$x'),
        ('S/N', '$y')
    ]

    p = figure(plot_width=800, plot_height=800,
               y_axis_type='log', x_axis_label='V',
               y_axis_label='S/N', title='Airmass = {:.2f}'.format(X))

    data = dict(
        x=mags,
        y1=sn10, y2=sn15,
        y3=sn20, y4=sn30,
        y5=sn45, y6=sn60,
        t1=sp.ones(10000) * 10 * 60,
        t2=sp.ones(10000) * 15 * 60,
        t3=sp.ones(10000) * 20 * 60,
        t4=sp.ones(10000) * 30 * 60,
        t5=sp.ones(10000) * 45 * 60,
        t6=sp.ones(10000) * 60 * 60
    )
    source = ColumnDataSource(data)

    p1 = p.line(x='x', y='y1', color='cornflowerblue',
                line_width=4, source=source)
    p2 = p.line(x='x', y='y2', color='chocolate',
                line_width=4, source=source)
    p3 = p.line(x='x', y='y3', color='lightseagreen',
                line_width=4, source=source)
    p4 = p.line(x='x', y='y4', color='mediumseagreen',
                line_width=4, source=source)
    p5 = p.line(x='x', y='y5', color='firebrick',
                line_width=4, source=source)
    p6 = p.line(x='x', y='y6', color='darkmagenta',
                line_width=4, source=source)
    p.add_tools(HoverTool(
        renderers=[p1],
        tooltips=[
            ('t', '@t1'),
            ('V', '$x'),
            ('S/N', '@y1')
        ], mode='mouse')
    )
    p.add_tools(HoverTool(
        renderers=[p2],
        tooltips=[
            ('t', '@t2'),
            ('V', '$x'),
            ('S/N', '@y2')
        ], mode='mouse')
    )
    p.add_tools(HoverTool(
        renderers=[p3],
        tooltips=[
            ('t', '@t3'),
            ('V', '$x'),
            ('S/N', '@y3')
        ], mode='mouse')
    )
    p.add_tools(HoverTool(
        renderers=[p4],
        tooltips=[
                ('t', '@t4'),
                ('V', '$x'),
                ('S/N', '@y4')
                ], mode='mouse')
                )
    p.add_tools(HoverTool(
        renderers=[p5],
        tooltips=[
                ('t', '@t5'),
                ('V', '$x'),
                ('S/N', '@y5')
                ], mode='mouse')
                )
    p.add_tools(HoverTool(
        renderers=[p6],
        tooltips=[
                ('t', '@t6'),
                ('V', '$x'),
                ('S/N', '@y6')
                ], mode='mouse')
                )
    p.x_range = Range1d(21, 5)
    show(p)
