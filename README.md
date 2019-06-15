## Crevasse detection with OIB/ATM laser altimeter data

Airborne Topographic Mapper (ATM) laser altimetry data collected by NASA's Operation IceBridge. 

The procedure follows the analysis used to detect the morphology of sea ice features (e.g. pressure ridges) across the Arctic using OIB data (Petty et al., 2016) but applied to crevasses over specific glaciers in Greenland. Some important differences that need to be noted:

1. For sea ice we assume a relative flat reference surface (sea level) over the 500 m - 1000 m segments we analyze, whereas for these glaciers the underlying reference surface can vary considerably. We thus experiment with fitting planes to the elevations (e.g. a quadratic plane) and assess heigh anomalies relative to this surface.

2. Crevasse depths can be a lot bigger than pressure ridge heights!

3. Uncertainity regarding the ability of the off-nadir conical scanning laser to accureately pentrate down to the crevasse base.


### Setup

Note that individual script descriptions should be included at the top of each script.

The plan is for this to all be in Python 3.7. The original code was in 2.7 so I still need to check this porting hasn't messed anything up.

The 'calc' scripts are all used to process IceBridge ATM data into 'topography' datasets used by the 'plot' plotting scripts. 

Provide some conda info

### References

Petty, A. A., M. C. Tsamados, N. T. Kurtz, S. L. Farrell, T. Newman, J. P. Harbeck, D. L. Feltham, and J. A. Richter-Menge (2016), Characterizing Arctic sea ice topography using high-resolution IceBridge data</a>, The Cryosphere, 10(3), 1161–1179, doi:10.5194/tc-10-1161-2016.



