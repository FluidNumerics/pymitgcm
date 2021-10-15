# pymitgcm
Copyright 2020-2021 Fluid Numerics LLC


## Installation
The pymitgcm package can be installed via pip using this public git repository as the source.
```
python -m pip install git+https://github.com/FluidNumerics/pymitgcm.git
```



## Usage

```
#!/usr/bin/env python3

from pymitgcm import pymitgcm
import matplotlib.pyplot as plt


# The path should be a path to a mitgcm simulation directory that has
# subdirectories : build/, input/, and run/
#
model = pymitgcm(directory="/path/to/mitgcm/simulation",
                 iterate=0,
                 dateTimeStart=[2000,1,1,0,0,0])

# Query monitoring statistics from stdout
model.setMonitorStatistics()
plt.plot(model.monstats['time_secondsf'],model.monstats['advcfl_wvel_max'])
```
