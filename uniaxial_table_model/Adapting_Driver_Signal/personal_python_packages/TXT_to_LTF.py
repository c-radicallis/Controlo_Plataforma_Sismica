# Python standard library
import sys

from pathlib import Path
# pxc modules
if (ppfolder := str(Path('..') / 'personal_python_packages')) not in sys.path:
    sys.path.append(ppfolder)

import numpy as np
from pathlib import Path
from scipy.constants import g

from LNECSPA import LTFdb

import pandas as pd

def txt_to_ltf(file_path , out_dir):
    df = pd.read_csv(file_path, delim_whitespace=True)
    time = df['time']
    PosT = df['PosT']
    PosL = df['PosL']
    PosV = df['PosV']
    accT = df['accT']
    accL = df['accL']
    accV = df['accV']
    number_columns = len(df.columns)-1
    ltfA = LTFdb()

    ltfA.update(t0=np.repeat(ltfA.t0, number_columns),
                dt=np.repeat(ltfA.dt*time[1], number_columns),
                dataformat=np.repeat(ltfA.dataformat, number_columns),
                IDstring=ltfA.IDstring*number_columns,
                scalefactor=np.repeat(ltfA.scalefactor, number_columns),
                offset=np.repeat(ltfA.offset, number_columns),
                data = [PosT*1e3, PosL*1e3, PosV*1e3, accT/g, accL/g, accV/g], 
                names=['DispT', 'DispL', 'DispV', 'AccT', 'AccL', 'AccV'],
                units=['mm','mm','mm', 'g','g','g'],
                types=['Displacement','Displacement','Displacement', 'Acceleration','Acceleration','Acceleration'],
                info=ltfA.info*6)
    # Make sure output directory exists
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Determine output filename
    file_path = Path(file_path)
    filename = file_path.name + ".acq"
    output_path = out_dir / filename
    # Write the .ltf file into the specified folder
    ltfA.write(str(output_path))

    return str(output_path)  # return path for confirmation