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
from scipy.integrate import cumulative_trapezoid

def txt_to_ltf(file_path, out_dir):
    """
    This function writes a .acq ot .tgt file from .txt file

    Reads a whitespace-delimited text file without headers, where:
      - column 0 is time
      - column 1 is PosT
      - column 2 is PosL
    Creates PosV as zeros, then writes a .drv via LTFdb.
    """
    df = pd.read_csv(file_path, delim_whitespace=True)
    number_columns = len(df.columns) - 1

    if number_columns == 6: # this is the case for which we'll want to write a .acq
        time = df['time']
        PosT = df['PosT']
        PosL = df['PosL']
        PosV = df['PosV']
        accT = df['accT']
        accL = df['accL']
        accV = df['accV']
        ltfA = LTFdb()

        ltfA.update(t0=np.repeat(ltfA.t0, number_columns),
                    dt=np.repeat(ltfA.dt * time[1], number_columns),
                    dataformat=np.repeat(ltfA.dataformat, number_columns),
                    IDstring=ltfA.IDstring * number_columns,
                    scalefactor=np.repeat(ltfA.scalefactor, number_columns),
                    offset=np.repeat(ltfA.offset, number_columns),
                    data=[PosT * 1e3, PosL * 1e3, PosV * 1e3, accT / g, accL / g, accV / g],
                    names=['DispT', 'DispL', 'DispV', 'AccT', 'AccL', 'AccV'],
                    units=['mm', 'mm', 'mm', 'g', 'g', 'g'],
                    types=['Displacement', 'Displacement', 'Displacement',
                        'Acceleration', 'Acceleration', 'Acceleration'],
                    info=ltfA.info * 6)

        # Make sure output directory exists
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Determine output filename
        file_path = Path(file_path)
        raw_name = file_path.name  # e.g. "LAquilaReducedScale_0.ACQ.txt" or "data.txt" etc.

        # If it ends in ".ACQ.txt" (case-insensitive), strip only the ".txt"
        if raw_name.lower().endswith('.acq.txt'):
            new_name = raw_name[:-4]  # remove the last 4 characters ".txt"
        else:
            # otherwise append ".acq" as before
            new_name = raw_name + ".acq"

        output_path = out_dir / new_name

        # Write the .ltf/.acq file into the specified folder
        ltfA.write(str(output_path))

    else:
        df = pd.read_csv(
        file_path,
        delim_whitespace=True,
        header=None,
        names=['time', 'accT', 'accL'])

        # Create PosV column of zeros, same length as the data
        time = df['time']
        accT = df['accT']
        accL = df['accL']

        df['PosV'] = 0
        df['accV'] = 0

        velT = cumulative_trapezoid(accT, time, initial=0.0)
        PosT = cumulative_trapezoid(velT, time, initial=0.0)

        velL = cumulative_trapezoid(accL, time, initial=0.0)
        PosL = cumulative_trapezoid(velL, time, initial=0.0)

        number_columns = 6

        # Initialize LTFdb instance
        ltfA = LTFdb()
        ltfA.update(t0=np.repeat(ltfA.t0, number_columns),
                    dt=np.repeat(ltfA.dt * time[1], number_columns),
                    dataformat=np.repeat(ltfA.dataformat, number_columns),
                    IDstring=ltfA.IDstring * number_columns,
                    scalefactor=np.repeat(ltfA.scalefactor, number_columns),
                    offset=np.repeat(ltfA.offset, number_columns),
                    data=[PosT * 1e3, PosL * 1e3, PosV * 1e3, accT / g, accL / g, accV / g],
                    names=['DispT', 'DispL', 'DispV', 'AccT', 'AccL', 'AccV'],
                    units=['mm', 'mm', 'mm', 'g', 'g', 'g'],
                    types=['Displacement', 'Displacement', 'Displacement',
                        'Acceleration', 'Acceleration', 'Acceleration'],
                    info=ltfA.info * 6)

        # Ensure output directory exists
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        file_path = Path(file_path)
        raw_name = file_path.name
        new_name = raw_name[:-4]+ ".tgt"  # remove the last 4 characters ".txt"

        output_path = out_dir / new_name

        # Write out the file
        ltfA.write(str(output_path))


    return str(output_path)




def txt_to_drv(file_path, out_dir):
    """
    Reads a whitespace-delimited text file without headers, where:
      - column 0 is time
      - column 1 is PosT
      - column 2 is PosL
    Creates PosV as zeros, then writes a .drv via LTFdb.
    """
    # Read file without headers, assign column names
    df = pd.read_csv(
        file_path,
        delim_whitespace=True,
        header=None,
        names=['time', 'PosT', 'PosL'],
    )

    # Create PosV column of zeros, same length as the data
    df['PosV'] = 0

    # Extract series
    time = df['time']
    PosT = df['PosT']
    PosL = df['PosL']
    PosV = df['PosV']

    # Number of channels/columns to write (time is not written directly; you keep 3 data columns)
    number_columns = 3

    # Initialize LTFdb instance
    ltfA = LTFdb()

    # Update with the three channels:
    # Note: the original used dt=np.repeat(ltfA.dt * time[1], ...)
    #       You may want to double-check that logic: time[1] is the second time value.
    #       If instead you need the sample interval, you could compute dt from time.diff(),
    #       but here we keep the original pattern.
    ltfA.update(
        t0=np.repeat(ltfA.t0, number_columns),
        dt=np.repeat(ltfA.dt * time.iloc[1], number_columns),
        dataformat=np.repeat(ltfA.dataformat, number_columns),
        IDstring=ltfA.IDstring * number_columns,
        scalefactor=np.repeat(ltfA.scalefactor, number_columns),
        offset=np.repeat(ltfA.offset, number_columns),
        data=[PosT * 1e3, PosL * 1e3, PosV * 1e3],
        names=['DispT', 'DispL', 'DispV'],
        units=['mm', 'mm', 'mm'],
        types=['Displacement', 'Displacement', 'Displacement'],
        info=ltfA.info * number_columns
    )

    # Ensure output directory exists
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Determine output filename: strip last 4 chars of original name (assumed ".txt"), append ".acq"
    file_path = Path(file_path)
    raw_name = file_path.name
    new_name = raw_name[:-4]+ ".DRV"  # remove the last 4 characters ".txt"

    output_path = out_dir / new_name

    # Write out the file
    ltfA.write(str(output_path))

    return str(output_path)
