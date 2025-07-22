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

import matplotlib.pyplot as plt


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

    else: #this is the case for when we want to write .tgt from the text files with fault normal and parallel acceleration from Fernando. 
            #This also corrects displacement drift by computing the slope m = (xF - xi) / (tF - ti), and subtrating the drift at every time instant disp_drift(t) = m*t
        df = pd.read_csv(
        file_path,
        delim_whitespace=True,
        header=None,
        names=['time', 'accT', 'accL'])

        # Create PosV column of zeros, same length as the data
        time = df['time']
        accT = df['accT']#0.4*
        accL = df['accL']#0.4*

        df['PosV'] = 0
        df['accV'] = 0
        PosV = df['PosV']
        accV = df['accV']

        velT = cumulative_trapezoid(accT, time, initial=0.0)
        PosT = cumulative_trapezoid(velT, time, initial=0.0)
        m_T = (PosT[-1] - PosT[0])/time.values[-1]
        PosT = PosT - m_T*time

        velL = cumulative_trapezoid(accL, time, initial=0.0)
        PosL = cumulative_trapezoid(velL, time, initial=0.0)
        m_L = (PosL[-1] - PosL[0])/time.values[-1]
        PosL = PosL - m_L*time

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
    

        # # # assume `time` is a 1‑D numpy array (or pandas Index) of monotonically increasing times,
        # # # and `PosT` is your detrended position array of the same length.
        # # # 1) compute velocity as dPos/dt
        # # d_PosT = np.gradient(PosT, time)
        # # # 2) compute acceleration as dVel/dt
        # # dd_PosT = np.gradient(d_PosT, time)

        # def second_derivative_time(y, dt):
        #     """
        #     Compute the second time‐derivative of a 1D signal.

        #     Parameters
        #     ----------
        #     y : array_like, shape (N,)
        #         Input signal (list or 1D NumPy array) with at least 3 samples.
        #     dt : float
        #         Constant time step between samples.

        #     Returns
        #     -------
        #     d2y : ndarray, shape (N,)
        #         Second‐derivative approximation at each sample.
        #     """
        #     y = np.asarray(y, dtype=float)
        #     if y.ndim != 1:
        #         raise ValueError("y must be a 1D array or list")
        #     N = y.size
        #     if N < 3:
        #         raise ValueError("Input signal must have at least 3 samples.")

        #     # Preallocate output
        #     d2y = np.empty_like(y)

        #     # Forward‐difference (second order) at index 0
        #     d2y[0] = (y[2] - 2 * y[1] + y[0]) / dt**2

        #     # Central‐difference (second order) for interior points
        #     # d2y[i] = (y[i+1] - 2*y[i] + y[i-1]) / dt^2
        #     d2y[1:-1] = (y[2:] - 2 * y[1:-1] + y[:-2]) / dt**2

        #     # Backward‐difference (second order) at last index
        #     d2y[-1] = (y[-1] - 2 * y[-2] + y[-3]) / dt**2

        #     return d2y

        # dd_PosT = second_derivative_time(PosT, 0.005)

        # # 3) plot position, velocity, and acceleration
        # fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
        # axs[0].plot(time, PosT*1e3, label='∫∫original accel-m*t * 1e3', marker='.')
        # axs[0].set_ylabel('Position (mm)')
        # axs[0].grid(True)

        # # #axs[1].plot(time, d_PosT, label='2nd order central diff of ∫∫accell-m*t', marker='.')
        # # axs[1].plot(time, velT, label='∫original accel', marker='.')
        # # axs[1].set_ylabel('Velocity')
        # # axs[1].legend(loc='best')
        # # axs[1].grid(True)

        # axs[1].plot(time, dd_PosT, label='2nd order central diff of ∫∫original accel-m*t * 1e3', marker='.')
        # axs[1].plot(time, accT, label='Original Acceleration', marker='.')
        # axs[1].set_ylabel('Acceleration (m/s^2)')
        # axs[1].set_xlabel('Time')
        # axs[1].grid(True)

        # tgtA = LTFdb()
        # tgtA.read(str(output_path))
        # arr = np.array(tgtA._data, dtype=float)
        


        # axs[0].plot(time, arr[0], label='x_tgt_T (discretized signal in the .tgt file)', marker='.')
        # x_tgt = arr[0]
        # ddx_tgt = arr[3]*g
        # #axs[1].plot(time, d_PosT, label='d_∫∫accell-m*t', marker='.')
        # axs[1].plot(time, arr[3]*g, label='ddx_tgt_T (discretized original signal in the .tgt file)', marker='.', linestyle='--')
        # d2_x_tgt = second_derivative_time(arr[0]*1e-3, 0.005)
        # axs[1].plot(time, d2_x_tgt, label='2nd order central diff from x_tgt_T', marker='.')

        # print('Média accel original= ',np.mean(accT))
        # print('Média ddx_tgt= ',np.mean(ddx_tgt))
        # print('Média d2_x_tgt= ',np.mean(d2_x_tgt))

        # axs[0].legend(loc='best')
        # axs[1].legend(loc='best')   

        # # axs[0].set_xlim(0, 0.05)       # x-axis limits
        # # axs[0].set_ylim(-0.02, 0.01)   # y-axis limits
        # # # axs[1].set_xlim(0, 0.05)       # x-axis limits
        # # # axs[1].set_ylim(-0.001, 0.001)   # y-axis limits
        # # axs[1].set_xlim(0, 0.05)       # x-axis limits
        # # axs[1].set_ylim(-0.15, 0.3)   # y-axis limits

        # plt.title('The errors resulting from the discretization of the Position when writing it to the .tgt file (see plot above), further propagate when computing the 2nd derivative of that signal (see red line in plot below)')
        # plt.tight_layout()
        # plt.show()


    return str(output_path)




def txt_to_drv(file_path, out_dir):
    # Read file without headers, assign column names
    df = pd.read_csv(
        file_path,
        delim_whitespace=True,
        header=0,
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
        dt=np.repeat(ltfA.dt * time[1], number_columns),
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

    # assume file_path is a string or Path to your input file
    file_path = Path(file_path)
    raw_name = file_path.name.lower()

    # figure out the base name (without extension)
    if raw_name.endswith('.drv.txt'):
        base = file_path.stem[:-4]    # .stem removes .txt, so strip ".drv" off that
    elif raw_name.endswith('.txt'):
        base = file_path.stem         # just .stem to drop the .txt
    else:
        # no recognized extension—fall back to the full name without suffix
        base = file_path.stem

    # build the new name
    new_name = base + '.DRV'

    # assemble the output path
    output_path = out_dir / new_name

    # write it out
    ltfA.write(str(output_path))

    return str(output_path)


## Example / Test

# out_dir = Path(r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\pink_noise_test')
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\pink_noise\pink_noise_40Hz_T3mm_0.drv.txt'
# txt_to_drv( filepath, out_dir)

# Pink Noise
# out_dir = Path(r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\pink_noise')
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\pink_noise\pink_noise_40Hz_T3mm_0.drv.txt'
# txt_to_ltf(filepath, out_dir)

# # El Centro
# out_dir = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_ElCentro_scaled'
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\elcentro.txt'
# txt_to_ltf(filepath, out_dir)

# # Erzikan
# out_dir = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Erzikan'
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\erzikan.txt'
# txt_to_ltf(filepath, out_dir)


 # Jiji
out_dir = Path(r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Jiji_corrected')
filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\jiji.txt'
txt_to_ltf(filepath, out_dir)

# # Kobe
# out_dir = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Kobe'
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\kobe.txt'
# txt_to_ltf(filepath, out_dir)

# # Newhall
# out_dir = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Newhall_corrected'
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\newhall.txt'
# txt_to_ltf(filepath, out_dir)

# # Rinaldi
# out_dir = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Rinaldi'
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\rinaldi.txt'
# txt_to_ltf(filepath, out_dir)

# # sylmar
# out_dir = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Sylmar'
# filepath = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\AcaoSismica\Sismos\sylmar.txt'
# txt_to_ltf(filepath, out_dir)