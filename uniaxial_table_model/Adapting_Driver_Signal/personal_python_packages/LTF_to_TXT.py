from pathlib import Path
import numpy as np
from scipy.constants import g
from LNECSPA import LTFdb

def ltf_to_txt(file_path, out_dir):
    """
    Read LNEC .ltf file export to TXT, automatically classifying
    channels as displacement or acceleration based on types or units.
    Displacement units: 'mm' -> convert to meters.
    Acceleration units: 'g'  -> convert to m/s^2.
    After conversion, update units to SI and fill missing types.
    Directly modifies tgtA._data so that write_txt sees the converted arrays.
    """
    file_path = Path(file_path)
    tgtA = LTFdb()
    tgtA.read(file_path)

    disp_inds = [0 , 1 ,2]
    acc_inds  = [ 3, 4, 5]

    # # Determine channel indices by type or unit
    # for idx, (typ, unit) in enumerate(zip(tgtA.types, tgtA.units)):
    #     t = (typ or '').strip().lower()
    #     u = (unit or '').strip().lower()
    #     if t == 'acceleration' or (not t and u == 'g'):
    #         acc_inds.append(idx)
    #         if not tgtA.types[idx].strip():
    #             tgtA.types[idx] = 'acceleration'
    #     elif t == 'displacement' or (not t and u == 'mm'):
    #         disp_inds.append(idx)
    #         if not tgtA.types[idx].strip():
    #             tgtA.types[idx] = 'displacement'
    #     else:
    #         print(f"Channel {idx}: Unknown type '{typ}' and/or unit '{unit}'. Skipping conversion for this channel.")
    #         continue

    # Convert displacement channels
    for ch_idx in disp_inds:
        arr = np.array(tgtA._data[ch_idx], dtype=float)
        original_unit = (tgtA.units[ch_idx] or '').strip().lower()
        if original_unit == 'mm':
            arr = arr * 1e-3
            tgtA.units[ch_idx] = 'm'
        elif original_unit == 'm':
            pass
        else:
            print(f"Displacement channel {ch_idx}: Unrecognized unit '{original_unit}'. No conversion applied.")
        tgtA._data[ch_idx] = arr

    # Convert acceleration channels
    if file_path.suffix.upper() != ".DRV":
        for ch_idx in acc_inds:
            print("ch_idx=", ch_idx)
            arr = np.array(tgtA._data[ch_idx], dtype=float)
            original_unit = (tgtA.units[ch_idx] or '').strip().lower()
            if original_unit in ['g', 'G' ]:
                arr = arr * g
                tgtA.units[ch_idx] = 'm/s^2'
            elif original_unit in ['m/s^2', 'm/s²']:
                pass
            else:
                print(f"Acceleration channel {ch_idx}: Unrecognized unit '{original_unit}'. No conversion applied.")
            tgtA._data[ch_idx] = arr

    # Ensure output directory exists
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    filename = file_path.name + ".txt"
    output_path = out_dir / filename

    tgtA.write_txt(str(output_path))
    return str(output_path)  # return path for confirmation

# file_path = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Jiji\jiji_0.DRV'
# folder = r'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Jiji'
# ltf_to_txt(file_path , folder)