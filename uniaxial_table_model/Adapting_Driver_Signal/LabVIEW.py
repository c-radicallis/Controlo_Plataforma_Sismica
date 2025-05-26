#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""LabVIEW.py - LabVIEW function library

Created on Sat Jun  9 17:50:13 2018

@author: Paulo Xavier Candeias
"""

# Developed for Python 3.11.5
__author__     = 'Paulo Xavier Candeias'
__copyright__  = 'Copyright 2013, Paulo Xavier Candeias'
__credits__    = ['Paulo Xavier Candeias']
__license__    = 'GPL-3.0'
__version__    = '2024-01-27'
__maintainer__ = 'Paulo Xavier Candeias'
__email__      = 'pxcandeias@gmail.com'
__status__     = 'Production' # one of "Prototype", "Development", or "Production"

# Python standard library
# Nothing here

# 3rd party modules
import numpy as np

# pxc modules
from utilities import nextpow2


class LV_waveform():
    """LabVIEW waveform.

    A LabVIEW waveform is composed of the following elements:
    t0 - initial time stamp
    dt - time sampling
    y - signal

    Additional attributes can be defined by the user.
    
    References
    ----------
    https://zone.ni.com/reference/en-XX/help/371361R-01/lvwave/waveform_functions_and_vis/
    """
    def __init__(self, t0=0., dt=1., y=np.zeros(1)):
        """Initialize class.
        
        Parameters
        ----------
        t0 : float, optional
            Initial time stamp,
        dt : float, optional
            Time sampling.
        y : ndarray, optional
            Signal.
        """
        self.t0 = t0
        self.dt = dt
        self.y = y

    def time_stamps(self):
        """Generate time stamps.
        
        Returns
        -------
        out : nadarray
            Time stamps.
        """
        return self.t0 + np.arange(self.y.size) * self.dt


def LV_adc(analog, resolution: int, full_scale: float, data_format: int = 1):
    """LabVIEW Analog to Digital conversion.
    
    Parameters
    ----------
    analog : array_like
        Analog waveform to be converted to digital waveform or digital data.

    resolution : int
        Number of bits represented in the digital waveform or digital data.

    full_scale : float
        Total peak-to-peak range, or the difference between the minimum and maximum, for the digital waveform or digital data.

    data_format : int, optional
        Binary representation of the digital waveform or digital data (0: unsigned binary, 1: offset binary, 2: 2's complement).
        Default is 1
    
    Returns
    -------
    out : tuple
        Tuple containing (digital, resolution, full_scale).
    
    References
    ----------
    https://zone.ni.com/reference/en-XX/help/371361R-01/lvwave/analog_to_digital_wf/
    """
    # validate parameters
    analog = np.asarray(analog)

    assert 0 < resolution <= 64 , f'`resolution` must be in the range ]0, 64] ({resolution}).'

    assert full_scale > 0, f'`full_scale` must be greater than zero ({full_scale}).'

    assert data_format in range(3), f'`data_format` must be one of [0, 1, 2] ({data_format}).'

    resolution = max(nextpow2(resolution), 8) # minimum resolution is 8

    if data_format == 0:
        dtype = np.dtype(f'>u{resolution // 8}')
    elif data_format == 1:
        dtype = np.dtype(f'>i{resolution // 8}')
    elif data_format == 2:
        dtype = np.dtype(f'<i{resolution // 8}')
    else:
        raise ValueError(f'Wrong `data_format` value ({data_format}).')

    digital = np.asarray(analog / full_scale * 2**resolution, dtype=dtype)

    return digital, resolution, full_scale


def LV_dac(digital, full_scale: float, data_format: int = 1):
    """LabVIEW Digital to Analog conversion.
    
    Parameters
    ----------
    digital : array_like
        Digital waveform.

    full_scale : float
        Total peak-to-peak range, or the difference between the minimum and maximum, for the digital waveform or digital data.

    data_format : int, optional
        Binary representation of the digital waveform or digital data (0: unsigned binary, 1: offset binary, 2: 2's complement).
        Default is 1
    
    Returns
    -------
    out : tuple
        Tuple containing (analog, resolution, full_scale).
    
    References
    ----------
    https://zone.ni.com/reference/en-XX/help/371361R-01/lvwave/digital_to_analog_wf/
    """
    # validate parameters
    digital = np.asarray(digital)

    assert full_scale > 0, f'`full_scale` must be greater than zero ({full_scale}).'

    assert data_format in range(3), f'`data_format` must be one of [0, 1, 2] ({data_format}).'

    resolution = 8 * digital.itemsize

    if data_format == 0:
        dtype = np.dtype(f'>u{resolution // 8}')
    elif data_format == 1:
        dtype = np.dtype(f'>i{resolution // 8}')
    elif data_format == 2:
        dtype = np.dtype(f'<i{resolution // 8}')
    
    if dtype != digital.dtype:
        raise ValueError(f'`data_format` not compatible with `digital.dtype` ({dtype}, {digital.dtype}).')

    analog = np.asarray(digital * full_scale / 2**resolution)

    return analog, resolution, full_scale


def LV_to_integer(values, dtype):
    """LabVIEW To____Integer function family.
    
    Correspondence between LabVIEW types and numpy.dtype:
        Byte - int8
        Word - int16
        Long - int32
        Quad - int64
    
    This function rounds all floating-point and fixed-point numeric values to
    the nearest EVEN integer.
    http://zone.ni.com/reference/en-XX/help/371361J-01/glang/to_word_integer/
    
    LabVIEW coerces out-of-range values to the minimum or maximum value of the
    integer.
    http://zone.ni.com/reference/en-XX/help/371361J-01/lvhowto/numeric_conversion/
    
    Parameters
    ----------
    values : ndarray
        Numpy array of values.

    dtype - numpy.dtype
        Numpy integer dtype.
    
    Returns
    -------
    out : ndarray
        Numpy array of coerced `dtype` values.
    
    See also
    --------
    `numpy.dtype.itemsize`
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361J-01/glang/conversion_functions/
    """
    values = np.around(values)
    Drange = 2**(8*dtype.itemsize-1)
    values[values<-Drange] = -Drange
    values[values>=Drange] = Drange-1
    
    return values.astype(dtype)


def LV_to_unsigned_integer(values, dtype):
    """LabVIEW ToUnsigned____Integer function family.
    
    Correspondence between LabVIEW types and numpy.dtype:
        UnsignedByte - uint8
        UnsignedWord - uint16
        UnsignedLong - uint32
        UnsignedQuad - uint64
    
    This function rounds all floating-point and fixed-point numeric values to
    the nearest EVEN unsigned integer.
    http://zone.ni.com/reference/en-XX/help/371361J-01/glang/to_word_integer/
    
    LabVIEW coerces out-of-range values to the minimum or maximum value of the
    unsigned integer.
    http://zone.ni.com/reference/en-XX/help/371361J-01/lvhowto/numeric_conversion/
    
    Parameters
    ----------
    values : ndarray
        Numpy array of values.

    dtype : ndarray
        Numpy unsigned integer dtype.
    
    Returns
    -------
    out : ndarray
        Numpy array of coerced dtype values.
    
    See also
    --------
    `numpy.dtype.itemsize`
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361J-01/glang/conversion_functions/
    """
    values = np.around(values)
    Drange = 2**8 * dtype.itemsize
    values[values<0] = 0
    values[values>=Drange] = Drange-1
    
    return values.astype(dtype)


def LV_decimate(x, decimating_factor=1, averaging=False, start_index=0):
    """Signal decimation, LabVIEW version.
    
    Parameters
    ----------
    x : ndarray
        Signal.

    decimating_factor : int, optional
        Decimating factor.
        Default is 1

    averaging : bool, optional
        Controls whether the signal is averaged or not.
        Default is False

    start_index : int, optional
        Default is 0
    
    Returns
    -------
    y : ndarray
        Decimated signal.
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361H-01/lvanls/decimatecont/
    """
    if averaging:
        y = 0
        
        for v in range(decimating_factor):
            y += x[start_index+v:v-decimating_factor:decimating_factor]
        
        y /= float(decimating_factor)
    else:
        y = x[start_index::decimating_factor]
    
    return y


def LV_resample():
    """Signal resample, LabVIEW version.
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361J-01/lvwave/resample_waveforms_cont/
    """
    pass


def LV_cosine_tapered_window(N, r=0.2):
    """Cosine tapered window, LabVIEW version.
    
    Parameters
    ----------
    N : int
        Window length.

    r : float, optional
        Ratio between the cosine length and the window length.
        Default is 0.2
    
    Returns
    -------
    out : ndarray
        Cosine tapered window.
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361P-01/TOC166.htm
    http://zone.ni.com/reference/en-XX/help/371361P-01/lvanls/cosine_tapered_window/
    """
    m = np.floor(N*r/2.)
    i = np.arange(m)
    out = np.ones(N)
    out[i] = 0.5 * (1.-np.cos(2.*np.pi*i/(2.*m)))
    out[N-i-1] = out[i]
    
    return out


def LV_force_window(N, duty_cycle=0.5):
    """Force window, LabVIEW version.
    
    The force window is used to analyse transients.
    
    Parameters
    ----------
    N : int
        Window length.

    duty_cycle : float, optional
        Ratio between the force length and the window length.
        Default is 0.5
    
    Returns
    -------
    out : ndarray
        Force window.
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361P-01/TOC166.htm
    http://zone.ni.com/reference/en-XX/help/371361P-01/lvanls/force_window/
    """
    out = np.zeros(N, dtype=float)
    n = int(N * duty_cycle)
    out[:n] = 1.
    
    return out


def LV_exponential_window(N, final_value=0.1):
    """Exponential window, LabVIEW version.
    
    The exponential window is used to analyse transients.
    
    Parameters
    ----------
    N : int
        Window length.

    final_value : float, optional
        Coefficient of the last point of the window.
        Default is 0.1
    
    Returns
    -------
    out : ndarray
        Exponential window.
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361P-01/TOC166.htm
    http://zone.ni.com/reference/en-XX/help/371361P-01/lvanls/exponential_window/
    """
    out = np.geomspace(1., final_value, num=N, endpoint=True)
    
    return out


def LV_force_exponential_window(N, force_window=0.5, exp_window=0.1):
    """Force - Exponential window, LabVIEW version.
    
    The combined force and exponential windows are used to analyse transients.
    
    Parameters
    ----------
    N : int
        Window length.

    force_window : float
        Length of the force window.
        Default is 0.5
        force_window specifies the duration of the force window as a fraction
        of the total duration of the signal.
        Setting `force_window` at 1 has the effect of not applying any window
        on the stimulus signal.

    exp_window : float
        Decay rate of the exponential window.
        Default is 0.1
        `exp_window` specifies the remaining level of the applied exponential
        window at the end of the signal as a fraction.
    
    Returns
    -------
    out : ndarray
        Force - Exponential window.
    
    References
    ----------
    http://zone.ni.com/reference/en-XX/help/371361P-01/TOC166.htm
    http://zone.ni.com/reference/en-XX/help/372416K-01/sndvibtk/freqncy_response_mag_phase/
    """
    Fw = LV_force_window(N, duty_cycle=force_window)
    Ew = LV_exponential_window(N, final_value=exp_window)
    out = Fw * Ew
    
    return out


def LV_low_side_lobe_window(N):
    """Low side lobe window.

    The Low Sidelobe window reduces the level of the side lobe at the cost of broadening the main lobe. 
    
    Parameters
    ----------
    N : int
        Window length.
    
    Returns
    -------
    out : ndarray
        Low side lobe window.
    
    Reference
    ---------
    https://www.ni.com/docs/en-US/bundle/labwindows-cvi/page/advancedanalysisconcepts/lvac_low_sidelobe.html
    """
    n = np.arange(N)
    w = 2. * np.pi * n / N
    a = np.array((0.471492057, 0., 0.17553428, 0.028497078, 0.001261367))
    k = np.arange(a.size)

    out = np.sum((-1) ** k * a * np.cos(k * w[:, np.newaxis]), axis=1)

    return out


def main():
    """Test suite."""
    print(f'{nextpow2(8)=}')
    print(f'{max(nextpow2(5), 8)=}')
    print(f'{LV_low_side_lobe_window(8)=}')
    pass


if __name__ == '__main__':    # run tests if called from command-line
    main()
