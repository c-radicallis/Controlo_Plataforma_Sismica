#! /usr/bin/env python
# -*- coding: utf-8 -*-

""" utilities.py - Utilities

Created on Mon Jun 23 23:20:08 2014

@author: Paulo Xavier Candeias
"""

# Developed for Python 3.11.7
__author__     = 'Paulo Xavier Candeias'
__copyright__  = 'Copyright 2014, Paulo Xavier Candeias'
__credits__    = ['Paulo Xavier Candeias']
__license__    = 'GPL-3.0'
__version__    = '2024-03-22'
__maintainer__ = 'Paulo Xavier Candeias'
__email__      = 'pxcandeias@gmail.com'
__status__     = 'Production' # one of "Prototype", "Development", or "Production"

# Python standard library
import os
import math

from io import StringIO
from numbers import Number

# 3rd party modules
import numpy as np


def correct_filename(fullpath: str, prefix: str ='', suffix: str ='_CORR') -> str:
    """Correct filename.
    
    Parameters
    ----------
    fullpath: str
        File full path.

    prefix: str, optional
        Default is ''

    suffix: str, optional
        Default is '_CORR'
    
    Returns
    -------
    corrected_fullpath: str
        Corrected file full path
    """
    folder, filename = os.path.split(fullpath)
    name, extension = os.path.splitext(filename)
    corrected_fullpath = os.path.join(folder, f'{prefix}{name}{suffix}{extension}')
    
    return corrected_fullpath


class Trace(object):
    """Trace class to be used as decorator.
    
    References
    ----------
    https://en.wikibooks.org/wiki/Python_Programming/Decorators
    """
    def __init__(self, func):
        """..."""
        self.func = func
    
    def __call__(self, *args, **kwargs):
        """..."""
        print('Entering function {}'.format(self.func.__name__))
        
        for i, arg in enumerate(args):
            print('arg {0}: {1}'.format(i, arg))
        
        return self.func(*args, **kwargs)


def dump_args(func):
    """This decorator dumps out the arguments passed to a function before calling it.
    
    Parameters
    ----------
    func : callable
    
    Notes
    -----
    Defines a closure (https://www.learnpython.org/en/Closures).

    References
    ----------
    https://wiki.python.org/moin/PythonDecoratorLibrary#Easy_Dump_of_Function_Arguments
    """
    argnames = func.func_code.co_varnames[:func.func_code.co_argcount]
    fname = func.func_name

    def echo_func(*args,**kwargs):
        print(fname, ":", ', '.join(
            '%s=%r' % entry for entry in zip(argnames,args) + kwargs.items()))
        
        return func(*args, **kwargs)

    return echo_func


def hex_dump(f_in, start=0, length=None, linesize=16):
    """Prints length bytes from a binary file including start byte.
    
    Parameters
    ----------
    f_in : file object

    start : int, optional

    length : int, optional

    linesize : int, optional
    
    Returns
    -------
    out : str
        hexadecimal dump
    """
    if length is None:
        f_in.seek(0, 2)
        length = f_in.tell()
    
    first = start - (start % linesize)
    f_in.seek(first, 0) # begining of a linesize byte sequence
    length -= f_in.tell()
    lines = ['+' + ('-' * 8) + '+' + ('-' * (3 * linesize + 1)) +
             '+' + ('-' * (linesize + 2)) + '+']
    
    for line in range(length // linesize + 1):
        bin_string = f_in.read(linesize)
        text = list('| {0:06x} | '.format(first + line * linesize))

#        if first != start:
#            bin_string = bin_string[start-first:]
#            text.extend('   '* (start-first))
        
        text.extend(['{0:02x} '.format(char) for char in bin_string])
        
        if len(bin_string) % linesize > 0:
            text.extend('   ' * (linesize - (len(bin_string) % linesize)))
        
        text.extend('| ')
        text.extend([chr(char) if char>=0x20 else '.' for char in bin_string])
        
        if len(bin_string) % linesize > 0:
            text.extend(' ' * (linesize - (len(bin_string) % linesize)))
        
        text.extend(' |')
        lines.append(''.join(text))
    
    lines.append(lines[0])
    out = '\n'.join(lines)

    return out


def BCD(hexvalue, packed=False):
    """Binary coded decimal.
    
    Parameters
    ----------
    hexvalue : str
        Hexadecimal string.

    packed : bool, optional
        If `True`, packed BCD, if `False`, unpacked BCD.
    
    Returns
    -------
    out : str
        Formated time string.
    
    References
    ----------
    http://en.wikipedia.org/wiki/Binary-coded_decimal
    """
    if packed:
        print('TBD!')
        decvalue = None
    else:
        decvalue = [int(hex(v)[2:]) for v in hexvalue]
    
    return decvalue


def BCDtime(hexvalue, sep=','):
    """Binary coded decimal time.
    
    Parameters
    ----------
    hexvalue : str
        Hexadecimal string.

    sep : char, optional
        Fractional seconds separator.
    
    Returns
    -------
    out : str
        Formated time string
    """
    decvalue = BCD(hexvalue)
    formatstring = sep.join(('{:02d}:{:02d}:{:02d}','{:02d}'))
    
    return formatstring.format(*decvalue)


def BCDdate(hexvalue, sep='-'):
    """Binary coded decimal date.
    
    Parameters
    ----------
    hexvalue : str
        Hexadecimal string.

    sep : char, optional
        Day-month-year separator.
    
    Returns
    -------
    date : str
        Formated date string
    """
    decvalue = BCD(hexvalue)
    formatstring = sep.join(('{:02d}','{:02d}','{:02d}{:02d}'))
    
    return formatstring.format(*decvalue)


def CoG(mxyz):
    """Compute centre of gravity.
    
    Parameters
    ----------
    mxyz : mass and X, Y and Z coordinates (columns of an array)
    
    Returns
    -------
    mass : float
        Total mass.

    cogx : float
        Centre of gravity X coordinate.

    cogy : float
        Centre of gravity Y coordinate.

    cogz : float
        Centre of gravity Z coordinate.
    """
    mass = np.sum(mxyz[:,0])
    cogx = np.sum(mxyz[:,0] * mxyz[:,1]) / mass
    cogy = np.sum(mxyz[:,0] * mxyz[:,2]) / mass
    cogz = np.sum(mxyz[:,0] * mxyz[:,3]) / mass
    
    return mass, cogx, cogy, cogz


def Ipolygon(mass, points):
    """Mass moment of inertia for a plane polygon.
    
    Parameters
    ----------
    mass : float
        Total mass (uniformly distributed).

    points : array_like
        Counterclockwise sequence of points (closed contour).
    
    Returns
    -------
    out : float
        Mass moment of inertia
    
    References
    ----------
    http://en.wikipedia.org/wiki/List_of_moments_of_inertia
    """
    num = den = 0.
    
    for n in range(len(points)-1):
        temp0 = np.inner(points[n], points[n])
        temp1 = np.inner(points[n+1], points[n])
        temp2 = np.inner(points[n+1],points[n+1])
        num += np.cross(points[n+1], points[n]) * (temp0 + temp1 + temp2)
        den += np.cross(points[n+1], points[n])
    
    out = mass / 6. * num / den

    return out


def csv_en2pt(csv_en, csv_pt):
    """Translate english to portuguese decimal separator.
    
    Parameters
    ----------
    csv_en - str

    csv_pt - str
    
    Returns
    -------
    None
    """
    en2pt = str.maketrans('.,', ',;')
    
    with open(csv_en, mode='rb') as file_en:
        
        with open(csv_pt, mode='wb') as file_pt:
            
            for line in file_en:
                file_pt.write(line.translate(en2pt))


def csv_pt2en(csv_pt, csv_en):
    """Translate portuguese to english decimal separator.
    
    Parameters
    ----------
    csv_pt - str

    csv_en - str
    
    Returns
    -------
    None
    """
    pt2en = str.maketrans(',;', '.,')
    
    with open(csv_pt, mode='rb') as file_pt:
        
        with open(csv_en, mode='wb') as file_en:
            
            for line in file_pt:
                file_en.write(line.translate(pt2en))


def decimalmark(filename, old=',', new='.'):
    """Replace decimal mark.
    
    Parameters
    ----------
    name : str
        File name.

    old : char, optional
        Old decimal mark.

    new : char, optional
        New decimal mark.
    
    Returns
    -------
    None
    """
    stream = StringIO()
    
    with open(filename, 'r') as ifile:
        
        for line in ifile:
            stream.write(line.replace(old,new))
        
    with open(filename, 'w') as ofile:
        ofile.writelines(stream.getvalue())
    
    stream.close()


def nextpow2(n: int) -> float:
    """Next power of two of number `n`.
    
    Parameters
    ----------
    n : int or float
        Positive integer number.
    
    Returns
    -------
    out : int
        Next power of two of number `n`.
    
    Raises
    ------
    ValueError
        If `n` is not positive.
    
    See also
    --------
    math.ceil
    math.log2
    
    References
    ----------
    http://docs.obspy.org/_modules/obspy/signal/util.html
    
    Examples
    --------
    >>> nextpow2(5)
    8
    >>> nextpow2(250)
    256
    """
    # parameter validation
    if n <= 0:
        raise ValueError(f'n must be positive ({n}).')
    
    out = 2 ** int(math.ceil(math.log2(n)))

    return out


def sma(a, window) -> np.array:
    """Compute a simple moving average (sma) of `a`.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose `sma` is desired.

    window : int
        Size of moving window.
    
    Returns
    -------
    out : ndarray
        Array of `sma` values.
    
    See also
    --------
    `numpy.convolve`
    
    Notes
    -----
    Based on convolution.
    
    References
    ----------
    https://en.wikipedia.org/wiki/Moving_average
    """
    weights = np.ones(window)/window # normalise weights
    out = np.convolve(a, weights, mode='valid')

    return out


def ema(a, window) -> np.array:
    """Compute an exponential moving average (ema) of `a`.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose `ema` is desired.

    window : int
        Size of moving window.
    
    Returns
    -------
    out : ndarray
        Array of `ema` values.
    
    Warnings
    --------
    NOT TESTED YET!
    """
    weights = np.exp(np.linspace(-1., 0., window))
    weights /= weights.sum()
    out = np.convolve(a, weights, mode='full')[:len(a)]
    out[:window] = out[window]
    
    return out


def maxabs(a, axis=None, out=None, keepdims=np._NoValue) -> float:
    """Return the absolute maximum of an array or along an axis.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose absolute maximum is desired.
        If `a` is not an array, a conversion is attempted.
    
    axis : None or int or tuple of ints, optional
        Axis or axes along which the absolute maximum are computed.
        The default is to compute the mean of the flattened array.
        If this is a tuple of ints, a mean is performed over multiple axes,
        instead of a single axis or all the axes as before.
    
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.
        See `doc.ufuncs` for details.
    
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the input array.
        If the default value is passed, then `keepdims` will not be
        passed through to the `mean` method of sub-classes of
        `ndarray`, however any non-default value will be.  If the
        sub-classes `sum` method does not implement `keepdims` any
        exceptions will be raised.
    
    Returns
    -------
    maxabs : ndarray or scalar
        Maximum absolute value of array `a`.
    
    See also
    --------
    numpy.amax
    numpy.absolute
    
    Examples
    --------
    >>> maxabs([-2, -1, 0, 1])
    2
    """
    # parameter validation
    a = np.asarray(a)

    # compute result
    out = np.amax(np.absolute(a), axis=axis, out=out, keepdims=keepdims)
    
    return out


def minabs(a, axis=None, out=None, keepdims=np._NoValue) -> float:
    """Return the absolute minimum of an array or along an axis.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose absolute minimum is desired.
        If `a` is not an array, a conversion is attempted.
    axis : None or int or tuple of ints, optional
        Axis or axes along which the absolute maximum are computed.
        The default is to compute the mean of the flattened array.
        If this is a tuple of ints, a mean is performed over multiple axes,
        instead of a single axis or all the axes as before.
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.
        See `doc.ufuncs` for details.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the input array.
        If the default value is passed, then `keepdims` will not be
        passed through to the `mean` method of sub-classes of
        `ndarray`, however any non-default value will be.  If the
        sub-classes `sum` method does not implement `keepdims` any
        exceptions will be raised.
    
    Returns
    -------
    minabs : ndarray or scalar
        Maximum absolute value of array `a`.
    
    See also
    --------
    numpy.amin
    numpy.absolute
    
    Examples
    --------
    >>> minabs([-2, -1, 0, 1])
    0
    """
    # parameter validation
    a = np.asarray(a)

    # compute result
    out = np.amin(np.absolute(a), axis=axis, out=out, keepdims=keepdims)
    
    return out


def srss(a, axis=None, dtype=None, out=None, keepdims=np._NoValue) -> np.array:
    """Compute the square root of sum of squares (srss) of `a` along the specified axis.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose `srss` is desired.
        If `a` is not an array, a conversion is attempted.
    
    axis : None or int or tuple of ints, optional
        Axis or axes along which the means are computed.
        The default is to compute the mean of the flattened array.
        If this is a tuple of ints, a mean is performed over multiple axes,
        instead of a single axis or all the axes as before.

    dtype : data-type, optional
        Type to use in computing the mean.  For integer inputs, the default
        is `float64`; for floating point inputs, it is the same as the
        input dtype.

    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.
        See `doc.ufuncs` for details.

    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the input array.
        If the default value is passed, then `keepdims` will not be
        passed through to the `sum` method of sub-classes of
        `ndarray`, however any non-default value will be.  If the
        sub-class' method does not implement `keepdims` any exceptions
        will be raised.

    Returns
    -------
    out : ndarray, see dtype parameter above
        If `out=None`, returns a new array containing the `srss` values,
        otherwise a reference to the output array is returned.

    See also
    --------
    numpy.sqrt
    numpy.sum
    numpy.square
    
    Notes
    -----
    The `srss` function is implemented based on `numpy` functions `sqrt`, `sum` and `square`.
    `float64` intermediate and return values are used for integer inputs.

    Examples
    --------
    >>> srss(np.arange(10))
    16.881943016134134
    """
    # parameter validation
    a = np.asarray(a)

    # compute result
    out = np.sqrt(np.sum(np.square(a), axis=axis, dtype=dtype, out=out, keepdims=keepdims))

    return out


def rms(a, axis=None, dtype=None, out=None, keepdims=np._NoValue) -> np.array:
    """Compute the root mean square (rms) of `a` along the specified axis.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose `rms` is desired.
        If `a` is not an array, a conversion is attempted.

    axis : None or int or tuple of ints, optional
        Axis or axes along which the means are computed.
        The default is to compute the mean of the flattened array.
        If this is a tuple of ints, a mean is performed over multiple axes,
        instead of a single axis or all the axes as before.

    dtype : data-type, optional
        Type to use in computing the mean. For integer inputs, the default
        is `float64`; for floating point inputs, it is the same as the
        input dtype.

    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.
        See `doc.ufuncs` for details.

    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the input array.
        If the default value is passed, then `keepdims` will not be
        passed through to the `mean` method of sub-classes of
        `ndarray`, however any non-default value will be.  If the
        sub-class' method does not implement `keepdims` any exceptions
        will be raised.
    
    Returns
    -------
    out : ndarray, see dtype parameter above
        If `out=None`, returns a new array containing the `rms` values,
        otherwise a reference to the output array is returned.
    
    See also
    --------
    numpy.sqrt
    numpy.mean
    numpy.square
    
    Notes
    -----
    The `rms` function is implemented based on `numpy` functions `sqrt`, `mean` and `square`.
    `float64` intermediate and return values are used for integer inputs.
    
    References
    ----------
    http://en.wikipedia.org/wiki/Root_mean_square
    
    Examples
    --------
    >>> rms(np.arange(10))
    5.3385391260156556
    """
    # parameter validation
    a = np.asarray(a)

    # compute result
    out = np.sqrt(np.mean(np.square(a), axis=axis, dtype=dtype, out=out, keepdims=keepdims))
    
    return out


def crest_factor(a, axis=None, dtype=None, out=None, keepdims=np._NoValue) -> np.array:
    """Compute the crest factor of `a`.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose crest factor is desired.
        If `a` is not an array, a conversion is attempted.
    
    axis : None or int or tuple of ints, optional
        Axis or axes along which the absolute maximum are computed.
        The default is to compute the mean of the flattened array.
        If this is a tuple of ints, a mean is performed over multiple axes,
        instead of a single axis or all the axes as before.

    dtype : data-type, optional
        Type to use in computing the `rms`.  For integer inputs, the default
        is `float64`; for floating point inputs, it is the same as the
        input dtype.
    
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.
        See `doc.ufuncs` for details.
    
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the input array.
        If the default value is passed, then `keepdims` will not be
        passed through to the `mean` method of sub-classes of
        `ndarray`, however any non-default value will be.  If the
        sub-classes `sum` method does not implement `keepdims` any
        exceptions will be raised.
    
    Returns
    -------
    out : float
        Crest factor of `a`.
    """
    # parameter validation
    a = np.asarray(a)

    # compute result
    out = maxabs(a, axis=axis, out=out, keepdims=keepdims) / rms(a, axis=axis, dtype=dtype, out=out, keepdims=keepdims)

    return out


def cumulative_mean(a, axis=None, dtype=None, out=None) -> np.array:
    """Calculate the cumulative mean of `a`.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose cumulative mean is desired.
    
    axis : int, optional
        Axis along which the cumulative sum is computed.
        The default (None) is to compute the cumsum over the flattened array.
    
    dtype : dtype, optional
        Type of the returned array and of the accumulator in which the elements are summed.
        If dtype is not specified, it defaults to the dtype of a, unless a has an integer dtype with a precision less than that of the default platform integer.
        In that case, the default platform integer is used.
    
    out : ndarray, optional
        Alternative output array in which to place the result.
        It must have the same shape and buffer length as the expected output but the type will be cast if necessary.
        See Output type determination for more details.
    
    Returns
    -------
    out : ndarray
        Array of cumulative mean values.
    
    See also
    --------
    https://numpy.org/doc/stable/reference/generated/numpy.cumsum.html
    
    Notes
    -----
    Based on `numpy.cumsum`.
    
    Examples
    --------
    """
    # parameter validation
    a = np.asarray(a)

    if axis is None:
        N = a.size
    else:
        N = a.shape[axis]
    
    # compute result
    out = np.cumsum(a, axis=axis, dtype=dtype, out=out) / (1 + np.arange(N))
    
    return out


def cumulative_rms(a, axis=None, dtype=None, out=None) -> np.array:
    """Calculate a cumulative rms.
    
    Parameters
    ----------
    a : array_like
        Array containing numbers whose cumulatve rms is desired.
    
    axis : int, optional
        Axis along which the cumulative sum is computed.
        The default (None) is to compute the cumsum over the flattened array.
    
    dtype : dtype, optional
        Type of the returned array and of the accumulator in which the elements are summed.
        If dtype is not specified, it defaults to the dtype of a, unless a has an integer dtype with a precision less than that of the default platform integer.
        In that case, the default platform integer is used.
    
    out : ndarray, optional
        Alternative output array in which to place the result.
        It must have the same shape and buffer length as the expected output but the type will be cast if necessary.
        See Output type determination for more details.
    
    Returns
    -------
    out : ndarray
        Array of cumulative rms values.
    
    See also
    --------
    https://numpy.org/doc/stable/reference/generated/numpy.cumsum.html
    
    Notes
    -----
    Based on `numpy.cumsum`.
    
    Examples
    --------
    """
    # parameter validation
    a = np.asarray(a)

    if axis is None:
        N = a.size
    else:
        N = a.shape[axis]
    
    # compute result
    out = np.sqrt(np.cumsum(np.square(a), axis=axis, dtype=dtype, out=out) / (1 + np.arange(N)))
    
    return out


def rms_filter(a, window_length, mode: str = 'same') -> np.array:
    """Calculate a root mean square filter.
    
    Parameters
    ----------
    a : array_like
        Array of numbers whose rms filter is desired.
    window_length : float
        Size of the window.
    mode : str
        The mode parameter determines how the input array is extended beyond its boundaries.
    
    Returns
    -------
    out : ndarray
        Array of rms filtered values.
    
    See also
    --------
    https://numpy.org/doc/stable/reference/generated/numpy.convolve.html
    
    Reference
    ---------
    https://stackoverflow.com/questions/8245687/numpy-root-mean-squared-rms-smoothing-of-a-signal
    """
    # parameter validation
    a = np.asarray(a)
    
    if np.iscomplexobj(a):
        raise TypeError(f'Complex type not supported ({a.dtype}).')
    
    if window_length < 1:
        raise RuntimeError(f'Incorrect filter size ({window_length}).')
    
    # compute result
    a2 = np.square(a)
    w1 = np.ones(window_length) / float(window_length)
    out = np.sqrt(np.convolve(a2, w1, mode=mode))
    
    return out


def interpxy(x, xp, yp, left=None, right=None, xscale='lin', yscale='lin'):
    """Interpolation using different scales in the x and y axes.
    
    Parameters
    ----------
    x : array_like
        The x-coordinates of the interpolated values.
    xp : array_like
        The x-coordinates of the data points.
    yp : array_like
        The y-coordinates of the data points, same length as xp.
    left : float, optional
        Value to return for x < xp[0], default is yp[0].
    right : float, optional
        Value to return for x > xp[-1], default is yp[-1].
    xscale : str, optional
        Valid options are 'lin', 'log'. Default is 'lin'.
    yscale : str, optional
        Valid options are 'lin', 'log'. Default is 'lin'.
    
    Returns
    -------
    y : ndarray
        The y-coordinates of the interpolated series.
    
    See also
    --------
    numpy.interp

    Notes
    -----
    Based on numpy.interp, it allows different scales on both axes.
    
    Examples
    --------
    >>> 
    """
    # parameter validation
    xp = np.asarray(xp)
    yp = np.asarray(yp)
    
    assert np.all(np.diff(xp) > 0), 'xp is not increasing'
    
    if xp.size != yp.size:
        raise ValueError(f'`xp` and `yp` have different sizes ({xp.size, yp.size})')
    
    if xscale == 'log':
        x = np.log10(x)
        xp = np.log10(xp)
    elif xscale != 'lin':
        raise ValueError(f'Wrong `xscale` value ({xscale})')

    if yscale == 'log':
        yp = np.log10(yp)
    elif yscale != 'lin':
        raise ValueError(f'Wrong `yscale` value ({yscale})')
    
    # compute result
    y = np.interp(x, xp, yp, left, right)
    
    if yscale == 'log':
        y = np.power(10, y)
    
    return y


def rangenum(xp, xscale='lin', num=50, endpoint=True):
    """Interpolated range.
    
    Parameters
    ----------
    xp : array_like
        The x-coordinates of the data points.
    xscale : str, optional
        Valid options are 'lin', 'log'. Default is 'lin'.
    num : int, array_like, optional
        The number of interpolated values in each range. Defaults to 50
    endpoint : boolean, optional
        If true, `xp[-1]` is the last sample. Otherwise, it is not included.
        Default is True.
    
    Returns
    -------
    x : ndarray
        The x-coordinates of the interpolated range.
    
    See also
    --------
    numpy.linspace
    numpy.geomspace

    Notes
    -----
    Based on `numpy.linspace` and `numpy.geomspace`.
    
    Examples
    --------
    >>> rangenum([1, 2], num=10)
    array([1. , 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. ])
    """
    # parameter validation
    xp = np.asarray(xp)
    
    assert np.all(np.diff(xp) > 0), 'xp is not increasing'
    
    if isinstance(num, Number):
        num *= np.ones(xp.size-1, dtype='int')
    elif len(num) != xp.size-1:
        raise ValueError(f'Wrong `num` array_like length ({len(num)})')
    
    if xscale == 'lin':
        fun = np.linspace
    elif xscale == 'log':
        fun = np.geomspace
    else:
        raise ValueError(f'Wrong `xscale` value ({xscale})')
    
    # compute result
    x = []
    
    for x0, x1, n in zip(xp[:-1], xp[1:], num):
        x.extend(fun(x0, x1, num=n, endpoint=False))
    
    if endpoint == True:
        x.append(xp[-1])
    
    return np.asarray(x)


def rangexy(xp, yp, left=None, right=None, xscale='lin', yscale='lin', num=50,
            endpoint=True):
    """Interpolated ranges using different scales in the x and y axes.
    
    Parameters
    ----------
    xp : array_like
        The x-coordinates of the data points.
    yp : array_like
        The y-coordinates of the data points, same length as xp.
    left : float, optional
        Value to return for x < xp[0], default is yp[0].
    right : float, optional
        Value to return for x > xp[-1], default is yp[-1].
    xscale : str, optional
        Valid options are 'lin', 'log'. Default is 'lin'.
    yscale : str, optional
        Valid options are 'lin', 'log'. Default is 'lin'.
    num : int, array_like, optional
        The number of interpolated values in each range. Defaults to 50
    endpoint : boolean, optional
        If true, `xp[-1]` is the last sample. Otherwise, it is not included.
        Default is True.
    
    Returns
    -------
    x : ndarray
        The x-coordinates of the interpolated range.
    y : ndarray
        The y-coordinates of the interpolated ranges series.
    
    See also
    --------
    numpy.linspace
    numpy.geomspace

    Notes
    -----
    Based on `rangenum` and `interpxy`, it allows different scales on both axes.
    
    Examples
    --------
    >>> 
    """
    # parameter validation
    xp = np.asarray(xp)
    yp = np.asarray(yp)
    
    # compute result
    x = rangenum(xp, xscale=xscale, num=num, endpoint=endpoint)
    y = interpxy(x, xp, yp, left=left, right=right,
                 xscale=xscale, yscale=yscale)
    
    return x, y


def ULenvelope(data, ref_level=0):
    """Return upper and lower envelope of data above and below ref_level.

    The envelopes are in a peak and hold format.
    
    Parameters
    ----------
    data : array_like
        Data array.
    ref_level : float, optional
        Reference level. Defaults to 0.
    
    Returns
    -------
    upper_env : ndarray
        Array of upper envelope steps.
    lower_env : ndarray
        Array of lower envelope steps.
    
    See also
    --------
    scipy.signal.find_peaks
    
    Notes
    -----
    Algorithm is based on a forward-backward traversal of data.
    It is able to deal with a single wide peak (data increases and then decreases)
    but not with a succession of peaks and valleys.
    
    References
    ----------
    https://www.groundai.com/project/an-heuristic-approach-to-obtain-signal-envelope-with-a-simple-software-implementation/2
    https://ryukau.github.io/filter_notes/peak_hold_envelope/peak_hold_envelope.html
    https://gist.github.com/aanastasiou/480d81361abcdc794783
    
    Examples
    --------
    >>>> ULenvelope([0., 1., 0.5, 2., -0.5, 0., -1., 0.])
    (array([0., 1., 1., 2., 0., 0., 0., 0.]),
     array([ 0. ,  0. ,  0. ,  0. , -0.5, -0.5, -1. ,  0. ]))
    """
    # parameter validation
    data = np.asarray(data)

    # initial assignements
    data_above = np.maximum(data, ref_level)
    data_below = np.minimum(data, ref_level)
    upper_env = np.maximum(data, ref_level)
    lower_env = np.minimum(data, ref_level)

    for n in range(data.size-1): # forward traversal

        if data_above[n+1] < data_above[n]:
            data_above[n+1] = data_above[n]

        if data_below[n+1] > data_below[n]:
            data_below[n+1] = data_below[n]

    for n in range(data.size-1): # backward traversal

        if upper_env[-n-2] < upper_env[-n-1]:
            upper_env[-n-2] = upper_env[-n-1]

        if lower_env[-n-2] > lower_env[-n-1]:
            lower_env[-n-2] = lower_env[-n-1]

    # compute result
    upper_env = np.minimum(data_above, upper_env)
    lower_env = np.maximum(data_below, lower_env)

    return upper_env, lower_env


class Shapes():
    """Build continuous time signal from simple shapes."""
    def __init__(self, dt):
        """..."""
        self.dt = dt
        self.height = 0
        self._signal = []

    def __call__(self):
        """..."""
        signal = np.hstack(self._signal)
        time = np.arange(signal.size) * self.dt
        
        return (time, signal)
    
    def step(self, duration, height=0):
        """..."""
        signal = np.ones(duration/self.dt) * height
        self._signal.append(signal)
        self.height = height
    
    def ramp(self, duration, slope=None, height=None):
        """..."""
        if slope is None:
            slope = height / duration
        elif height is None:
            height = slope * duration
        else:
            Exception('')
        
        signal = np.linspace(0., 1., num=duration/self.dt, endpoint=False)
        self._signal.append(signal)
        self.height = height


def main():
    """Function to run tests."""
    pass


if __name__ == '__main__': # run tests if called from command-line
    main()
