#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""LNECSPA.py - LNEC Signal Processing and Analysis binary database format

Created on Wed Mar 19 18:15:54 2014

@author: Paulo Xavier Candeias
"""

# Developed for Python 3.12.7
__author__     = 'Paulo Xavier Candeias'
__copyright__  = 'Copyright 2013, Paulo Xavier Candeias'
__credits__    = ['Paulo Xavier Candeias']
__license__    = 'GPL-3.0'
__version__    = '2024-10-19'
__maintainer__ = 'Paulo Xavier Candeias'
__email__      = 'pxcandeias@gmail.com'
__status__     = 'Production' # one of "Prototype", "Development", or "Production"

# Python standard library
import os.path
from datetime import datetime

# 3rd party modules
import numpy as np

# pxc modules
from LVBF import LV_fd
from LabVIEW import LV_to_integer


class Signal():
    """..."""
    def __init__(self, super_self):
        """..."""
        self.super_self = super_self
    
    def __getitem__(self, key):
        """..."""
        return {'dt': self.super_self.t0[key], 'name': self.super_self.names[key]}


class LTFdb(LV_fd):
    """LTF database.
    
    LNEC Test File attributes
    -------------------------
    AQTD : str
        LabVIEW string.
    stage : str
        LabVIEW string.
    obs : str
        LabVIEW string.
    date : str
        LabVIEW string.
    t0 : ndarray
        LabVIEW array of timestamp.
    dt : ndarray
        LabVIEW array of float64.
    cluster : list of tuples [(scale, offset, digital),]
        scale : float
            LabVIEW float64.
        offset : float
            LabVIEW float64.
        digital : ndarray
            LabVIEW 1D array of int16.
    dataformat : ndarray
        LabVIEW array of int32.
    IDstring : list of strings
        LabVIEW array of string.
    scalefactor : ndarray
        LabVIEW array of float64.
    offset : ndarray
        LabVIEW array of float64.
    names : list of strings
        LabVIEW array of string.
    units : list of strings
        LabVIEW array of string.
    types : list of strings
        LabVIEW array of string.
    info : list of strings
        LabVIEW array of string.
    
    Data attributes
    ---------------
    _filename : str
        Name of file to read/write.
    _eof : bool
        End-of-file indicator (True/False).
    _data : list of ndarray
        Analog data as list of 1D `ndarray`.
    data : ndarray
        Analog data as 2D `ndarray`. Deprecation warning: to be replaced by `_data`!
    
    Properties
    ----------
    filename : LTF filename
    eof : LTF end-of-file has been reached?
    metadata : `LTFdb` metadata
    signals : `LTFdb` number of channels
    samples : `LTFdb` number of samples per channel
    
    Methods
    -------
    __init__ :
    __str__ :
    read : read LTF file
    write : write LTF file
    keys : channels keys
    values : channels values
    items : channels items
    copy : 
    update : update `LTFdb` in memory
    append : append a new channel at the end of `LTFdb`
    index : index channel names
    extend : extend with data coming from a new LTF file or from another `LTFdb` instance
    __len__ :
    __getitem__ :
    __setitem__ :
    __deltitem__ :
    __missing__ :
    __iter__ :
    __contains__ :
    read_cluster : read LTF cluster
    write_cluster : write LTF cluster
    time : returns 1D time array
    print_metadata : print metadata
    print_data : print data
    write_txt : write text file
    waveforms : return LabVIEW waveforms (t0, dt, array)
    dac : digital-to-analog conversion
    adc : analog-to-digital conversion
    """
    def __init__(self, endian='>', encoding='cp1252'):
        """Initialize class instance.
        
        Parameters
        ----------
        endian : str, optional
            Defaults to '>'
        encoding : str, optional
            Defaults to 'cp1252'
        
        Notes
        -----
        For endianness test:
            On Windows 10
                CLI ... returns ...
                Python console 'sys.byteorder' returns 'little'-
            On Linux
                CLI 'lscpu | grep Endian' returns 'Byte Order: Little Endian'
                Python console 'sys.byteorder' returns 'little'.
        For encoding test:
            On Windows 10
                CLI 'chcp' returns 'cp1252'.
                Python console 'sys.stdin.encoding' returns ...
            On Linux
                CLI 'locale | grep LANG' returns 'LANG=pt_PT.UTF-8'
                Python console 'sys.stdin.encoding' returns 'UTF-8'.
        """
        super().__init__(endian, encoding)
        
        # LTF data attributes

        # AQ terminal data
        self.AQTD = ['AQTD']
        self.stage = ['Stage']
        self.obs = ['Obs']
        self.date = ['Date']

        # Initial time (counting from 01-01-1904, 00:00:00.00, LabVIEW date/time origin)
        self.t0 = np.zeros(1, dtype=self.LVtimestamp)
        self.dt = np.ones(1, dtype=self.LVfloat64)

        # Digital time series
        self.cluster = [(1., 0., np.zeros(1, dtype=self.LVint16))]

        # Channel info
        self.dataformat = np.zeros(1, dtype=self.LVint32)
        self.IDstring = ['IDstring']
        self.scalefactor = np.ones(1, dtype=self.LVfloat64)
        self.offset = np.array([0.], dtype=self.LVfloat64)
        self.names = ['names']
        self.units = ['units']
        self.types = ['types']
        self.info = ['info']
        
        # private data attributes
        self._filename = None
        self._eof = None
        self._data = None

        # EXPERIMENTAL!!!
        self._signal = Signal(self)

    # `LTFdb` properties
    filename = property(lambda self: self._filename)
    eof = property(lambda self: self._eof)
    metadata = property(lambda self: [self.AQTD,
                                      self.stage,
                                      self.obs,
                                      self.date,
                                      self.t0,
                                      self.dt,
                                      self.cluster,
                                      self.dataformat,
                                      self.IDstring,
                                      self.scalefactor,
                                      self.offset,
                                      self.names,
                                      self.units,
                                      self.types,
                                      self.info,
                                     ])
    signals = property(lambda self: self.__len__())
    samples = property(lambda self: np.array([v.size for v in self.data]))
    
    def __str__(self):
        """LTF string representation."""
        return f'Filename: {self.filename}\nAQTD: {self.AQTD}\nStage: {self.stage}\nObs: {self.obs}\nDate: {self.date}\nSignals: {self.__len__()}'
    
    def read(self, filename):
        """Read LTF file.
        
        Parameters
        ----------
        filename : str
            LTF file name to read.
            The value is stored in the attribute `self._filename`.
        
        Returns
        -------
        self

        Notes
        -----
        The digital to analog conversion (`dac`) of the `cluster` to obtain `_data` is performed after reading.
        """
        self._filename = filename
        
        with open(filename, 'rb') as self.fobj:
            self.AQTD = self.read_string()
            self.stage = self.read_string()
            self.obs = self.read_string()
            self.date = self.read_string()
            self.t0 = self.read_array(self.read_timestamp, self.LVtimestamp)
            self.dt = self.read_array(self.read_numeric, self.LVfloat64)
            self.cluster = self.read_cluster(self.LVint16)
            self.dataformat = self.read_array(self.read_numeric, self.LVint32)
            self.IDstring = self.read_array(self.read_string)
            self.scalefactor = self.read_array(self.read_numeric, self.LVfloat64)
            self.offset = self.read_array(self.read_numeric, self.LVfloat64)
            self.names = self.read_array(self.read_string)
            self.units = self.read_array(self.read_string)
            self.types = self.read_array(self.read_string)
            self.info = self.read_array(self.read_string)
            self._eof = self.EOD()
        
        # Analog time series
        self._data = self.dac()
        self.data = np.asarray(self._data) # Deprecation warning: to be replaced by `_data`!
        
        return self

    def write(self, filename):
        """Write LTF file.
        
        Parameters
        ----------
        filename : str
            LTF file name to write.
            The `self._filename` attribute is updated.

        Notes
        -----
        The analog to digital conversion (`adc`) of the `_data` to obtain `cluster` is performed before writing.
        """
        self._filename = filename

        # Update digital cluster
        self.cluster = self.adc(self.LVint16)
        
        with open(filename, 'wb') as self.fobj:
            self.write_string(self.AQTD)
            self.write_string(self.stage)
            self.write_string(self.obs)
            self.write_string(self.date)
            self.write_array(self.write_timestamp, self.t0)
            self.write_array(self.write_numeric, self.dt)
            self.write_cluster(self.cluster, self.LVint16)
            self.write_array(self.write_numeric, self.dataformat)
            self.write_array(self.write_string, self.IDstring)
            self.write_array(self.write_numeric, self.scalefactor)
            self.write_array(self.write_numeric, self.offset)
            self.write_array(self.write_string, self.names)
            self.write_array(self.write_string, self.units)
            self.write_array(self.write_string, self.types)
            self.write_array(self.write_string, self.info)

    def keys(self):
        """`LTFdb` channel keys.
        
        Returns
        -------
        out : list
            List of keys.
        """
        return ['t0',
                'dt',
                'cluster',
                'dataformat',
                'IDstring',
                'scalefactor',
                'offset',
                'names',
                'units',
                'types',
                'info',
                '_data',
                'data', # Deprecation warning: to be replaced by `_data`!
               ]
    
    def values(self):
        """`LTFdb` channel values.
        
        Returns
        -------
        out : list
            List of values.
            See `keys()` for existing keys.
        """
        return [self.__dict__[key] for key in self.keys()]
    
    def items(self):
        """`LTFdb` channel items.
        
        Returns
        -------
        out : list
            List of (key, value) pairs.
            See `keys()` for existing keys.
        """
        return [(key, self.__dict__[key]) for key in self.keys()]    
    
    def copy(self, selection=None):
        """Make a copy of `LTFdb` with selected channels.
        
        Parameters
        ----------
        selection : array_like or None, optional
            Channels to copy to the new `LTFdb` instance.
            If `None`, all channels will be copied.
            Default is `None`.

        Returns
        -------
        ltf : `LTFdb`
            New `LTFdb` instance.
        """
        # Update digital cluster
        self.cluster = self.adc(self.LVint16)

        if selection is None:
            selection = range(self.__len__())

        ltf = LTFdb(endian=self.endian, encoding=self.encoding)
        ltf._filename = self._filename
        ltf.update(AQTD=self.AQTD,
                    stage=self.stage,
                    obs=self.obs,
                    date=self.date,
                    t0=[self.t0[v] for v in selection],
                    dt=[self.dt[v] for v in selection],
                    cluster=[self.cluster[v] for v in selection],
                    dataformat=[self.dataformat[v] for v in selection],
                    IDstring=[self.IDstring[v] for v in selection],
                    scalefactor=[self.scalefactor[v] for v in selection],
                    offset=[self.offset[v] for v in selection],
                    names=[self.names[v] for v in selection],
                    units=[self.units[v] for v in selection],
                    types=[self.types[v] for v in selection],
                    info=[self.info[v] for v in selection],
                  )
        
        return ltf

    def update(self, **kwargs):
        """Update `LTFdb` in memory.
        
        This method is not suitable to reshape the data, only to update it.
        
        Parameters
        ----------
        **kwargs : dict
            LNEC Test File data atributes.
        
        Raises
        ------
        ValueError
            When both `cluster`, `data` are supplied.

        KeyError
            When key is not valid.

        Notes
        -----
        At most one of `cluster`, `data` can be supplied.
        The other one, either `data` or `cluster`, will be computed.
        """
        if ('cluster' in kwargs) and ('data' in kwargs):
            raise ValueError('At most one of `cluster`, `data` can be supplied.')
        
        for key in kwargs:
            
            if key in ('AQTD', 'stage', 'obs', 'date'):
                value = kwargs[key]
            elif key == 't0':
                value = np.asarray(kwargs[key], dtype=self.LVtimestamp)
            elif key == 'dt':
                value = np.asarray(kwargs[key], dtype=self.LVfloat64)
            elif key == 'cluster':
                value = kwargs[key]
            elif key == 'dataformat':
                value = np.asarray(kwargs[key], dtype=self.LVint32)
            elif key == 'IDstring':
                value = kwargs[key]
            elif key in ('scalefactor', 'offset'):
                value = np.asarray(kwargs[key], dtype=self.LVfloat64)
            elif key in ('names', 'units', 'types', 'info'):
                value = kwargs[key]
            elif key == 'data':
                value = np.asarray(kwargs[key])
            else:
                raise KeyError(f'Invalid key ({key}).')
            
            self.__dict__[key] = value
            
        if 'cluster' in kwargs:
            self._data = self.dac()
            self.data = np.asarray(self._data) # Deprecation warning: to be replaced by `_data`!
        elif 'data' in kwargs:
            self._data = self.data.tolist()
            self.cluster = self.adc(self.LVint16)

    def append(self, t0=None, dt=None, cluster=None, dataformat=None,
               IDstring='IDstring', scalefactor=1., offset=0.,
               names='name', units='unit', types='type', info='info',
               data=None):
        """Append a new channel at the end of the LTF data.
        
        This method reshapes the LTF data.
        
        Parameters
        ----------
        t0 : tuple or `None`, optional
            If None then use `t0` from last channel, otherwise a timestamp
            tuple (0., 0.) should be supplied. Default is `None`.

        dt : float or None, optional
            If None then use `dt` from last channel, otherwise a float value
            should be supplied. Default is `None`.

        cluster : tuple of (`scale`, `offset`, `digital`) or `None`, optional
            Digital data. Default is `None`.

        dataformat : int or None, optional
            Signal data format. Default is `None`.

        IDstring : str, optional
            Signal ID string.

        scalefactor : float, optional
            Signal scale factor.

        offset : float, optional
            Signal offset.

        names : str, optional
            Signal name.

        units : str, optional
            Signal units.

        types : str, optional
            Signal type.

        info : str, optional
            Signal info.

        data : array_like or None, optional
            Analog data. Default is `None`.
        
        Raises
        ------
        ValueError
            When none of `cluster`, `_data` are supplied.

        Notes
        -----
        At least one of `cluster`, `data` must be supplied.
        When both are supplied, `cluster` takes precedence over `data`.
        The other one, either `data` or `cluster`, will be computed.
        """
        kwargs = {}
        kwargs['t0'] = self.t0[-1] if t0 is None else t0
        kwargs['dt'] = self.dt[-1] if dt is None else dt
        if dataformat is not None: kwargs['dataformat'] = dataformat
        if IDstring is not None: kwargs['IDstring'] = [IDstring]
        if scalefactor is not None: kwargs['scalefactor'] = scalefactor
        if offset is not None: kwargs['offset'] = offset
        if names is not None: kwargs['names'] = [names]
        if units is not None: kwargs['units'] = [units]
        if types is not None: kwargs['types'] = [types]
        if info is not None: kwargs['info'] = [info]
        
        if cluster is not None:
            kwargs['cluster'] = cluster
        elif data is not None:
#            kwargs['_data'] = data
            kwargs['data'] = data[np.newaxis,:] # Deprecation warning: to be replaced by `_data`!
        else:
            raise ValueError('At least one of `cluster`, `data` must be supplied.')
        
        ltf = LTFdb(endian=self.endian, encoding=self.encoding)
        ltf.update(**kwargs)
        self.extend(ltf=ltf)
        del ltf

    def index(self, x, start=None, end=None):
        """Index channel name.
        
        Parameters
        ----------
        x : str

        start : int or `None`, optional

        end : int or `None`, optional
        
        Returns
        -------
        """
        if start is None:
            start = 0
        
        if end is None:
            end = self.__len__()
        
        return self.names.index(x, start, end)

    def extend(self, filename=None, ltf=None, update_header=False):
        """Extend LTF data with another LTF file or `LTFdb` instance.
        
        This method reshapes the data.
        
        Parameters
        ----------
        filename : str or `None`, optional
            LTF data file name to extend.
            Default is `None`.

        ltf : `LTFdb` instance or None, optional
            `LTFdb` data to extend.
            Default is `None`.

        update_header : bool, optional
            If `True` update the AQ terminal data.
            Default is `False`.
        
        Raises
        ------
        ValueError
            When none of `filename`, `ltf` are supplied.
        """
        if filename is not None:
            ltf = LTFdb(endian=self.endian, encoding=self.encoding)
            ltf.read(filename)
        elif ltf is None:
            raise ValueError('Either `filename` or `ltf` must be given')
        
        if update_header:
            self.AQTD = ltf.AQTD
            self.stage = ltf.stage
            self.obs = ltf.obs
            self.date = ltf.date
        
        self.t0 = np.hstack((self.t0, ltf.t0)).astype(self.t0.dtype)
        self.dt = np.hstack((self.dt, ltf.dt)).astype(self.dt.dtype)
        self.cluster.extend(ltf.cluster)
        self.dataformat = np.hstack((self.dataformat, ltf.dataformat)).astype(self.dataformat.dtype)
        self.IDstring.extend(ltf.IDstring)
        self.scalefactor = np.hstack((self.scalefactor, ltf.scalefactor)).astype(self.scalefactor.dtype)
        self.offset = np.hstack((self.offset, ltf.offset)).astype(self.offset.dtype)
        self.names.extend(ltf.names)
        self.units.extend(ltf.units)
        self.types.extend(ltf.types)
        self.info.extend(ltf.info)
        self._data.extend(ltf._data)
        self.data = np.vstack((self.data, ltf.data)).astype(self.data.dtype) # Deprecation warning: to be replaced by `_data`!
        
        if filename is not None:
            del ltf

    def __len__(self):
        """Number of channels."""
        return len(self.names)

    def __getitem__(self, channel):
        """Get channel data.

        A full dictionary is given according to `self.keys()` plus `data`.
        
        Parameters
        ----------
        channel : int or str
            Number (or name) of the channel to get.
        
        Returns
        -------
        LTF_dict : dict
            Channel data.
            See `keys()` method for existing keys.
        
        Raises
        ------
        KeyError
            When a slice or an invalid channel are used.
            (see https://medium.com/@beef_and_rice/mastering-pythons-getitem-and-slicing-c94f85415e1c)
        """
        if isinstance(channel, str):
            channel = self.index(channel)
        elif isinstance(channel, slice):
            print(channel.start, channel.stop, channel.step)
            raise KeyError(f'Slice objects are not accepted ({channel}).')
        elif channel >= self.__len__():
            raise KeyError(f'Invalid channel ({channel}).')

        LTF_dict = {key: self.__dict__[key][channel] for key in self.keys()} # does not include `_data`
        LTF_dict['data'] = self._data[channel]

        return LTF_dict

    def __setitem__(self, channel, LTF_dict):
        """Set channel data.

        A full dictionary should be given according to `self.keys()` plus `data`.

        Parameters
        ----------
        channel : int or str
            Number (or name) of the channel to set.

        LTF_dict : dict
            Channel data.
            See `keys()` method for required keys.
        
        Raises
        ------
        KeyError
            When a slice or an invalid channel are used.
            (see https://medium.com/@beef_and_rice/mastering-pythons-getitem-and-slicing-c94f85415e1c)

        Notes
        -----
        Use the `append()` method when only one of `cluster`, `_data` are supplied.
        """
        if isinstance(channel , str):
            channel = self.index(channel)
        elif isinstance(channel, slice):
            print(channel.start, channel.stop, channel.step)
            raise KeyError(f'Slice objects are not accepted ({channel}).')
        elif channel >= self.__len__():
            raise KeyError(f'Invalid channel ({channel}).')

        if len(self.keys()) != len(LTF_dict.keys()):
            raise KeyError('Wrong number of keys.')

        for key, value in LTF_dict.items():

            if key == 'data':
                self._data[channel] = value
            elif key in self.keys():
                self.__dict__[key][channel] = value
            else:
                raise KeyError(f'Keys ({LTF_dict.keys()}) do not match LTF keys.')

    def __delitem__(self, channel):
        """Delete channel data.
        
        Parameters
        ----------
        channel : int or str
            Number (or name) of the channel to delete.
        
        Raises
        ------
        KeyError
            When a slice or an invalid channel are used.
        """
        if isinstance(channel , str):
            channel = self.index(channel)
        elif isinstance(channel, slice):
            print(channel.start, channel.stop, channel.step)
            raise KeyError(f'Slice objects are not accepted ({channel}).')
        elif channel >= self.__len__():
            raise KeyError(f'Invalid channel ({channel}).')

        self.t0 = np.delete(self.t0, channel, 0)
        self.dt = np.delete(self.dt, channel, 0)
        del self.cluster[channel]
        self.dataformat = np.delete(self.dataformat, channel, 0)
        del self.IDstring[channel]
        self.scalefactor = np.delete(self.scalefactor, channel, 0)
        self.offset = np.delete(self.offset, channel, 0)
        del self.names[channel]
        del self.units[channel]
        del self.types[channel]
        del self.info[channel]
        del self._data[channel]
        self.data = np.delete(self.data, channel, 0) # Deprecation warning: to be replaced by `_data`!

    def __missing__(self, channel):
        """Missing channel.
        
        Parameters
        ----------
        channel : int or str
            Number (or name) of the missing channel.
        
        Raises
        ------
        ValueError
        """
        raise ValueError(f'Name ({channel}) not in `LTFdb`.')

    def __iter__(self):
        """Iterate through `LTFdb` channels.
        
        Yields
        ------
        out : dict
            Channel data.
            See `keys()` method for existing keys.
        """
        for channel in range(self.__len__()):
            yield {key: self.__dict__[key][channel] for key in self.keys()}

    def __contains__(self, name):
        """Membership test operator.
        
        Parameters
        ----------
        name : str
            Channel name to lookup.
        
        Returns
        -------
        out : bool
        """
        return name in self.names

    def read_cluster(self, dtype):
        """Read LTF digital array cluster.
        
        Parameters
        ----------
        dtype : numpy.dtype
            Digital data type.
        
        Returns
        -------
        cluster : list of tuples [(scale, offset, digital),]
            scale : LabVIEW float64
            offset : LabVIEW float64
            digital : LabVIEW 1D array of int16
        """
        n = self.read_numeric(self.LVint32)
        cluster = []
        
        for v in range(n):
            (scale, offset) = self.read_numeric(self.LVfloat64, 2)
            digital = self.read_array(self.read_numeric, dtype)
            cluster.append((scale, offset, digital))
        
        return cluster

    def write_cluster(self, cluster, dtype):
        """Write LTF digital array cluster.
        
        Parameters
        ----------
        cluster : list of tuples [(scale, offset, digital),]
            scale : LabVIEW float64
            offset : LabVIEW float64
            digital : LabVIEW 1D array of int16

        dtype : numpy.dtype
            Digital data type.
        """
        self.write_numeric(np.array(len(cluster), self.LVint32))
        
        for (scale, offset, digital) in cluster:
            self.write_numeric(np.array(scale, self.LVfloat64))
            self.write_numeric(np.array(offset, self.LVfloat64))
            self.write_array(self.write_numeric, np.array(digital, dtype=dtype))

    def time(self, column):
        """Time array corresponding to data column.
        
        Parameters
        ----------
        column : int
            Data column to retrieve.
        
        Returns
        -------
        out : numpy array
            time array
        """
        return np.arange(self.samples[column])*self.dt[column]

    def print_metadata(self):
        """Print metadata."""
        print(f'AQTD: {self.AQTD}')
        print(f'Stage: {self.stage}')
        print(f'Obs: {self.obs}')
        print(f'Date: {self.date}')
        print(f't0: {self.t0}')
        print(f'dt: {self.dt}')
        print(f'Cluster: {self.cluster}')
        print(f'Dataformat: {self.dataformat}')
        print(f'IDstring: {self.IDstring}')
        print(f'Scalefactor: {self.scalefactor}')
        print(f'Offset: {self.offset}')
        print(f'Names: {self.names}')
        print(f'Units: {self.units}')
        print(f'Types: {self.types}')
        print(f'Info: {self.info}')

    def print_data(self):
        """Print data."""
        print(f'Data: {self._data}')

    def write_txt(self, filename, sep='\t'):
        """Write text file.
        
        Convenience method for LTFX.
        
        Parameters
        ----------
        filename : str
            File name.

        sep : str, optional
            Column separator.
            Defaults to '\t'.
        """
        with open(filename, 'w') as fobj:
            fobj.write('LNEC Test File eXtractor\n')
            fobj.write('\n')
            fobj.write(f'File: {os.path.basename(self._filename)}\n')
            fobj.write('\n')
            fobj.write(f'Created in: {datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")}\n')
            fobj.write('\n')
            fobj.write(f'Number of channels: {np.size(self.dt)}\n')
            fobj.write('\n')
            fobj.write(f'Name:\t{sep.join(self.names)}\n')
            fobj.write(f'Type:\t{sep.join(self.types)}\n')
            fobj.write(f'Unit:\t{sep.join(self.units)}\n')
            fobj.write(f'No. of Samples:\t{sep.join(str(v) for v in self.samples)}\n')
            fobj.write(f'Time step [s]:\t{sep.join(str(v) for v in self.dt)}\n')
            temp = [datetime.fromtimestamp(self.Unix_time(v)) for v in self.t0]
            fobj.write(f'Trigger date:\t{sep.join([v.strftime("%Y-%m-%d") for v in temp])}\n')
            fobj.write(f'Trigger time:\t{sep.join([v.strftime("%H:%M:%S.%f") for v in temp])}\n')
            fobj.write(f'Duration [s]:\t{sep.join(str((np.size(n)-1)*d) for (n,d) in zip(self._data,self.dt))}\n')
            fobj.write('Data:\n')
            fobj.write(f'{sep.join(self.names)}\n')
            
            for v in range(np.size(self._data[0])):
                fobj.write(f'{sep.join([str(n[v]) for n in self._data])}\n')

    def waveforms(self, channels=None):
        """LabVIEW waveform.
        
        Parameters
        ----------
        channels : iterable or None (default)
            List of channel numbers.
        
        Returns
        -------
        out : tuples of (t0, dt, data)
            LabVIEW waveforms.
        
        References
        ----------
        """
        if channels == None:
            channels = range(self.t0.size)
        
        t0 = self.t0[channels]
        dt = self.dt[channels]
        data = [self._data[v] for v in channels]
        
        return zip(t0, dt, data)

    def dac(self, Arange=20.):
        """Convert LTF digital cluster to analog array.
        
        Parameters
        ----------
        Arange : float, optional
            Full scale analog range (peak to peak amplitude) [V].
        
        Returns
        -------
        analog : list
            List of 1D `ndarray` [physical unit].
        
        Notes
        -----
        THIS CODE HAS BEEN VALIDATED, DO NOT CHANGE IT!

        References
        ----------
        http://zone.ni.com/reference/en-XX/help/371361F-01/lvwave/digital_to_analog_wf/
        """
        analog = []
        
        for (scale, offset, digital) in self.cluster:
            Drange = 2 ** (8 * digital.itemsize) # digital range in bits
            analog.append((digital * Arange / Drange - offset) * scale)
        
        return analog
    
    def adc(self, dtype, Arange=20.):
        """Convert analog array to LTF digital cluster.
        
        Parameters
        ----------
        dtype : numpy.dtype
            Digital data type.

        Arange : float, optional
            Full scale analog range (peak to peak amplitude) [V].
        
        Returns
        -------
        cluster : list of tuples [(scale, offset, digital),]
            scale : LabVIEW float64
            offset : LabVIEW float64
            digital : LabVIEW 1D `ndarray` of `int16`
        
        Notes
        -----
        THIS CODE HAS BEEN VALIDATED, DO NOT CHANGE IT!

        References
        ----------
        http://zone.ni.com/reference/en-XX/help/371361F-01/lvwave/analog_to_digital_wf/
        """
        Drange = 2 ** (8 * dtype.itemsize) # full scale digital range
        cluster = []
        
        for analog in self.data:
            amplitude = abs(analog.max() - analog.min()) # peak to peak amplitude
            
            if amplitude == 0.: # to prevent division by zero
                scale = 1.
                offset = -analog.max()
                digital = np.zeros_like(analog).astype(dtype)
            else:
                scale = amplitude/Arange
                offset = -(analog.max() + analog.min()) / 2. / scale
                digital = LV_to_integer(Drange / Arange * (analog / scale + offset), dtype)
            
            cluster.append((scale, offset, digital))
            
        return cluster


class cc2db(LV_fd):
    """cc2 database.
    
    LNEC-SPA channels configuration file.
    """
    def __init__(self, endian='>', encoding='cp1252'):
        """Initialize instance.
        
        Parameters
        ----------
        endian : str, optional
            Defaults to '>'

        encoding : str, optional
            Defaults to 'cp1252'

        Notes
        -----
        Test for endianness:
            On Windows 10
                CLI ... returns ...
                Python console 'sys.byteorder' returns 'little'-
            On Linux
                CLI 'lscpu | grep Endian' returns 'Byte Order: Little Endian'
                Python console 'sys.byteorder' returns 'little'.
        Test for encoding:
            On Windows 10
                CLI 'chcp' returns 'cp1252'.
                Python console 'sys.stdin.encoding'.
            On Linux
                CLI 'locale | grep LANG' returns 'LANG=pt_PT.UTF-8'
                Python console 'sys.stdin.encoding' returns 'UTF-8'.
        """
        super().__init__(endian, encoding)

    metadata = property(lambda self: ['Name',
                                      'Physical Ch.',
                                      'Sensor',
                                      'Cable(s)',
                                      'Is Active',
                                      'Offset [V]',
                                      'SF Sensor [EGU/V]',
                                      'SF Board Gain Effect',
                                      'SF External Gain',
                                      'SF Final [EGU/V]',
                                      'Max. Value [EGU]',
                                      'Min. Value [EGU]',
                                      'Type',
                                      'Unit',
                                      'Terminal Cfg.',
                                      'Device Type',
                                      'Basic-Res1',
                                      'Basic-Res2',
                                      'Mesh Node',
                                      'Mesh Direc',
                                      'Mesh-Res1',
                                      'Mesh-Res2',
                                      'V1140-Gain',
                                      'V1140-Res1',
                                      'V1140-Res2',
                                      'A153x-ExcitSource',
                                      'A153x-ExcitValue',
                                      'A153x-Gain',
                                      'A153x-LowpassFreq',
                                      'A153x-Res1',
                                      'A153x-Res2',
                                      'V1121-Gain',
                                      'Res.32',
                                      'Res.33',
                                      'Res.34',
                                      'Res.35',
                                      'Res.36',
                                      'Res.37',
                                      'Res.38',
                                      'Res.39',
                                      'Res.40'])
    signals = property(lambda self: self.__len__())
    samples = property(lambda self: np.asarray([v.size for v in self.data]))

    def read(self, filename):
        """Read cc2 cluster from file.
        
        Parameters
        ----------
        filename : str
            cc2 file name to read.
        
        Returns
        -------
        self
        """
        self.filename = filename
        
        with open(filename, 'rb') as self.fobj:
            self.data = self.read_array(self.read_string, ndims=2)
            self.eof = self.EOD()
        
        return self

    def write(self, filename):
        """Write cc2 cluster to file.
        
        Parameters
        ----------
        filename : str
            cc2 file name to write.
        """
        self.filename = filename
        
        with open(filename, 'wb') as self.fobj:
            self.write_array(self.write_string, self.data)

    def write_txt(self, filename, sep='\t'):
        """Write data to text file.
        
        Parameters
        ----------
        filename : str
            txt file name to write.
        """
        with open(filename, 'w') as fobj:
            fobj.write(f'{sep.join(self.metadata)}\n')
            
            for row in self.data:
                fobj.write(f'{sep.join(row)}\n')


class acpdb(LV_fd):
    """acp database.
    
    LNEC-SPA acquisition parameters file.
    """
    def __init__(self, endian='>', encoding='cp1252'):
        """Initialize instance.
        
        Parameters
        ----------
        endian : str, optional
            Defaults to '>'

        encoding : str, optional
            Defaults to 'cp1252'
        
        Notes
        -----
        Test for endianness:
            On Windows 10
                CLI ... returns ...
                Python console 'sys.byteorder' returns 'little'-
            On Linux
                CLI 'lscpu | grep Endian' returns 'Byte Order: Little Endian'
                Python console 'sys.byteorder' returns 'little'.
        Test for encoding:
            On Windows 10
                CLI 'chcp' returns 'cp1252'.
                Python console 'sys.stdin.encoding'.
            On Linux
                CLI 'locale | grep LANG' returns 'LANG=pt_PT.UTF-8'
                Python console 'sys.stdin.encoding' returns 'UTF-8'.
        """
        super().__init__(endian, encoding)

    def read(self, filename):
        """Read acp cluster from file.
        
        Parameters
        ----------
        filename : str
            acp file name to read.
        
        Returns
        -------
        self
        """
        self.filename = filename
        
        with open(filename, 'rb') as self.fobj:
            # Timing parameters
            self.TMGsource = self.read_string() # Timing source
            self.TMGsamples = self.read_numeric(self.LVint32) # Number of samples
            self.TMGsampling = self.read_numeric(self.LVfloat64) # Sampling frequency
            self.TMGbuffer = self.read_numeric(self.LVuint32) # Buffer size
            self.TMGscans = self.read_numeric(self.LVint32) # Samples per bin
            self.TMGreading = self.read_numeric(self.LVint32) # Reading relative to
            # Trigger parameters
            self.TRGsource = self.read_string() # Trigger source
            self.TRGtype = self.read_numeric(self.LVint32) # Start trig type
            self.TRGedge = self.read_numeric(self.LVint32) # Trigger edge/slope
            self.TRGlevel = self.read_numeric(self.LVfloat64) # Trigger level
            # Auto-processing parameters
            self.APenabled = self.read_boolean()
            self.APbins = self.read_numeric(self.LVint32)
            self.APfilterdesign = self.read_numeric(self.LVint32)
            self.APfilterorder = self.read_numeric(self.LVint32)
            self.APfiltertype = self.read_numeric(self.LVint32)
            self.APcutoff1 = self.read_numeric(self.LVfloat64)
            self.APcutoff2 = self.read_numeric(self.LVfloat64)
            self.APdecimate = self.read_numeric(self.LVuint16)
            # End of file
            self.eof = self.EOD()
        
        return self

    def write_txt(self, filename, sep='\t'):
        """Write data to text file.
        
        Parameters
        ----------
        filename : str
            txt file name to write.
            
        sep : str, optional
            Defaults to '\t'
        """
        with open(filename, 'w') as fobj:
            fobj.write('Timing parameters\n')
            fobj.write(f'Timing source: {self.TMGsource}\n')
            fobj.write(f'Number of samples: {self.TMGsamples}\n')
            fobj.write(f'Sampling frequency: {self.TMGsampling}\n')
            fobj.write(f'Buffer size: {self.TMGbuffer}\n')
            fobj.write(f'Samples per bin: {self.TMGscans}\n')
            fobj.write(f'Reading relative to: {self.TMGreading}\n')
            fobj.write('Trigger parameters\n')
            fobj.write(f'Trigger source: {self.TRGsource}\n')
            fobj.write(f'Start trigger type: {self.TRGtype}\n')
            fobj.write(f'Trigger edge/slope: {self.TRGedge}\n')
            fobj.write(f'Trigger level: {self.TRGlevel}\n')
            fobj.write('Auto-processing parameters\n')
            fobj.write(f'Enabled: {self.APenabled}\n')
            fobj.write(f'Apply to each bins: {self.APbins}\n')
            fobj.write(f'Filter design: {self.APfilterdesign}\n')
            fobj.write(f'Filter order: {self.APfilterorder}\n')
            fobj.write(f'Filter type: {self.APfiltertype}\n')
            fobj.write(f'Cut-off frequencies: {self.APcutoff1}, {self.APcutoff2}\n')
            fobj.write(f'Decimate factor: {self.APdecimate}\n')


def main():
    """Test suite."""
    pass


if __name__ == '__main__':    # run tests if called from command-line
    main()
