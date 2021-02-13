import netCDF4 as nc4
import numpy as np
from dateutil.parser import parse
import pandas as pd

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def read_csv(fname):
    with open(fname) as f:
        # first row indicates size
        # first col indicates time
        data_raw = np.loadtxt(fname, delimiter=',')

        size = data_raw[0,1:]
        time_hours = data_raw[1:,0]

        time_min = np.zeros_like(time_hours)
        time_min[0] = time_hours[0]*60

        iday = 0
        for i in range(1,len(time_hours)):
            if time_hours[i-1] > time_hours[i]:
                iday += 1
            time_min[i] = time_hours[i]*60 + iday*1440

        v_edges = np.empty(shape=len(size)+1)
        v_edges[:-1] = size
        v_edges[-1] = v_edges[-2] + (v_edges[-2] - v_edges[-3])
        print(v_edges[-10:])
        data = {'t_min':time_min, 'data':data_raw[1:,1:], 'v_edges':v_edges}

    return data

data = read_csv('size_distribution/zinser_psd.csv')
# new: set initial time to zero
data['t_min'] -= data['t_min'][0]

data_fig3 = pd.read_csv('ground_truth/zinser_figure3.csv')
data['par'] = (data_fig3['Estar'].values*data_fig3['Ek2'].values).astype(float)

zinser_data = pd.read_csv('ground_truth/zinser_figure2a.csv')
zinser_abundances = (0.5 * (zinser_data['cells A'].values + zinser_data['cells B'].values)).astype(int) # mean of both columns

create_plots = True


#
# specify v
#

v_min = 16
m = 27 
delta_v_inv = 8 

import sys
if len(sys.argv) > 1:
    m = int(sys.argv[1])
    if len(sys.argv) > 2:
        delta_v_inv = int(sys.argv[2])
        if len(sys.argv) > 3:
            v_min = float(sys.argv[3])

delta_v = 1.0/delta_v_inv
v = v_min * 2**(np.arange(m+1)*delta_v) # to get m intervals, we need m+1 edge

print(data['data'].shape)

# data['data'] is num_t x num_v

def regrid(v0, w0, v1, permit_smallgrid=False):
    # assert that old grid is contained in new one
    #assert v1[0] <= v0[0]
    if not permit_smallgrid:
        assert v1[-1] >= v0[-1]
    # check sizes
    if v0.size != w0.shape[0]+1:
        raise ValueError('Size of v0 ({}) does not match shape of w0 {}.'.format(v0.size, w0.shape))

    w1 = np.zeros((v1.size-1,w0.shape[1]))
    #v_data = data['v_edges']

    i1 = 1 # first right edge
    for i0 in range(1,v0.size):
        if v0[i0] < v1[i1]:
            if v0[i0-1] >= v1[i1-1]:
                w1[i1-1,:] += w0[i0-1,:]
            else:
                pass
        elif i1+1 == v1.size:
            # can only happen for permit_smallgrid
            a = (v1[i1]-v0[i0-1])/(v0[i0]-v0[i0-1])
            # left side
            w1[i1-1,:] += w0[i0-1,:]*a
            # right side is not covered by v1
        elif v0[i0] < v1[i1+1]:
            a = (v1[i1]-v0[i0-1])/(v0[i0]-v0[i0-1])
            # left side
            w1[i1-1,:] += w0[i0-1,:]*a
            # right side
            w1[i1,:] += w0[i0-1,:]*(1-a)
            i1 += 1
        else:
            raise NotImplementedError('Case not yet covered')

    if v1[-1] >= v0[-1] and v1[0] <= v0[0]:
        assert np.all(np.abs(np.sum(w0,axis=0)-np.sum(w1,axis=0))<1e-9)
    else:
        print('maximum loss of data due to reduction in grid coverage: {}'.format(np.max(np.abs(np.sum(w0,axis=0)-np.sum(w1,axis=0)))))
    return w1

w = regrid(data['v_edges'], data['data'].T, v, permit_smallgrid=True)

#
# write to netCDF
#

fname = 'size_distribution/zinser_test.nc'.format(m, delta_v_inv)
with nc4.Dataset(fname,'w') as nc:
    nc.createDimension('time', w.shape[1])
    nc.createDimension('size', w.shape[0])
    nc.createDimension('size_bounds', w.shape[0]+1)

    nc.createVariable('time', int, ('time',), fill_value=False)
    nc.variables['time'][:] = data['t_min']
    nc.variables['time'].units = 'minutes since start of experiment'
    nc.variables['time'].size == len(zinser_abundances)

    nc.createVariable('size_bounds', float, ('size_bounds',), fill_value=False)
    nc.variables['size_bounds'][:] = v
    nc.variables['size_bounds'].units = 'um3'

    nc.createVariable('w_obs', float, ('size','time'), fill_value=False)
    nc.variables['w_obs'][:] = w/np.sum(w,axis=0)

    nc.createVariable('abundance', int, ('time',), fill_value=False) 
    nc.variables['abundance'][:] = zinser_abundances
    nc.variables['abundance'].units = 'cells ml-1'
    nc.variables['abundance'].long_name = 'cell abundance'
       
    nc.createVariable('count', int, ('time',), fill_value=False)
    nc.variables['count'][:] = (np.round(zinser_abundances * 0.0002)).astype(int) # 0.02 mL of 100X diluted samples was analyzed
    nc.variables['count'].units = 'cells'
    nc.variables['count'].long_name = 'cell count'

    nc.createVariable('PAR', float, ('time',), fill_value=False)
    nc.variables['PAR'][:] = data['par']
    nc.variables['PAR'].units = 'umol photons/m2/s'

    nc.createVariable('m', int, fill_value=False)
    nc.variables['m'][:] = m

    nc.createVariable('delta_v_inv', int, fill_value=False)
    nc.variables['delta_v_inv'][:] = delta_v_inv

    nc.createVariable('v_min', float, fill_value=False)
    nc.variables['v_min'][:] = v_min