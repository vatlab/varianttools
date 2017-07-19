from scipy.sparse import csr_matrix,csc_matrix,rand,coo_matrix
import tables as tb
import numpy as np
import pandas as pd

import time
# from .utils import env
import math
import os
import numpy as np



def store_full_genotype_into_HDF5(file,store):
    start=time.time()
    storeHDF5=pd.HDFStore(store)
    count=0
    for chunk in pd.read_table(file,index_col=0,chunksize=1000,dtype=np.int8):
         storeHDF5.append("dataTable",chunk,complevel=9,min_itemsize={'index':100}, complib='blosc')
         print(storeHDF5.get_storer("dataTable").nrows)
    # storeHDF5["dataTable"]=pd.read_table(file,index_col=0)
    storeHDF5.close()
    print('load time: {0}'.format(time.time()-start))





def store_csr_arrays_into_earray_HDF5(data,indices,indptr,shape,ids,groupName,chr,store="check.h5"):
    rownames=np.array(ids)
    filters = tb.Filters(complevel=9, complib='blosc')
    with tb.open_file(store,'a') as f:
        node="/chr"+chr
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName
        if node in f:
            group=f.get_node(node)
        else:
            if "/chr"+chr not in f:
                group=f.create_group("/","chr"+chr,"chromosome")
            if len(groupName)>0:
                group=f.create_group("/chr"+chr,groupName)
        for par in ('data', 'indices', 'indptr', 'shape'):
            arr = None
            atom = None
            if (par=='data'):
                arr=np.array(data)
                atom=tb.Atom.from_dtype(np.dtype(np.float64))
            elif (par=='indices'):
                arr=np.array(indices)
                atom=tb.Atom.from_dtype(np.dtype(np.int32))
            elif (par=='indptr'):
                arr=np.array(indptr)
                atom=tb.Atom.from_dtype(np.dtype(np.int32))
            elif (par=='shape'):
                arr=np.array(shape)
                atom=tb.Atom.from_dtype(np.dtype(np.int32))
            # if (arr.shape[0]!=0):
            ds = f.create_earray(group, par, atom, (0,),filters=filters)
            ds.append(arr)
       
        ds = f.create_earray(group, "rownames",  tb.StringAtom(itemsize=200), (0,),filters=filters)
        ds.append(rownames)
        # ds = f.create_earray(f.root, "rownames",  tb.StringAtom(itemsize=200), (0,),filters=filters)
        # ds.append(rownames)
    f.close()




def append_csr_arrays_into_earray_HDF5(data,indices,indptr,shape,ids,groupName,chr,store="check.h5"):
    rownames=np.array(ids)

    with tb.open_file(store,'a') as f:
        node="/chr"+chr
        group=None
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName
        if node in f:
            group=f.get_node(node)
        for par in ('data', 'indices', 'indptr'):
            arr = None
            atom = None
            ds=None
            if (par=='data' and data is not None):
                arr=np.array(data,dtype=np.dtype(np.float64))
                ds=group.data
                ds.append(arr)
            elif (par=='indices' and indices is not None):
                arr=np.array(indices,dtype=np.dtype(np.int32))
                ds=group.indices
                ds.append(arr)
            elif (par=='indptr' and indptr is not None):
                arr=np.array(indptr,dtype=np.dtype(np.int32))
                ds=group.indptr
                ds.append(arr)
        # print(len(group.rownames[:]))
        group.rownames.append(rownames)
        group.shape[1]=shape[1]
        group.shape[0]=len(group.rownames[:])
    f.close()



def store_csr_genotype_into_earray_HDF5(data, name, ids,store='check.h5'):
    #data: a pandas data frame
    #name: the name of the tables
    if (data.__class__ != csr_matrix):
        m=csr_matrix(data.as_matrix())
        rownames=data.index.values
        colnames=data.columns.values
    else:
        m=data
        rownames=np.array(ids)
    msg = "This code only works for csr matrices"
    assert(m.__class__ == csr_matrix), msg
    filters = tb.Filters(complevel=9, complib='blosc')
    with tb.open_file(store,'a') as f:
        for par in ('data', 'indices', 'indptr', 'shape'):
            full_name = '%s_%s' % (name, par)
            try:
                n = getattr(f.root, full_name)
                n._f_remove()
            except AttributeError:
                pass
            arr = np.array(getattr(m, par))
            atom = tb.Atom.from_dtype(arr.dtype)
            if (par=='data'):
                atom=tb.Atom.from_dtype(np.dtype(np.int8))
            if (arr.shape[0]!=0):
                ds = f.create_earray(f.root, full_name, atom, (0,),filters=filters)
                ds.append(arr)
        ds = f.create_earray(f.root, name+"_rownames",  tb.StringAtom(itemsize=200), (0,),filters=filters)
        ds.append(rownames)
        # ds = f.create_carray(f.root, name+"_colnames", tb.StringAtom(itemsize=100), colnames.shape,filters=filters)
        # ds[:] = colnames
    f.close()


#passs in a data matrix, save into HDF5
def store_csr_genotype_into_one_HDF5(data, groupName, ids,index,chr,store='check.h5'):

    if (data.__class__ != csr_matrix):
        m=csr_matrix(data.as_matrix())
        rownames=data.index.values
        colnames=data.columns.values
    else:
        m=data
        rownames=np.array(ids)
    msg = "This code only works for csr matrices"
    assert(m.__class__ == csr_matrix), msg
    filters = tb.Filters(complevel=9, complib='blosc')
    with tb.open_file(store,'a') as f:
        node="/chr"+chr
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName
        if node in f:
            group=f.get_node(node)
        else:
            if "/chr"+chr not in f:
                group=f.create_group("/","chr"+chr,"chromosome")
            if len(groupName)>0:
                group=f.create_group("/chr"+chr,groupName)
        for par in ('data', 'indices', 'indptr', 'shape'):
            try:
                n = getattr(group, par)
                n._f_remove()
            except AttributeError:
                pass
            arr = np.array(getattr(m, par))
            atom = tb.Atom.from_dtype(arr.dtype)
            if (arr is not None):
                ds = f.create_carray(group, par, atom, arr.shape,filters=filters)
                ds[:] = arr
        ds = f.create_carray(group, "rownames",  tb.StringAtom(itemsize=200), rownames.shape,filters=filters)
        ds[:] = rownames
        # ds = f.create_carray(f.root, name+"_colnames", tb.StringAtom(itemsize=100), colnames.shape,filters=filters)
        # ds[:] = colnames
    f.close()




if __name__ == '__main__':
    directory="/Users/jma7/Development/VAT/chr22_50t/"
    file=directory+"output.txt"
    store_full_genotype_into_HDF5(file,store=directory+"full_chr22_test.h5")
    # store_fake_genotype_into_HDF5(file,store=directory+"fake_chr22_test.h5")

