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



def store_csc_genotype_into_one_HDF5(data, name, store='check.h5'):
    #data: a pandas data frame
    #name: the name of the tables
    m=csc_matrix(data.as_matrix())
    # rownames=data.index.values
    # colnames=data.columns.values
    msg = "This code only works for csc matrices"
    assert(m.__class__ == csc_matrix), msg
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
            ds = f.create_carray(f.root, full_name, atom, arr.shape,filters=filters)
            ds[:] = arr
        # ds = f.create_carray(f.root, name+"_rownames",  tb.StringAtom(itemsize=200), rownames.shape,filters=filters)
        # ds[:] = rownames
        # ds = f.create_carray(f.root, name+"_colnames", tb.StringAtom(itemsize=100), colnames.shape,filters=filters)
        # ds[:] = colnames
    f.close()


def store_csr_arrays_into_earray_HDF5(data,indices,indptr,shape,ids,name,chr,store="check.h5"):
    rownames=np.array(ids)
    filters = tb.Filters(complevel=9, complib='blosc')
    full_name=None
    with tb.open_file(store,'a') as f:
        group=f.create_group("/","chr"+chr,"chromosome")
        for par in ('data', 'indices', 'indptr', 'shape'):
            if len(name)>0:
                full_name = '%s_%s' % (name, par) 
            else:
                full_name = '%s' % (par)
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
            if (arr.shape[0]!=0):
                ds = f.create_earray(group, full_name, atom, (0,),filters=filters)
                ds.append(arr)
        if(len(name)>0): 
            ds = f.create_earray(group, name+"_rownames",  tb.StringAtom(itemsize=200), (0,),filters=filters)
        else:
            ds = f.create_earray(group, "rownames",  tb.StringAtom(itemsize=200), (0,),filters=filters)
        ds.append(rownames)
        # ds = f.create_earray(f.root, "rownames",  tb.StringAtom(itemsize=200), (0,),filters=filters)
        # ds.append(rownames)
    f.close()


def append_csr_arrays_into_earray_HDF5(data,indices,indptr,shape,ids,chr,store="check.h5"):
    rownames=np.array(ids)
    with tb.open_file(store,'a') as f:
        group=f.get_node("/chr"+chr)

        for par in ('data', 'indices', 'indptr'):
            full_name = '%s' % (par)
            arr = None
            atom = None
            ds=None
            if (par=='data'):
                arr=np.array(data,dtype=np.dtype(np.float64))
                ds=group.data
            elif (par=='indices'):
                arr=np.array(indices,dtype=np.dtype(np.int32))
                ds=group.indices
            elif (par=='indptr'):
                arr=np.array(indptr,dtype=np.dtype(np.int32))
                ds=group.indptr
            if (arr.shape[0]!=0):
                ds.append(arr)
        group.shape[0]=shape[0]
        group.rownames.append(rownames)
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


def store_csr_genotype_into_one_HDF5(data, name, ids,index,chr,store='check.h5'):
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
        node="/chr"+chr
        if node in f:
            group=f.get_node(node)
        else:
            group=f.create_group("/","chr"+chr,"chromosome")
        for par in ('data', 'indices', 'indptr', 'shape'):
            full_name = '%s_%s' % (name, par)
            try:
                n = getattr(group, full_name)
                n._f_remove()
            except AttributeError:
                pass
            arr = np.array(getattr(m, par))
            atom = tb.Atom.from_dtype(arr.dtype)

            if (arr.shape[0]!=0):
                ds = f.create_carray(group, full_name, atom, arr.shape,filters=filters)
                ds[:] = arr
        ds = f.create_carray(group, name+"_rownames",  tb.StringAtom(itemsize=200), rownames.shape,filters=filters)
        ds[:] = rownames
        # ds = f.create_carray(f.root, name+"_colnames", tb.StringAtom(itemsize=100), colnames.shape,filters=filters)
        # ds[:] = colnames
    f.close()


def store_genotype_into_nested_HDF5(data, name, store='check'):
    colnames=data.columns.values
    for index,colname in enumerate(colnames):
        print(colname)
        filename=str(index+1)
        # store_genotype_into_one_HDF5(data[[colname]],name,store+filename+".h5")
        store_genotype_into_one_HDF5(data[[colname]],name+filename,store+"multiple.h5")




if __name__ == '__main__':
    directory="/Users/jma7/Development/VAT/chr22_50t/"
    #store in spase 
    #file=pd.read_table(directory+"output.txt",index_col=0,dtype=np.int8)
    #store_csr_genotype_into_one_HDF5(file,"dataTable",store=directory+"chr_50t_genotype.h5")
    # store_csc_genotype_into_one_HDF5(file,"dataTable",store=directory+"csc_chr22_test.h5")
    #store in full
    file=directory+"output.txt"
    store_full_genotype_into_HDF5(file,store=directory+"full_chr22_test.h5")
    # store_fake_genotype_into_HDF5(file,store=directory+"fake_chr22_test.h5")

