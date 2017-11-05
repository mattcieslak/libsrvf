# distutils: language=c++ 
import cython 

import numpy as np
cimport numpy as np
np.import_array()

from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.deque cimport deque
from libcpp cimport bool

cdef extern from "include/srvf/matrix.h" namespace "srvf":
    cdef cppclass Matrix:
        enum PackingMethod:
            POINT_PER_ROW=0
            POINT_PER_COLUMN=1
        Matrix()
        Matrix(size_t, size_t)
        Matrix(size_t, size_t, double *data)
        size_t rows()
        size_t cols()
        vector[double] data()
        double norm()


cdef extern from "include/srvf/point.h" namespace "srvf":
    cdef cppclass Point:
        Point() except +
        Point(size_t) except +
        Point(size_t, double) except +
        size_t dim()
        double norm()
        double operator[](size_t)
        Point operator+(Point)
        Point operator-(Point)
        
cdef class ndPoint:
    cdef Point _thisptr
    def __cinit__(self, size_t ndims=1, double v=0.0):
        self._thisptr = Point(ndims, v)
    
    cpdef double norm(self):
        return self._thisptr.norm()
        
    cpdef size_t dim(self):
        return self._thisptr.dim()
        
            

cdef extern from "include/srvf/pointset.h" namespace "srvf":
        
    cdef cppclass Pointset:
        enum PackingMethod:
            POINT_PER_ROW=0
            POINT_PER_COLUMN=1
        Pointset() except +
        Pointset(size_t, size_t, double *data) except +
        Pointset(size_t, size_t, double *data, PackingMethod) except +
        Pointset(size_t, size_t, vector[double]) except +
        Pointset(size_t, size_t, vector[double], PackingMethod) except +
        Pointset(size_t, size_t) except +
        Pointset(size_t, size_t, double v) except +
        size_t dim()
        size_t npts()
        Point operator[](size_t)
        double distance(size_t, size_t)
        Point centroid()
        vector[double] to_vector()
        vector[double] to_vector(PackingMethod)
        void translate(size_t, Point)
        void scale(double)
        
    cdef double distance(Pointset, Pointset, int, int)
    cdef double dot_product(Pointset, Pointset, size_t, size_t)
    
cdef class ndPointset:
    cdef Pointset _thisptr
    def __cinit__(self, np.ndarray[np.float_t, ndim=1] data ):
        cdef size_t npts = data.shape[0]
        cdef size_t ndims = 1
        self._thisptr = Pointset(ndims, npts, &data[0])

    def get_data(self):
        """ Returns a copy of the internal data"""
        cdef vector[double] vec_data = self._thisptr.to_vector()
        cdef size_t num_pts = self._thisptr.npts()        
        cdef np.ndarray[double, ndim=1] data = np.empty(num_pts)
        cdef size_t i
        for i in range(num_pts):
            data[i] = vec_data[i]
        return data
        
    def scale(self, double scale_factor):
        self._thisptr.scale(scale_factor)
    def dim(self):
        return self._thisptr.dim()
    def npts(self):
        return self._thisptr.npts()
        
    
        
cdef extern from "include/srvf/util.h" namespace "srvf::util":
    cdef vector[double] linspace(double, double, size_t)    
    cdef vector[double] random_vector(size_t)
        

#cdef class ndMatrix:
    #cdef Matrix _thisptr
    #def __cinit__(self, np.ndarray[np.float_t, ndim=2] data ):
    #    cdef size_t npts = data.shape[0]
    #    cdef size_t ndims = 1
    #    self._thisptr = Pointset(ndims, npts, &data[0])

    #def get_data(self):
    #    """ Returns a copy of the internal data"""
    #    cdef vector[double] vec_data = self._thisptr.to_vector()
    #    cdef size_t num_pts = self._thisptr.npts()        
    #    cdef np.ndarray[double, ndim=1] data = np.empty(num_pts)
    #    cdef size_t i
    #    for i in range(num_pts):
    #        data[i] = vec_data[i]
    #    return data
        


cdef extern from "include/srvf/plf.h" namespace "srvf":
    cdef cppclass Plf:
        Plf()
        Plf(Pointset)
        double domain_ub()
        double domain_lb()
        double arc_length()
        Point centroid()
        void scale_to_unit_arc_length()
        Pointset samps()
        vector[double] params()
        size_t ncp()
        
cdef class ndPlf:
    cdef Plf _thisptr
    def __cinit__(self, np.ndarray[double, ndim=1] data):
        cdef size_t npts = data.shape[0]
        cdef size_t ndim = 1
        cdef Pointset temp_pointset = Pointset(ndim, npts, &data[0])
        self._thisptr = Plf(temp_pointset)
        
cdef extern from "include/srvf/srvf.h" namespace "srvf":
    cdef cppclass Srvf:
        Srvf()
        Srvf(Pointset)
        size_t dim()
        size_t ncp()
        double domain_lb()
        double domain_ub()
        void scale(double)
        void scale_to_unit_norm()
        Pointset samps()
        vector[double] params()
    

#cdef extern from "include/srvf/opencurves.h" namespace "srvf::opencurves":
#    cdef Srvf karcher_mean(vector[Srvf], bool, bool, size_t, size_t, size_t)
#    cdef Srvf shooting_vector(Srvf, Srvf)


cdef extern from "include/srvf/functions.h" namespace "srvf::functions":
    ctypedef pair[size_t,size_t] match_vertex_t
    cdef Srvf karcher_mean(vector[Srvf] Qs, double tol, size_t max_iters)
    cdef vector[Plf] groupwise_build_gammas(Srvf, vector[Srvf])
    cdef vector[Plf] groupwise_optimal_reparam(Srvf, vector[Srvf])
    cdef vector[Plf] build_gammas(Srvf, Srvf, deque[match_vertex_t])


    
cdef extern from "include/srvf/qmap.h" namespace "srvf":
    cdef Srvf plf_to_srvf(Plf)
    cdef Plf srvf_to_plf(Srvf)
    
# Utility functions    
cdef np.ndarray[double, ndim=1] get_srvf_samples(Srvf input_srvf):
    cdef size_t num_points = input_srvf.ncp()
    cdef vector[double] vec_data = input_srvf.samps().to_vector()
    cdef np.ndarray[double, ndim=1] output = np.empty(num_points)
    cdef size_t index
    for index in range(num_points):
        output[index] = vec_data[index]
    return output
    
cdef np.ndarray[double, ndim=1] get_plf_samples(Srvf input_srvf):
    cdef size_t num_points = input_srvf.ncp()
    cdef vector[double] vec_data = input_srvf.samps().to_vector()
    cdef np.ndarray[double, ndim=1] output = np.empty(num_points)
    cdef size_t index
    for index in range(num_points):
        output[index] = vec_data[index]
    return output




cdef np.ndarray[double, ndim=1] get_srvf_times(Srvf input_srvf):
    cdef size_t num_points = input_srvf.ncp()
    cdef np.ndarray[double, ndim=1] output = np.empty(num_points)
    cdef size_t index
    cdef vector[double] time_vec = input_srvf.params()
    for index in range(num_points):
        output[index] = time_vec[index]
    return output
    
    
cdef np.ndarray[double, ndim=1] get_plf_times(Plf input_plf):
    cdef size_t num_points = input_plf.ncp()
    cdef np.ndarray[double, ndim=1] output = np.empty(num_points)
    cdef size_t index
    cdef vector[double] time_vec = input_plf.params()
    for index in range(num_points):
        output[index] = time_vec[index]
    return output


# Python functions !    

cpdef np.ndarray[double, ndim=2] functions_to_srvfs(np.ndarray[double, ndim=2] orig_functions):
    """ Returns the SRVFs of the input functions"""
    cdef size_t num_functions = orig_functions.shape[0]
    cdef size_t num_samples = orig_functions.shape[1]
    cdef size_t function_num = 0
    cdef size_t ndims = 1
    cdef size_t sample_num
    cdef np.ndarray[double, ndim=2] srvf_mat = np.empty((num_functions,num_samples))
    cdef np.ndarray[double, ndim=1] temp_data
    cdef Pointset orig_pointset, q_pointset
    cdef Plf temp_plf
    cdef Srvf temp_srvf
    cdef vector[double] q_vector = q_pointset.to_vector()
    
    for function_num in range(num_functions):
        temp_data = orig_functions[function_num]
        orig_pointset = Pointset(ndims, num_samples, &temp_data[0]) 
        
        # Make a Plf
        temp_plf = Plf(orig_pointset)
        temp_srvf = plf_to_srvf(temp_plf)
        q_pointset = temp_srvf.samps()
        q_vector = q_pointset.to_vector()
        for sample_num in range(temp_srvf.ncp()):
            srvf_mat[function_num,sample_num] =  q_vector[sample_num]
            
    return srvf_mat
        
        
        


def calculate_karcher_mean(np.ndarray[double, ndim=2] orig_functions,
                double tol=1e-3, size_t max_iters = 0, 
                bool do_rots=True, bool do_reparams=True):
    """
    Calculates the Karcher Mean of a set of input functions.
    
    Parameters:
    ===========
    
    orig_functions:np.ndarray
      Each row is a pointwise-linear function, each column 
    """
    cdef size_t num_functions = orig_functions.shape[0]
    cdef size_t num_samples = orig_functions.shape[1]
    cdef size_t function_num = 0
    cdef size_t ndims = 1
    cdef size_t sample_num
    cdef np.ndarray[double, ndim=1] temp_data
    cdef Pointset orig_pointset, q_pointset
    cdef Plf temp_plf
    
    cdef Srvf q_mean
    cdef vector[Srvf] Qs
    
    # Create the SRVFs 
    for function_num in range(num_functions):
        temp_data = orig_functions[function_num]
        orig_pointset = Pointset(ndims, num_samples, &temp_data[0]) 
        temp_plf = Plf(orig_pointset)    
        Qs.push_back( plf_to_srvf(temp_plf))
        Qs[function_num].scale_to_unit_norm() # Required
        #Qs[function_num] = constant_speed_param(Qs[function_num]) # Required
    
    #  actually do it
    q_mean = karcher_mean(Qs, tol, max_iters)
    return get_srvf_samples(q_mean), get_srvf_times(q_mean)
    