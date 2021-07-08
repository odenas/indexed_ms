# distutils: language = c++

from InputSpec cimport InputSpec
from StreeOhleb cimport StreeOhleb


cdef class PyInputSpec:
    cdef InputSpec c_ispec  # the C++ instance we're wrapping

    def __cinit__(self, spath, tpath):
        self.c_ispec = InputSpec(spath.encode(), spath.encode())

cdef class Cst:
    cdef StreeOhleb c_cst

    def __cinit__(self, ispec: PyInputSpec):
        StreeOhleb.load_or_build(self.c_cst, ispec.c_ispec, False, False)


def slow(s, t):
    for l in range(len(t), 0, -1):
        if s.find(t[0: l]) >= 0:
            return l
    return 0

def fast(s, t):
    ispec = PyInputSpec(s.encode(), t.encode())
    return ispec

