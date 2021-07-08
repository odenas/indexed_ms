from libcpp.string cimport string
from libcpp cimport bool
from InputSpec cimport InputSpec


cdef extern from "fd_ms/stree_sct3.hpp" namespace "fdms":
    cdef cppclass StreeOhleb[t_csa=*, t_lsp=*, t_bp_support=*, t_bv=*, t_rank=*, t_sel=*]:

        StreeOhleb() except +

        @staticmethod
        void load_or_build(StreeOhleb[]&, InputSpec&, const bool, const bool)

        const size_t size()

