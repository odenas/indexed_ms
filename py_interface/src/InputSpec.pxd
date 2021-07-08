from libcpp.string cimport string
from libcpp cimport bool



cdef extern from "fd_ms/input_spec.hpp" namespace "fdms":
    cdef cppclass InputSpec:
        InputSpec() except +
        InputSpec(InputSpec&) except +
        InputSpec(string, string) except +

        string s_fname, t_fname;
        string fwd_cst_fname, rev_cst_fname;
        string fwd_maxrep_fname, rev_maxrep_fname;
        string rev_elst_fname;
        string fwd_nwdlst_fname;
        string runs_fname, ms_fname;

        InputSpec& operator=(InputSpec)

        string load_s(bool reverse)

        @staticmethod
        void reverse_in_place(string& s)

        const size_t s_size()
        const size_t t_size()

        @staticmethod
        string rdix_fname(const string& ms_fname, const size_t block_size)

        @staticmethod
        const size_t ridx_block_size(const string& ms_fname, const string ridx_fname)

