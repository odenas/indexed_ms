//
//  input_spec.hpp
//  fast_ms
//
//  Created by denas on 10/21/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef input_spec_h
#define input_spec_h

#include "query.hpp"

using namespace std;

namespace fdms {

    class InputSpec {
    public:
        string s_fname, t_fname;
        string fwd_cst_fname, rev_cst_fname;
        string fwd_maxrep_fname, rev_maxrep_fname;
        string rev_elst_fname;
        string fwd_nwdlst_fname;
        string runs_fname, ms_fname;

        InputSpec() {
        }

        InputSpec(const InputSpec &is) :
        s_fname{is.s_fname}, t_fname{is.t_fname},
        fwd_cst_fname{is.fwd_cst_fname},
        rev_cst_fname{is.rev_cst_fname},
        rev_maxrep_fname{is.rev_maxrep_fname},
        rev_elst_fname{is.rev_elst_fname},
        fwd_nwdlst_fname{is.fwd_nwdlst_fname},
        runs_fname{is.runs_fname}, ms_fname{is.ms_fname}
        {
        }

        InputSpec(string s_fn, string t_fn) : s_fname(s_fn), t_fname(t_fn) {
            fwd_cst_fname = s_fname + ".fwd.stree";
            rev_cst_fname = s_fname + ".rev.stree";
            fwd_maxrep_fname = s_fname + ".fwd.maxrep";
            rev_maxrep_fname = s_fname + ".rev.maxrep";
            rev_elst_fname = s_fname + ".rev.elst";
            fwd_nwdlst_fname = s_fname + ".fwd.nwdlst";
            runs_fname = t_fname + ".runs";
            ms_fname = t_fname + ".ms";
        }

        InputSpec& operator=(const InputSpec& other) {
            if (this != &other) {
                s_fname = string(other.s_fname);
                t_fname = string(other.t_fname);
                fwd_cst_fname = string(other.fwd_cst_fname);
                rev_cst_fname = string(other.rev_cst_fname);
                fwd_maxrep_fname = string(other.fwd_maxrep_fname);
                rev_maxrep_fname = string(other.rev_maxrep_fname);
                rev_elst_fname = string(other.rev_elst_fname);
                fwd_nwdlst_fname = string(other.fwd_nwdlst_fname);
                runs_fname = string(other.runs_fname);
                ms_fname = string(other.ms_fname);
            }
            return *this;
        }

        string load_s(bool reverse = false) const {
            string s;
            char ch{};
            std::ifstream s_file{s_fname, ios::binary | ios::in};
            while (s_file.get(ch))
                s += ch;
            if (reverse)
                reverse_in_place(s);
            return s;
        }

        static void reverse_in_place(string& s) {
            unsigned long long n = s.size();

            for (int i = 0; i < n / 2; i++) {
                char c = s[i];
                s[i] = s[n - 1 - i];
                s[n - 1 - i] = c;
            }
        }

        const size_t s_size() const {
            return (size_t) Query::query_length(s_fname);
        }

        const size_t t_size() const {
            return (size_t) Query::query_length(t_fname);
        }

        static string rdix_fname(const string ms_fname, const size_t block_size) {
            return ms_fname + "." + std::to_string(block_size) + ".ridx";
        }

        static const size_t ridx_block_size(const string ms_fname, const string ridx_fname) {
            size_t suff_len = 5; // '.ridx'.size()
            size_t l = ms_fname.size();

            string s = ridx_fname.substr(l + 1,
                    ridx_fname.size() - (l + 1) - suff_len);
            return static_cast<size_t> (std::stoi(s));
        }
    };
};

#endif /* input_spec_h */
