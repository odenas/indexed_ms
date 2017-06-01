//
//  utils.h
//  fast_ms
//
//  Created by denas on 4/3/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include <iostream>
#include <string>
#include <vector>

#include <sdsl/suffix_trees.hpp>
#include <sdsl/bit_vectors.hpp>

#include "stree_sct3.hpp"
#include "maxrep_construction.hpp"

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))


using namespace std;


namespace fdms{
    typedef uint8_t char_type;
    typedef unsigned long long size_type;
    typedef std::pair<size_type, size_type> Interval;
    typedef map<std::string, size_type> Counter;


    using timer = std::chrono::high_resolution_clock;

    class InputSpec{
    private:
        sdsl::bit_vector parse_bitstr(string& s){
            sdsl::bit_vector b(s.size());

            for(size_type i = 0; i < s.size(); i++)
                b[i] = ((unsigned char)s[i] - 48);
            return b;
        }

    public:
        string s_fname, fwd_cst_fname, rev_cst_fname, rev_maxrep_fname;

        InputSpec(string s_fn) : s_fname(s_fn){
            fwd_cst_fname = s_fname + ".fwd.stree";
            rev_cst_fname = s_fname + ".rev.stree";
            rev_maxrep_fname = s_fname + ".rev.maxrep";
        }

        string load_s() const {
            string s;
            std::ifstream s_file {s_fname};
            while(s_file >> s)
                ;
            return s;
        }

        sdsl::bit_vector bps(string &s){
            sdsl::cst_sada<> st_of_s;
            sdsl::construct_im(st_of_s, s, 1);
            return st_of_s.bp;
        }

    };


    // copied from http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
    class OptParser{
    public:
        std::string empty = "0";
        OptParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            return empty;
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
            != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
    };


    class InputFlags{
    private:
        void check(){
            if (use_maxrep && lazy){
                cerr << "lazy and use_maxrep cannot be active at the same time" << endl;
                exit(1);
            }
            if(nthreads == 0){
                cerr << "nr. of threads (parallelism) should be a positive number (got " << nthreads << ")" << endl;
                exit(1);
            }
        }

    public:
        bool lazy, rank_fail, use_maxrep, lca_parents;
        bool space_usage, time_usage;
        bool answer;
        bool verbose;
        bool load_stree,load_maxrep;
        size_type runs_progress, ms_progress;
        size_type nthreads;

        InputFlags(bool lazy_wl, bool use_rank_fail, bool use_maxrep, bool lca_parents,
                   bool space, bool time_,
                   bool ans, bool v,
                   size_type runs_prgs, size_type ms_prgs,
                   bool load_stree, bool load_maxrep,
                   size_type nthreads) :
        lazy{lazy_wl}, rank_fail{use_rank_fail}, use_maxrep{use_maxrep}, lca_parents{lca_parents},
        space_usage {space},
        time_usage {time_},
        answer {ans},
        verbose{v},
        load_stree{load_stree}, load_maxrep{load_maxrep},
        runs_progress{runs_prgs}, ms_progress{ms_prgs},
        nthreads{nthreads}
        {
            check();
        }

        InputFlags (OptParser input) :
        lazy {input.getCmdOption("-lazy_wl") == "1"},             // lazy winer links
        rank_fail {input.getCmdOption("-rank_fail") == "1"},      // use the rank-and-fail strategy
        use_maxrep {input.getCmdOption("-use_maxrep") == "1"},    // use the maxrep vector
        lca_parents {input.getCmdOption("-lca_parents") == "1"},  // use lca insted of conscutive parent calls
        space_usage {input.getCmdOption("-space_usage") == "1"},  // space usage
        time_usage {input.getCmdOption("-time_usage") == "1"},    // time usage
        answer {input.getCmdOption("-answer") == "1"},            // answer
        verbose{input.getCmdOption("-verbose") == "1"},           // verbose
        load_stree{input.getCmdOption("-load_cst") == "1"},       // load CST of S and S'
        load_maxrep{input.getCmdOption("-load_maxrep") == "1"},   // load MAXREP of S'
        runs_progress{static_cast<size_type>(std::stoi(input.getCmdOption("-runs_progress")))},
        ms_progress{static_cast<size_type>(std::stoi(input.getCmdOption("-ms_progress")))},
        nthreads{static_cast<size_type>(std::stoi(input.getCmdOption("-nthreads")))}
        {
            nthreads = (nthreads <= 0 ? 1 : nthreads);
            check();
        }

        parent_seq_method get_parent_seq_method() const {
            return (lca_parents ? &StreeOhleb<>::maxrep_ancestor : &StreeOhleb<>::parent_sequence);
        }

        double_rank_method get_rank_method() const {
            return (rank_fail ?
                    &sdsl::bwt_of_csa_wt<sdsl::csa_wt<>>::double_rank_and_fail :
                    &sdsl::bwt_of_csa_wt<sdsl::csa_wt<>>::double_rank);
        }

        wl_method_t1 get_wl_method() const {
            if(lazy)
                return (rank_fail ? &StreeOhleb<>::lazy_double_rank_fail_wl : &StreeOhleb<>::lazy_double_rank_wl);
            return (rank_fail ? &StreeOhleb<>::double_rank_fail_wl : &StreeOhleb<>::double_rank_nofail_wl);
        }



        void show() const {
            cerr << "**********" << endl;
            cerr << "[wl]     lazy: " << lazy << ". rank_fail: " << rank_fail << "." << endl;
            cerr << "[parent] lca_parents: " << lca_parents << "." << endl;
            cerr << "[maxrep] use: " << use_maxrep << ". load: " << load_maxrep << "." << endl;
            cerr << "[cst]    load: " << load_stree << "." << endl;
            cerr << "[report] space_usage: " << space_usage << ". time_usage: " << time_usage << ". answer: " << answer << "." << endl;
            cerr << "nthreads: " << nthreads << endl;
            cerr << "[progress] ms: " << ms_progress << ". runs: " << runs_progress << "." << endl;
            cerr << "**********" << endl;
        }

        string ms_strategy_string(const size_type indent_level, const bool ellipsis) const {
            assert (indent_level < 10);
            string indent = {"**********", indent_level};

            return (" " + indent +
                    " computing with a " + (lazy ? "lazy" : "non-lazy") + ", " +
                    (rank_fail ? "double_rank_and_fail" : "double_rank_no_fail") + " strategy, " +
                    ("over " + std::to_string(nthreads) + " thread" + (nthreads > 1 ? "s" : "")) +
                    (use_maxrep ? "Using maxrep" : "") +
                    (ellipsis ? " ... " : ".")
                    );
        }

        string runs_strategy_string(const size_type indent_level, const bool ellipsis) const {
            assert (indent_level < 10);
            string indent = {"**********", indent_level};

            return (" " + indent +
                    " computing with a, " + (rank_fail ? "double_rank_and_fail" : "double_rank_no_fail") + " / " +
                    (lca_parents ? "lca_parent" : "consecutive_parents") + " strategy, " +
                    ("over " + std::to_string(nthreads) + " thread" + (nthreads > 1 ? "s" : "")) +
                    (ellipsis ? " ... " : "."));
        }

    };

    void dump_ms(sdsl::bit_vector& ms){
        size_type k = 0;
        for (size_type i = 0; i < ms.size(); i++){
            if(ms[i] == 1){
                cout << i - (2*k) << " ";
                k += 1;
            }
        }
        //cout << endl;
    }

    void reverse_in_place(string& s){
        size_type n = s.size();

        for(int i = 0; i < n / 2; i++){
            char c = s[i];
            s[i] = s[n - 1 - i];
            s[n - 1 - i] = c;
        }
    }


    std::vector<pair<size_type, size_type>> slice_input(const size_type input_size, const size_type nthreads){
        size_type chunk = input_size / nthreads;
        size_type extra = input_size % nthreads;
        size_type step = 0;

        std::vector<pair<size_type, size_type>> slices (nthreads);
        for(size_type i=0, from = 0; i<nthreads; i++){
            step = chunk + (i < extra ? 1 : 0);
            slices[i] = std::make_pair(from, from + step);
            from += step;
        }
        return slices;
    }

    template<class tree_type>
    size_type load_st(tree_type &st, const string &s, const string potential_stree_fname, const bool load){
        auto start = timer::now();
        if(load){
            cerr << " * loading the CST from " << potential_stree_fname << " ";
            sdsl::load_from_file(st, potential_stree_fname);
        } else {
            cerr << " * building the CST of length " << s.size() << " ";
            sdsl::construct_im(st, s, 1);
        }
        auto stop = timer::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }

    size_type load_maxrep(sdsl::bit_vector& vec, const StreeOhleb<>& st, const string& s, const string& load_fname, const bool load){
        auto start = timer::now();
        if(load){
            cerr << " * loading MAXREP from " << load_fname << " ";
            sdsl::load_from_file(vec, load_fname);
        } else {
            cerr << " * building MAXREP of length " << s.size() << " ";
            vec.resize(s.size() + 1);
            sdsl::util::set_to_value(vec, 0);
            build_maxrep_ohleb<StreeOhleb<>, sdsl::bit_vector, StreeOhleb<>::node_type>(st, vec);
        }
        auto stop = timer::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }


    /*
    size_type load_ohleb_st(StreeOhleb<> &st, const string &s, const string potential_stree_fname, const bool load){
        auto start = timer::now();
        if(load){
            cerr << " * loading the CST T(s') from " << potential_stree_fname << " ";
            sdsl::load_from_file(st, potential_stree_fname);
        } else {
            cerr << " * building the CST T(s') of length " << s.size() << " ";
            sdsl::construct_im(st, s, 1);
        }
        auto stop = timer::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    }
    */

    void report_progress(timer::time_point start_time, size_type curr_idx, size_type total){
        timer::time_point stop_time = timer::now();
        size_type elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() + 1;
        cerr << endl << "[" << elapsed / 1000 << " s] " << 100.0 * curr_idx / total << "% @ " << (1.0 * curr_idx / elapsed) << " KHz";
    }
    

}

#endif /* utils_h */
