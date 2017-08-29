//
//  Header.h
//  fast_ms
//
//  Created by denas on 8/28/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef Header_h
#define Header_h

#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>
#include <string.h>


using namespace std;

namespace fdms {
    
    
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
    

    class InputSpec{
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
        size_t runs_progress, ms_progress;
        size_t nthreads;
        
        InputFlags(bool lazy_wl, bool use_rank_fail, bool use_maxrep, bool lca_parents,
                   bool space, bool time_,
                   bool ans, bool v,
                   size_t runs_prgs, size_t ms_prgs,
                   bool load_stree, bool load_maxrep,
                   size_t nthreads) :
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
        runs_progress{static_cast<size_t>(std::stoi(input.getCmdOption("-runs_progress")))},
        ms_progress{static_cast<size_t>(std::stoi(input.getCmdOption("-ms_progress")))},
        nthreads{static_cast<size_t>(std::stoi(input.getCmdOption("-nthreads")))}
        {
            nthreads = (nthreads <= 0 ? 1 : nthreads);
            check();
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
        
        string ms_strategy_string(const size_t indent_level, const bool ellipsis) const {
            assert (indent_level < 10);
            string indent = {"**********", indent_level};
            
            return (" " + indent +
                    " computing with a " + (lazy ? "lazy" : "non-lazy") + ", " +
                    (rank_fail ? "double_rank_and_fail" : "double_rank_no_fail") + " / " +
                    (lca_parents ? "lca_parent" : "consecutive_parents") + " strategy, " +
                    ("over " + std::to_string(nthreads) + " thread" + (nthreads > 1 ? "s" : "")) +
                    (use_maxrep ? "Using maxrep" : "") +
                    (ellipsis ? " ... " : ".")
                    );
        }
        
        string runs_strategy_string(const size_t indent_level, const bool ellipsis) const {
            assert (indent_level < 10);
            string indent = {"**********", indent_level};
            
            return (" " + indent +
                    " computing with a, " + (rank_fail ? "double_rank_and_fail" : "double_rank_no_fail") + " / " +
                    (lca_parents ? "lca_parent" : "consecutive_parents") + " strategy, " +
                    ("over " + std::to_string(nthreads) + " thread" + (nthreads > 1 ? "s" : "")) +
                    (ellipsis ? " ... " : "."));
        }
        
    };
}


#endif /* Header_h */
