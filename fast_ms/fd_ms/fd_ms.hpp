//
//  fd_ms.hpp
//  fast_ms
//
//  Created by denas on 11/2/16.
//  Copyright © 2016 denas. All rights reserved.
//

#ifndef fd_ms_h
#define fd_ms_h

#include <iostream>
#include <string>
#include <vector>

#include <mach/mach.h>
#include <mach/task.h>

#include <sdsl/select_support.hpp>
#include "basic.hpp"
#include "fd_ms_algorithms.hpp"

using namespace std;


namespace fdms{
    namespace monitor{
        typedef std::map<std::string, std::string> str_dict;
        typedef std::map<std::string, size_type> size_dict;

        const char K_BWT[]      = "bwt";
        const char K_STREE[]    = "stree";
        const char K_TOTAL[]    = "total";
    };

    typedef std::pair<monitor::size_dict, monitor::size_dict>  performance_monitor;



    class InputSpec{
    private:
        sdsl::bit_vector parse_bitstr(string& s){
            sdsl::bit_vector b(s.size());

            for(size_type i = 0; i < s.size(); i++)
                b[i] = ((unsigned char)s[i] - 48);
            return b;
        }

    public:
        string s_fname;

        InputSpec(string s_fn) : s_fname(s_fn){}

        string load_s(){
            string s;
            std::ifstream s_file {s_fname};
            while(s_file >> s)
                ;
            return s;
        }
    };


    // copied from http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
    class InputParser{
    public:
        std::string empty = "0";
        InputParser (int &argc, char **argv){
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
    public:
        bool lazy, sada;
        bool space_usage, time_usage;
        bool space_or_time_usage;
        bool answer;
        bool verbose;
        bool load_stree;
        size_type runs_progress, ms_progress;

        InputFlags(bool lazy_wl, bool sada_st, bool space, bool time_, bool ans, bool v, size_type runs_prgs, size_type ms_prgs, bool load_stree) :
        lazy{lazy_wl}, sada{sada_st},
        space_usage {space},
        time_usage {time_},
        answer {ans},
        verbose{v},
        load_stree{load_stree},
        runs_progress{runs_prgs}, ms_progress{ms_prgs}
        {
            space_or_time_usage = (space_usage || time_usage);
        }

        InputFlags (InputParser input) :
        lazy {input.getCmdOption("-lazy_wl") == "1"},             // lazy winer links
        sada {input.getCmdOption("-sada") == "1"},                // sadakane's suffix tree (rather tha ohleb)
        space_usage {input.getCmdOption("-space_usage") == "1"},  // space usage
        time_usage {input.getCmdOption("-time_usage") == "1"},    // time usage
        answer {input.getCmdOption("-answer") == "1"},            // answer
        verbose{input.getCmdOption("-verbose") == "1"},           // verbose
        load_stree{input.getCmdOption("-load_cst") == "1"},           // load CST of S and S'
        runs_progress{static_cast<size_type>(std::stoi(input.getCmdOption("-runs_progress")))},
        ms_progress{static_cast<size_type>(std::stoi(input.getCmdOption("-ms_progress")))}
        {
            space_or_time_usage = (space_usage || time_usage);
        }
    };


    int getmem (unsigned long *rss, unsigned long *vs)
    {
        //task_t task = MACH_PORT_NULL;
        struct task_basic_info t_info;
        mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

        if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))
        {
            return -1;
        }
        *rss = t_info.resident_size;
        *vs  = t_info.virtual_size;
        return 0;
    }

    void _dump_ms(sdsl::bit_vector& ms){
        auto get_ms = [] (sdsl::bit_vector& __ms, size_type __k) -> size_type {
            if(__k == -1)
                return (size_type) 1;
            return sdsl::select_support_mcl<1,1> (&__ms)(__k + 1) - (2 * __k);
        };

        for (size_type i = 0; i < ms.size() / 2; i++)
            cout << get_ms(ms, i) << " ";
        cout << endl;
    }

    void dump_ms(sdsl::bit_vector& ms){
        size_type k = 0;
        for (size_type i = 0; i < ms.size(); i++){
            if(ms[i] == 1){
                cout << i - (2*k) << " ";
                k += 1;
            }
        }
        cout << endl;
    }

}

#endif /* fd_ms_h */
