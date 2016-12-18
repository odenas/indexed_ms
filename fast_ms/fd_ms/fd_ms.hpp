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
        std::string empty = "";
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
        bool space_usage, mach_space_usage;
        bool time_usage;
        bool answer;
        bool verbose;

        InputFlags(bool space, bool mach_space, bool time_, bool ans, bool v) :
        space_usage {space},
        mach_space_usage {mach_space},
        time_usage {time_},
        answer {ans},
        verbose{v}
        {}

        InputFlags (InputParser input) :
        space_usage {input.getCmdOption("-s") == "1"},      // space usage
        mach_space_usage {input.getCmdOption("-S") == "1"}, // resident/vm memory usage
        time_usage {input.getCmdOption("-t") == "1"},       // time usage
        answer {input.getCmdOption("-a") == "1"},           // answer
        verbose{input.getCmdOption("-v") == "1"}            // verbose
        {}
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

    void dump_ms(sdsl::bit_vector& ms){
        auto get_ms = [] (sdsl::bit_vector& __ms, size_type __k) -> size_type {
            if(__k == -1)
                return (size_type) 1;
            return sdsl::select_support_mcl<1,1> (&__ms)(__k + 1) - (2 * __k);
        };

        for (size_type i = 0; i < ms.size() / 2; i++)
            cout << get_ms(ms, i);
        cout << endl;
    }


}

#endif /* fd_ms_h */
