#ifndef help_h
#define help_h

#include <map>
#include <string>

using namespace std;

namespace fdms{
    string help__s_path{"\t-s_path <path to a text file>: input\n"};
    string help__t_path{"\t-t_path <path to a text file>: query\n"};
    string help__ms_path{"\t-ms_path <path to a binary>: the ms bit vector\n"};
    string help__load_cst{"\t-load_cst 1: load a previously computed cst. Will look for <s_path>.rev.stree\n"};
    string help__load_maxrep{"\t-load_cst 1: load a previously computed maxrep. Will look for <s_path>.rev.stree and <s_path>.fwd.stree\n"};
    string help__time_usage{"\t-time_usage 1: report time usage on the standard output\n"};
};
#endif