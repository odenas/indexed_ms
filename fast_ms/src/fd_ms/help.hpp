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
    string help__ridx_path{"\t-ridx_path <path to binary file>: the range query index; only used if block_size > 0\n"};
    string help__from_idx{"\t-from_idx <non-negative int>: start of a 0-based half-open interval [from_idx, to_idx)\n"};
    string help__to_idx{"\t-to_idx <non-negative int>: end of a 0-based half-open interval [from_idx, to_idx)\n"};
    string help_block_size{"\t-block_size <non-negative int>: range query index block size. If 0, do not use index.\n"};
};
#endif