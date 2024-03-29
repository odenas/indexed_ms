#ifndef help_h
#define help_h

#include <map>
#include <string>

using namespace std;

namespace fdms{
    string help__s_path{"\t-s_path <path to a text file>: input\n"};
    string help__t_path{"\t-t_path <path to a text file>: query\n"};
    string help__ms_path{"\t-ms_path <path to a binary>: the ms bit vector\n"};
    string help__freq_path{"\t-freq_path <path to a binary>: the MS frequency vector\n"};
    string help__load_cst{"\t-load_cst 1: load a previously computed cst. Will look for <s_path>.rev.stree\n"};
    string help__load_maxrep{"\t-load_maxrep 1: load a previously computed maxrep. Will look for <s_path>.rev.stree and <s_path>.fwd.stree\n"};
    string help__time_usage{"\t-time_usage 1: report time usage on the standard output\n"};
    string help__ridx_path{"\t-ridx_path <path to binary file>: the range query index; only used if block_size > 0\n"};
    string help__from_idx{"\t-from_idx <non-negative int>: start of a 0-based half-open interval [from_idx, to_idx)\n"};
    string help__to_idx{"\t-to_idx <non-negative int>: end of a 0-based half-open interval [from_idx, to_idx)\n"};
    string help__block_size{"\t-block_size <non-negative int>: range query index block size. If 0, do not use index.\n"};
    string help__compression{"\t-compression <compression type>: One of: rrr, hybrid, rle, delta, succint, nibble.\n"};
    string help__algo{"\t-algo <trivial|djamal> whether to use the trivial or djamal algorithm.\n"};
    string help__rangeop{"\t-op <sum|max> range operator sum or max of elements in a range.\n"};
    // lossy compression specific
    string help__threshold{"\t-threshold <1> loose ms entries below this value.\n"};
    string help__nzeros{"\t-nZeroes <non-negative int> max number of zeros pre-computed in a table.\n"};
    string help__nones{"\t-nOnes <non-negative int> max number of ones pre-computed in a table.\n"};
    string help__negative{"\t-negative <1> allow negative values.\n"};
    string help__greddy{"\t-greedy <1> use the greedy strategy for too large windows.\n"};
    string help__verbose{"\t-verbose <1> be verbose.\n"};
};
#endif
