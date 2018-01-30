
#include <sdsl/suffix_trees.hpp>
#include <iostream>

#include "input_spec.hpp"
#include "opt_parser.hpp"

using namespace sdsl;
using namespace std;

typedef cst_sct3<> t_cst;
typedef typename t_cst::node_type node_type;
typedef typename t_cst::size_type size_type;

class InputFlags{
public:
    bool load_cst;
    size_type min_node_depth, max_str_depth, dump_size;

    InputFlags() {}

    InputFlags(const bool load_cst_, const size_type min_ndepth, const size_type max_sdepth, const size_type dump_size) :
        load_cst{load_cst_}, min_node_depth{min_ndepth}, max_str_depth{max_sdepth}, dump_size{dump_size} {}

    InputFlags(const InputFlags& i) :
        load_cst{i.load_cst}, min_node_depth{i.min_node_depth}, max_str_depth{i.max_str_depth}, dump_size{i.dump_size} {}

    InputFlags(const OptParser& args) :
        load_cst {(args.getCmdOption("-load_cst") == "1")},
        min_node_depth{static_cast<size_type>(std::stoi(args.getCmdOption("-min_node_depth")))},
        max_str_depth{static_cast<size_type>(std::stoi(args.getCmdOption("-max_str_depth")))},
        dump_size{static_cast<size_type>(std::stoi(args.getCmdOption("-dump_size")))}
    {}
};


string node_label(const node_type& v, const t_cst& cst)
{
    if(v == cst.root())
        return "";
    return extract(cst, v);
}

size_type lowest_with_wl(const t_cst& st, node_type v, const char c){
	node_type u = st.wl(v, c);
	while(u == st.root() && v != st.root()){
		v = st.parent(v);
		u = st.wl(v, c);
	}
	return st.node_depth(v);
}

size_type output_node(const node_type& v, const t_cst& cst, const size_type min_ndepth, const size_type max_sdepth)
{
	char c = 'z';

    size_type v_ndepth = cst.node_depth(v);
    size_type v_sdepth = cst.depth(v);
    assert(v_ndepth <= v_sdepth);

    if(v_sdepth > max_sdepth)
        return 0;

    size_type u_ndepth = lowest_with_wl(cst, v, c);
    assert(u_ndepth <= v_ndepth);

    if((v_ndepth - u_ndepth) < min_ndepth)
        return 0;

    string s = node_label(v, cst);

    cerr << (v_ndepth - u_ndepth) << " , " << v_sdepth << ",";
    char end_of_s = s[s.size() - 1];

    if (s[s.size() - 1] == '\000'){
		string s_cpy (s.substr(0, s.size() - 2));
		reverse(s_cpy.begin(), s_cpy.end());
        cout << s_cpy << c;
	} else {
		string s_cpy (s);
		reverse(s_cpy.begin(), s_cpy.end());
        cout << s_cpy << c;
	}
    return s.size();
}

void run(const InputSpec& ispec, const InputFlags& flags)
{
    t_cst cst;
	string s = ispec.load_s(true);

    construct_im(cst, s, 1);
    cerr << "done constructing cst of size " << cst.size() << endl;
    size_type size = 0;
    for (auto v : cst) {
        size_type curr_size = output_node(v, cst, flags.min_node_depth, flags.max_str_depth);
        if (curr_size > 0){
            size += curr_size;
            cerr << size << endl;
        }
        if(size > flags.dump_size)
            break;
    }
    /*
    auto v = cst.select_leaf(2);
    for (auto it = cst.begin(v); it != cst.end(v); ++it) {
        output_node(*it, cst);
    }
    cout<<"--"<<endl;
    v = cst.parent(cst.select_leaf(4));
    for (auto it = cst.begin(v); it != cst.end(v); ++it) {
        output_node(*it, cst);
    }
    */
}

int main(int argc, char **argv)
{
    OptParser input(argc, argv);
    InputSpec sfwd_spec;
    InputFlags flags;

    if(argc == 1){
        sfwd_spec = InputSpec("/home/brt/code/matching_statistics/indexed_ms/tests/input_stats_data/exploration/index_based_query/rep_2.s");
        flags = InputFlags(true, 10, 5, 1);
    } else {
        sfwd_spec = InputSpec(input.getCmdOption("-s_path"));
        flags = InputFlags(input);
    }
    run(sfwd_spec, flags);
}
