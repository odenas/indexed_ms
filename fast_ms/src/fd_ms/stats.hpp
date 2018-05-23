
#ifndef STATS_HPP
#define STATS_HPP


#include <iostream>
#include <fstream>

using namespace std;

template<typename cst_t, typename maxrep_t>
class NodeProperty {
public:
    typedef typename cst_t::node_type node_type;
    typedef typename cst_t::char_type char_type;

    const cst_t& m_st;
    const maxrep_t& m_mrep;

    NodeProperty(const cst_t& cst, const maxrep_t& mr) : m_st{cst}, m_mrep{mr}
    {
    }

    static bool is_wide(const node_type v) {
        return (((v.i) >> 8) != ((v.j) >> 8));
    }

    bool has_wl(const node_type v, char c) const {
        node_type u = m_st.double_rank_nofail_wl(v, c);
        return !m_st.is_root(u);
    }

    bool is_max(const node_type v) const {
        return m_mrep.is_maximal(v);
    }

    string runs_node_label(const node_type v, const char_type c) const {
        string ch = {(char) c};
        if (c < 'a')
            ch = "?";
        string key = (ch + "_" + // char
                "na" + "_" + // maximality
                (has_wl(v, c) ? "wl" : "nowl") + "_" + // has wl(c)
                (is_wide(v) ? "wide" : "narrow")); // interval width
        return key;
    }

    string ms_node_label(const node_type v, const char_type c) const {
        string ch = {(char) c};
        if (c < 'a')
            ch = "?";
        string key = (ch + "_" + // char
                (is_max(v) ? "maxrep" : "nomaxrep") + "_" + // maximality
                (has_wl(v, c) ? "wl" : "nowl") + "_" + // has wl(c)
                (is_wide(v) ? "wide" : "narrow")); // interval width
        return key;
    }
};

template<typename cst_t, typename maxrep_t>
class Stats {
    typedef typename cst_t::node_type node_type;
    typedef typename cst_t::size_type size_type;
    typedef typename cst_t::char_type char_type;

    typedef map<size_type, size_type> hist_type;
    typedef map<string, size_type> counter_type;

private:
    size_type len_s, len_t;

    template<typename map_t>
    void dump_map(std::ostream& out, string measuring, string where, map_t& values) {
        for (auto item : values)
            out << len_s << "," << len_t << "," << measuring << "," << where << "," << item.first << "," << item.second << endl;
    }

    void _reg_pseq(const cst_t& st, const NodeProperty<cst_t, maxrep_t>& np,
            const node_type& from, const node_type& to, char_type c,
            const bool runs) {

        assert(from != to);
        size_type i = 0;
        node_type u = from, v = to;
        while (u != v) {
            i += 1;
            if (runs)
                runs_wl_calls[np.runs_node_label(u, c)] += 1;
            else
                ms_wl_calls[np.ms_node_label(u, c)] += 1;

            u = st.parent(u);
        }
        if (runs)
            runs_pcalls_seq[i] += 1;
        else
            ms_pcalls_seq[i] += 1;
    }

public:
    hist_type runs_pcalls_seq, ms_pcalls_seq, ms_wlcalls_seq;
    counter_type runs_wl_calls, ms_wl_calls;

    Stats(size_type sl, size_type tl) : len_s{sl}, len_t{tl}
    {
    };

    void register_runs_pseq(const cst_t& st, const NodeProperty<cst_t, maxrep_t>& np,
            const node_type& from, const node_type& to, char_type c) {
        _reg_pseq(st, np, from, to, c, true);
    }

    void register_ms_pseq(const cst_t& st, const NodeProperty<cst_t, maxrep_t>& np,
            const node_type& from, const node_type& to, char_type c) {
        _reg_pseq(st, np, from, to, c, false);
    }

    void dump_stats(std::ostream &out) {
        out << "len_s,len_t,measuring,where,key,value" << endl;
        dump_map<counter_type>(out, "wl_calls", "runs", runs_wl_calls);
        dump_map<counter_type>(out, "wl_calls", "ms", ms_wl_calls);

        dump_map<hist_type>(out, "pseq", "runs", runs_pcalls_seq);
        dump_map<hist_type>(out, "pseq", "ms", ms_pcalls_seq);
        dump_map<hist_type>(out, "wlseq", "ms", ms_wlcalls_seq);
    }
};


#endif /* stats_h */
