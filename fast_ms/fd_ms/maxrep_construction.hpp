//
//  maxrep_construction.hpp
//  fast_ms
//
//  Created by denas on 5/22/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef maxrep_construction_h
#define maxrep_construction_h

namespace fdms{
    template <typename cst_t, typename vector_t, typename node_type>
    std::pair<size_t, size_t> build_maxrep_ohleb_debug(const cst_t& st, vector_t& maxrep){ // vanilla
        node_type currnode = st.root(), nextnode = st.root();
        bool direction_down = true;
        char c = 0;
        std::pair<size_t, size_t> maximal_count(0, 0);

        do{
            if(direction_down){
                if(!st.is_leaf(currnode)){ // process currnode
                    c = st.csa.bwt[currnode.j];
                    std::pair<size_t, size_t> ni = st.csa.bwt.double_rank(currnode.i, currnode.j + 1, c);
                    size_t count = ni.second - ni.first;
                    if(count != currnode.j - currnode.i + 1){
                        maxrep[currnode.i] = maxrep[currnode.j] = 1;
                        maximal_count.first += 1;
                    } else {
                        maximal_count.second += 1;
                    }
                }

                nextnode = st.first_child(currnode);
                if(st.is_root(nextnode))
                    direction_down = false;
                else
                    currnode = nextnode;
            } else {
                nextnode = st.sibling(currnode);
                if(st.is_root(nextnode)) {
                    currnode = st.parent(currnode);
                } else {
                    currnode = nextnode;
                    direction_down = true;
                }
            }
        } while(!st.is_root(currnode));
        return maximal_count;
    }

    template <typename cst_t, typename vector_t, typename node_type>
    void build_maxrep_ohleb(const cst_t& st, vector_t& maxrep){
        node_type currnode = st.root(), nextnode = st.root();
        bool direction_down = true, currnode_maximal = false;
        char c = 0;

        do{
            if(direction_down){
                if(!st.is_leaf(currnode)){ // process currnode
                    //c = st.csa.bwt[currnode.j];
                    //auto ni = st.csa.bwt.double_rank(currnode.i, currnode.j + 1, c);
                    //size_t count = ni.second - ni.first;
                    //if(count != currnode.j - currnode.i + 1){ // maximal
                    //    maxrep[currnode.i] = maxrep[currnode.j] = currnode_maximal = 1;
                    int wl_count = 0;
                    for(int i=1; i<st.csa.sigma; i++)
                        wl_count += (st.is_root(st.single_rank_wl(currnode, st.csa.comp2char[i])) ? 0 : 1);
                    if(wl_count > 1){
                        maxrep[currnode.i] = maxrep[currnode.j] = currnode_maximal = 1;

                        // try going down the subtree
                        nextnode = st.first_child(currnode);
                        if(st.is_root(nextnode) or !currnode_maximal)
                            direction_down = false;
                        else
                            currnode = nextnode;
                    } else { // no node in subtree is maximal
                        direction_down = false;
                    }
                } else {
                    direction_down = false;
                }

            } else {
                nextnode = st.sibling(currnode);
                if(st.is_root(nextnode)) {
                    currnode = st.parent(currnode);
                } else {
                    currnode = nextnode;
                    direction_down = true;
                }
            }
        } while(!st.is_root(currnode));
    }
}

#endif /* maxrep_construction_h */
