//
//  cst_iterator.hpp
//  fast_ms
//
//  Created by denas on 10/16/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef cst_iterator_h
#define cst_iterator_h


namespace fdms {

    template<class Cst>
    class cst_dfslr_iterator : public std::iterator<std::output_iterator_tag, typename Cst::node_type> {
    public:
        typedef typename Cst::node_type value_type;
        typedef typename Cst::size_type size_type;
        typedef const value_type const_reference;
        typedef cst_dfslr_iterator<Cst> iterator;

    private:
        const Cst* m_cst; // Pointer to the cst.
        value_type nextnode;
        bool direction_down, node_processed;
        value_type currnode;

    public:
        size_type nodes_visited;

        //! Constructor

        /*!
         * \param cst   Pointer to the compressed suffix tree.
         */
        cst_dfslr_iterator(const Cst* cst) {
            m_cst = cst;
            currnode = m_cst->root();
            direction_down = true;
            node_processed = false;
            nodes_visited = 0;
            ++(*this);
        }

        //! Returns the current number of nodes in the queue.

        size_type size()const {
            return m_cst->size();
        }

        //! Method for dereferencing the iterator.

        const_reference operator*()const {
            return currnode;
        }

        //! Prefix increment of the iterator.

        iterator& operator++() {
            Cst st = *m_cst;

            do {
                if (direction_down) {
                    if (!st.is_leaf(currnode)) {
                        // process currnode
                        if (!node_processed) {
                            node_processed = true;
                            return *this;
                        } else {
                            nodes_visited += 1;
                            node_processed = false;
                        }
                    }

                    nextnode = st.first_child(currnode);
                    if (st.is_root(nextnode))
                        direction_down = false;
                    else
                        currnode = nextnode;
                } else {
                    nextnode = st.sibling(currnode);
                    if (st.is_root(nextnode)) {
                        currnode = st.parent(currnode);
                    } else {
                        currnode = nextnode;
                        direction_down = true;
                    }
                }
            } while (!st.is_root(currnode));
            return *this;
        }

        //! Postfix increment of the iterator.

        iterator operator++(int) {
            iterator it = *this;
            ++(*this);
            return it;
        }

        bool end() {
            return (m_cst->is_root(currnode) && nodes_visited > 0);
        }

        //! Equality operator.

        bool operator==(const iterator& it)const {
            if (m_cst->size() != it.m_cst->size()) { // if the queue size is different
                return false; // the state of the to iterator are different
            }
            return (it.m_cst == m_cst); // iterator belongs to the same cst
        }

        //! Inequality operator.

        bool operator!=(const iterator& it)const {
            return !(*this == it);
        }
    };
};
#endif /* cst_iterator_h */
