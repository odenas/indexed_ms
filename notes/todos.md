# TODOS from email thread 11/30-12/6

## Optimizations
- I don't think we need any additional data structure on `runs`, since we just use it to compute, for every position `k`, the position of the first zero to the right of `k`, and this probably cumulates to a single linear scan of
runs? **DONE**
- `runs` can be stored in the last `|t|-1` bits of bitvector `ms`, i.e. we don't need an additional bitvector for it?
-  should we try to come up with a realistic class of inputs for which compressed `stringDepth` (in sdsl and compressed ST) takes too much space in practice? maybe highly repetitive strings? (i haven't thought about it yet).

TODO:

> `runs` can be stored at the end of `ms`. find a (repetitive) for which `stringDepth` takes a lot of space

## Bp-support and Weiner links
### Problem
So, `parent` is not the only needed operation.

We also need some suffix-tree support for `wl`, although `wl` is just two rank queries on a bwt-interval. We need some mechanism  to convert between suffix-tree node identifiers and bwt-intervals in both directions. 

Given a suffix tree node identifier `n`, 

- `l` <- `rank(leftmost(n))` and `r` <- `rank(rightmost(n))` (suffix tree --> bwt)
- a `L <- leafSelect(backward(l))` and `R <- leafSelect(backward(r))` (bwt --> suffix tree)
- `lca(L, R)`

So, in total `wl`  needs 9 operations which is a lot.

- 4 opeations to convert a node to a bwt-interval
- 2 rank queries to do the wl operation on bwt intervals 
- 3 operations to get node id. 


> This can be seen as a bad thing from an angle,  but from another, it seems there is a lot of room for optimization. 

In [cst++](https://www.researchgate.net/publication/221580034_CST) they use a representation of suffix tree as intervals  of suffixes (intervals of leaves also called bwt-intervals or SA-intervals, LCP-intervals). Basically, the only thing needed is NSV/PSV operations (next smaller/larger value) on the LCP array (see page 11 in the paper above). 

This is implemented in `cst_sct3.hpp` and indeed has less costly `wl` but more costly `parent` operations.

### Lazy Suffix Tree
`wl` in `cst_sct3.hpp` can be optimized. In practice if we do `d` consecutive `wl` operations, then between the two operations we will not need to recompute all the fields of the node object. Essentially, except for the last `wl` operations, we will only need to compute the two rank queries (we can skip computation of `ipos`, `icpos` and `jp1pos` and only need to compute `i` and `j`). I think the same optimization can be applied to `cst_sada.hpp`. 
 
The second optimization for `wl` is that the two rank queries can often be solved  more efficiently into one combined operation, since when the two rank queries are on  small range, two independent rank queries would repeat many operations twice.  I think sdsl might have such an optimized `double_rank` query (I remember I saw it in a paper in which Simon Gog was co-author). 

> implement the "lazy suffix tree" and the `double_rank` operation

## maxrep pre-processing (ignore for now)

**F**: the bwt interval of every node of the ST that is not a
maximal repeat contains exactly one distinct character. so if we can
afford n additional bits we could avoid rank operations altogether at
such nodes? (not sure this is faster in practice...)

**D**: I do not think we can avoid the two rank operations, but we could avoid 
just one of the two rank operations on the bwt, but require on the other 
hand require a rank operation to determine whether the given bwt interval 
contains one single character. Overall I believe an optimized double rank 
operation will optimize well such cases. In general I think most of the intervals  will be small so that the computations for the two two rank queries will often  overlap by a lot. 

**F**: i meant that we precompute all maximal repeats and mark all maxrep
nodes in an auxiliary bitvector of length n. then, when we are in the
bwt interval of such maxrep, we know from the corresponding bit that
it contains exactly one distinct character.

**D**: Nice idea. The only potential problem is that this is incompatible with the idea  that we do not systematically recompute the node preorder after each `wl`, but  only recompute it after when we need to do parent. Notice that this optimization  is well motivated by the fact that we have more `wl` than `parent` operations, since a `parent` operation could potentially decrease the string depth by more than one, whereas a successful `wl` will increase it by exactly one and the amount of depth increase is in total at least the amount of decrease. 

... some other emails for future 


## Failing `wl` queries
In the case of matching statistics problem we will have many failing `wl` queries (about as many as parent queries). It thus makes sense to optimize also for negative queries to be fast. I think this is easy for the case of the `rank_and_get(B,c,i)` query, which I think we should rename to `rank_and_check(B,c,i)`. That is, when traversing down the WT, we stop as soon as we have determined that `B[i]` does not equal `c`,  which could happen even at first level (without even counting in the first bitvector). I think this `rank_and_check(B,c,i)` is even more useful when we are at a leaf, since there we just can have one single `rank_and_check(B,c,i)` instruction instead of two rank queries. 

> Optimize failing `wl` queries (see above how)

## Parallelization