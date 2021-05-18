

#ifndef RANGE_QUERY_HPP
#define RANGE_QUERY_HPP

enum class RangeAlgorithm {trivial, djamal, none};
enum class RangeOperation {r_sum, r_max};

template<typename size_type>
class rq_result {
    public:
        size_type index;
        size_type value;

        rq_result() : index{0}, value{0}{}

        rq_result(const rq_result& r) : index{r.index}, value{r.value}{}

        rq_result(const size_type index, const size_type value) : index{index}, value{value}{}
};

#endif /* RANGE_QUERY */



