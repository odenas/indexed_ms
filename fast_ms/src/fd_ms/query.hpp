#ifndef query_h
#define query_h

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

namespace fdms{
    class Query {
    protected:
        char* buffer; 
        size_t buff_idx, buff_size, q_size;
        ifstream query_fd;

        size_t _size(){
            query_fd.seekg(0, query_fd.end);
            return (size_t) (query_fd.tellg());
        }

        /* [from, from + buffer_isize) */
        void read_slice(const size_t from){
            size_t to = from + buff_size;
            to = (to > q_size ? q_size : to);

            //cerr << "reading [" << from << ", " << to << ")" << endl;
            query_fd.seekg(from);
            query_fd.read(buffer, buff_size);
            buff_idx = from;
        }

        size_t buffer_start() const { return buff_idx; }
        size_t buffer_end() const { return buff_idx + buff_size; }

        void index_in_bounds(size_t index){
            if(index >= size()){
                cerr <<  "index > size()" << endl;
                throw "index > size()";
            }
        }

    public:
        string fname;
        Query(const string fname, const size_t buffer_size) {
            query_fd.open(fname, ios::binary | ios::in);
            q_size = _size();
            buff_size = (buffer_size > q_size ? q_size : buffer_size);
            buffer = new char[buff_size];

            read_slice(q_size - buff_size);
            if(buffer[buff_size - 1] == '\n'){
                q_size -= 1;
                buff_size = (buffer_size > q_size ? q_size : buffer_size);
                buffer = new char[buff_size];
                read_slice(q_size - buff_size);
            }
        }

        ~Query() { query_fd.close(); }

        size_t size() const { return q_size; }

    };

    class Query_rev : public Query {
    public:
        Query_rev(const string fname, const size_t buffer_size) : Query(fname, buffer_size) {}

        char operator[](size_t index){
            index_in_bounds(index);
            if(buffer_start() <= index and index < buffer_end())
                return buffer[index - buff_idx];
            read_slice(index > buff_size - 1 ? index - buff_size + 1 : 0);
            return (*this)[index];
        }
    };

    class Query_fwd : public Query {
    public:
        Query_fwd(const string fname, const size_t buffer_size) : Query(fname, buffer_size) {}

        char operator[](size_t index){
            if(buffer_start() <= index and index < buffer_end())
                return buffer[index - buff_idx];
            read_slice(index);
            return (*this)[index];
        }
    };
}
#endif /* query_h */
