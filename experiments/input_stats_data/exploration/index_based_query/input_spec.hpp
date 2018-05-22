#ifndef input_spec_h
#define input_spec_h


using namespace std;

class InputSpec{
public:
    string s_fname, fwd_cst_fname, rev_cst_fname, rev_maxrep_fname, rev_elst_fname, fwd_nwdlst_fname;
    
    InputSpec(){}
    
    InputSpec(const InputSpec &is) :
    s_fname{is.s_fname},
    fwd_cst_fname{is.fwd_cst_fname},
    rev_cst_fname{is.rev_cst_fname},
    rev_maxrep_fname{is.rev_maxrep_fname},
    rev_elst_fname{is.rev_elst_fname},
    fwd_nwdlst_fname{is.fwd_nwdlst_fname} {}
    
    InputSpec(string s_fn) : s_fname(s_fn){
        fwd_cst_fname = s_fname + ".fwd.stree";
        rev_cst_fname = s_fname + ".rev.stree";
        rev_maxrep_fname = s_fname + ".rev.maxrep";
        rev_elst_fname = s_fname + ".rev.elst";
        fwd_nwdlst_fname = s_fname + ".fwd.nwdlst";
    }
    
    InputSpec& operator=(const InputSpec& other){
        if(this != &other){
            s_fname = string(other.s_fname);
            fwd_cst_fname = string(other.fwd_cst_fname);
            rev_cst_fname = string(other.rev_cst_fname);
            rev_maxrep_fname = string(other.rev_maxrep_fname);
            rev_elst_fname = string(other.rev_elst_fname);
            fwd_nwdlst_fname = string(other.fwd_nwdlst_fname);
        }
        return *this;
    }
    
    string load_s(bool reverse = false) const {
        string s;
        std::ifstream s_file {s_fname};
        while(s_file >> s)
            ;
        
        if(reverse)
            reverse_in_place(s);
        return s;
    }
    
    static void reverse_in_place(string& s){
        unsigned long long n = s.size();
        
        for(int i = 0; i < n / 2; i++){
            char c = s[i];
            s[i] = s[n - 1 - i];
            s[n - 1 - i] = c;
        }
    }
};
#endif
