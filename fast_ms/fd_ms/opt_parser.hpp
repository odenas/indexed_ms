//
//  opt_parser.hpp
//  fast_ms
//
//  Created by denas on 10/21/17.
//  Copyright © 2017 denas. All rights reserved.
//

#ifndef opt_parser_h
#define opt_parser_h


namespace fdms {

    // copied from http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
    class OptParser{
    public:
        std::string empty = "0";
        OptParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            return empty;
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
            != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
    };

};

#endif /* opt_parser_h */
