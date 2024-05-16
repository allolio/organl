#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <random>
#include <chrono>
#include <omp.h>

using namespace std;
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()
#define HAVE_CONFIG_NG
#define STR_MAXCHARS 2000
#define INVALID_INDEX -1
#define INVALID_EDGE 0xFFFFFFFF

    int Tokenize(const string& str,
                      vector<string>* tokens = NULL,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
    int ntokens=1;
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        if (tokens!=NULL) tokens->push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
	ntokens++;
    }
    return ntokens;
}

bool isNumber(string s)
{
    std::size_t char_pos(0);

    // skip the whilespaces
    char_pos = s.find_first_not_of(' ');
    if (char_pos == s.size()) return false;


    // check the significand
    if (s[char_pos] == '+' || s[char_pos] == '-') ++char_pos; // skip the sign if exist

    int n_nm, n_pt;
    for (n_nm = 0, n_pt = 0; std::isdigit(s[char_pos]) || s[char_pos] == '.'; ++char_pos) {
        s[char_pos] == '.' ? ++n_pt : ++n_nm;
    }
    if (n_pt>1 || n_nm<1) // no more than one point, at least one digit
        return false;

    // skip the trailing whitespaces
    while (s[char_pos] == ' ') {
        ++ char_pos;
    }

    return char_pos == s.size();  // must reach the ending 0 of the string
}


#include "minimath.hpp"
/*#ifndef HAVE_FT
//#include "filetool.hpp"
#endif*/
