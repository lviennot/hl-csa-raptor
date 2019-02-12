#ifndef STRING_UTIL_HH
#define STRING_UTIL_HH

#include "string.h"
#include <iostream>


// trim from right
static std::string rtrim(std::string s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
    return s;
}

static std::vector<std::string> split(const std::string &s,
                                      const char delim) {
    std::vector<std::string> v;
    /*
      std::string field{""};
      for (auto c : s) {
      if (c != delim) field += c;
      else { v.push_back(field); field = ""; }
      }
      v.push_back(field);
    */
    size_t pos = 0, pos_prev = 0;
    while ((pos = s.find(delim, pos_prev)) != std::string::npos) {
        v.push_back(s.substr(pos_prev, pos - pos_prev));
        pos_prev = pos + 1;
    }
    v.push_back(s.substr(pos_prev, s.size() - pos_prev));
    return v;
};

template<typename T>
static void rev_vector(std::vector<T> &v) {
    int l = 0, r = v.size() - 1;
    while (l < r) {
        std::swap(v[l],v[r]);
        ++l; --r;
    }
}

#endif // STRING_UTIL_HH
