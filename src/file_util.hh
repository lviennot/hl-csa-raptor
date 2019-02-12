#ifndef FILE_UTIL_HH
#define FILE_UTIL_HH

#include "string.h"
#include <iostream>

#include <zlib.h>
#include "string_util.hh"


class file_or_gz {
private:
    bool gzipped;
    union { FILE *in; gzFile gz_in; };
    char *line;
    const size_t max_line_size;
    
public:
    file_or_gz(std::string filename, size_t mls=10000) : max_line_size(mls) {
        gzipped = filename.size() > 3
            && filename.substr(filename.size() - 3) == ".gz";
        if(gzipped) {
            gz_in = gzopen(filename.c_str(), "r");
        } else {
            in = filename != "-" ? fopen(filename.c_str(), "r") : stdin;
        }
        assert(gzipped ? gz_in  != nullptr : in  != nullptr);
        line = new char[max_line_size];
    }

    ~file_or_gz() {
        delete[] line;
    }

    // returns "" at end of file
    std::string get_line() {
        bool eof = (gzipped ? gzgets(gz_in, line, max_line_size)
                            : fgets(line, max_line_size, in)) == NULL;
        if (eof) return std::string("");
        else return std::string(line);
    }

    void close() {
        if (gzipped) {
            gzclose(gz_in);
        } else {
            if(in != stdin) fclose(in);
        }
    }

};

static std::vector<std::vector<std::string> >
read_tuples(const std::string filename, const size_t ncols) {
    std::vector<std::vector<std::string> > rows;
    file_or_gz in(filename);
    while(true) {
        std::string line = rtrim(in.get_line());
        if (line == "") break;
        auto v = split(line, ' ');
        if (v.size() != ncols) {
            std::cerr <<"wrong ncols : '"<< line <<"'\n";
            assert(v.size() == ncols);
        }
        rows.push_back(v);
    }
    return rows;
}
    
static std::vector<std::vector<std::string> >
read_csv(const std::string filename, const size_t ncol, ...) { // ... = column names
    std::vector<std::vector<std::string> > rows;
    file_or_gz in(filename);
    bool first = true;

    std::vector<std::string> colnames(ncol);
    va_list args;
    va_start(args, ncol);
    for (int j = 0; j < ncol; ++j) {
        colnames[j] = va_arg(args, char *);
    }
    va_end(args);

    std::vector<int> cols(ncol, -1);
    while(true) {
        std::string line = rtrim(in.get_line());
        if (line == "") break;
        if (first) {
            first = false;
            int i = 0;
            for (auto s : split(line, ',')) {
                for (int j = 0; j < ncol; ++j) {
                    if (s == colnames[j]) cols[j] = i;
                }
                ++i;
            }
            for (int j = 0; j < ncol; ++j) {
                if (cols[j] < 0)
                    throw std::invalid_argument("missing column : "
                                                + colnames[j] + " in:\n'" + line + "'");
            }
        } else {
            auto v = split(line, ',');
            std::vector<std::string> r(ncol);
            for (int j = 0; j < ncol; ++j) {
                assert(cols[j] < v.size());
                r[j] = v[cols[j]];
            }
            rows.push_back(r);
        }
    }
    in.close();
        
    return rows;
}


#endif // FILE_UTIL_HH
