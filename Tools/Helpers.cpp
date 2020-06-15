#include "Helpers.h"


vector<scalar> string_to_numeric_vector(const string &s) {
    string line = s;
    size_t sz;
    vector<scalar> numbers;
    while (line.length() > 0) {
        numbers.push_back(stod(line, &sz));
        line = line.substr(sz);
    }
    return numbers;
}
