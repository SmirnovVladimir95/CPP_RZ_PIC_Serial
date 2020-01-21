#ifndef CPP_RZ_PIC_LOGGER_H
#define CPP_RZ_PIC_LOGGER_H

#include <string>
#include <iostream>
#include <fstream>
using namespace std;


template <class Type>
void element_logging(const Type& element, const string& path, string end = "\n", unsigned int mode = ios::app);

template <class Type>
void element_logging(const Type& element, const string& path, string end, unsigned int mode) {
    ofstream output(path, mode);
    if (output) {
        output << element << end;
    } else {
        cout << "Can't create/open file" << endl;
        throw;
    }
}

void clear_file(const string& path);

bool check_file(const string& path);


#endif //CPP_RZ_PIC_LOGGER_H
