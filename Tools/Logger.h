#ifndef CPP_RZ_PIC_LOGGER_H
#define CPP_RZ_PIC_LOGGER_H

#include <string>
#include <fstream>
using namespace std;


template <class T>
void element_logging(const T& element, const string& path, string end = "\n", unsigned int mode = ios::app);

template <class T>
void element_logging(const T& element, const string& path, string end, unsigned int mode) {
    ofstream output(path, mode);
    if (output) {
        output << element << end;
    }
}


#endif //CPP_RZ_PIC_LOGGER_H
