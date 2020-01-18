#ifndef CPP_RZ_PIC_TEST_LOGGER_H
#define CPP_RZ_PIC_TEST_LOGGER_H

#include "../Tools/Logger.h"


void test_Logger() {
    int a = 1;
    char b = '#';
    element_logging(a, "logger_test.txt", " ");
    element_logging(b, "logger_test.txt", " ");
}

#endif //CPP_RZ_PIC_TEST_LOGGER_H
