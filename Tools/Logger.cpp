#include "Logger.h"


void clear_file(const string& path) {
    ofstream output(path, ios::out);
    if (output) {
        return;
    } else {
        cout << "Can't clear file" << endl;
        throw;
    }
}

bool check_file(const string &path) {
    ifstream file(path);
    return file.is_open();
}