//
// Created by Luis Yanes (EI) on 23/11/2017.
//

#ifndef BSG_FILESYSTEM_H
#define BSG_FILESYSTEM_H

#include <string>
#include <cerrno>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>

namespace sglib {
    bool check_or_create_directory(std::string &output_prefix);
    void remove_directory(std::string path);
}
#endif //BSG_FILESYSTEM_H
