#ifndef PDBFILE_H
#define PDBFILE_H

#include "type.hpp"
#include "frame.h"
#include "vector.h"
#include <vector>
#include <string>

class PdbFile {
public:
    std::string                 filename;
    std::vector<std::string>    content;
    INTG                        line_index;

    PdbFile(std::string pdb_filename);
    ~PdbFile();

    Frame   get_one_frame();
};

#endif  // PDBFILE_H