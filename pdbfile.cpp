#include "pdbfile.h"
#include "stropt.h"
#include <string>
#include <iostream>
#include <ifstream>

using namespace std;

PdbFile::PdbFile(std::string pdb_filename)
    : filename(pdb_filename), line_index(0)
{
    ifstream    file(filename);

    string  line;
    while (getline(file, line)) {
        remove_whitespace(line);
        if (!line.empty())
            content.push_back(line);
    }
    
    file.close();
}

PdbFile::~PdbFile()
{}

Frame   PdbFile::get_one_frame()
{
    Frame   ret_frame;

    if (line == 0)  line = 1;

    while (line < content.size()) {
        
    }

    return ret_frame;
}