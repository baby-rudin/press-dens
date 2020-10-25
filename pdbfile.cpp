#include "pdbfile.h"
#include "molecule.h"
#include "stropt.h"
#include <string>
#include <iostream>
#include <ifstream>
#include <sstream>


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

    if (line_index == 0)  line_index = 1;

    INTG        mol_num = 0;
    Molecule    mol_tmp;

    while (line_index < content.size()) {
        string  line = content[line_index];
        if (line == "END")  break;

        string  atom;
        INTG    atom_idx;
        string  atom_name;
        string  tip3u;
        INTG    mol_idx;
        REAL    x, y, z;
        REAL    zero1, zero2;
        string  u;

        istringstream is(line);

        is >> atom >> atom_idx >> atom_name >> tip3u
            >> mol_idx >> x >> y >> z >> zero1 >> zero2 >> u;

        if (mol_idx != mol_num) {

            if (mol_num != 0) {
                ret_frame.nMol ++;
                ret_frame.mols.push_back(mol_tmp);
            }

            mol_num = mol_idx;
            mol_tmp = Molecule();
        }

        mol_tmp.nAtom ++;
        mol_tmp.zval.push_back(get_atom_id(atom_name.substr(0,1)));
        mol_tmp.geom.push_back(VecReal(x, y, z));

        line_index ++;
    }

    return ret_frame;
}