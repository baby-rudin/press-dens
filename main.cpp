#include <iostream>
#include <fstream>

#include "type.hpp"
#include "vector.h"
#include "atom.h"
#include "molecule.h"
#include "frame.h"
#include "constant.h"
#include "stropt.h"
#include "pdbflie.h"

using namespace std;

int main(int argc, char* argv[])
{
    string  filename;
    INTG    nFrame, nStart;
    REAL    temper;

    filename = string(argv[1]);
    nFrame  = atoi(argv[2]);
    nStart  = atoi(argv[3]);
    temper  = atof(argv[4]);

    // number of frames for calculating
    INTG    NumFrame = nFrame - nStart + 1;

    // record frames for calculating
    vector<Frame>   frames(NumFrame);

    // create PdbFile object
    PdbFile pdb(filename);

    // read in frames
    for (INTG i = 1; i < nStart; i++)
        pdb.get_one_frame();
    
    for (INTG i = 0; i < NumFrame; i++)
        frames[i] = pdb_read_one_frame();

    

    // deal with each frame
    for (INTG iFrame = 0; iFrame < NumFrame; iFrame++) {
        // init pair interaction
            // TODO

        // for each pair
            // TODO

    } // end for iFrame

    return 0;
}
