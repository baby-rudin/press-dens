#ifndef FRAME_H
#define FRAME_H

#include "type.hpp"
#include "molecule.h"
#include "vector.h"
#include <vector>

class Frame {
public:
    INTG                    nMol;
    std::vector<Molecule>   mols;
    std::vector<VecReal>    center;

    Frame(void);
    Frame(const Frame &fra);

    void        translate(VecReal R);

    VecReal     mass_center() const;

};

#endif  // FRAME_H