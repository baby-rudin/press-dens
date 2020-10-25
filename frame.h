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

    Frame(void);
    Frame(const Frame &frame);

    void        translate(VecReal R);

    REAL        frame_mass() const; 
    VecReal     mass_center() const;
};

#endif  // FRAME_H