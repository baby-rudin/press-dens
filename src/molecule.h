#ifndef MOLECULE_H
#define MOLECULE_H

#include "type.hpp"
#include "vector.h"
#include <string>
#include <vector>
#include <iostream>

class Molecule
{
public:
    size_t                  nAtom;
    INTG                    charge;
    std::vector<INTG>       zval;
    std::vector<VecReal>    geom;


    Molecule(void);
    Molecule(const Molecule &mol);

    Molecule&       operator= (const Molecule &mol);
    Molecule&       operator+=(const Molecule &mol);

    std::string     print_geom() const;

    void            translate(VecReal R);
    
    REAL            bond(size_t atom1, size_t atom2) const;
    REAL            angle(size_t atom1, size_t atom2, size_t atom3) const;
    REAL            mol_mass() const;
    VecReal         mass_center() const;
};

Molecule    operator+ (const Molecule &mol1, const Molecule &mol2);

#endif // MOLECULE_H
