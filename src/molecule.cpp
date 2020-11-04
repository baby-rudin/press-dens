#include "molecule.h"
#include "constant.h"
#include "atom.h"
#include <string>
#include <cmath>
#include <vector>


#define BUFF_LEN 1024

using namespace std;

//==================== Molecule ====================
Molecule::Molecule(void)
    : nAtom(0), charge(0)
{}

Molecule::Molecule(const Molecule &mol)
    : nAtom(mol.nAtom), charge(mol.charge), 
      zval(mol.zval), geom(mol.geom)
{}


// functions of Molecule
string  Molecule::print_geom() const
{
    string  geometry;
    char buff[BUFF_LEN];

    for (size_t i = 0; i < nAtom; i++) {
        sprintf(buff, " %-2s %14.8lf %14.8lf %14.8lf",
                get_atom_name(zval[i]).c_str(), geom[i].x, geom[i].y, geom[i].z);
        geometry += string(buff) + string("\n");
    }
    return  geometry;
}

void    Molecule::translate(VecReal R)
{
    for (size_t i = 0; i < nAtom; i++)
        geom[i] += R;
}

REAL    Molecule::bond(size_t a, size_t b) const
{
    return (geom[a] - geom[b]).len();
}

REAL    Molecule::angle(size_t a, size_t b, size_t c) const
{
    VecReal b_a = geom[b] - geom[a];
    VecReal b_c = geom[b] - geom[c];
    return  acos((b_a * b_c) / (b_a.len() * b_c.len()));
}

REAL    Molecule::mol_mass() const
{
    REAL    MassTotal = 0.0;

    for (size_t i = 0; i < nAtom; i++) {
        MassTotal += get_atom_mass(zval[i]);
    }

    return MassTotal;
}

VecReal Molecule::mass_center() const
{
    VecReal MassCent;
    REAL    MassTotal = 0.0;

    for (size_t i = 0; i < nAtom; i++) {
        MassTotal += get_atom_mass(zval[i]);
        MassCent += get_atom_mass(zval[i]) * geom[i];
    }
    MassCent /= MassTotal;

    return MassCent;
}

// operator =
Molecule&   Molecule::operator= (const Molecule &mol)
{
    nAtom   = mol.nAtom;
    charge  = mol.charge;
    zval    = mol.zval;
    geom    = mol.geom;

    return *this;
}

// operator +=
Molecule&   Molecule::operator+=(const Molecule &mol)
{
    nAtom   += mol.nAtom;
    charge  += mol.charge;
    zval.insert(zval.end(), mol.zval.begin(), mol.zval.end());
    geom.insert(geom.end(), mol.geom.begin(), mol.geom.end());

    return *this;
}

// operator +
Molecule    operator+ (const Molecule &mol1, const Molecule &mol2)
{
    Molecule    mol(mol1);
    mol += mol2;

    return mol;
}
