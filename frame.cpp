#include "frame.h"

using namespace std;

Frame::Frame(void)
    : nMol(0)
{}

Frame::Frame(const Frame &frame)
    : nMol(frame.nMol), mols(frame.mols)
{}

void    Frame::translate(VecReal R)
{
    for (INTG i = 0; i < nMol; i++)
        mols[i].translate(R);
}

REAL    Frame::frame_mass() const
{
    REAL    mass_total = 0.0;

    for (INTG i = 0; i < nMol; i++)
        mass_total += mols[i].mol_mass();

    return mass_total;
}

VecReal Frame::mass_center() const
{
    VecReal center;

    for (INTG i = 0; i < nMol; i++)
        center += mols[i].mass_center() * mols[i].mol_mass();
    
    center /= frame_mass();

    return center;
}