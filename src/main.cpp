#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>

#include "type.hpp"
#include "vector.h"
#include "atom.h"
#include "molecule.h"
#include "frame.h"
#include "constant.h"
#include "stropt.h"
#include "pdbfile.h"

using namespace std;

int main(int argc, char **argv)
{
    REAL    dR  = 0.03;
    INTG    NdR = 100;

    REAL    q_O = -0.834;
    REAL    q_H = 0.417;
    REAL    f   = 138.95;   // kj/mol nm

    REAL    epsla_O = 0.64;
    REAL    epsla_H = 0.192;
    REAL    sigma_O = 0.354;
    REAL    sigma_H = 0.0449;

    // REAL    epsla_O = 0.636;
    // REAL    epsla_H = 0.000;
    // REAL    sigma_O = 0.315;
    // REAL    sigma_H = 0.000;

    string  filename;
    cout << "Input file name: ";
    cin >> filename;

    // create PdbFile object
    PdbFile pdb(filename);

    // read in frames
    vector<Frame>   frames;
    while (true) {
        Frame frame_tmp = pdb.get_one_frame();
        if (frame_tmp.nMol == 0) break;
        frames.push_back(frame_tmp);
    }

    cout << "Number of frames: " << frames.size() << endl;

    INTG nStart;
    INTG nEnd;
    cout << "Start at frame: ";
    cin >> nStart;
    cout << "End at frame: ";
    cin >> nEnd;

    // pressure records
    vector<REAL>    pressU(NdR+1, 0.0);
    vector<REAL>    pressK(NdR+1, 0.0);
    vector<REAL>    pressN(NdR+1, 0.0);

    // deal with each frame
    for (INTG iFrame = nStart - 1; iFrame < nEnd; iFrame++) {
        // read in frame
        Frame   this_frame = frames[iFrame];

        // translate to mass center
        this_frame.translate(-1.0 * this_frame.mass_center());


        for (INTG ix = 0; ix < this_frame.nMol; ix++)
        for (INTG jx = 0; jx < this_frame.mols[ix].nAtom; jx++)
            this_frame.mols[ix].geom[jx] /= 10.0;


        for (INTG im = 0; im < this_frame.nMol; im++) {
            Molecule imol = this_frame.mols[im];
            VecReal ic = imol.mass_center();   // mass center

            INTG bin = (INTG) (ic.len() / dR);
            REAL dV  = M_PI * 4.0/3.0 * ( pow((bin+1)*dR,3) - pow(bin*dR,3) );
            if (bin < NdR) pressK[bin] += 1/dV;
        }

        // for each pair
        for (INTG im = 0; im < this_frame.nMol - 1; im++) {
            Molecule imol = this_frame.mols[im];
            VecReal ic = imol.mass_center();   // mass center


            for (INTG jm = im + 1; jm < this_frame.nMol; jm++) {
                Molecule jmol = this_frame.mols[jm];
                VecReal jc = jmol.mass_center();   // mass center
                
                // force vector
                VecReal fij_xyz;

                for (INTG ia = 0; ia < imol.nAtom; ia++) {
                    for (INTG ja = 0; ja < jmol.nAtom; ja++) {
                        
                        // Vector atom i --> atom j
                        VecReal raij = imol.geom[ia] - jmol.geom[ja];

                        // parameters for atom i and j
                        INTG izval  = imol.zval[ia];
                        INTG jzval  = jmol.zval[ja];

                        REAL iq     = (izval == 1 ? q_H : q_O);
                        REAL jq     = (jzval == 1 ? q_H : q_O);

                        REAL iepsla = (izval == 1 ? epsla_H : epsla_O);
                        REAL jepsla = (jzval == 1 ? epsla_H : epsla_O);

                        REAL isigma = (izval == 1 ? sigma_H : sigma_O);
                        REAL jsigma = (jzval == 1 ? sigma_H : sigma_O);

                        REAL epslaij = sqrt(iepsla * jepsla);
                        REAL sigmaij = (isigma + jsigma) / 2.0;


                        // Coulomb forces
                        REAL    fijq = f * iq * jq  / raij.len2();
                        fij_xyz += (fijq * raij) / raij.len();

                        // vdw forces, LJ potential
                        REAL    first_part  = pow(sigmaij,12)/pow(raij.len(),13);
                        REAL    second_part = pow(sigmaij,6)/pow(raij.len(),7);
                        REAL    fijv = 12.0 * epslaij * (first_part - second_part);
                        fij_xyz += (fijv * raij) / raij.len();

                        // REAL    sr2  = isigma * jsigma / raij.len2();
                        // REAL    sr6  = sr2 * sr2 * sr2;
                        // REAL    sr12 = sr6 * sr6;
                        // REAL    fijv = (sr12*2.0-sr6)/raij.len2() * sqrt(iepsla * jepsla) * 24.0;
                        // fij_xyz += (fijv * raij);
                    }
                }


                // interaction
                REAL    ffij = (fij_xyz * (ic-jc)) / (ic-jc).len();

                // cout << ffij << endl;

                // distance between line_ij and point O
                REAL d = sqrt(ic.len2() - (ic*(ic-jc))*(ic*(ic-jc))/(ic-jc).len2());
                // cout << d << endl;

                for (INTG idR = 1; idR <= NdR; idR++) {
                    REAL    R = idR * dR;

                    REAL cosine = R > d ? sqrt(R*R - d*d)/R : 0.0;

                    // cout << cosine << endl;

                    if ( (ic.len() - R) * (jc.len() - R) < 0.0 ) {
                        pressU[idR] += ffij * cosine * 1 / (4.0*M_PI*R*R);
                    }
                    else if (ic.len() > R && jc.len() > R
                            && fabs(ic.len2() - jc.len2()) < R * R
                            && d < R) {
                        pressU[idR] += ffij * cosine * 2 / (4.0*M_PI*R*R);
                    }
                } // end for idR

            } // end for jm
        } // end for im

    } // end for iFrame



    for (INTG i = 0; i < pressN.size(); i++) {
        pressK[i] *= 8.314 * 298 / 1000.0;
        pressK[i] /= (nEnd - nStart + 1);
        pressU[i] /= (nEnd - nStart + 1);

        pressN[i] = pressU[i] + pressK[i];
    }

    for (INTG i = 1; i <= NdR; i++) {
        printf("%8.2lf ", i * dR);
    }

    cout << endl;
    cout << endl;

    for (INTG i = 1; i <= NdR; i++) {
        printf("%8.2lf ", pressU[i]);
    }

    cout << endl;
    cout << endl;

    for (INTG i = 1; i <= NdR; i++) {
        printf("%8.2lf ", pressN[i]);
    }

    return 0;
}
