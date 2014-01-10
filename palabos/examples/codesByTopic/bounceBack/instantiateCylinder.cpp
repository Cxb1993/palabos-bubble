/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "palabos2D.h"
#include "palabos2D.hh"

#include "poiseuille.h"
#include "poiseuille.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

template<typename T>
class CylinderShapeDomain2D : public DomainFunctional2D {
public:
    CylinderShapeDomain2D(plint cx_, plint cy_, plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(util::sqr(radius))
    { }
    virtual bool operator() (plint iX, plint iY) const {
        return util::sqr(iX-cx) + util::sqr(iY-cy) <= radiusSqr;
    }
    virtual CylinderShapeDomain2D<T>* clone() const {
        return new CylinderShapeDomain2D<T>(*this);
    }
private:
    plint cx;
    plint cy;
    plint radiusSqr;
};

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    // Define numeric parameters.
    IncomprFlowParam<T> parameters (
        (T) 1e-2,  // uMax
        (T) 300.,  // Re
        40,        // N
        6.,        // lx
        1.         // ly 
    );

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              parameters.getNx(), parameters.getNy(),
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    createPoiseuilleBoundaries(lattice, parameters, *boundaryCondition);
    plint N = parameters.getResolution();
    defineDynamics(lattice, lattice.getBoundingBox(),
                   new CylinderShapeDomain2D<T>(N, N/2+2, N/4),
                   new BounceBack<T,DESCRIPTOR>);

    lattice.initialize();

    // Main loop over time iterations.
    for (plint iT=0; iT<80000; ++iT) {
        if (iT%1000==0) {
            pcout << "Writing image at dimensionless time " << iT*parameters.getDeltaT() << endl;
            ImageWriter<T> imageWriter("leeloo");
            imageWriter.writeScaledGif (
                    createFileName("velocity", iT, 6),
                    *computeVelocityNorm(lattice) );
            pcout << computeAverageEnergy(lattice) << endl;
        }
        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
