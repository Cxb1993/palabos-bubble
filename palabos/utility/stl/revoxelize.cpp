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

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;

const int typeOfVoxelization = voxelFlag::inside;
const plint extraLayer = 0;
const plint blockSize = 20;
const plint extendedEnvelopeWidth = 2;
const plint borderWidth = 1;
const plint margin = 4; // Extra margin of allocated cells around the obstacle. 

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(false);
    plint resolution=0;
    plint referenceDirection=0;

    string stlFileName, outFileName;
    try {
        global::argv(1).read(stlFileName);
        global::argv(2).read(outFileName);
        global::argv(3).read(resolution);
        global::argv(4).read(referenceDirection);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " inputSTL.stl outputSTL.stl resolution referenceDirection[0,1,2]" << std::endl;
        exit(-1);
    }

    TriangleSet<T>* triangleSet = 0;
    try {
        triangleSet = new TriangleSet<T>(stlFileName, FLT);
    }
    catch (PlbIOException& exception) {
        pcout << "Error, could not read STL file " << stlFileName
              << ": " << exception.what() << std::endl;
        exit(-1);
    }
    DEFscaledMesh<T> mesh(*triangleSet, resolution, referenceDirection, margin, extraLayer);
    TriangleBoundary3D<T> boundary(mesh);
    boundary.getMesh().inflate();
    T dx = boundary.getDx();
    Array<T,3> location(boundary.getPhysicalLocation());
    VoxelizedDomain3D<T> voxelizedDomain (
            boundary, typeOfVoxelization, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize );

    MultiScalarField3D<T> flagMatrix((MultiBlock3D&)voxelizedDomain.getVoxelMatrix());
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::inside, flagMatrix.getBoundingBox(), 1.0);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::innerBorder, flagMatrix.getBoundingBox(), 1.0);
    pcout << "Number of inside cells: " << computeSum(flagMatrix) << std::endl;

   //TriangleSet<T> triangles (
           //vofToTriangles<T,descriptors::D3Q19Descriptor>( 
               //flagMatrix, 0.5, flagMatrix.getBoundingBox().enlarge(-2) ) );

    std::vector<TriangleSet<T>::Triangle> triangles;
    isoSurfaceMarchingCube(triangles, voxelizedDomain, flagMatrix.getBoundingBox());
    TriangleSet<T> newTriangleSet(triangles, FLT);
    newTriangleSet.writeBinarySTL(outFileName);


    DEFscaledMesh<T> newMesh(newTriangleSet);
    TriangleBoundary3D<T> newBoundary(newMesh);
    newBoundary.getMesh().scale(dx);
    newBoundary.getMesh().translate(location);
    newBoundary.getMesh().writeBinarySTL(outFileName);

    return 0;
}

