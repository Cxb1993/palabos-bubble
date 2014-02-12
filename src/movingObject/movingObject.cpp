/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2014 FlowKit Sarl
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

#include <cstdlib>
#include <vector>
#include <string>

using namespace plb;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DESCRIPTOR descriptors::ForcedD3Q27Descriptor

#define PADDING 8

typedef double T;

std::string outDir("./tmp/");

struct SimulationParameters {
    /*
     * Parameters set by the user.
     */

    T lx, ly, lz;                                   // Geometric parameters.
    T fluidPoolHeight;
    std::string fileName;
    T rMax;

    T rho;                                          // Material parameters.
    T Re;
    T Ga;
    T La;
    T contactAngle;

    plint resolution;                               // Numerical parameters.
    T u_LB;
    plint maxIter, ibIter, startIter;
    T cSmago;
    Array<T,3> rotationAxis, rotationAxisPoint;
    bool useBubblePressureModel;
    bool freezeLargestBubble;
    T bubbleVolumeRatio;
    T alpha, beta;

    plint statIter;                                 // Output parameters.
    plint outIter;

    /*
     * Parameters NOT set by the user.
     */

    T dx;
    T dt;
    plint nx, ny, nz;
    T fluidPoolHeight_LB;
    T rMax_LB;
    Array<T,3> rotationAxisPoint_LB;
    T angularVelocity;
    T angularVelocity_LB;
    Array<T,3> g_LB;
    T rho_LB;
    T nu_LB;
    T sigma_LB;
    T omega;
};

SimulationParameters param;

void readUserDefinedSimulationParameters(std::string xmlInputFileName, SimulationParameters& param)
{
    XMLreader document(xmlInputFileName);

    document["geometry"]["simulationDomain"]["lx"].read(param.lx);
    document["geometry"]["simulationDomain"]["ly"].read(param.ly);
    document["geometry"]["simulationDomain"]["lz"].read(param.lz);
    document["geometry"]["fluidPoolHeight"].read(param.fluidPoolHeight);
    document["geometry"]["fileName"].read(param.fileName);
    document["geometry"]["rMax"].read(param.rMax);

    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["Re"].read(param.Re);
    document["fluid"]["Ga"].read(param.Ga);
    document["fluid"]["La"].read(param.La);
    document["fluid"]["contactAngle"].read(param.contactAngle);
    T pi = acos((T) -1);
    param.contactAngle *= pi/180.0;

    document["numerics"]["resolution"].read(param.resolution);
    document["numerics"]["uLB"].read(param.u_LB);
    document["numerics"]["maxIter"].read(param.maxIter);
    document["numerics"]["ibIter"].read(param.ibIter);
    document["numerics"]["startIter"].read(param.startIter);
    document["numerics"]["cSmago"].read(param.cSmago);
    document["numerics"]["rotationAxis"].read<T,3>(param.rotationAxis);
    param.rotationAxis /= norm<T,3>(param.rotationAxis);
    document["numerics"]["rotationAxisPoint"].read<T,3>(param.rotationAxisPoint);
    document["numerics"]["useBubblePressureModel"].read(param.useBubblePressureModel);
    document["numerics"]["freezeLargestBubble"].read(param.freezeLargestBubble);
    document["numerics"]["bubbleVolumeRatio"].read(param.bubbleVolumeRatio);
    document["numerics"]["alpha"].read(param.alpha);
    document["numerics"]["beta"].read(param.beta);

    document["output"]["statIter"].read(param.statIter);
    document["output"]["outIter"].read(param.outIter);
}

void calculateDerivedSimulationParameters(SimulationParameters& param)
{
    // Derived quantities.

    param.dx = 1.0 / (param.resolution - 1.0);
    param.dt = param.u_LB * param.dx;

    T l_LB = param.resolution - 1.0;

    param.nx = util::roundToInt(param.lx / param.dx);
    param.ny = util::roundToInt(param.ly / param.dx);
    param.nz = util::roundToInt(param.lz / param.dx);

    param.fluidPoolHeight_LB = param.fluidPoolHeight / param.dx;
    param.rMax_LB = param.rMax / param.dx;
    param.rotationAxisPoint_LB = (1.0 / param.dx) * param.rotationAxisPoint;
    param.angularVelocity = 1.0 / param.rMax;
    param.angularVelocity_LB = param.angularVelocity * param.dt;

    param.rho_LB = 1.0;
    param.nu_LB = param.u_LB * l_LB / param.Re;
    T norm_g_LB = param.Ga * param.nu_LB * param.nu_LB / (l_LB * l_LB * l_LB);
    param.g_LB = Array<T,3>((T) 0, (T) 0, -norm_g_LB);
    T mu_LB = param.nu_LB * param.rho_LB;
    param.sigma_LB = param.La * mu_LB * mu_LB / (param.rho_LB * l_LB);
    param.omega = 1.0 / (3.0 * param.nu_LB + 0.5);
}

void printSimulationParameters(SimulationParameters const& param)
{
    pcout << "lx = " << param.lx << std::endl;
    pcout << "ly = " << param.ly << std::endl;
    pcout << "lz = " << param.lz << std::endl;
    pcout << "fluidPoolHeight = " << param.fluidPoolHeight << std::endl;
    pcout << "fileName = " << param.fileName << std::endl;
    pcout << "rMax = " << param.rMax << std::endl;

    pcout << "rho = " << param.rho << std::endl;
    pcout << "Re = " << param.Re << std::endl;
    pcout << "Ga = " << param.Ga << std::endl;
    pcout << "La = " << param.La << std::endl;
    pcout << "contractAngle = " << param.contactAngle * 180.0 / acos((T) -1) << std::endl;

    pcout << "resolution = " << param.resolution << std::endl;
    pcout << "uLB = " << param.u_LB << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;
    pcout << "ibIter = " << param.ibIter << std::endl;
    pcout << "startIter = " << param.startIter << std::endl;
    pcout << "cSmago = " << param.cSmago << std::endl;
    pcout << "rotationAxis = (" << param.rotationAxis[0] << ", "
                                << param.rotationAxis[1] << ", "
                                << param.rotationAxis[2] << ")" << std::endl;
    pcout << "rotationAxisPoint = (" << param.rotationAxisPoint[0] << ", "
                                     << param.rotationAxisPoint[1] << ", "
                                     << param.rotationAxisPoint[2] << ")" << std::endl;
    pcout << "useBubblePressureModel = " << (param.useBubblePressureModel ? "true" : "false") << std::endl;
    pcout << "freezeLargestBubble = " << (param.freezeLargestBubble ? "true" : "false") << std::endl;
    pcout << "bubbleVolumeRatio = " << param.bubbleVolumeRatio << std::endl;
    pcout << "alpha = " << param.alpha << std::endl;
    pcout << "beta = " << param.beta << std::endl;

    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;

    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "dt / (dx * dx) = " << param.dt / (param.dx * param.dx) << std::endl;
    pcout << "nx = " << param.nx << std::endl;
    pcout << "ny = " << param.ny << std::endl;
    pcout << "nz = " << param.nz << std::endl;
    pcout << "fluidPoolHeight_LB = " << param.fluidPoolHeight_LB << std::endl;
    pcout << "rMax_LB = " << param.rMax_LB << std::endl;
    pcout << "rotationAxisPoint_LB = (" << param.rotationAxisPoint_LB[0] << ", "
                                        << param.rotationAxisPoint_LB[1] << ", "
                                        << param.rotationAxisPoint_LB[2] << ")" << std::endl;
    pcout << "angularVelocity = " << param.angularVelocity << std::endl;
    pcout << "angularVelocity_LB = " << param.angularVelocity_LB << std::endl;
    pcout << "g_LB = (" << param.g_LB[0] << ", " << param.g_LB[1] << ", " << param.g_LB[2] << ")" << std::endl;
    pcout << "rho_LB = " << param.rho_LB << std::endl;
    pcout << "nu_LB = " << param.nu_LB << std::endl;
    pcout << "sigma_LB = " << param.sigma_LB << std::endl;
    pcout << "omega = " << param.omega << std::endl;
    pcout << "tau = " << 1.0 / param.omega << std::endl;
    pcout << std::endl;
}

bool insideFluid(T x, T y, T z)
{
    if (z <= param.fluidPoolHeight_LB) {
        return true;
    }
    return false;
}

// Specifies the initial condition for the fluid (each cell is assigned the
// flag "fluid", "empty", or "wall").
int initialFluidFlags(plint iX, plint iY, plint iZ)
{
    if (insideFluid(iX, iY, iZ)) {
        return twoPhaseFlag::fluid;
    }
    return twoPhaseFlag::empty;
}

T getAngularVelocity(T targetValue, plint iIter, plint startIter)
{
    static T pi = acos(-1.0);

    if (iIter >= startIter) {
        return targetValue;
    }

    if (iIter < 0) {
        iIter = 0;
    }

    return (targetValue * sin(pi * iIter / (2.0 * startIter)));
}

Array<T,3> getVelocity(Array<T,3> const& pos, Array<T,3> const& rotationAxis, Array<T,3> const& rotationAxisPoint,
        T angularVelocity)
{
    Array<T,3> omega = angularVelocity * rotationAxis;
    Array<T,3> p = pos - rotationAxisPoint;
    Array<T,3> r = p - dot<T,3>(p, rotationAxis) * rotationAxis;
    return crossProduct<T>(omega, r);
}

class RotationalVelocity {
public:
    RotationalVelocity(Array<T,3> const& rotationAxis_, Array<T,3> const& rotationAxisPoint_, T angularVelocity_)
        : rotationAxis(rotationAxis_),
          rotationAxisPoint(rotationAxisPoint_),
          angularVelocity(angularVelocity_)
    { }
    Array<T,3> operator()(Array<T,3> const& pos)
    {
        return getVelocity(pos, rotationAxis, rotationAxisPoint, angularVelocity);
    }
private:
    Array<T,3> rotationAxis, rotationAxisPoint;
    T angularVelocity;
};

void writeResults(FreeSurfaceFields3D<T,DESCRIPTOR> *fields, MultiScalarField3D<plint> *tagMatrix, plint iT)
{
    std::auto_ptr<MultiScalarField3D<T> > smoothVF(lbmSmoothen<T,DESCRIPTOR>(fields->volumeFraction,
                fields->volumeFraction.getBoundingBox()));

    std::vector<T> isoLevels;
    isoLevels.push_back(0.5);

    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, smoothVF->getBoundingBox().enlarge(-1));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(param.dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "smoothedInterface_", iT, PADDING)+".stl");
    }
    triangles.clear();
    isoSurfaceMarchingCube(triangles, fields->volumeFraction, isoLevels, fields->volumeFraction.getBoundingBox().enlarge(-1));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(param.dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "interface_", iT, PADDING)+".stl");
    }

    //T coef = 1.0 / 3.0;
    VtkImageOutput3D<T> vtkOut(createFileName(outDir + "volumeData_", iT, PADDING), param.dx);
    std::auto_ptr<MultiTensorField3D<T,3> > v = computeVelocity(fields->lattice);
    std::auto_ptr<MultiScalarField3D<T> > rho = computeDensity(fields->lattice);
    vtkOut.writeData<3,float>(*v, "velocity", param.dx / param.dt);
    //vtkOut.writeData<float>(*rho, "pressure", param.rho * coef * (param.dx * param.dx) / (param.dt * param.dt));
    vtkOut.writeData<float>(fields->volumeFraction, "volumeFraction", 1.0);
    vtkOut.writeData<float>(*smoothVF, "smoothedVolumeFraction", 1.0);
    //vtkOut.writeData<float>(fields->outsideDensity, "outsidePressure",
    //        param.rho * coef * (param.dx * param.dx) / (param.dt * param.dt));
    if (tagMatrix != 0) {
        vtkOut.writeData<float>(*copyConvert<plint,double>(*tagMatrix), "bubbleTags", 1.0);
    }
}


int main(int argc, char **argv)
{
    plbInit(&argc, &argv);

    std::cout.precision(10);
    std::scientific(std::cout);

    // Command-line arguments

    if (argc != 2) {
        pcout << "Usage: " << argv[0] << " xml-input-file-name" << std::endl;
        exit(1);
    }

    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);

    // Set the simulation parameters.

    readUserDefinedSimulationParameters(xmlInputFileName, param);
    calculateDerivedSimulationParameters(param);
    printSimulationParameters(param);

    SparseBlockStructure3D blockStructure(createRegularDistribution3D(param.nx, param.ny, param.nz));

    Dynamics<T,DESCRIPTOR>* dynamics = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);

    bool useImmersedWalls = true;
    FreeSurfaceFields3D<T,DESCRIPTOR> fields(blockStructure, dynamics->clone(), param.rho_LB,
            param.sigma_LB, param.contactAngle, param.g_LB, useImmersedWalls);

    // Initialization

    pcout << "Setting up initial condition." << std::endl;

    analyticalIniVolumeFraction(fields.volumeFraction, fields.flag, insideFluid, 32);

    Box3D bottom(0, param.nx-1, 0, param.ny-1, 0, 0);
    Box3D top(0, param.nx-1, 0, param.ny-1, param.nz-1, param.nz-1);
    Box3D lateral1(0, 0, 0, param.ny-1, 0, param.nz-1);
    Box3D lateral2(param.nx-1, param.nx-1, 0, param.ny-1, 0, param.nz-1);
    Box3D lateral3(0, param.nx-1, 0, 0, 0, param.nz-1);
    Box3D lateral4(0, param.nx-1, param.ny-1, param.ny-1, 0, param.nz-1);

    setToConstant(fields.flag, bottom, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, top, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral1, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral2, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral3, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral4, (int) twoPhaseFlag::wall);

    fields.periodicityToggleAll(false);
    fields.partiallyDefaultInitialize();

    // Moving Object.

    TriangleSet<T> surfaceTriangleSet(param.fileName, DBL);
    plint maxRefinements = 64;
    bool succeeded = surfaceTriangleSet.refineRecursively(param.dx, maxRefinements);
    if (!succeeded) {
        pcout << std::endl;
        pcout << "WARNING: The target maximum triangle edge length " << param.dx
              << " for the immersed surface was not reached after " << maxRefinements
              << " refinement iterations." << std::endl;
        pcout << std::endl;
    }
    FileName newSurfaceFileName(param.fileName);
    newSurfaceFileName.setName(outDir + newSurfaceFileName.getName() + "_refined");
    surfaceTriangleSet.writeBinarySTL(newSurfaceFileName.get());

    surfaceTriangleSet.scale(1.0 / param.dx);

    ConnectedTriangleSet<T> *connectedTriangleSet = new ConnectedTriangleSet<T>(surfaceTriangleSet);
    surfaceTriangleSet.scale(param.dx);
    pcout << "The refined immersed surface has " << connectedTriangleSet->getNumVertices()
          << " vertices and " << connectedTriangleSet->getNumTriangles() << " triangles." << std::endl;

    plint numVertices = connectedTriangleSet->getNumVertices();
    std::vector<Array<T,3> > vertices(numVertices);
    std::vector<T> areas(numVertices);
    for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
        vertices[iVertex] = connectedTriangleSet->getVertex(iVertex);
        T area;
        Array<T,3> normal;
        connectedTriangleSet->computeVertexAreaAndUnitNormal(iVertex, area, normal);
        areas[iVertex] = area;
    }
    delete connectedTriangleSet;

    MultiContainerBlock3D container(fields.rhoBar);

    plint iniIter = 0;

    BubbleMatch3D *bubbleMatch = 0;
    BubbleHistory3D<T> *bubbleHistory = 0;
    FILE *fp = 0;

    if (param.useBubblePressureModel) {
        bubbleMatch = new BubbleMatch3D(fields.flag);
        bubbleHistory = new BubbleHistory3D<T>(fields.flag);
        std::string fname = outDir + "bubbles.log";
        fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != 0);
    }

    pcout << std::endl;

    // Main iteration loop.

    T angle = 0.0; // Rotation angle.
    T previousAngle = 0.0; // Rotation angle for output purposes.
    for (plint iT = iniIter; iT < param.maxIter; iT++) {
        if (iT % param.statIter == 0 || iT == param.maxIter-1) {
            pcout << "At iteration " << iT << ", t = " << iT * param.dt << std::endl;
            T avE = computeAverageEnergy(fields.lattice);
            avE *= param.rho * param.dx * param.dx / (param.dt * param.dt);
            pcout << "Average kinetic energy: " << avE << std::endl;
            plint numIntCells = fields.lattice.getInternalStatistics().getIntSum(0);
            pcout << "Number of interface cells: " << numIntCells << std::endl;
            if (iT != iniIter) {
                pcout << "Time spent for each iteration: "
                      << global::timer("iteration").getTime() / (T) param.statIter << std::endl;
                global::timer("iteration").reset();
            }
            pcout << std::endl;
        }

        if (iT % param.outIter == 0 || iT == param.maxIter-1) {
            pcout << "Writing results at iteration " << iT << ", t = " << iT * param.dt << std::endl;
            global::timer("images").restart();
            if (param.useBubblePressureModel) {
                writeResults(&fields, bubbleMatch->getTagMatrix(), iT);
            } else {
                writeResults(&fields, 0, iT);
            }
            surfaceTriangleSet.writeBinarySTL(outDir + createFileName("object_", iT, PADDING) + ".stl");
            global::timer("images").stop();
            pcout << "Time spent for writing results: " << global::timer("images").getTime() << std::endl;
            pcout << std::endl;
        }

        global::timer("iteration").start();
        if (param.useBubblePressureModel) {
            T bubbleVolumeRatio = iT == 0 ? 1.0 : param.bubbleVolumeRatio;
            bubbleMatch->execute(fields.flag, fields.volumeFraction);
            bubbleHistory->transition(*bubbleMatch, iT, bubbleVolumeRatio);
            bubbleHistory->updateBubblePressure(fields.outsideDensity, param.rho_LB, param.alpha, param.beta);
        }

        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();
        global::timer("iteration").stop();

        if (param.useBubblePressureModel) {
            if (iT == 0 && param.freezeLargestBubble) {
                bubbleHistory->freezeLargestBubble();
            }
            if (iT % param.statIter == 0 || iT == param.maxIter-1) {
                bubbleHistory->timeHistoryLog(outDir + "bubbleTimeHistory.log");
                bubbleHistory->fullBubbleLog(outDir + "fullBubbleRecord.log");

                // We do not log frozen bubbles.
                std::map<plint,BubbleInfo3D>::const_iterator it = bubbleHistory->getBubbles().begin();
                T totalBubbleVolume = T();
                plint numBubbles = 0;
                TriangleSet<T> bubbles;
                T pi = acos((T) -1);
                for (; it != bubbleHistory->getBubbles().end(); ++it) {
                    if (it->second.isFrozen()) {
                        continue;
                    }
                    numBubbles++;
                    T v = it->second.getVolume() * param.dx * param.dx * param.dx;
                    T r = pow(v * 3.0 / 4.0 / pi, 1.0 / 3.0);
                    Array<T,3> center = it->second.getCenter() * param.dx;
                    pcout << "Bubble " << it->first
                          << " has volume " << v
                          << ", radius " << r
                          << " and center (" << center[0] << "," << center[1] << "," << center[2] << ")" << std::endl;
                    totalBubbleVolume += v;
                    fprintf(fp, "Time: % .8e, Bubble id: %8ld, Volume: % .8e, Radius: % .8e, Center: (% .8e, % .8e, % .8e)\n",
                            iT * param.dt, it->first, v, r, center[0], center[1], center[2]);

                    // Save STL files of bubbles.

                    if (!util::fpequal(r, (T) 0.0, std::numeric_limits<T>::epsilon())) {
                        plint minNumOfTriangles = 64;
                        TriangleSet<T> triangleSet = constructSphere(center, r, minNumOfTriangles);
                        bubbles.append(triangleSet);
                    }
                }
                std::string fname = createFileName(outDir + "bubbles_", iT, PADDING) + ".stl";
                bubbles.writeBinarySTL(fname);
                pcout << "At iteration " << iT << ", the number of bubbles is " << numBubbles << std::endl;
                pcout << "The total volume of bubbles is: " << totalBubbleVolume << std::endl;
                pcout << std::endl;
                fflush(fp);
            }
        }

        // Immersed Boundary Method.

        plint nextIt = iT + 1;

        T nextAngularVelocity_LB = getAngularVelocity(param.angularVelocity_LB, nextIt, param.startIter);
        T currentAngularVelocity_LB = getAngularVelocity(param.angularVelocity_LB, iT, param.startIter);

        T meanAngularVelocity_LB = 0.5 * (nextAngularVelocity_LB + currentAngularVelocity_LB);

        for (plint iVertex = 0; iVertex < (plint) vertices.size(); iVertex++) {
            vertices[iVertex] += getVelocity(vertices[iVertex], param.rotationAxis, param.rotationAxisPoint_LB,
                    meanAngularVelocity_LB); // Trapezoidal rule with dt_LB = 1.
        }
        angle += atan(meanAngularVelocity_LB);
        if (nextIt % param.outIter == 0 || nextIt == param.maxIter-1) {
            // Rotate the geometry in physical units.
            T theta = angle - previousAngle;
            surfaceTriangleSet.translate(-param.rotationAxisPoint);
            surfaceTriangleSet.rotateAtOrigin(param.rotationAxis, theta);
            surfaceTriangleSet.translate(param.rotationAxisPoint);
            previousAngle = angle;
        }

        instantiateImmersedWallData(vertices, areas, container);
        for (plint i = 0; i < param.ibIter; i++) {
            inamuroIteration(
                    RotationalVelocity(param.rotationAxis, param.rotationAxisPoint_LB, nextAngularVelocity_LB),
                    fields.rhoBar, fields.j, container, 1.0 / param.omega);
        }
    }

    if (param.useBubblePressureModel) {
        delete bubbleHistory;
        delete bubbleMatch;
        fclose(fp);
    }

    delete dynamics;

    exit(0);
}

