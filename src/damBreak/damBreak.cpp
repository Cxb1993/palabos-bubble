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

#include <cmath>
#include <cstdlib>

using namespace plb;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DESCRIPTOR descriptors::ForcedD3Q27Descriptor

#define PADDING 8

typedef double T;


/*
 * Simulation parameters.
 */

static const T lx = 3.22;
static const T ly = 1.0;
static const T lz = 1.0;

static const T rho_LB = 1.0;
    
static plint maxIter;
static plint outIter;
static plint statIter;

static plint resolution;
static plint nx, ny, nz;
static T dt, dx, u_LB;
static Array<T,3> externalForce;
static T rho;
static T Je, Ga;
static T nu_LB, omega, sigma_LB, contactAngle;
static T cSmago;
static bool useBubblePressureModel, freezeLargestBubble;
static T bubbleVolumeRatio;
static T alpha, beta;

static plint obstacleCenterXYplane, obstacleLength, obstacleWidth, obstacleHeight;
static plint beginWaterReservoir, waterReservoirHeight;

static std::string outDir("./tmp/");


void readUserDefinedSimulationParameters(std::string xmlInputFileName)
{
    XMLreader document(xmlInputFileName);

    document["fluid"]["rho"].read(rho);
    document["fluid"]["Je"].read(Je);
    document["fluid"]["Ga"].read(Ga);
    document["fluid"]["contactAngle"].read(contactAngle);
    T pi = acos((T) -1);
    contactAngle *= pi/180.0;

    document["numerics"]["resolution"].read(resolution);
    document["numerics"]["uLB"].read(u_LB);
    document["numerics"]["maxIter"].read(maxIter);
    document["numerics"]["cSmago"].read(cSmago);
    document["numerics"]["useBubblePressureModel"].read(useBubblePressureModel);
    document["numerics"]["freezeLargestBubble"].read(freezeLargestBubble);
    document["numerics"]["bubbleVolumeRatio"].read(bubbleVolumeRatio);
    document["numerics"]["alpha"].read(alpha);
    document["numerics"]["beta"].read(beta);

    document["output"]["statIter"].read(statIter);
    document["output"]["outIter"].read(outIter);
}

void calculateDerivedSimulationParameters()
{
    dx = lz / (resolution - 1.0);
    T l_LB = resolution - 1.0;

    nx = util::roundToInt(lx / dx);
    ny = util::roundToInt(ly / dx);
    nz = util::roundToInt(lz / dx);

    dt = u_LB * dx;

    // Gravity in lattice units (gravity is 1.0 in dimensionless units).
    T g_LB = dt * dt / dx;
    externalForce  = Array<T,3>((T) 0, (T) 0, -g_LB);

    nu_LB = sqrt(g_LB * l_LB * l_LB * l_LB / Ga);
    omega = 1.0 / (3.0 * nu_LB + 0.5);
    
    sigma_LB = Je * rho_LB * g_LB * l_LB * l_LB;

    obstacleCenterXYplane = util::roundToInt(0.744*l_LB);
    obstacleLength        = util::roundToInt(0.403*l_LB);
    obstacleWidth         = util::roundToInt(0.161*l_LB);
    obstacleHeight        = util::roundToInt(0.161*l_LB);
    beginWaterReservoir   = util::roundToInt((0.744+1.248)*l_LB);
    waterReservoirHeight  = util::roundToInt(0.55*l_LB);
}

void printSimulationParameters()
{
    pcout << "lx = " << lx << std::endl;
    pcout << "ly = " << ly << std::endl;
    pcout << "lz = " << lz << std::endl;

    pcout << "rho = " << rho << std::endl;
    pcout << "Je = " << Je << std::endl;
    pcout << "Ga = " << Ga << std::endl;
    pcout << "contractAngle = " << contactAngle * 180.0 / acos((T) -1) << std::endl;

    pcout << "resolution = " << resolution << std::endl;
    pcout << "uLB = " << u_LB << std::endl;
    pcout << "maxIter = " << maxIter << std::endl;
    pcout << "cSmago = " << cSmago << std::endl;
    pcout << "useBubblePressureModel = " << (useBubblePressureModel ? "true" : "false") << std::endl;
    pcout << "freezeLargestBubble = " << (freezeLargestBubble ? "true" : "false") << std::endl;
    pcout << "bubbleVolumeRatio = " << bubbleVolumeRatio << std::endl;
    pcout << "alpha = " << alpha << std::endl;
    pcout << "beta = " << beta << std::endl;

    pcout << "statIter = " << statIter << std::endl;
    pcout << "outIter = " << outIter << std::endl;

    pcout << "dx = " << dx << std::endl;
    pcout << "dt = " << dt << std::endl;
    pcout << "dt / (dx * dx) = " << dt / (dx * dx) << std::endl;
    pcout << "nx = " << nx << std::endl;
    pcout << "ny = " << ny << std::endl;
    pcout << "nz = " << nz << std::endl;
    pcout << "g_LB = (" << externalForce[0] << ", " << externalForce[1] << ", " << externalForce[2] << ")" << std::endl;
    pcout << "rho_LB = " << rho_LB << std::endl;
    pcout << "nu_LB = " << nu_LB << std::endl;
    pcout << "sigma_LB = " << sigma_LB << std::endl;
    pcout << "omega = " << omega << std::endl;
    pcout << "tau = " << 1.0 / omega << std::endl;
    pcout << std::endl;
}

// Specifies the initial condition for the fluid (each cell is assigned the
// flag "fluid", "empty", or "wall").
int initialFluidFlags(plint iX, plint iY, plint iZ)
{
    // Place an obstacle on the left end, which is hit by the fluid.
    bool insideObstacle =
        iX >= obstacleCenterXYplane-obstacleWidth/2 &&
        iX <= obstacleCenterXYplane+obstacleWidth/2 &&
        iY >= ny/2-obstacleLength/2 &&
        iY <= ny/2+obstacleLength/2 &&
        iZ <= obstacleHeight+1;
    
    if (insideObstacle) {
        return twoPhaseFlag::wall;
    }
    else if (iX >= beginWaterReservoir && iZ <= waterReservoirHeight) {
        return twoPhaseFlag::fluid;
    }
    else {
        return twoPhaseFlag::empty;
    }
}

void writeResults(FreeSurfaceFields3D<T,DESCRIPTOR> *fields, MultiScalarField3D<plint> *tagMatrix, plint iT)
{
    std::auto_ptr<MultiScalarField3D<T> > smoothVF(lbmSmoothen<T,DESCRIPTOR>(fields->volumeFraction,
                fields->volumeFraction.getBoundingBox()));

    std::vector<T> isoLevels;
    isoLevels.push_back(0.5);

    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, smoothVF->getBoundingBox().enlarge(-2));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "smoothedInterface_", iT, PADDING)+".stl");
    }
    triangles.clear();
    isoSurfaceMarchingCube(triangles, fields->volumeFraction, isoLevels, fields->volumeFraction.getBoundingBox().enlarge(-2));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "interface_", iT, PADDING)+".stl");
    }

    //T coef = 1.0 / 3.0;
    VtkImageOutput3D<T> vtkOut(createFileName(outDir + "volumeData_", iT, PADDING), dx);
    std::auto_ptr<MultiTensorField3D<T,3> > v = computeVelocity(fields->lattice);
    std::auto_ptr<MultiScalarField3D<T> > density = computeDensity(fields->lattice);
    vtkOut.writeData<3,float>(*v, "velocity", dx / dt);
    //vtkOut.writeData<float>(*density, "pressure", rho * coef * (dx * dx) / (dt * dt));
    vtkOut.writeData<float>(fields->volumeFraction, "volumeFraction", 1.0);
    vtkOut.writeData<float>(*smoothVF, "smoothedVolumeFraction", 1.0);
    //vtkOut.writeData<float>(fields->outsideDensity, "outsidePressure",
    //        rho * coef * (dx * dx) / (dt * dt));
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

    readUserDefinedSimulationParameters(xmlInputFileName);
    calculateDerivedSimulationParameters();
    printSimulationParameters();

    SparseBlockStructure3D blockStructure(createRegularDistribution3D(nx, ny, nz));
    Dynamics<T,DESCRIPTOR>* dynamics = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega, cSmago);
    FreeSurfaceFields3D<T,DESCRIPTOR> fields(blockStructure, dynamics->clone(), rho_LB, sigma_LB,
            contactAngle, externalForce, false);

    // Initialization.

    // Set all outer-wall cells to "wall" (here, bulk-cells are also set to "wall", but it
    // doesn't matter, because they are overwritten on the next line).
    setToConstant(fields.flag, fields.flag.getBoundingBox(), (int)twoPhaseFlag::wall);
    // In the bulk (all except outer wall layer), initialize the flags as specified by
    // the function "initialFluidFlags".
    setToFunction(fields.flag, fields.flag.getBoundingBox().enlarge(-1), initialFluidFlags);
    
    fields.defaultInitialize();

    BubbleMatch3D *bubbleMatch = 0;
    BubbleHistory3D<T> *bubbleHistory = 0;
    FILE *fp = 0;

    if (useBubblePressureModel) {
        bubbleMatch = new BubbleMatch3D(fields.flag);
        bubbleHistory = new BubbleHistory3D<T>(fields.flag);
        std::string fname = outDir + "bubbles.log";
        fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != 0);
    }

    // Main iteration loop.

    for (plint iT = 0; iT < maxIter; iT++) {
        if (iT % statIter == 0 || iT == maxIter-1) {
            pcout << "At iteration " << iT << ", t = " << iT * dt << std::endl;
            T avE = computeAverageEnergy(fields.lattice);
            avE *= rho * dx * dx / (dt * dt);
            pcout << "Average kinetic energy: " << avE << std::endl;
            plint numIntCells = fields.lattice.getInternalStatistics().getIntSum(0);
            pcout << "Number of interface cells: " << numIntCells << std::endl;
            if (iT != 0) {
                pcout << "Time spent for each iteration: "
                      << global::timer("iteration").getTime() / (T) statIter << std::endl;
                global::timer("iteration").reset();
            }
            pcout << std::endl;
        }

        if (iT % outIter == 0 || iT == maxIter-1) {
            pcout << "Writing results at iteration " << iT << ", t = " << iT * dt << std::endl;
            global::timer("images").restart();
            if (useBubblePressureModel) {
                writeResults(&fields, bubbleMatch->getTagMatrix(), iT);
            } else {
                writeResults(&fields, 0, iT);
            }
            global::timer("images").stop();
            pcout << "Time spent for writing results: " << global::timer("images").getTime() << std::endl;
            pcout << std::endl;
        }

        global::timer("iteration").start();
        if (useBubblePressureModel) {
            T bubbleVolumeRatio = iT == 0 ? 1.0 : bubbleVolumeRatio;
            bubbleMatch->execute(fields.flag, fields.volumeFraction);
            bubbleHistory->transition(*bubbleMatch, iT, bubbleVolumeRatio);
            bubbleHistory->updateBubblePressure(fields.outsideDensity, rho_LB, alpha, beta);
        }

        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();
        global::timer("iteration").stop();

        if (useBubblePressureModel) {
            if (iT == 0 && freezeLargestBubble) {
                bubbleHistory->freezeLargestBubble();
            }
            if (iT % statIter == 0 || iT == maxIter-1) {
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
                    T v = it->second.getVolume() * dx * dx * dx;
                    T r = pow(v * 3.0 / 4.0 / pi, 1.0 / 3.0);
                    Array<T,3> center = it->second.getCenter() * dx;
                    pcout << "Bubble " << it->first
                          << " has volume " << v
                          << ", radius " << r
                          << " and center (" << center[0] << "," << center[1] << "," << center[2] << ")" << std::endl;
                    totalBubbleVolume += v;
                    fprintf(fp, "Time: % .8e, Bubble id: %8ld, Volume: % .8e, Radius: % .8e, Center: (% .8e, % .8e, % .8e)\n",
                            iT * dt, it->first, v, r, center[0], center[1], center[2]);

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
    }

    if (useBubblePressureModel) {
        delete bubbleHistory;
        delete bubbleMatch;
        fclose(fp);
    }

    delete dynamics;

    exit(0);
}

