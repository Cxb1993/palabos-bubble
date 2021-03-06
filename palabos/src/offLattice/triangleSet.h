/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
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

/* Main author: Dimitrios Kontaxakis */

#ifndef TRIANGLE_SET_H
#define TRIANGLE_SET_H

#include "core/globalDefs.h"
#include "core/geometry3D.h"
#include "core/array.h"
#include <string>
#include <vector>
#include <cstdio>

namespace plb {

template<typename T>
class TriangleSet {
public:
    typedef Array<Array<T,3>,3> Triangle;
public:
    TriangleSet(Precision precision_ = FLT);
    TriangleSet(std::vector<Triangle> const& triangles_, Precision precision_ = FLT);
    // Currently STL and OFF files are supported by this class.
    TriangleSet(std::string fname, Precision precision_ = FLT, SurfaceGeometryFileFormat fformat = STL);
    std::vector<Triangle> const& getTriangles() const;
    Precision getPrecision() const { return precision; }
    void setPrecision(Precision precision_);

    /// Translate the triangle set surface mesh.
    void translate(Array<T,3> const& vector);
    /// Scale the triangle set surface mesh.
    void scale(T alpha);
    /// Rotate the triangle set surface mesh.
    ///   The arguments of this function are the Euler angles in radians.
    ///   The so-called "x-convention" is used, in which the rotation is
    ///   given by the three angles (phi, theta, psi), where:
    ///   1.] The first rotation is by an angle phi about the z-axis,
    ///   2.] The second rotation is by an angle theta in [0, pi] about
    ///       the new x-axis,
    ///   3.] The third rotation is by an angle psi about the new z-axis.
    void rotate(T phi, T theta, T psi);
    void rotateAtOrigin(Array<T,3> const& normedAxis, T theta);

    /// Clear the triangle mesh.
    void clear();

    /// Erase the current triangle set surface mesh, and merge into it the new meshes.
    ///   This function currently does not check for duplicate
    ///   triangles in the new meshes, and does not handle
    ///   cases when the resulting merged surface mesh is
    ///   non-manifold. It is in the user's jurisdiction to
    ///   make sure that the resulting mesh is well defined.
    ///   Practically, this can be achieved if the new triangle set
    ///   meshes given as arguments are mutually disjoint.
    void merge(std::vector<TriangleSet<T>*> meshes);

    /// Append to the current triangle mesh, the mesh that is passed as an argument.
    ///   This function currently does not check for duplicate
    ///   triangles in the two meshes, and does not handle
    ///   cases when the resulting merged surface mesh is
    ///   non-manifold. It is in the user's jurisdiction to
    ///   make sure that the resulting mesh is well defined.
    ///   Practically, this can be achieved if the two triangle set
    ///   meshes are mutually disjoint.
    void append(TriangleSet<T> const& mesh);

    /// Refine the current triangle set surface mesh by splitting each
    ///   original triangle into four smaller triangles constructed by
    ///   joining the middle points of the original edges. The old mesh
    ///   is deleted.
    void refine();
    /// Refine the current triangle set surface mesh by considering
    ///   edge lengths. If the length of an edge is greater than or
    ///   equal to the edgeLengthThreshold, then this edge is splitted
    ///   in two. The old mesh is deleted.
    void refine(T edgeLengthThreshold);
    /// Refine the current triangle set surface mesh by considering
    ///   edge lengths. If the length of an edge is greater than or
    ///   equal to the targetMaxEdgeLength, then this edge is splitted
    ///   in two. The procedure is repeated until the maximum edge length
    ///   is less than the one provided, or until the maxNumIterations
    ///   has been reached. The old mesh is deleted. The function returns
    ///   true upon success or false otherwise.
    bool refineRecursively(T targetMaxEdgeLength, plint maxNumIterations);

    /// A very simple orientation reversing function.
    void reverseOrientation();

    /// Cast all coordinates to float and then back to type T.
    void toFloatAndBack();

    /// Check if the triangle set has triangles with zero area.
    bool hasDegenerateTriangles() const;

    /// Export the mesh as an ASCII STL file.
    void writeAsciiSTL(std::string fname) const;
    /// Export the mesh as an binary STL file.
    void writeBinarySTL(std::string fname) const;

    /// Cut the current triangle set mesh by a plane "plane" which is
    ///   defined by a point and a normal. This cutting operation will
    ///   produce a new triangle set "newTriangleSet" which is
    ///   the one that belongs to the half-space that the normal of the
    ///   cutting plane points outside of. The function returns 1 if
    ///   the cut was successful, 0 if there was no intersection between
    ///   the original triangle set and the plane provided, and -1 if an
    ///   error occured.
    int cutWithPlane(Plane<T> const& plane, TriangleSet<T>& newTriangleSet) const;

    /// Cut the current triangle set mesh by a plane "plane" which is
    ///   defined by a point and a normal and a cuboid "cuboid" which is
    ///   defined by a lower left corner and a right upper corner. This
    ///   cutting operation will cut off all triangles, or parts of triangles
    ///   that are fully contained in the cuboid and are positioned in
    ///   the half-space that the normal of the cutting plane points into.
    ///   Caution is needed, because this function does not use the cuboid
    ///   to cut, but only to select parts of the original triangle set,
    ///   to be then cut by the plane. Obviously, at least a part of the
    ///   cutting plane must be contained in the cuboid, for the cutting
    ///   to make any sense. If not used wisely, this function can lead to
    ///   broken STL files. The function returns 1 if the cut was successful,
    ///   0 if there was no intersection between the original triangle set
    ///   the plane and the cuboid provided, and -1 if an error occured.
    int cutWithPlane(Plane<T> const& plane, Cuboid<T> const& cuboid,
            TriangleSet<T>& newTriangleSet) const;

    T getMinEdgeLength() const { return minEdgeLength; }
    T getMaxEdgeLength() const { return maxEdgeLength; }

    Cuboid<T> getBoundingCuboid() const { return boundingCuboid; }
    Array<T,3> getCentroid() const;

private:
    void readSTL(std::string fname);
    void readAsciiSTL(FILE* fp);
    void readBinarySTL(FILE* fp);
    void readOFF(std::string fname);
    void readAsciiOFF(FILE* fp);
    void checkForDegenerateTriangles(Triangle const& triangle, Array<T,3>& computedNormal) const;
    bool checkForDegenerateTrianglesNoAbort(Triangle const& triangle, Array<T,3>& computedNormal) const;
    void checkForDegenerateTrianglesAndFixOrientation(Triangle& triangle, Array<T,3> const& n) const;
    bool checkForDegenerateTrianglesAndFixOrientationNoAbort(Triangle& triangle, Array<T,3> const& n) const;
    void computeMinMaxEdge(pluint iTriangle, T& minEdge, T& maxEdge) const;
    void computeMinMaxEdges();
    void computeBoundingCuboid();
    /// Cut the current triangle by a plane "plane" which is
    ///   defined by a point and a normal. This cutting operation will
    ///   add or not one or more triangles to the triangle set "newTriangleSet".
    ///   The triangles added belong to the half-space that the normal of
    ///   the plane given points outside of. The function returns 1 if the
    ///   cut was successful, 0 if there was no intersection between the
    ///   original triangle set and the plane and -1 if an error occured.
    int cutTriangleWithPlane(Plane<T> const& plane, Triangle const& triangle,
            TriangleSet<T>& newTriangleSet) const;
    /// Skip nLines number of lines in the file.
    void skipLines(plint nLines, FILE* fp) const;
    /// Read the file character-by-character. This function returns 0 if a non-white-space
    /// character is found, or EOF if the end-of-file is found. It "consumes" the
    /// rest of the line if the "commentCharacter" is found. Obviously, it works only
    /// with text files.
    int readAhead(FILE *fp, char commentCharacter) const;
    /// Check if the buffer which holds an input line of text is full.
    bool checkForBufferOverflow(char* buf) const;

private:
    std::vector<Triangle> triangles;
    T minEdgeLength, maxEdgeLength;
    Cuboid<T> boundingCuboid;
    /// Enumeration constant that sets the precision for triangle set checks.
    /// Single precision (float): FLT
    /// Double precision (double): DBL
    /// Extended precision (long double): LDBL
    /// "Infinite" precision: INF
    Precision precision;
};

} // namespace plb

#endif  // TRIANGLE_SET_H
