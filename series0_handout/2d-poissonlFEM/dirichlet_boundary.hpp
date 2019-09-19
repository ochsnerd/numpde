#pragma once
#include <Eigen/Core>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/setdiff.h>
#include <igl/colon.h>


//! Finds the boundary and interior vertices, and evaluates 
//! the function g on the boundary vertices.
//!
//! @param[out] u should be of size <number of vertices>. At the end:
//!             u(index) = g(vertices(index))  for each boundary index.
//!
//! @param[out] interiorVertexIndices the list of vertices that are not on 
//!                                   the boundary.
//!
//! @param[in] vertices a list of vertices for the mesh
//!
//! @param[in] triangles the elements represented as triangles
//!
//! @param[in] g the boundary value function g.
void setDirichletBoundary(Eigen::VectorXd& u,
                          Eigen::VectorXi& interiorVertexIndices,
                          const Eigen::MatrixXd& vertices,
                          const Eigen::MatrixXi& triangles,
                          const std::function<double(double, double)>& g)
{
    // Find boundary edges
    Eigen::MatrixXi boundaryEdgeIndices;
    igl::boundary_facets(triangles, boundaryEdgeIndices);

    // Find boundary vertices
    Eigen::VectorXi boundaryVertexIndices, IA, IC;
    igl::unique(boundaryEdgeIndices, boundaryVertexIndices, IA, IC);

    for ( int i = 0; i < boundaryVertexIndices.size(); ++i) {
        const int index = boundaryVertexIndices(i);
        const auto& x = vertices.row(index);
        u (index) = g(x(0), x(1));
    }

    Eigen::VectorXi allIndices;
    igl::colon<int>(0, vertices.rows() - 1, allIndices);
    igl::setdiff(allIndices, boundaryVertexIndices, interiorVertexIndices, IA);
}
