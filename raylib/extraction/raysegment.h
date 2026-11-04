// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYSEGMENT_H
#define RAYLIB_RAYSEGMENT_H

#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"


namespace ray
{

/// Vertex structure used in shortest path algorithm
struct RAYLIB_EXPORT Vertex
{
  Vertex(const Eigen::Vector3d &pos, const Eigen::Vector3d &start)
    : pos(pos), start(start)
    , parent(-1)
    , root(-1)
    , distance_to_ground(std::numeric_limits<double>::max())
    , distance_to_end(0.0)
    , score(std::numeric_limits<double>::max())
    , visited(false)
  {}
  Eigen::Vector3d pos;
  Eigen::Vector3d start;
  int parent, root;
  double distance_to_ground;
  double distance_to_end;
  double score;
  bool visited;
};

/// Converts a ray cloud to a set of points @c points connected by the shortest path to the ground @c mesh
/// the returned vector of index sets provides the root points for each separated tree
/// @c height_min minimum height that counts as a tree
/// @c max_diameter maximum diameter of a tree trunk
/// @c distance_limit maximum distance between points that can be connected
/// @c gravity_factor controls how far laterally the shortest paths can travel
/// @c point_labels optional tree label per ray in @c cloud. When supplied, one tree is
/// generated per distinct label, rather than segmenting the trees from the point distribution. Negative
/// labels mean unlabelled, these points are attached to whichever labelled tree reaches them
/// @c graph_cache optional file for the shortest-path graph (each point's parent, root and distance to ground):
/// loaded when it exists and matches the point count, otherwise computed and saved there. It is only valid for the
/// same cloud, mesh, labels, @c max_diameter, @c distance_limit and @c gravity_factor, which the caller must key.
std::vector<std::vector<int>> RAYLIB_EXPORT getRootsAndSegment(std::vector<Vertex> &points, const Cloud &cloud, const Mesh &mesh,
                                                               double max_diameter, double distance_limit, double height_min,
                                                               double gravity_factor,
                                                               const std::vector<int> *point_labels = nullptr,
                                                               const std::string &graph_cache = std::string());

}  // namespace ray
#endif  // RAYLIB_RAYSEGMENT_H