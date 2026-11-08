// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_CLUSTERS_H
#define RAYLIB_CLUSTERS_H

#include "../rayutils.h"
#include "raylib/raylibconfig.h"

namespace ray
{
/// An optional record of one agglomeration (rayextract trees --export_sections): every neighbour link between two
/// different clusters, in the order tried, with the box size (max side) the merge would give and whether it merged
struct RAYLIB_EXPORT ClusterLink
{
  int a, b;       // the linked points
  int na, nb;     // the sizes of their clusters when tried
  double dist, box;
  bool merged;    // box < max_diameter
};

/// Generate clusters from a set of points based on min and max diameter criteria
/// This is an agglomerative clustering method, whereby clusters iteratively grow and fuse
void RAYLIB_EXPORT generateClusters(std::vector<std::vector<int>> &point_clusters, const std::vector<Eigen::Vector3d> &points,
                                    double min_diameter, double max_diameter,
                                    std::vector<ClusterLink> *links = nullptr);
}  // namespace ray
#endif  // RAYLIB_RAYEXTRACT_WOODS_H