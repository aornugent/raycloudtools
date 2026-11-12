// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#ifndef RAYLIB_RAYEXTRACT_TREES_H
#define RAYLIB_RAYEXTRACT_TREES_H

#include "../raycloud.h"
#include "../raymesh.h"
#include "../rayutils.h"
#include "raylib/raylibconfig.h"
#include "raysegment.h"

namespace ray
{
/// structure containing the parameters used in tree reconstruction
struct RAYLIB_EXPORT TreesParams
{
  TreesParams();
  double max_diameter;              // maximum tree diameter. Trees wider than this may be segmented into multiple trees
  double crop_length;               // distance to end of branch where it crops the branch, and doesn't generate further geometry
  double distance_limit;            // maximum distance between points that can count as connected
  double height_min;                // minimum height for a tree. Lower values are considered undergrowth and excluded
  double girth_height_ratio;        // how far up tree to measure girth
  double cylinder_length_to_width;  // the slenderness of the branch segment cylinders
  double gap_ratio;                 // max gap per branch length
  double span_ratio;                // points that span a larger width determine that a branch has become two
  double gravity_factor;   // preferences branches that are less lateral, so penalises implausable horizontal branches
  double grid_width;       // used on a grid cell with overlap, to remove trees with a base in the overlap zone
  bool segment_branches;   // flag to output the ray cloud coloured by branch segment index rather than by tree index
  double global_taper;     // forced global taper, uses global_taper_factor to define how much it is applied
  double global_taper_factor; // 0 estimates per-tree tapering, 1 uses per-scan tapering, 0.5 is mid-way on mid-weight trees
  double taper_weight;  // if > 0, the weight of the global taper prior, instead of the run's running sum(w^2)/sum(w)
  double taper_prior_sd;  // if > 0, a per-tree conjugate taper blend in log space: prior log T ~ N(log T0(H), sd^2)
  double taper_exponent;  // the prior mean T0(H) = global_taper * (H / 10 m)^taper_exponent
  double taper_vote_cv;   // coefficient of variation of one taper vote from a well-measured section
  double radius_floor;    // radius measurement floor delta (m): votes below 2 delta are dropped, others gain (delta / r)^2
  double branch_vote_sd;  // extra log sd of a branch vote, for the share-rule assumption it rests on
  double stem_taper;      // radius also scales by (d / d_root)^stem_taper along each path, d = distance to the tips
  double trunk_gap;       // if > 0, the trunk split's gap (m) instead of gap_ratio x length, so stems split apart
  double trunk_span_ratio; // if > 0, the trunk split's span ratio instead of span_ratio (branches keep theirs)
  int trunk_splits;       // how many times one root group may split into trees (stems)
  bool use_rays; // use the full rays in order to estimate a smaller radius when points are not all on the real branch
  std::string graph_cache; // if set, load the shortest-path graph from this file, or compute it and save it there
  bool geometric_centre; // after the algebraic centre step, iterate a robust geometric circle fit with a radius prior
  bool slab_sections; // fit sections made by a split on their share of the slab of points, not on the path chains
  bool export_sections; // record what the fit measured and decided, for saveSections (changes no output)
};

struct BranchSection;  // forwards declaration

/// The class for a set of trees, stored as a list of (connected) branch sections
/// together with the function for their extrsction from a ray cloud
/// This is the class used in the tool rayextrsact trees, which reconstructs the branch structures
/// of the specified ray cloud, and saves them as a text file
class RAYLIB_EXPORT Trees
{
public:
  /// Constructs the piecewise cylindrical tree structures from the input ray cloud @c cloud
  /// The ground @c mesh defines the ground and @params are used to control the reconstruction.
  /// When @c point_labels is supplied (one tree label per ray in @c cloud, -1 for unlabelled) the
  /// cloud is treated as pre-segmented, so only the tree models are estimated, one per label
  Trees(Cloud &cloud, const Eigen::Vector3d &offset, const Mesh &mesh, const TreesParams &params, bool verbose,
        const std::vector<int> *point_labels = nullptr);

  /// save the trees representation to a text file
  bool save(const std::string &filename, const Eigen::Vector3d &offset, bool verbose) const;
  /// write one line per tree (in _trees.txt order) with its taper blend and measured trunk radius
  bool saveTaper(const std::string &filename, const Eigen::Vector3d &offset) const;
  /// with params.export_sections: write <stub>_sections.txt, _sections.bin and _modules.txt
  bool saveSections(const std::string &stub, const Eigen::Vector3d &offset) const;

private:
  /// a tip's module (export_sections): the points calculateSectionIds assigns past its ends, summarised
  struct Module
  {
    int n{0}, n_ends{0};  // the points, and of them the tip's own end points
    Eigen::Vector3d mean{0, 0, 0};
    Eigen::Matrix3d cov{Eigen::Matrix3d::Zero()};  // population covariance of their positions
    double top{0}, reach{0};  // the highest z; graph distance from the owned ends to the farthest tip
  };
  /// summarise the module of each tip from the final ownership and each tip's points past its ends
  void summariseModules(const std::vector<int> &section_ids, const std::vector<std::vector<int>> &beyond);
  std::vector<Module> modules_;        // per section, set for tips (export_sections)
  std::vector<int> ray_section_ids_;   // calculateSectionIds' ownership per ray of the cloud (contiguous ids), -1 for none

  /// The piecewise cylindrical represenation of all of the trees
  std::vector<BranchSection> sections_;

  /// calculate the distance to farthest connected branch tip, for each point in the cloud
  void calculatePointDistancesToEnd();
  /// create the start branch segments at the root positions
  void generateRootSections(const std::vector<std::vector<int>> &roots_list);
  /// finalise the attributes of an end (tip) of a branch
  void setBranchTip();
  /// get the root position for the current section
  Eigen::Vector3d getRootPosition() const;
  /// find the points and end points within this branch section
  void extractNodesAndEndsFromRoots(std::vector<int> &nodes, const Eigen::Vector3d &base,
                                    const std::vector<std::vector<int>> &children, double min_dist, double max_dist);
  /// find separate clusters of points within the branch section
  std::vector<std::vector<int>> findPointClusters(const Eigen::Vector3d &base, bool &points_removed,
                                                  double thickness, double span, double gap);
  bool reachesGround(const std::vector<int> &ends, const std::vector<int> &wood) const;
  /// split the branch section to one branch for each cluster
  void bifurcate(const std::vector<std::vector<int>> &clusters, double thickness, std::vector<std::vector<int>> &children, bool clip_tree, bool add_offshoots,
                 const std::vector<int> *slab = nullptr, double prior_radius = 0.0);
  /// share the slab of points among the clusters' path chains (nearest chain point), one list per cluster
  std::vector<std::vector<int>> partitionSlab(const std::vector<std::vector<int>> &clusters, const std::vector<int> &slab,
                                              double max_dist) const;
  /// find the points within the branch section from its end points
  void extractNodesFromEnds(std::vector<int> &nodes);
  /// set the branch section tip position from the supplied list of Vertex IDs
  Eigen::Vector3d calculateTipFromVertices(const std::vector<int> &nodes) const;
  /// estimate the vector to the cylinder centre from the set of nodes
  Eigen::Vector3d vectorToCylinderCentre(const std::vector<int> &nodes, const Eigen::Vector3d &dir, double prior_radius = 0.0) const;
  /// estimate the cylinder's radius
  double estimateCylinderRadius(const std::vector<int> &nodes, const Eigen::Vector3d &dir, double &accuracy);
  /// estimate the cylinder's taper rate from its centre, @c dir and set of nodes
  void estimateCylinderTaper(double radius, double accuracy, bool extract_from_ends);
  /// add a new section to continue reconstructing the branch
  void addChildSection();
  /// calculate the ownership, what branch section does each point belong to
  void calculateSectionIds(std::vector<int> &section_ids, const std::vector<std::vector<int>> &children,
                           std::vector<std::vector<int>> *beyond = nullptr);
  /// set ids that are locel (0-based) per tree
  void generateLocalSectionIds();
  /// if using an overlapping grid, then remove trees with base outside the non-overlapping cell bounds
  void removeOutOfBoundSections(const Cloud &cloud, Eigen::Vector3d &min_bound, Eigen::Vector3d &max_bound, const Eigen::Vector3d &offset);
  /// colour the cloud based on the section id for each point
  void segmentCloud(Cloud &cloud, std::vector<int> &root_segs, const std::vector<int> &section_ids);
  /// remove points from the ray cloud if outside of the non-overlapping grid cell bounds
  void removeOutOfBoundRays(Cloud &cloud, const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound,
                            const std::vector<int> &root_segs);
  /// estimate the mean taper for the specified section
  double meanTaper(const BranchSection &section) const;
  double priorTaper(int root) const;  // the prior taper alone, before any vote
  /// the per-tree taper blend's parts: prior mean, the tree's own mean vote, data share and result
  struct TaperBlend { double prior, data, data_var, data_share, result; };
  TaperBlend taperBlend(int root) const;
  /// estimate the radius for the specified section
  double radius(const BranchSection &section) const;
  /// the along-path taper factor (d / d_root)^stem_taper, 1 at the tree base
  double stemTaper(const BranchSection &section) const;
  /// remove elements of nodes that are too distant to the set of end points
  bool removeDistantPoints(std::vector<int> &nodes);

  // cached data that is used throughout the processing method
  int sec_;
  const TreesParams *params_;
  const std::vector<int> *point_labels_; // per-ray tree labels when the cloud is pre-segmented, otherwise null
  std::vector<Vertex> points_;
  double forest_taper_{0};
  double forest_weight_{0};
  double forest_weight_squared_{0};
  std::vector<int> contiguous_section_ids_; // converts section ids to contiguous (empty trees removed) sectino ids for output

};

/// The structure for a single (cylindrical) branch section
struct RAYLIB_EXPORT BranchSection
{
  BranchSection()
    : tip(0, 0, 0)
    , taper(0.0)
    , weight(0.0)
    , parent(-1)
    , root(-1)
    , id(-1)
    , max_distance_to_end(0.0)
    , total_taper(0)
    , total_weight(0)
    , trunk_precision(0)
    , trunk_log_precision(0)
    , branch_precision(0)
    , branch_log_precision(0)
    , trunk_radius(0)
    , trunk_z(0)
    , len(0)
    , accuracy(0)
    , junction_weight(0)
    , radius_scale(1)
    , split_count(0)
    , tree_height(0)
  {}
  Eigen::Vector3d tip;

  double taper;
  double weight;
  int parent;
  int root;
  int id;  // 0 based per tree
  double max_distance_to_end;
  double total_taper; // a weighted taper estimate
  double total_weight; // the weighting
  double trunk_precision, trunk_log_precision;    // sums of 1 / var and log(vote) / var over trunk votes (roots only)
  double branch_precision, branch_log_precision;  // the same over branch votes, before their shared share-rule error
  double trunk_radius, trunk_z; // the measured trunk radius and the height of its girth slab (root sections only)
  double len, accuracy, junction_weight;
  double radius_scale;
  int split_count;
  double tree_height;
  std::vector<int> roots;  // root points
  std::vector<int> ends;
  std::vector<int> children;
  std::vector<int> slab;  // after a split: this section's share of the slab of points (slab_sections)
  double measured_radius{0};  // the radius vote's measured radius, for diagnostics
  int n_nodes{0};  // points the section was fitted on, for diagnostics
  int from_ends{0};  // 1 if fitted after a split
  Eigen::Vector3d centroid{0, 0, 0};  // the fitted points' centroid, before the centre step (for diagnostics)
  struct SplitTest  // the last split test at this section's top, for diagnostics (a split's new sections start empty)
  {
    int tests{0};             // split tests run here
    double gap{0}, span{0};   // the last one's thresholds
    int ends{0}, removed{0};  // the end points it clustered, and the clusters it dropped as too short (trunks)
    std::vector<int> clusters;  // each cluster's size in end points, before dropping any
  } split;
};

/// Converts an index in to a unique colour
/// only for ints up to 255*255*255-1   (as it leaves black as a special colour)
inline void RAYLIB_EXPORT convertIntToColour(int x, RGBA &colour)
{
  colour.red = colour.green = colour.blue = 0;
  x++;
  for (int i = 0; i < 24; ++i)
  {
    if (x & (1 << i))
    {
      int channel = i % 3;
      int offset = i / 3;
      uint8_t &ch = channel == 0 ? colour.red : (channel == 1 ? colour.green : colour.blue);
      ch |= (uint8_t)(1 << (7 - offset));
    }
  }
}

/// Converts a colour to a unique integer (index) value
/// returns -1 for the special case of the colour black, otherwise the integer index that the colour represents
inline int RAYLIB_EXPORT convertColourToInt(const RGBA &colour)
{
  int result = 0;
  for (int i = 0; i < 24; ++i)
  {
    int channel = i % 3;
    int offset = i / 3;
    const uint8_t &ch = channel == 0 ? colour.red : (channel == 1 ? colour.green : colour.blue);
    if (ch & (1 << (7 - offset)))
      result |= 1 << i;
  }
  return result - 1;
}
}  // namespace ray
#endif  // RAYLIB_RAYEXTRACT_TREES_H
