// Copyright (c) 2022
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Thomas Lowe
#include "raytrees.h"
#include <nabo/nabo.h>
#include "rayclusters.h"
#include <cstdlib>

namespace ray
{
TreesParams::TreesParams()
  : max_diameter(0.9)
  , crop_length(1.0)
  , distance_limit(1.0)
  , height_min(2.0)
  , girth_height_ratio(0.12)
  , cylinder_length_to_width(4.0)
  , gap_ratio(0.016)
  , span_ratio(4.5)
  , gravity_factor(0.3)
  , grid_width(0.0)
  , segment_branches(false)
  , global_taper(0.012)
  , global_taper_factor(0.3)
  , taper_weight(0.0)
  , taper_prior_sd(0.0)
  , taper_exponent(0.0)
  , taper_vote_cv(0.1)
  , radius_floor(0.01)
  , branch_vote_sd(1.0)
  , stem_taper(0.0)
  , trunk_gap(0.0)
  , trunk_span_ratio(0.0)
  , trunk_splits(2)
  , geometric_centre(false)
  , slab_sections(false)
  , export_sections(false)
{}

/// The main reconstruction algorithm
/// It is based on finding the shortest paths using Djikstra's algorithm, followed
/// by an agglomeration of paths, with repeated splitting from root to tips
Trees::Trees(Cloud &cloud, const Eigen::Vector3d &offset, const Mesh &mesh, const TreesParams &params, bool verbose,
             const std::vector<int> *point_labels)
{
  // firstly, get the full set of shortest paths from ground to tips, and the set of roots
  params_ = &params;
  point_labels_ = point_labels;

  std::vector<std::vector<int>> roots_list =
    getRootsAndSegment(points_, cloud, mesh, params_->max_diameter, params_->distance_limit, params_->height_min,
                       params_->gravity_factor, point_labels_, params_->graph_cache);

  // Now we want to convert these paths into a set of branch sections, from root to tips
  // splitting as we go up...

  // calculate the maximum distance to tip for every point in the cloud
  calculatePointDistancesToEnd();

  // then generate all the root sections
  generateRootSections(roots_list);

  // the set of child points will be useful later, so generate them now
  std::vector<std::vector<int>> children(points_.size());
  for (size_t i = 0; i < points_.size(); i++)
  {
    if (points_[i].parent != -1)
    {
      children[points_[i].parent].push_back(static_cast<int>(i));
    }
  }

  // first do a special case for all the trunks. This is where we estimate mean taper
  ray::Cloud debug_cloud;
  for (sec_ = 0; sec_ < (int)sections_.size(); sec_++)
  {
    double best_accuracy = -1.0;
    std::vector<int> nodes;  // all the points in the section

    // find height to tip by following the children:
    std::vector<int> list = sections_[sec_].roots;
    double max_height = -1e10;
    for (size_t j = 0; j < list.size(); j++)
    {
      max_height = std::max(max_height, points_[list[j]].pos[2]);
      list.insert(list.end(), children[list[j]].begin(), children[list[j]].end());
    }
    Eigen::Vector3d base = getRootPosition();
    double tree_height = std::max(0.01, max_height - base[2]);
    sections_[sec_].tree_height = tree_height;

    // Use a usr defined taper to control the height up the trunk to calculate the radius at
    double girth_height = params_->girth_height_ratio * tree_height; // sections_[sec_].max_distance_to_end;
    double estimated_radius = 1e10;
    double best_dist = 0.0;
    Eigen::Vector3d best_tip;
    
    std::vector<int> best_nodes;
    std::vector<int> best_ends;
    for (int j = 1; j<=3; j++)
    {
      double max_dist = girth_height * (double)j / 2.0; // range from 0.5 to 1.5 times the specified height
      nodes.clear();
      sections_[sec_].ends.clear();
      extractNodesAndEndsFromRoots(nodes, base, children, max_dist * 2.0/3.0, max_dist);
      if (nodes.size() < 2)
      {
        continue;
      }
      sections_[sec_].tip = calculateTipFromVertices(nodes);
      if (verbose)
      {
        for (auto &node: nodes)
        {
          debug_cloud.addRay(Eigen::Vector3d(0,0,0), points_[node].pos, 0.0, ray::RGBA(j==1 ? 255 : 0, j==2 ? 255:0, j==3 ? 255:0, 255));
        }
      }         
      if (removeDistantPoints(nodes))
      {
        sections_[sec_].tip = calculateTipFromVertices(nodes);
      }
      if (verbose)
      {
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), sections_[sec_].tip + Eigen::Vector3d(0,0,0.03), 0.0, ray::RGBA(255,255,0, 255));
      }
      // shift to cylinder's centre
      Eigen::Vector3d up(0,0,1);
      const Eigen::Vector3d centroid = sections_[sec_].tip;
      sections_[sec_].tip += vectorToCylinderCentre(nodes, up);
      // now find the segment radius
      double accuracy;
      double radius = estimateCylinderRadius(nodes, up, accuracy);

      ray::RGBA col(j==1 ? 255 : 127, j==2 ? 255:127, j==3 ? 255:127, 255);
      if (verbose)
      {
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), sections_[sec_].tip, 0.0, col);
        for (double ang = 0; ang < 2.0*ray::kPi; ang += 0.1)
        {
          debug_cloud.addRay(Eigen::Vector3d(0,0,0), sections_[sec_].tip + radius * Eigen::Vector3d(std::sin(ang), std::cos(ang),0), 0.0, col);
        }
      }
      if (radius < estimated_radius)
      {
        best_accuracy = accuracy;
        estimated_radius = radius;
        best_dist = max_dist;
        best_tip = sections_[sec_].tip;
        sections_[sec_].centroid = centroid;
        best_nodes = nodes;
        best_ends = sections_[sec_].ends;
      }
    }
    if (best_dist == 0.0)
    {
      std::cout << "warning: could not find any points on trunk " << sec_ << " at " << base.transpose() << " so removing the whole section" << std::endl;
      sections_[sec_].tip = base + Eigen::Vector3d(0,0,0.01);
      sections_[sec_].total_weight = 1e-10;
      sections_[sec_].ends.clear(); // so this trunk is not ever used
      continue; 
    }    
    sections_[sec_].tip = best_tip;
    sections_[sec_].ends = best_ends;
    nodes = best_nodes;
    // when the trees are pre-labelled we trust the labels, so we don't split a trunk into several trees
    if (!point_labels_ && sections_[sec_].split_count < params_->trunk_splits)
    {
      double thickness = best_dist; 
      bool points_removed = false;
      // the trunk split can have its own thresholds, so stems of a clump separate without changing how
      // branches split (papers/frontier W6)
      double gap = params_->trunk_gap > 0.0 ? params_->trunk_gap
                                            : params_->gap_ratio * sections_[sec_].max_distance_to_end;
      double span = (params_->trunk_span_ratio > 0.0 ? params_->trunk_span_ratio : params_->span_ratio) * estimated_radius;
      std::vector<std::vector<int>> clusters = findPointClusters(base, points_removed, thickness, span, gap);

      if (clusters.size() > 1 || (points_removed && clusters.size() > 0))  // a bifurcation (or an alteration)
      {
        sections_[sec_].split_count++;
        bifurcate(clusters, thickness, children, true, true);
        sec_--;
        continue;
      }
    }
    if (verbose)
    {
      for (auto &node: sections_[sec_].ends)
      {
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), points_[node].pos + Eigen::Vector3d(0,0,0.02), 0.0, ray::RGBA(255, 0, 255, 255));
      }
      for (double ang = 0; ang < 2.0*ray::kPi; ang += 0.1)
      {
        uint8_t shade = 255; // (uint8_t)(best_accuracy * 255.0);
        debug_cloud.addRay(Eigen::Vector3d(0,0,0), best_tip + (estimated_radius + 0.01) * Eigen::Vector3d(std::sin(ang), std::cos(ang),0), 0.0, ray::RGBA(shade,shade,shade,255));
      }
    }

    nodes.clear();
    sections_[sec_].ends.clear();
    extractNodesAndEndsFromRoots(nodes, base, children, 0.0, best_dist/2.0); // make it lower
    // The ring is measured before the trunk's ends are clustered, so on a clump of stems it reads their envelope, and
    // as the tree's heaviest vote it would set every section's length and span (and so the fans at the base). An
    // envelope can only over-read, so a ring wider than the prior's radius doesn't vote; the taper is then learnt
    // from the sections, each measured on its own share of the wood.
    const double prior = priorTaper(sec_);
    if (prior <= 0.0 || estimated_radius <= prior * sections_[sec_].tree_height)
    {
      estimateCylinderTaper(estimated_radius, best_accuracy, false); // update the expected taper
    }
    sections_[sec_].trunk_radius = estimated_radius;
    sections_[sec_].trunk_z = best_tip[2];
    sections_[sec_].measured_radius = estimated_radius;
    sections_[sec_].n_nodes = (int)best_nodes.size();
  }
  if (verbose)
  {
    debug_cloud.translate(offset);
    debug_cloud.save("debug.ply");
  }

  // now trace from root tree nodes upwards, getting node centroids
  // create new BranchSections as we go
  for (sec_ = 0; sec_ < (int)sections_.size(); sec_++)
  {
    const int par = sections_[sec_].parent;
    if (par == -1)
    {
      // now add the single child for this particular tree node, assuming there are still ends
      if (sections_[sec_].ends.size() > 0)
      {
        addChildSection();
        for (const auto &i: sections_[sec_].roots)
        {
          sections_[sec_].tip[2] = std::min(sections_[sec_].tip[2], points_[i].pos[2]);
        }
      }      
      else
      {
        std::cout << "weird, a trunk without end points! " << sec_ << std::endl;
      }
      continue;
    }
    if (!(sec_ % 10000))
    {
      std::cout << "generating segment " << sec_ << std::endl;
    }

    // this branch can come to an end as it is now too small
    if (sections_[sec_].max_distance_to_end < params_->crop_length)
    {
      setBranchTip();
      continue;
    }

    std::vector<int> nodes;  // all the points in the section
    bool extract_from_ends = sections_[sec_].ends.size() > 0;
    // if the branch section has no end points recorded, then we need to examine this branch to
    // find end points and potentially branch (bifurcate)
    if (!extract_from_ends)
    {
      Eigen::Vector3d base = getRootPosition();
      double span_rad = radius(sections_[sec_]);
      // at cylinder_length_to_width x the prior radius a crown section is 2-5 cm: a few points along one shortest path,
      // so its centroid follows the path's zigzag rather than the branch. A floor keeps sections long enough to average
      const double min_section_length = 0.1;  // m
      double thickness = std::max(params_->cylinder_length_to_width * span_rad, min_section_length);
      extractNodesAndEndsFromRoots(nodes, base, children, 0.0, thickness);
      
      bool points_removed = false;
      double gap = params_->gap_ratio * sections_[sec_].max_distance_to_end; // gap threshold for splitting
      double span = params_->span_ratio * span_rad; // span thershold for splitting
      std::vector<std::vector<int>> clusters = findPointClusters(base, points_removed, thickness, span, gap);

      if (clusters.size() > 1 || (points_removed && clusters.size() > 0))  // a bifurcation (or an alteration)
      {
        extract_from_ends = true; // don't trust the found nodes as it is now two separate tree nodes
        bool add_offshoots = true; // par != -1 && sections_[par].parent != -1;  // when this is false it can lead to whole branches missing.
        // if points have been removed then this only resets the current section's points
        // otherwise it creates new branch sections_ for each cluster and adds to the end of the sections_ list
        bifurcate(clusters, thickness, children, false, add_offshoots, &nodes, span_rad);
      }
    }

    if (extract_from_ends) // we have split the ends, so we need to extract the set of nodes in a backwards manner
    {
      extractNodesFromEnds(nodes);
      // the path chains are a thin strip along one side of the branch: fit on this section's share of the slab
      if (params_->slab_sections && sections_[sec_].slab.size() > 5)
      {
        nodes.swap(sections_[sec_].slab);
      }
      std::vector<int>().swap(sections_[sec_].slab);
    }
    // estimate the section's tip (the centre of the cylinder of points)
    sections_[sec_].tip = calculateTipFromVertices(nodes);
    sections_[sec_].centroid = sections_[sec_].tip;
    // get section's direction
    Eigen::Vector3d dir = par >= 0 ? (sections_[sec_].tip - sections_[par].tip).normalized() : Eigen::Vector3d(0, 0, 1);
    // shift to cylinder's centre
    sections_[sec_].tip += vectorToCylinderCentre(nodes, dir, radius(sections_[sec_]));
    // re-estimate direction
    dir = par >= 0 ? (sections_[sec_].tip - sections_[par].tip).normalized() : Eigen::Vector3d(0, 0, 1);
    // now find the segment radius
    double accuracy = 0.0;
    double rad = estimateCylinderRadius(nodes, dir, accuracy);
    // and estimate taper
    estimateCylinderTaper(rad / (sections_[sec_].radius_scale * stemTaper(sections_[sec_])), accuracy, extract_from_ends);
    sections_[sec_].measured_radius = rad;
    sections_[sec_].n_nodes = (int)nodes.size();
    sections_[sec_].from_ends = extract_from_ends ? 1 : 0;

    // now add the single child for this particular tree node, assuming there are still ends
    if (sections_[sec_].ends.size() > 0)
    {
      addChildSection();
    }
  }  // end of loop. We now have created all of the BranchSections

  // Now calculate the section ids for all of the points, for the segmented cloud
  std::vector<int> section_ids(points_.size(), -1);
  std::vector<std::vector<int>> beyond;  // per tip, its points past its ends (export_sections)
  calculateSectionIds(section_ids, children, params_->export_sections ? &beyond : nullptr);
  if (params_->export_sections)
  {
    summariseModules(section_ids, beyond);
  }

  generateLocalSectionIds();

  Eigen::Vector3d min_bound(0, 0, 0), max_bound(0, 0, 0);
  // remove all sections with a root out of bounds, if we have gridded the cloud with an overlap
  if (params_->grid_width)
  {
    removeOutOfBoundSections(cloud, min_bound, max_bound, offset);
  }

  std::vector<int> root_segs(cloud.ends.size(), -1);
  // now colour the ray cloud based on the segmentation
  segmentCloud(cloud, root_segs, section_ids);

  if (params_->grid_width)  // also remove rays from the segmented cloud
  {
    removeOutOfBoundRays(cloud, min_bound, max_bound, root_segs);
  }

  std::cout << "cloud's estimated mean taper ratio (diameter / length): " << 2.0 * forest_taper_ / forest_weight_ << std::endl;
  std::cout << "cloud's taper weight sum(w^2)/sum(w): " << forest_weight_squared_ / forest_weight_ << " (sum w " << forest_weight_
            << ", sum w^2 " << forest_weight_squared_ << ")" << std::endl;
}

// If 1 tree plus some low lying foliage is in the node, then it won't be split if the foliage doesn't reach to the 
// top of the section (so doesn't have an end point), but it can create a very large radius estimate. 
// here we remove points that are too far from any end position. 
// this allows us to use most nodes to get an accurate radius estimate, but still uses only the end nodes to decide whether
// to split
bool Trees::removeDistantPoints(std::vector<int> &nodes)
{
  Eigen::Vector3d up(0,0,1);
  double max_rad = params_->gap_ratio * sections_[sec_].max_distance_to_end; // gap threshold for splitting
  double max_rad_sqr = max_rad * max_rad;
  bool found = false;
  for (int k = 0; k<(int)nodes.size(); k++)
  {
    if (nodes.size() <= 6)
    {
      break;
    }
    bool all_distant = true;
    for (auto &end: sections_[sec_].ends)
    {
      Eigen::Vector3d dif = points_[nodes[k]].pos - points_[end].pos;
      dif[2] = 0.0;
      if (dif.squaredNorm() < max_rad_sqr)
      {
        all_distant = false;
        break;
      }
    }
    if (all_distant)
    {
      nodes[k] = nodes.back();
      nodes.pop_back();
      k--;
      found = true;
    }
  }
  return found;
}

Trees::TaperBlend Trees::taperBlend(int root) const
{
  // Prior: log T ~ N(log T0(H), sd^2), T0(H) = global_taper * (H / 10)^taper_exponent.
  // Data: two pooled measurements of log T, combined by inverse variance (estimateCylinderTaper):
  //  - the trunk votes, each with log variance (taper_vote_cv^2 + (radius_floor / r)^2) / accuracy;
  //  - the branch votes, pooled the same way, plus branch_vote_sd^2 once: their error from the share rule they
  //    assume is shared by the whole tree, so it doesn't average away over many twigs.
  // Votes from radii under 2 radius_floor are dropped (the trust floor: below it a radius is mostly floor).
  // Nothing here reads another tree, so a tree's taper depends only on its own votes (research/20, frontier H 2.3).
  const BranchSection &r = sections_[root];
  TaperBlend b{};
  b.prior = params_->global_taper * std::pow(std::max(r.tree_height, 0.01) / 10.0, params_->taper_exponent);
  double precision = 0.0, log_precision = 0.0;
  if (r.trunk_precision > 0.0)
  {
    precision += r.trunk_precision;
    log_precision += r.trunk_log_precision;
  }
  if (r.branch_precision > 0.0)
  {
    double var = 1.0 / r.branch_precision + params_->branch_vote_sd * params_->branch_vote_sd;
    precision += 1.0 / var;
    log_precision += (r.branch_log_precision / r.branch_precision) / var;
  }
  if (precision <= 0.0)
  {
    b.result = b.prior;
    return b;
  }
  b.data = std::exp(log_precision / precision);
  b.data_var = 1.0 / precision;
  double var_prior = params_->taper_prior_sd * params_->taper_prior_sd;
  b.data_share = var_prior / (var_prior + b.data_var);
  b.result = std::exp(b.data_share * std::log(b.data) + (1.0 - b.data_share) * std::log(b.prior));
  return b;
}

double Trees::priorTaper(int root) const
{
  if (params_->taper_prior_sd > 0.0)
  {
    return taperBlend(root).prior;
  }
  return params_->global_taper ? params_->global_taper : (forest_weight_ > 0.0 ? forest_taper_ / forest_weight_ : 0.0);
}

double Trees::meanTaper(const BranchSection &section) const
{
  int root = section.root;
  if (params_->taper_prior_sd > 0.0)
  {
    return taperBlend(root).result;
  }
  double mean_taper = priorTaper(root);
  // the prior's weight: by default the running sum(w^2)/sum(w) over every taper vote so far in this run, which
  // couples each tree to the others in the cloud and to the processing order; a fixed taper_weight decouples them
  double mean_weight = params_->taper_weight > 0.0 ? params_->taper_weight : forest_weight_squared_ / forest_weight_;
  // a tree with no votes yet (its ring left out as an envelope, no section measured) takes the prior
  double taper = sections_[root].total_weight > 0.0 ? sections_[root].total_taper / sections_[root].total_weight
                                                    : mean_taper;
  double weight = sections_[root].total_weight;
  double blend = params_->global_taper_factor * params_->global_taper_factor * params_->global_taper_factor;
  mean_weight *= blend;
  weight *= 1.0 - blend;

  double result = (mean_taper * mean_weight + taper*weight) / (mean_weight + weight);
  if (!(result == result))
  {
    std::cout << "bad taper estimate" << std::endl;
    std::cout << "root: " << root << " estimated taper: " << result << ", mt: " << mean_taper << ", mw: " << mean_weight << ", t: " << taper << ", w: " << weight << ", blend: " << blend << std::endl;
  }
  return result;
}

double Trees::radius(const BranchSection &section) const 
{ 
  return sections_[section.root].tree_height * meanTaper(section) * section.radius_scale * stemTaper(section); // scaled down from root's estimated radius
}

double Trees::stemTaper(const BranchSection &section) const
{
  if (params_->stem_taper == 0.0)
  {
    return 1.0;
  }
  double d_root = std::max(sections_[section.root].max_distance_to_end, 1e-3);
  double ratio = std::min(1.0, std::max(section.max_distance_to_end / d_root, 0.02));
  return std::pow(ratio, params_->stem_taper);
}

void Trees::calculatePointDistancesToEnd()
{
  for (size_t i = 0; i < points_.size(); i++)
  {
    int id = static_cast<int>(i);
    int parent = points_[i].parent;
    points_[i].visited = false;  // reset the visited flag while we're here
    while (parent != -1)
    {
      const double end_dist =
        points_[id].distance_to_end + (points_[id].distance_to_ground - points_[parent].distance_to_ground);
      if (points_[parent].distance_to_end > end_dist)
      {
        break;
      }
      points_[parent].distance_to_end = end_dist;
      id = parent;
      parent = points_[id].parent;
    }
  }
}

// generate the initial branch sections, which are the root locations, defined by @c roots_list
// we also fill in the @c section_ids, which record which section each point belongs to
void Trees::generateRootSections(const std::vector<std::vector<int>> &roots_list)
{
  // first we initialise to one branch section per root
  for (size_t i = 0; i < roots_list.size(); i++)
  {
    BranchSection base;
    base.roots = roots_list[i];
    base.root = static_cast<int>(sections_.size());
    base.radius_scale = 1.0;
    sections_.push_back(base);
  }

  // we can now find the maximum distance to tip per root section, which allows us to estimate the
  // section radius for all the root sections.
  for (int i = static_cast<int>(sections_.size()) - 1; i >= 0; i--)
  {
    for (auto &root : sections_[i].roots)
    {
      sections_[i].max_distance_to_end = std::max(sections_[i].max_distance_to_end, points_[root].distance_to_end);
    }
  }
}

void Trees::setBranchTip()
{
  // before we finish, we need to calculate the tip of this branch, from its end points
  sections_[sec_].tip.setZero();
  for (auto &i : sections_[sec_].ends)
  {
    sections_[sec_].tip += points_[i].pos;
  }
  size_t sum = sections_[sec_].ends.size();
  // if it doesnt have end points then we'll have to use the root points
  if (sections_[sec_].ends.empty())
  {
    for (auto &i : sections_[sec_].roots)
    {
      sections_[sec_].tip += points_[i].pos;
    }
    sum = sections_[sec_].roots.size();
  }
  if (sum > 0)
  {
    sections_[sec_].tip /= static_cast<double>(sum);
  }
}

Eigen::Vector3d Trees::getRootPosition() const
{
  // get a base position for the section
  Eigen::Vector3d base(0, 0, 0);
  for (auto &root : sections_[sec_].roots)
  {
    base += points_[root].pos;
  }
  base /= static_cast<double>(sections_[sec_].roots.size());
  return base;
}


// find the node and end points for this section, from the root points
// note that the nodes contain everything between min_dist and max_dist and also the adjacent nodes either side
// i.e. a root node and an end node. This means the minimum number of nodes is 2
void Trees::extractNodesAndEndsFromRoots(std::vector<int> &nodes, const Eigen::Vector3d &base,
                                         const std::vector<std::vector<int>> &children, double min_dist, double max_dist)
{
  const int par = sections_[sec_].parent;
  std::vector<int> iteration_nodes = sections_[sec_].roots;

  // 2. find all the points (and the end points) for this section:
  for (unsigned int index = 0; index < iteration_nodes.size(); index++)
  {
    // we cannot use the for (int i: nodes) syntax above, as we are pushing into nodes as we go
    const int i = iteration_nodes[index];
    bool use_parent = false;
    for (auto &child : children[i])
    {
      const double dist = par == -1 ? points_[child].pos[2] - base[2] : (points_[child].pos - base).norm();
      if (dist > min_dist)
      {
        if (index < sections_[sec_].roots.size())
        {
          if (!use_parent)
          {
            nodes.push_back(i);
          }
          use_parent = true;
        }
        nodes.push_back(child);
      }
      if (dist < max_dist)  // in same slot, so accumulate
      {
        iteration_nodes.push_back(child);  // so we recurse on this child too
      }
      else
      {
        sections_[sec_].ends.push_back(child);
      }
    }
  }
}

// find clusters of points from the root points up the shortest paths, up to the cylinder length
// The algebraic (Kasa) circle through 2D points: its centre, radius and the rms of |q - c| - r. False if degenerate.
static bool fitCircle(const std::vector<Eigen::Vector2d> &qs, Eigen::Vector2d &centre, double &radius, double &rms)
{
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  Eigen::Vector3d b = Eigen::Vector3d::Zero();
  for (const auto &q : qs)
  {
    const Eigen::Vector3d row(q[0], q[1], 1.0);
    A += row * row.transpose();
    b += row * -q.squaredNorm();
  }
  const Eigen::Vector3d x = A.ldlt().solve(b);  // q.q + D x + E y + F = 0
  centre = Eigen::Vector2d(-x[0] / 2.0, -x[1] / 2.0);
  const double r2 = centre.squaredNorm() - x[2];
  if (!x.allFinite() || r2 <= 0.0)
  {
    return false;
  }
  radius = std::sqrt(r2);
  double sum = 0.0;
  for (const auto &q : qs)
  {
    sum += ((q - centre).norm() - radius) * ((q - centre).norm() - radius);
  }
  rms = std::sqrt(sum / static_cast<double>(qs.size()));
  return true;
}

// Split `cluster` (ids into ps) into two lobes when its points, seen along `axis`, are two stems side by side (a
// figure 8) rather than one: two circles, one per half of a 2-means split, fit them with under half the residual of
// one circle, each smaller than the one, side by side rather than overlapping. A ring or a one-sided arc of a single
// stem fits one circle, so it stays whole. The second lobe goes to `other`.
static bool splitLobes(const std::vector<Eigen::Vector3d> &ps, const Eigen::Vector3d &axis, std::vector<int> &cluster,
                       std::vector<int> &other)
{
  const size_t min_lobe = 4;           // points: a circle needs a few beyond its three
  const double residual_gain = 0.5;    // two circles must fit with under half one circle's residual
  const double side_by_side = 0.8;     // their centres over 0.8 x the sum of their radii apart
  if (cluster.size() < 2 * min_lobe)
  {
    return false;
  }
  const Eigen::Vector3d u = axis.unitOrthogonal(), v = axis.cross(u);
  std::vector<Eigen::Vector2d> qs;
  for (int i : cluster)
  {
    qs.emplace_back(ps[i].dot(u), ps[i].dot(v));
  }
  Eigen::Vector2d c1, cs[2];
  double r1, rms1, rs[2], rmss[2];
  if (!fitCircle(qs, c1, r1, rms1))
  {
    return false;
  }
  // 2-means, seeded with the two points farthest apart
  size_t a = 0, b = 0;
  for (size_t i = 0; i < qs.size(); i++)
  {
    if ((qs[i] - qs[0]).norm() > (qs[a] - qs[0]).norm()) a = i;
  }
  for (size_t i = 0; i < qs.size(); i++)
  {
    if ((qs[i] - qs[a]).norm() > (qs[b] - qs[a]).norm()) b = i;
  }
  Eigen::Vector2d m[2] = {qs[a], qs[b]};
  std::vector<int> side(qs.size(), 0);
  for (int it = 0; it < 10; it++)
  {
    Eigen::Vector2d sum[2] = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};
    int n[2] = {0, 0};
    for (size_t i = 0; i < qs.size(); i++)
    {
      side[i] = (qs[i] - m[1]).norm() < (qs[i] - m[0]).norm();
      sum[side[i]] += qs[i];
      n[side[i]]++;
    }
    if (n[0] == 0 || n[1] == 0)
    {
      return false;
    }
    m[0] = sum[0] / n[0];
    m[1] = sum[1] / n[1];
  }
  std::vector<Eigen::Vector2d> lobe[2];
  for (size_t i = 0; i < qs.size(); i++)
  {
    lobe[side[i]].push_back(qs[i]);
  }
  for (int k = 0; k < 2; k++)
  {
    if (lobe[k].size() < min_lobe || !fitCircle(lobe[k], cs[k], rs[k], rmss[k]) || rs[k] >= r1)
    {
      return false;
    }
  }
  const double rms2 = std::sqrt((rmss[0] * rmss[0] * static_cast<double>(lobe[0].size()) +
                                 rmss[1] * rmss[1] * static_cast<double>(lobe[1].size())) /
                                static_cast<double>(qs.size()));
  if (rms2 >= residual_gain * rms1 || (cs[0] - cs[1]).norm() <= side_by_side * (rs[0] + rs[1]))
  {
    return false;
  }
  std::vector<int> keep;
  other.clear();
  for (size_t i = 0; i < qs.size(); i++)
  {
    (side[i] ? other : keep).push_back(cluster[i]);
  }
  cluster = keep;
  return true;
}

std::vector<std::vector<int>> Trees::findPointClusters(const Eigen::Vector3d &base, bool &points_removed,
                              double thickness, double span, double gap)
{
  const int par = sections_[sec_].parent;
  std::vector<int> &all_ends = sections_[sec_].ends;

  // convert to a structure that is better for the cluster function
  std::vector<Eigen::Vector3d> ps;
  std::vector<int> v_indices;
  ps.reserve(all_ends.size());
  v_indices.reserve(all_ends.size());

  // cluster end points to find if we have separate branches
  // first, get interpolated edge points_. i.e. interpolate between two connected points inside and outside
  // the branch section's cylinder, so that the set edge_pos are right on the top boundary of the section
  for (auto &j : all_ends)
  {
    const double dist1j = par == -1 ? points_[j].pos[2] - base[2] : (points_[j].pos - base).norm();
    const double dist0j =
      par == -1 ? points_[points_[j].parent].pos[2] - base[2] : (points_[points_[j].parent].pos - base).norm();
    const double blendj = dist1j == dist0j ? 0.0 : (thickness - dist0j) / (dist1j - dist0j);
    Eigen::Vector3d edge_pos = points_[points_[j].parent].pos * (1.0 - blendj) + points_[j].pos * blendj;
    ps.push_back(edge_pos);
    v_indices.push_back(j);
  }

  // cluster these end points based on two separation criteria (gap_ratio and span_ratio)
  std::vector<std::vector<int>> clusters;
  generateClusters(clusters, ps, gap, span);
  // a cluster can still hold two stems side by side (a figure 8): touching stems link within the gap, and the span is
  // set by the prior radius of both together. Split it where two circles fit its points much better than one
  const Eigen::Vector3d axis = par == -1 || ps.empty() ? Eigen::Vector3d(0, 0, 1) : (mean(ps) - base).normalized();
  for (size_t c = 0; c < clusters.size(); c++)
  {
    std::vector<int> other;
    if (splitLobes(ps, axis, clusters[c], other))
    {
      clusters.push_back(other);
    }
  }
  BranchSection::SplitTest &test = sections_[sec_].split;  // record the test, for diagnostics (export_sections)
  test.tests++;
  test.gap = gap;
  test.span = span;
  test.ends = static_cast<int>(ps.size());
  test.removed = 0;
  test.clusters.clear();
  for (const auto &cluster : clusters)
  {
    test.clusters.push_back(static_cast<int>(cluster.size()));
  }
  // adjust back to global ids
  for (auto &cluster : clusters)
  {
    for (auto &id : cluster)
    {
      id = v_indices[id];
    }
  }

  points_removed = false;
  if (par == -1)  // if this is the root section
  {
    // then remove children that are smaller than the minimum tree height
    for (int i = static_cast<int>(clusters.size()) - 1; i >= 0; i--)
    {
      double max_dist = 0.0;
      for (auto &end : clusters[i])
      {
        max_dist = std::max(max_dist, points_[end].distance_to_end);
      }
      if (max_dist < params_->height_min)
      {
        test.removed++;
        clusters[i] = clusters.back();
        clusters.pop_back();
        points_removed = true;
      }
    }
  }
  return clusters;
}

// split into multiple branches and add as new branch sections to the end of the sections_ list
// that is being iterated through.
std::vector<std::vector<int>> Trees::partitionSlab(const std::vector<std::vector<int>> &clusters,
                                                   const std::vector<int> &slab, double max_dist) const
{
  // the path chains of each cluster: from each end point up its shortest-path parents to the section's roots
  std::vector<Eigen::Vector3d> chain_pos;
  std::vector<int> chain_cluster;
  std::vector<int> visited;
  for (size_t k = 0; k < clusters.size(); k++)
  {
    for (auto &end : clusters[k])
    {
      int node = end;
      while (node != -1)
      {
        if (std::find(visited.begin(), visited.end(), node) != visited.end())
        {
          break;
        }
        visited.push_back(node);
        chain_pos.push_back(points_[node].pos);
        chain_cluster.push_back(static_cast<int>(k));
        if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) != sections_[sec_].roots.end())
        {
          break;
        }
        node = points_[node].parent;
      }
    }
  }
  std::vector<std::vector<int>> shares(clusters.size());
  const double max_dist_sqr = max_dist * max_dist;
  for (auto &i : slab)
  {
    double best = max_dist_sqr;
    int best_k = -1;
    for (size_t j = 0; j < chain_pos.size(); j++)
    {
      const double d2 = (points_[i].pos - chain_pos[j]).squaredNorm();
      if (d2 < best)
      {
        best = d2;
        best_k = chain_cluster[j];
      }
    }
    if (best_k >= 0)
    {
      shares[best_k].push_back(i);
    }
  }
  return shares;
}

// whether any path down from `ends` reaches the current section's roots before meeting `wood` (another tree's nodes)
bool Trees::reachesGround(const std::vector<int> &ends, const std::vector<int> &wood) const
{
  const std::vector<int> &roots = sections_[sec_].roots;
  for (int end : ends)
  {
    for (int node = points_[end].parent; node != -1; node = points_[node].parent)
    {
      if (std::find(wood.begin(), wood.end(), node) != wood.end())
      {
        break;
      }
      if (std::find(roots.begin(), roots.end(), node) != roots.end())
      {
        return true;
      }
    }
  }
  return false;
}

void Trees::bifurcate(const std::vector<std::vector<int>> &clusters, double thickness, std::vector<std::vector<int>> &children, bool clip_tree, bool add_offshoots,
                      const std::vector<int> *slab, double prior_radius)
{
  const int par = sections_[sec_].parent;
  std::vector<std::vector<int>> shares;
  if (slab && params_->slab_sections && !clip_tree)
  {
    // a slab point is no further from its own chain than the far side of the branch (the prior under-reads)
    shares = partitionSlab(clusters, *slab, 3.0 * prior_radius + 0.01);
  }
  // find the maximum distance (to tip) for each cluster
  std::vector<double> max_distances(clusters.size(), 0.0);
  double maxmax = -1;
  int maxi = -1;
  for (size_t i = 0; i < clusters.size(); i++)
  {
    for (auto &end : clusters[i])
    {
      max_distances[i] = std::max(max_distances[i], points_[end].distance_to_end);
    }
    max_distances[i] += thickness;
    if (max_distances[i] > maxmax)
    {
      maxmax = max_distances[i];
      maxi = static_cast<int>(i);
    }
  }

  double total_area = 0.0;
  for (auto &dist: max_distances)
  {
    total_area += dist * dist;
  }

  // if clip_tree then we need to change the root nodes to be only those ones that
  // each path goes to... 
  
  // set the current branch section to the cluster with the
  // maximum distance to tip (i.e. the longest branch).
  sections_[sec_].ends = clusters[maxi];
  std::vector<int> all_nodes;
  std::vector<int> main_roots;

  if (clip_tree)
  {
    for (auto &end : sections_[sec_].ends)
    {
      int node = points_[end].parent;
      while (node != -1)
      {
        if (std::find(all_nodes.begin(), all_nodes.end(), node) != all_nodes.end())
        {
          break;
        }
        all_nodes.push_back(node);  // fill in the nodes in this branch section
        if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) != sections_[sec_].roots.end())
        {
          main_roots.push_back(node);
          break;
        }
        node = points_[node].parent;
      }
    }
  }

  sections_[sec_].max_distance_to_end = max_distances[maxi];

  if (add_offshoots)
  {
    // for all other clusters, add new sections to the list...
    for (size_t i = 0; i < clusters.size(); i++)
    {
      if (static_cast<int>(i) == maxi)
      {
        continue;
      }
      BranchSection new_node = sections_[sec_];
      new_node.split = BranchSection::SplitTest();  // its own tests are yet to run
      new_node.max_distance_to_end = max_distances[i];
      new_node.radius_scale *= max_distances[i] / std::sqrt(total_area);
      if (par == -1)
      {
        new_node.radius_scale = 1.0;
      }
      if (new_node.max_distance_to_end > params_->crop_length)  // but only add if they are large enough
      {
        // we only specify the end points at this stage. They will therefore enter the
        // extract_from_ends block below when the sec_ for loop gets to their section
        new_node.ends = clusters[i];
        new_node.slab.clear();
        if (!shares.empty())
        {
          new_node.slab = shares[i];
        }

        if (clip_tree && !reachesGround(clusters[i], all_nodes))
        {
          // every path down from this cluster meets the main tree's wood before the ground: it is a branch of that
          // tree, not a tree of its own. Left uncut, it stays among the trunk's ends and becomes a child section
          continue;
        }
        if (clip_tree)
        {
          std::vector<int> my_nodes;
          std::vector<int> my_roots;
          for (auto &end : new_node.ends)
          {
            int node = points_[end].parent;
            int child = end;
            while (node != -1)
            {
              if (std::find(my_nodes.begin(), my_nodes.end(), node) != my_nodes.end()) 
              {
                break; // just hit my own tree
              }
              if (std::find(all_nodes.begin(), all_nodes.end(), node) != all_nodes.end())
              {
                // hit another tree, so make a root here
                my_roots.push_back(child);
                // more than this, we need to unlink the children bit.
                auto &list = children[node];

                int find_count = 0;
                for (int j = 0; j<(int)list.size(); j++)
                {
                  if (list[j] == child)
                  {
                    find_count++;
                    list[j] = list.back();
                    list.pop_back();
                    j--;
                  }
                }
                points_[child].parent = -1;

                // we need to tell all of the children what the new root is
                std::vector<int> child_list = {child};
                for (size_t j = 0; j<child_list.size(); j++)
                {
                  points_[child_list[j]].root = child;
                  child_list.insert(child_list.end(), children[child_list[j]].begin(), children[child_list[j]].end());
                }
                break;
              }
              my_nodes.push_back(node);  // fill in the nodes in this branch section
              if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) != sections_[sec_].roots.end())
              {
                my_roots.push_back(node);
                break;
              }
              child = node;
              node = points_[node].parent;
            }
          }
          all_nodes.insert(all_nodes.end(), my_nodes.begin(), my_nodes.end());
          new_node.roots = my_roots;
        }

        if (par != -1)
        {
          sections_[par].children.push_back(static_cast<int>(sections_.size()));
        }
        if (par == -1)
        {
          new_node.root = (int)sections_.size();
        }
        sections_.push_back(new_node);
      }
    }
  }
  if (clip_tree)
  {
    sections_[sec_].roots = main_roots;
  }
  if (!shares.empty())
  {
    sections_[sec_].slab = shares[maxi];
  }
  if (par > -1)
  {
    sections_[sec_].radius_scale *= max_distances[maxi] / std::sqrt(total_area); // reduce branch radius due to the split
  }
}

// find the nodes between the section end points and the section root points
void Trees::extractNodesFromEnds(std::vector<int> &nodes)
{
  nodes = sections_[sec_].ends;
  for (auto &end : sections_[sec_].ends)
  {
    int node = points_[end].parent;
    if (node == -1)
    {
      std::cout << "shouldn't be parentless here" << std::endl;
      continue;
    }
    while (node != -1)
    {
      if (std::find(nodes.begin(), nodes.end(), node) != nodes.end())
      {
        break;
      }
      nodes.push_back(node);  // fill in the nodes in this branch section
      if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) != sections_[sec_].roots.end())
      {
        break;
      }
      node = points_[node].parent;
    }
  }
}

// use the list of Vertex IDs to estimate a tip location
Eigen::Vector3d Trees::calculateTipFromVertices(const std::vector<int> &vertices) const
{
  Eigen::Vector3d tip(0, 0, 0);
  if (vertices.empty())
  {
    std::cerr << "error: there shouldn't be empty nodes at this point in the processing" << std::endl;
  }
  for (auto &i : vertices)
  {
    tip += points_[i].pos;
  }
  if (!vertices.empty())
  {
    tip /= static_cast<double>(vertices.size());
  }

  return tip;
}

// Gauss-Newton on the geometric circle fit sum_i w_i (|q_i - c| - rho)^2 + w_p (rho - r0)^2, from the algebraic
// (Kasa) centre. The algebraic fit is pulled toward the visible arc when the noise is a sizeable part of the
// radius (one-sided, thin branches); the geometric one isn't. Huber weights shed foliage; the radius prior
// (worth a tenth of the points, and bounded to 0.5-3 r0) stops a short noisy arc's radius running away;
// each step is bounded by one radius, as the algebraic step is.
static Eigen::Vector2d geometricCircleCentre(const std::vector<Eigen::Vector3d> &ps, Eigen::Vector2d c,
                                             double approx_rad_sqr, double r0)
{
  const double n = static_cast<double>(ps.size());
  const double max_step = std::sqrt(approx_rad_sqr);
  const double max_shift = 2.0 * std::max(max_step, r0);
  double rho = 0.0;
  for (const auto &p : ps)
  {
    rho += (Eigen::Vector2d(p[0], p[1]) - c).norm();
  }
  rho /= n;
  const double w_p = r0 > 0.0 ? 0.1 * n : 0.0;
  if (r0 > 0.0)
  {
    rho = std::max(0.5 * r0, std::min(rho, 3.0 * r0));
  }
  for (int it = 0; it < 8; it++)
  {
    const double k = std::max(0.003, 0.25 * rho);  // Huber threshold
    Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
    Eigen::Vector3d g = Eigen::Vector3d::Zero();
    for (const auto &p : ps)
    {
      const Eigen::Vector2d v = Eigen::Vector2d(p[0], p[1]) - c;
      const double d = std::max(v.norm(), 1e-9);
      const double res = d - rho;
      const double w = std::abs(res) <= k ? 1.0 : k / std::abs(res);
      const Eigen::Vector3d J(-v[0] / d, -v[1] / d, -1.0);
      H += w * J * J.transpose();
      g += w * J * res;
    }
    H(2, 2) += w_p;
    g[2] += w_p * (rho - r0);
    H += 1e-12 * Eigen::Matrix3d::Identity();
    const Eigen::Vector3d delta = -H.ldlt().solve(g);
    if (!delta.allFinite())
    {
      break;
    }
    Eigen::Vector2d dc(delta[0], delta[1]);
    if (dc.norm() > max_step)
    {
      dc *= max_step / dc.norm();
    }
    c += dc;
    if (c.norm() > max_shift)
    {
      c *= max_shift / c.norm();
    }
    rho = std::max(rho + delta[2], 0.25 * max_step);
    if (r0 > 0.0)
    {
      rho = std::max(0.5 * r0, std::min(rho, 3.0 * r0));
    }
    if (dc.norm() < 1e-5)
    {
      break;
    }
  }
  return c;
}

Eigen::Vector3d Trees::vectorToCylinderCentre(const std::vector<int> &nodes, const Eigen::Vector3d &dir, double prior_radius) const
{
#define REAL_CENTROID  // finds a new centroid that is robust to branches scanned from a single side
#if defined REAL_CENTROID
  Eigen::Vector3d mean_p(0, 0, 0);
  std::vector<Eigen::Vector3d> ps;
  const Eigen::Vector3d vec(1, 2, 3);
  // obtain two orthogonal planes to the section's direction vector
  const Eigen::Vector3d ax1 = dir.cross(vec).normalized();
  const Eigen::Vector3d ax2 = dir.cross(ax1).normalized();
  double n = 0;
  // project points into a local space parabola
  for (auto &i : nodes)
  {
    const Eigen::Vector3d pos = points_[i].pos - sections_[sec_].tip;
    const Eigen::Vector2d offset(ax1.dot(pos), ax2.dot(pos));
    const Eigen::Vector2d xy = offset;
    const double l2 = xy.squaredNorm();
    const Eigen::Vector3d point(xy[0], xy[1], 0.5 * l2);  // a paraboloid that has gradient 1 at 1
    ps.push_back(point);
    mean_p += point;
    n++;
  }
  mean_p /= n;
  double approx_rad_sqr = 2.0 * mean_p[2];
  if (n > 5)  // assuming there are sufficient points for a resonable guess
  {
    // accumulation structure for plane least squares fitting
    struct Acc
    {
      Acc() { x2 = y2 = xy = xz = yz = 0; }
      double x2, y2, xy, xz, yz;
    };
    Acc plane;
    for (const auto &p : ps)
    {
      // fill in the parameters to estimate the plane of best fit to the paraboloid
      Eigen::Vector3d q = p - mean_p;
      plane.x2 += q[0] * q[0];
      plane.y2 += q[1] * q[1];
      plane.xy += q[0] * q[1];
      plane.xz += q[0] * q[2];
      plane.yz += q[1] * q[2];
    }
    const double eps = 1e-10;
    // is the plane is well-determined
    if (std::abs(plane.x2 * plane.y2 - plane.xy * plane.xy) > eps && std::abs(plane.y2) > eps)
    {
      // extract the local plane parameters
      const double A = (plane.xz * plane.y2 - plane.yz * plane.xy) / (plane.x2 * plane.y2 - plane.xy * plane.xy);
      const double B = (plane.yz - A * plane.xy) / plane.y2;

      Eigen::Vector2d shift(A, B);
      const double shift2 = shift.squaredNorm();
      if (shift2 > approx_rad_sqr)  // don't shift more than one radius each iteration, for safety
      {
        shift *= std::sqrt(approx_rad_sqr / shift2);
      }

      // trunk girth bands (no radius prior) stay on the algebraic fit: without a prior the geometric fit gave
      // flared or fused trunks bigger girth votes, raising the whole tree's taper (3dtrees research 54)
      if (params_->geometric_centre && prior_radius > 0.0)
      {
        shift = geometricCircleCentre(ps, shift, approx_rad_sqr, prior_radius);
      }
      // apply the plane parameters as a world-space shift in the branch section tip position
      return ax1 * shift[0] + ax2 * shift[1];
    }
  }
  else if (n > 0)  // if only a few points then use the mean instead
  {
    return ax1 * mean_p[0] + ax2 * mean_p[1];
  }
#endif
  return Eigen::Vector3d(0, 0, 0);
}

double Trees::estimateCylinderRadius(const std::vector<int> &nodes, const Eigen::Vector3d &dir, double &accuracy)
{
  double rad = 0.0;
  // get the mean radius
  double power = 0.25; // when 1 this is usual radius, but lower powers reduce outliers due to foliage
  std::vector<double> norms;
  norms.reserve(nodes.size());
  if (params_->use_rays)
  {
    for (auto &node : nodes)
    {
      Eigen::Vector3d offset = points_[node].pos - sections_[sec_].tip;
      offset -= dir * offset.dot(dir); // flatten    
      Eigen::Vector3d ray = points_[node].start - points_[node].pos;
      ray -= dir * ray.dot(dir); // flatten
      offset += ray*std::max(0.0, std::min(-offset.dot(ray)/ray.dot(ray), 1.0)); // move to closest point
      norms.push_back(offset.norm());
    }
  }
  else
  {
    for (auto &node : nodes)
    {
      Eigen::Vector3d offset = points_[node].pos - sections_[sec_].tip;
      offset -= dir * offset.dot(dir); // flatten    
      norms.push_back(offset.norm());
    }
  }

  for (auto &norm: norms)
  {
    rad += std::pow(norm, power);
  }
  const double eps = 1e-5; // prevent division by 0
  rad /= (double)nodes.size() + eps;
  rad = std::pow(rad, 1.0/power);
  double e = 0.0;
  for (auto &norm : norms)
  {
    e += std::abs(norm - rad);
  }
  e /= (double)nodes.size() + eps;
#define NEW_ACCURACY
#if defined NEW_ACCURACY
  const double sensor_noise = 0.02; // this prevents cylinders with a small number of very accurate points from totally dominating
  accuracy = rad / (e + sensor_noise);
#else // this one goes down to 0, which could be bad if all cylinder points were low accuracy
  accuracy = std::max(eps, rad - 2.0*e) / std::max(eps, rad); // so even distribution is accuracy 0, and cylinder shell is accuracy 1 
  accuracy *= std::min((double)nodes.size() / 3.0, 1.0);
#endif
  accuracy = std::max(accuracy, eps);
  rad = std::max(rad, eps);

  return rad;
}

void Trees::estimateCylinderTaper(double radius, double accuracy, bool extract_from_ends)
{
  int par = sections_[sec_].parent;
  int root = sections_[sec_].root;

  double l = sections_[sec_].radius_scale * sections_[root].tree_height;
  double L = sections_[root].tree_height;

  double junction_weight = 1.0;
  if (par > -1)
  {
    junction_weight = extract_from_ends ? 0.25 : 1.0; // extracting from ends means we are less confident about the quality
  }

  double weight = l*l*l * accuracy * junction_weight;
 // weight *= weight; // preference the strongest weight sections
  double taper = (radius/L) * weight;

  forest_taper_ += taper;
  forest_weight_ += weight;
  forest_weight_squared_ += weight*weight;

  sections_[root].total_taper += taper;
  sections_[root].total_weight += weight;
  // the per-tree blend's vote (taperBlend): log variance from the vote's own measured radius
  double measured = radius * sections_[sec_].radius_scale * stemTaper(sections_[sec_]);  // callers pass rad / (scale * stem taper)
  if (radius > 0.0 && measured >= 2.0 * params_->radius_floor)
  {
    double var = (params_->taper_vote_cv * params_->taper_vote_cv +
                  (params_->radius_floor / measured) * (params_->radius_floor / measured)) / std::max(accuracy, 0.05);
    double &prec = par > -1 ? sections_[root].branch_precision : sections_[root].trunk_precision;
    double &logprec = par > -1 ? sections_[root].branch_log_precision : sections_[root].trunk_log_precision;
    prec += 1.0 / var;
    logprec += std::log(radius / L) / var;
  }
  if (sections_[root].total_weight == 0.0)
  {
    std::cout << "bad section weight: " << l << "," << accuracy << "," << junction_weight << std::endl;
    std::cout << "sec: " << sec_ << ", root: " << root << std::endl;
  }

  sections_[sec_].taper = taper;
  sections_[sec_].weight = weight;
  sections_[sec_].len = l*l*l;
  sections_[sec_].accuracy = accuracy;
  sections_[sec_].junction_weight = junction_weight;
}

// add a child section to continue reconstructing the tree segments
void Trees::addChildSection()
{
  BranchSection new_node;
  new_node.parent = static_cast<int>(sec_);
  new_node.root = sections_[sec_].root;
  new_node.roots = sections_[sec_].ends;
  new_node.taper = sections_[sec_].taper;
  new_node.weight = sections_[sec_].weight;
  new_node.radius_scale = sections_[sec_].radius_scale;
  
  for (auto &root : new_node.roots)
  {
    new_node.max_distance_to_end = std::max(new_node.max_distance_to_end, points_[root].distance_to_end);
  }
  if (new_node.max_distance_to_end > params_->crop_length)  
  {
    sections_[sec_].children.push_back(static_cast<int>(sections_.size()));
    new_node.tip = calculateTipFromVertices(new_node.roots);
    sections_.push_back(new_node);
  }
}

// calculate what section every point belongs to
void Trees::calculateSectionIds(std::vector<int> &section_ids, const std::vector<std::vector<int>> &children,
                                std::vector<std::vector<int>> *beyond)
{
  // first we initialise to one branch section per root
  for (size_t i = 0; i < sections_.size(); i++)
  {
    if (sections_[i].parent == -1)
    {
      for (auto &root : sections_[i].roots)
      {
        section_ids[root] = static_cast<int>(i);
      }
    }
  }
  for (sec_ = 0; sec_ < (int)sections_.size(); sec_++)
  {
    std::vector<int> nodes;
    for (auto &end : sections_[sec_].ends)
    {
      int node = points_[end].parent;
      while (node != -1)
      {
        if (std::find(nodes.begin(), nodes.end(), node) != nodes.end())
        {
          break;
        }
        nodes.push_back(node);
        if (std::find(sections_[sec_].roots.begin(), sections_[sec_].roots.end(), node) !=
            sections_[sec_].roots.end())
        {
          break;
        }
        node = points_[node].parent;
      }
    }
    std::vector<int> ends = sections_[sec_].ends;
    if (sections_[sec_].children.size() == 0)  // then also add any offshoots...
    {
      for (size_t i = 0; i < ends.size(); i++)
      {
        for (auto &c : children[ends[i]])
        {
          ends.push_back(c);
        }
      }
      if (beyond)
      {
        beyond->resize(sections_.size());
        (*beyond)[sec_] = ends;
      }
    }
    for (auto &node : nodes)
    {
      section_ids[node] = static_cast<int>(sec_);
    }
    for (auto &end : ends)
    {
      section_ids[end] = static_cast<int>(sec_);
    }
  }

  // for points without a section id (e.g. end points)
  // look through their parents
  for (size_t i = 0; i < points_.size(); i++)
  {
    if (section_ids[i] == -1)  // an unfound point
    {
      int j = points_[i].parent;
      while (j != -1 && section_ids[j] == -1)
      {
        j = points_[j].parent;
      }
      if (j != -1)
      {
        section_ids[i] = section_ids[j];
      }
    }
  }
}

// generate local parent links. These are parent indices that are
// independent per tree
void Trees::generateLocalSectionIds()
{
  int num = 0;
  for (auto &section : sections_)
  {
    if (section.parent >= 0 || section.children.empty())
    {
      continue;
    }
    num++;
    int child_id = 0;
    section.id = child_id++;
    std::vector<int> children = section.children;
    for (unsigned int c = 0; c < children.size(); c++)
    {
      sections_[children[c]].id = child_id++;
      for (auto i : sections_[children[c]].children)
      {
        children.push_back(i);
      }
    }
  }
  std::cout << num << " trees saved" << std::endl;
}

void Trees::removeOutOfBoundSections(const Cloud &cloud, Eigen::Vector3d &min_bound, Eigen::Vector3d &max_bound, const Eigen::Vector3d &offset)
{
  const double width = params_->grid_width;
  cloud.calcBounds(&min_bound, &max_bound);
  const Eigen::Vector3d mid = (min_bound + max_bound)/2.0 + offset;
  const Eigen::Vector2d inds(std::round(mid[0] / width), std::round(mid[1] / width));
  min_bound[0] = width * (inds[0] - 0.5) - offset[0];
  min_bound[1] = width * (inds[1] - 0.5) - offset[1];
  max_bound[0] = width * (inds[0] + 0.5) - offset[0];
  max_bound[1] = width * (inds[1] + 0.5) - offset[1];
  std::cout << "min bound: " << (min_bound+offset).transpose() << ", max bound: " << (max_bound+offset).transpose() << std::endl;

  // disable trees out of bounds
  for (auto &section : sections_)
  {
    if (section.parent >= 0 || section.children.empty())
    {
      continue;
    }
    const Eigen::Vector3d pos = section.tip;
    if (pos[0] < min_bound[0] || pos[0] > max_bound[0] || pos[1] < min_bound[1] || pos[1] > max_bound[1])
    {
      section.children.clear();  // make it a non-tree
    }
  }
}

// colour the cloud by tree id, or by branch segment id
void Trees::segmentCloud(Cloud &cloud, std::vector<int> &root_segs, const std::vector<int> &section_ids)
{
  contiguous_section_ids_.resize(sections_.size(), -1); // these are different to the root section IDs as they exclude empty trees
  int num_trees = 0;
  bool first_non_root = false;
  for (size_t sec = 0; sec < sections_.size(); sec++)
  {
    const auto &section = sections_[sec];
    if (section.parent < 0 && first_non_root)
    {
      std::cout << "Error: root sections after some non-roots. ID: " << sec << std::endl;
    }
    if (section.parent >= 0)
      first_non_root = true;
    else if (section.children.empty())  // not a root section, so move on
    {
      continue;
    }
    contiguous_section_ids_[sec] = num_trees++;
  }

  if (params_->export_sections)
  {
    ray_section_ids_.assign(cloud.ends.size(), -1);
  }
  int j = -1;
  for (size_t i = 0; i < cloud.ends.size(); i++)
  {
    RGBA &colour = cloud.colours[i];
    if (cloud.rayBounded(i))
    {
      j++;
      int seg = section_ids[j];
      // the raw ownership, before the trunk cone below; only sections of written trees
      if (params_->export_sections && seg != -1 && !sections_[sections_[seg].root].children.empty())
      {
        ray_section_ids_[i] = contiguous_section_ids_[seg];
      }
      const int root_id = points_[j].root;
      root_segs[i] = root_id == -1 ? -1 : section_ids[root_id];

      // This block places a cone of gradient 2 around the trunk direction, to remove the square of initial ground points
      if (seg != -1 && sections_[seg].parent == -1)
      {
        Eigen::Vector3d dir(0,0,1e-10);
        for (auto &child : sections_[seg].children)
        {
          dir += sections_[child].tip - sections_[seg].tip;
        }
        dir.normalize();
        const double grad = 2.0; // larger cuts out a steeper (narrower) cone
        Eigen::Vector3d base = sections_[seg].tip - grad*dir*radius(sections_[seg]);
        Eigen::Vector3d dif = cloud.ends[i] - base;
        double h = dif.dot(dir);
        double w = (dif - dir*h).norm();
        if (grad*w > h)
        {
          colour.red = colour.green = colour.blue = 0;
          continue;
        }
      }

      if (!params_->segment_branches)
      {
        if (root_id == -1)
        {
          colour.red = colour.green = colour.blue = 0;
          continue;
        }
        seg = root_segs[i];
      }
      if (seg == -1)
      {
        colour.red = colour.green = colour.blue = 0;
        continue;
      }
      convertIntToColour(contiguous_section_ids_[seg], colour);
    }
    else
    {
      colour.red = colour.green = colour.blue = 0;
    }
  }
}

// remove rays from the ray cloud where the end points are out of bounds
void Trees::removeOutOfBoundRays(Cloud &cloud, const Eigen::Vector3d &min_bound, const Eigen::Vector3d &max_bound,
                                 const std::vector<int> &root_segs)
{
  for (int i = static_cast<int>(cloud.ends.size()) - 1; i >= 0; i--)
  {
    if (!cloud.rayBounded(i))
    {
      continue;
    }
    const Eigen::Vector3d pos = root_segs[i] == -1 ? cloud.ends[i] : sections_[root_segs[i]].tip;

    if (pos[0] < min_bound[0] || pos[0] > max_bound[0] || pos[1] < min_bound[1] ||
        pos[1] > max_bound[1])  // nope, can't do this here!
    {
      cloud.starts[i] = cloud.starts.back();
      cloud.starts.pop_back();
      cloud.ends[i] = cloud.ends.back();
      cloud.ends.pop_back();
      cloud.colours[i] = cloud.colours.back();
      cloud.colours.pop_back();
      cloud.times[i] = cloud.times.back();
      cloud.times.pop_back();
      if (!ray_section_ids_.empty())
      {
        ray_section_ids_[i] = ray_section_ids_.back();
        ray_section_ids_.pop_back();
      }
    }
  }
}

void Trees::summariseModules(const std::vector<int> &section_ids, const std::vector<std::vector<int>> &beyond)
{
  modules_.assign(sections_.size(), Module());
  for (size_t sec = 0; sec < beyond.size(); sec++)
  {
    Module &m = modules_[sec];
    std::vector<int> pts;  // its points past its ends that the final ownership keeps
    for (const int p : beyond[sec])
    {
      if (section_ids[p] == static_cast<int>(sec))
      {
        pts.push_back(p);
      }
    }
    for (const int e : sections_[sec].ends)
    {
      if (section_ids[e] == static_cast<int>(sec))
      {
        m.n_ends++;
        m.reach = std::max(m.reach, points_[e].distance_to_end);
      }
    }
    m.n = static_cast<int>(pts.size());
    if (pts.empty())
    {
      continue;
    }
    m.top = -1e10;
    for (const int p : pts)
    {
      m.mean += points_[p].pos;
      m.top = std::max(m.top, points_[p].pos[2]);
    }
    m.mean /= static_cast<double>(pts.size());
    for (const int p : pts)
    {
      const Eigen::Vector3d d = points_[p].pos - m.mean;
      m.cov += d * d.transpose();
    }
    m.cov /= static_cast<double>(pts.size());
  }
}

// --export_sections: what the fit measured and decided, beside _trees.txt and joined to it by (tree, seg)
bool Trees::saveSections(const std::string &stub, const Eigen::Vector3d &offset) const
{
  if (!params_->export_sections)
  {
    return false;
  }
  std::ofstream ofs((stub + "_sections.txt").c_str(), std::ios::out);
  std::ofstream bin((stub + "_sections.bin").c_str(), std::ios::binary);
  std::ofstream mfs((stub + "_modules.txt").c_str(), std::ios::out);
  auto opened = [&](const std::ofstream &f, const std::string &ext) {
    if (!f.is_open())
    {
      std::cerr << "Error: cannot open " << stub + ext << " for writing." << std::endl;
    }
    return f.is_open();
  };
  if (!opened(ofs, "_sections.txt") || !opened(bin, "_sections.bin") || !opened(mfs, "_modules.txt"))
  {
    return false;
  }
  auto key = [&](int sec) {  // tree, seg (the segment's index in its _trees.txt line), section_id
    return std::to_string(contiguous_section_ids_[sections_[sec].root]) + "," + std::to_string(sections_[sec].id) +
           "," + std::to_string(contiguous_section_ids_[sec]);
  };
  auto pos = [&](std::ostream &o, const Eigen::Vector3d &p) {
    o << std::fixed << std::setprecision(4) << p[0] + offset[0] << "," << p[1] + offset[1] << "," << p[2] + offset[2]
      << std::defaultfloat << std::setprecision(6);
  };
  // the sections of each written tree in _trees.txt order: the root, then its descendants breadth first
  std::vector<int> order;
  for (size_t sec = 0; sec < sections_.size(); sec++)
  {
    if (sections_[sec].parent >= 0 || sections_[sec].children.empty())
    {
      continue;
    }
    size_t first = order.size();
    order.push_back(static_cast<int>(sec));
    for (size_t c = first; c < order.size(); c++)
    {
      order.insert(order.end(), sections_[order[c]].children.begin(), sections_[order[c]].children.end());
    }
  }

  ofs << "# Sections, one per _trees.txt segment in the same order: tree,seg join to _trees.txt (seg = index in the "
         "tree's line), section_id = its section_id column = the ids in _sections.bin. rad: the measured radius (power "
         "mean) and its accuracy, from n_nodes points (from_ends: after a split); centroid: those points' centroid "
         "before the centre step (for a root, the girth band's, above the base _trees.txt gives); reach: graph distance "
         "to the farthest tip; split_*: the last of split_tests split tests at the section's top (clusters: end points "
         "per cluster, before removing those too short for a tree)"
      << std::endl;
  ofs << "tree,seg,section_id,parent_seg,rad,accuracy,n_nodes,from_ends,centroid_x,centroid_y,centroid_z,reach,"
         "radius_scale,split_tests,split_gap,split_span,split_ends,split_removed,split_clusters"
      << std::endl;
  ofs << std::setprecision(6);
  for (const int sec : order)
  {
    const BranchSection &s = sections_[sec];
    ofs << key(sec) << "," << (s.parent >= 0 ? sections_[s.parent].id : -1) << "," << s.measured_radius << ","
        << s.accuracy << "," << s.n_nodes << "," << s.from_ends << ",";
    pos(ofs, s.centroid);
    ofs << "," << s.max_distance_to_end << "," << s.radius_scale << "," << s.split.tests << "," << s.split.gap << ","
        << s.split.span << "," << s.split.ends << "," << s.split.removed << ",";
    for (size_t k = 0; k < s.split.clusters.size(); k++)
    {
      ofs << (k ? " " : "") << s.split.clusters[k];
    }
    ofs << std::endl;
  }

  // per ray of the cloud, in cloud order: the section_id that owns it (calculateSectionIds), -1 for none
  for (const int id : ray_section_ids_)
  {
    const int32_t v = id;
    bin.write(reinterpret_cast<const char *>(&v), sizeof(v));
  }

  mfs << "# Per tip (childless section): the points the ownership assigns past its end points, i.e. the ends and "
         "everything above them on the shortest-path graph. n_ends: of those, the tip's own end points; mean and cov "
         "(population, m^2) of their positions; top_z: the highest; reach: graph distance from those end points to the "
         "farthest tip"
      << std::endl;
  mfs << "tree,seg,section_id,n,n_ends,mean_x,mean_y,mean_z,cov_xx,cov_xy,cov_xz,cov_yy,cov_yz,cov_zz,top_z,reach"
      << std::endl;
  mfs << std::setprecision(6);
  for (const int sec : order)
  {
    if (!sections_[sec].children.empty())
    {
      continue;
    }
    const Module &m = modules_[sec];
    mfs << key(sec) << "," << m.n << "," << m.n_ends << ",";
    if (m.n == 0)
    {
      mfs << "nan,nan,nan,nan,nan,nan,nan,nan,nan,nan," << m.reach << std::endl;
      continue;
    }
    pos(mfs, m.mean);
    mfs << "," << m.cov(0, 0) << "," << m.cov(0, 1) << "," << m.cov(0, 2) << "," << m.cov(1, 1) << "," << m.cov(1, 2)
        << "," << m.cov(2, 2) << "," << std::fixed << std::setprecision(4) << m.top + offset[2] << std::defaultfloat
        << std::setprecision(6) << "," << m.reach << std::endl;
  }
  return ofs.good() && bin.good() && mfs.good();
}

// save the structure to a text file
bool Trees::save(const std::string &filename, const Eigen::Vector3d &offset, bool verbose) const
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  ofs << std::setprecision(4) << std::fixed;
  ofs << "# Tree file. Optional per-tree attributes (e.g. 'height,crown_radius, ') followed by 'x,y,z,radius' and any additional per-segment attributes:" << std::endl;
  ofs << "x,y,z,radius,parent_id,section_id"; // simple format
  if (verbose)
  {
    ofs << ",weight,len,accuracy,junction_weight";
  }
  ofs << std::endl;  
  for (size_t sec = 0; sec < sections_.size(); sec++)
  {
    const auto &section = sections_[sec];
    if (section.parent >= 0 || section.children.empty())  // not a root section, so move on
    {
      continue;
    }
    ofs << section.tip[0]+offset[0] << "," << section.tip[1]+offset[1] << "," << section.tip[2]+offset[2] << "," << radius(section) << ",-1," << contiguous_section_ids_[sec];
    if (verbose)
    {
      ofs << "," << section.weight << "," << section.len << "," << section.accuracy << "," << section.junction_weight;
    }
    int root = section.root;
    std::vector<int> children = section.children;
    for (unsigned int c = 0; c < children.size(); c++)
    {
      const BranchSection &node = sections_[children[c]];
      if (node.root != root)
      {
        std::cout << "bad format: " << node.root << " != " << root << std::endl;
      }
      ofs << ", " << node.tip[0]+offset[0] << "," << node.tip[1]+offset[1] << "," << node.tip[2]+offset[2] << "," << radius(node) << "," << sections_[node.parent].id << "," << contiguous_section_ids_[children[c]];
      if (verbose)
      {
        ofs << "," << node.weight << "," << node.len << "," << node.accuracy << "," << node.junction_weight;
      }
      for (auto i : sections_[children[c]].children)
      {
        children.push_back(i);
      }
    }
    ofs << std::endl;
  }
  return true;
}

bool Trees::saveTaper(const std::string &filename, const Eigen::Vector3d &offset) const
{
  std::ofstream ofs(filename.c_str(), std::ios::out);
  if (!ofs.is_open())
  {
    std::cerr << "Error: cannot open " << filename << " for writing." << std::endl;
    return false;
  }
  // one line per tree, in _trees.txt order; tapers are radius per length (half the CLI's diameter per length)
  ofs << "# per-tree taper blend: height (m), base z, measured trunk radius (m) and its slab height z, prior T0(H), "
         "the tree's inverse-variance mean vote and its log variance, data share, the taper used. With "
         "--taper_prior_sd the data share is the conjugate weight; otherwise data is RCT's weighted mean vote and "
         "data_share the default blend's data fraction" << std::endl;
  ofs << "height,base_z,trunk_radius,trunk_z,prior,data,data_var,data_share,taper" << std::endl;
  ofs << std::setprecision(6);
  for (size_t sec = 0; sec < sections_.size(); sec++)
  {
    const auto &section = sections_[sec];
    if (section.parent >= 0 || section.children.empty())  // the same filter as save()
    {
      continue;
    }
    int root = section.root;
    TaperBlend b = taperBlend(root);
    double used = meanTaper(section);
    if (params_->taper_prior_sd <= 0.0)  // report the default blend's data fraction for comparison
    {
      double prior = params_->global_taper ? params_->global_taper : (forest_taper_ / forest_weight_);
      b.prior = prior;
      b.data = sections_[root].total_weight > 0.0 ? sections_[root].total_taper / sections_[root].total_weight : 0.0;
      b.data_share = std::abs(b.data - prior) > 1e-12 ? (used - prior) / (b.data - prior) : 1.0;
    }
    ofs << section.tree_height << "," << section.tip[2] + offset[2] << ","
        << section.trunk_radius << "," << section.trunk_z + offset[2] << "," << b.prior << "," << b.data << "," << b.data_var
        << "," << b.data_share << "," << used << std::endl;
  }
  return true;
}

}  // namespace ray
