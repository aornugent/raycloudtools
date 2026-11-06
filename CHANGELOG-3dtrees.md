# 3dtrees changes to raycloudtools

What the [3dtrees](https://github.com/aornugent/3dtrees) fork adds to
[csiro-robotics/raycloudtools](https://github.com/csiro-robotics/raycloudtools). We pull from
upstream and never push to it.

- `main` mirrors upstream; `develop` is our trunk and merges upstream in.
- `git log main..develop` is everything that is ours.
- Each version below is a tag `3dtrees-N` on `develop`, newest first.
- New versions go on top. Older entries are rechunked now and then.

All flags are on `rayextract trees`. Unless a flag is set, output matches upstream.

## 3dtrees-1: 2026-11-05, upstream main `69a1e1d`, at `7307d86`

| Change | Changes fits? |
|---|---|
| Lower kNN memory | no: same output, 5.5 GB less peak memory |
| `--graph_cache FILE` | no: same output, reruns skip ~75 % of the work |
| `--labels_only` | no: writes `_labels.bin` instead of `_segmented.ply` |
| `--no_mesh` | no: skips `_trees_mesh.ply` |
| `--taper_weight W` | **yes, when set** (3dtrees always sets it, 195 for Alice Mulga) |

- **Lower kNN memory.** `connectPointsShortestPath` no longer keeps the full kNN distance matrix
  (5.5 GB of a 15.7 GB peak on 34.7M rays). It queries in chunks and recomputes distances where
  used.
- **`--graph_cache FILE`** saves the shortest-path graph (parent, root and distance to ground per
  point, 16 B each), and loads it again when the point count matches. Reruns that change only later
  flags, or that patch the tree fitting, skip the graph. The caller keys the file on the graph's
  inputs.
- **`--labels_only`** writes `<stub>_labels.bin` (one int32 tree id per ray, -1 for none) instead of
  rewriting every ray into `_segmented.ply` (1.67 GB for 34.7M rays).
- **`--no_mesh`** skips `_trees_mesh.ply`, which is only for viewing.
- **`--taper_weight W`: behaviour change.** Upstream shrinks every section radius toward a taper
  prior weighted by a running Σw²/Σw over all trees fitted so far. So each tree's radii depend on the
  other trees in the cloud and on processing order: crop the cloud or change one flag, and almost
  every tree moves. A fixed `W` (with the default nonzero `--global_taper`) makes each tree depend
  only on its own points. The run's Σw²/Σw is printed so a project value can be pooled. Default 0
  keeps upstream's behaviour.

### Not in `develop` yet

These are feature branches off `develop`, local only. Both reproduce the previous output at their
defaults, but change fits once their flags are set.

- `3dtrees/per-tree-taper` (`3f7c255`):
  - `--taper_prior_sd s` blends each tree's own taper with the prior in log space, weighted by the
    tree's evidence (trunk votes and pooled branch votes, by inverse variance).
  - `--taper_exponent e` makes the prior height-dependent: T0(H) = global_taper (H/10)^e.
  - `--stem_taper k` tapers radii along each path, so a stem is no longer one radius between forks.
  - Every run writes `<cloud>_taper.txt`, the per-tree taper diagnostics.
- `3dtrees/trunk-split` (`9cb9a61`): `--trunk_gap`, `--trunk_span_ratio` and `--trunk_splits` let a
  root group split into its stems without changing how branches split.
