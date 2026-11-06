# 3dtrees changes to raycloudtools

`main` mirrors upstream; `develop` is ours (`git log main..develop`). Newest first. All on `rayextract trees`; defaults match upstream.

- `--taper_weight W`: **changes fits.** Fixes the taper prior's weight so each tree depends only on its own points, not on the other trees or their order.
- `--graph_cache FILE`: save/reuse the shortest-path graph; same output, reruns ~4× faster.
- `--no_mesh`, `--labels_only`: skip `_trees_mesh.ply`; write per-ray tree ids instead of `_segmented.ply`.
- Lower kNN memory (−5.5 GB peak), same output.
