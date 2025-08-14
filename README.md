# Algebraic Multigrid Method for the Schwinger model.
Implementation of Algebric Multigrid Method for the Schwinger model in the context of Lattice simulations. 

Work in progress ...

- [x] Make aggregates and lattice blocks
- [x] Implement restriction and interpolation
- [ ] Implement coarse gauge fields
    - [x] Test implementation for the two-level case
    - [x] Implement boundary conditions for all the levels
    - [x] Implement Level::makeCoarseLinks(Level& next_level) in the general scenario, without calling the $U_\mu$ variables.
    - [x] Orthonormalize test vectors
    - [ ] Check that $P^\dagger D P$ coincides with my coarse gauge implementation for all levels 
            This currently works for the first two levels, i.e. level0 and level1. For level2 the result does not coincide. This
            could be related to the interpolators or the gauge links implementation, which I should check carefully. 
    - [ ] Test the implementation for several levels
- [x] Implement coarse grid matrix 
- [ ] Implement SAP for each level
- [ ] Integrate everything in a V-cycle
- [ ] Improve performance (if needed)
- [ ] Implement an K-cycle
