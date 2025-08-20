# Algebraic Multigrid Method for the Schwinger model.
Implementation of Algebric Multigrid Method for the Schwinger model in the context of Lattice simulations. 

Work in progress ...

- [x] Make aggregates and lattice blocks
- [x] Implement restriction and interpolation
- [x] Implement coarse gauge fields
    - [x] Test implementation for the two-level case
    - [x] Implement boundary conditions for all the levels
    - [x] Implement Level::makeCoarseLinks(Level& next_level) in the general scenario, without calling the $U_\mu$ variables.
    - [x] Orthonormalize test vectors
    - [x] Check that $P^\dagger D P$ coincides with my coarse gauge implementation for all levels 
    - [x] Test the implementation for several levels (Tested up to four levels)
- [x] Implement coarse grid matrix 
- [ ] Implement SAP for each level
    -[ ] Rewrite method for the Level0 and check that it coincide with the previous implementation
- [ ] Integrate everything in a V-cycle
- [ ] Improve performance (if needed)
- [ ] Implement an K-cycle
