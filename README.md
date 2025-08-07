# Algebraic Multigrid Method for the Schwinger model.
Implementation of Algebric Multigrid Method for the Schwinger model in the context of Lattice simulations. 

Work in progress ...

- [x] Make aggregates and lattice blocks
- [x] Implement restriction and interpolation
- [ ] Implement coarse gauge fields
    -[x] Test implementation for the two-level case
    - [ ] Implement boundary conditions for all the levels
    - [ ] Implement Level::makeCoarseLinks(Level& next_level) in the general scenario, without calling the $U_\mu$ variables.
- [ ] Implement coarse grid matrix 
- [ ] Implement SAP for each level
- [ ] Integrate everything in a V-cycle
- [ ] Improve performance (if needed)
- [ ] Implement an K-cycle