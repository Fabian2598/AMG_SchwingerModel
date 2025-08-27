# Algebraic Multigrid Method for the Schwinger model.
Implementation of Algebric Multigrid Method for the Schwinger model in the context of Lattice simulations. 

Work in progress ...

- [x] Make aggregates and lattice blocks.
- [x] Implement restriction and interpolation.
- [x] Implement coarse gauge fields.
    - [x] Test implementation for the two-level case.
    - [x] Implement boundary conditions for all the levels.
    - [x] Implement Level::makeCoarseLinks(Level& next_level) in the general scenario, without calling the $U_\mu$ variables.
    - [x] Orthonormalize test vectors.
    - [x] Check that $P^\dagger D P$ coincides with my coarse gauge implementation for all levels.
    - [x] Test the implementation for several levels (Tested up to four levels).
- [x] Implement coarse grid matrix.
- [x] Implement SAP for each level
    - [x] Rewrite method for the Level0 and check that it coincide with the previous implementation
    - [x] Test the method for the coarser levels. Compare the convergence of the solution with GMRES and check that the inversion is working properly.
- [x] Integrate everything in a V-cycle.
- [X] Implement a K-cycle.
- [ ] Write the SetUp phase properly.
- [ ] Use the method as a preconditioner for FGMRES.
- [ ] Perform simple convergence tests for $V = 32^2, 64^2, 128^2$.
- [ ] Improve the performance and clean the code. 
- [ ] Compare the two-grid, the K-cycle and the V-cycle used as preconditioners for the FGMRES with the confs that I already have.
- [ ] Extensively document the code.
