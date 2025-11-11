#include "Poisson1D.hpp"

void Poisson1D::setup() {
  // Create the mesh (through dealII functions for simple domains)
  // colorize = true assignis
  GridGenerator::subdivided_hyper_cube(mesh, N_el, 0.0, 1.0,
                                       /* colorize = */ true);

  std::cout << " Number of elements: " << mesh.n_active_cells() << std::endl;

  // Initialize the finite element space
  fe = std::make_unique<FE_SimplexP<dim>>(r);
  quadrature = std::make_unique<QGaussSimplex<dim>>(r + 1);

  std::cout << " DoF per cell: " << fe->dofs_per_cell << std::endl;

  dof_handler.reinit(mesh);
  dof_handler.distribute_dofs(*fe);

  std::cout << " Number of DoFs: " << dof_handler.n_dofs << std::endl;

  // Initialize linear algebra
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern_copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  system_rhs.reinit(dof_handler.n_dofs());
  solution.reinit(dof_handler.n_dofs());
}

void Poisson1D::assemble() {}

void Poisson1D::solve() {}

void Poisson1D::output() const {}