/*
  Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/boundary_fluid_pressure/density.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/melt.h>
#include <aspect/simulator_signals.h>
#include <utility>
#include <limits>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/mpi.h>
#include <algorithm>
#include <numeric>
#include <aspect/geometry_model/box.h>


namespace aspect
{
  namespace BoundaryFluidPressure
  {

    template <int dim>
    void
    Density<dim>::
    fluid_pressure_gradient (const types::boundary_id boundary_indicator,
                             const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                             const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                             const std::vector<Tensor<1,dim>> &normal_vectors,
                             std::vector<double> &fluid_pressure_gradient_outputs) const
    {
      const std::shared_ptr<const MaterialModel::MeltOutputs<dim>> melt_outputs
        = material_model_outputs.template get_additional_output_object<MaterialModel::MeltOutputs<dim>>();
      Assert(melt_outputs != nullptr, ExcMessage("Error, MeltOutputs are missing in fluid_pressure_gradient()"));
      for (unsigned int q=0; q<fluid_pressure_gradient_outputs.size(); ++q)
        {
          const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(material_model_inputs.position[q]);

          switch (density_formulation)
            {
              case DensityFormulation::solid_density:
              {
                fluid_pressure_gradient_outputs[q] = (material_model_outputs.densities[q] * gravity) * normal_vectors[q];
                break;
              }

              case DensityFormulation::fluid_density:
              {
                fluid_pressure_gradient_outputs[q] = (melt_outputs->fluid_densities[q] * gravity) * normal_vectors[q];
                break;
              }

              case DensityFormulation::average_density:
              {
                const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
                const double phi =  material_model_inputs.composition[q][porosity_idx];
                fluid_pressure_gradient_outputs[q] = ((1.0 - phi) * material_model_outputs.densities[q] * gravity
                                                      + phi * melt_outputs->fluid_densities[q] * gravity)
                                                     * normal_vectors[q];
                break;
              }

                case DensityFormulation::weighted_density:
              {
                fluid_pressure_gradient_outputs[q] = ((1.0 - fluid_density_weight) * material_model_outputs.densities[q] * gravity
                                                      + fluid_density_weight * melt_outputs->fluid_densities[q] * gravity)
                                                     * normal_vectors[q];
                break;
              }

              case DensityFormulation::side_boundary_magma_extraction:
              {
                const double y = material_model_inputs.position[q][dim-1];
                if (boundary_indicator == extraction_boundary_id &&
                    y >= extraction_y_min && y <= extraction_y_max)
                  {
                    const double depth  = this->get_geometry_model().depth(material_model_inputs.position[q]);
                    const double P_lith = interpolate_lith_pressure(depth);
                    const double rho_l  = melt_outputs->fluid_densities[q];
                    // (grad p_f)·n_outward should be negative for outward flux 
                    fluid_pressure_gradient_outputs[q] =-side_pressure_gradient_weight* (P_lith - rho_l * gravity.norm() * depth) / extraction_dx;
                  }
                else
                  {
                    fluid_pressure_gradient_outputs[q] = (melt_outputs->fluid_densities[q] * gravity) * normal_vectors[q];
                  }

                break;
              }


              default:
                Assert (false, ExcNotImplemented());
            }
        }
    }

 
    template <int dim>
    void
    Density<dim>::initialize ()
    {
      if (density_formulation == DensityFormulation::side_boundary_magma_extraction)
      {
        extraction_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id(extraction_boundary_name);

      // Connect to the start_timestep signal so recompute_lith_pressure_profile()
      // is called collectively on ALL MPI processes at the start of each timestep,
      // before Stokes assembly. 
        this->get_signals().start_timestep.connect([this](const SimulatorAccess<dim> &)
        {
          this->recompute_lith_pressure_profile();
        });
      }
    }

    template <int dim>
    void
    Density<dim>::recompute_lith_pressure_profile ()
    {
      // Verify that the geometry is a Box (required to get x-coordinates of boundaries
      // and the y-coordinate of the top surface).
      const GeometryModel::Box<dim> *box =
        dynamic_cast<const GeometryModel::Box<dim>*>(&this->get_geometry_model());
      AssertThrow(box != nullptr,
                  ExcMessage("The 'side boundary magma extraction' formulation "
                             "requires a box geometry."));

      // y-coordinate of the top surface (surface = largest y in a Box geometry)
      const double y_top     = box->get_origin()[dim-1] + box->get_extents()[dim-1];
      const unsigned int n_field = this->n_compositional_fields();
      const auto &intr = this->introspection();

      // FEFaceValues with one quadrature point at the center of each face (QMidpoint).
      // This lets us read the FE solution (T, composition) directly at boundary face centers
      // rather than evaluating at off-mesh points.
      const QMidpoint<dim-1> face_quadrature;
      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                        this->get_dof_handler().get_fe(),
                                        face_quadrature,
                                        update_values | update_quadrature_points);

      // Loop over locally owned cells and collect face-center data for all faces that
      // lie on the extraction boundary within the extraction y-range.
      // Each face contributes one entry in local_packed with format:
      //   [y_face, T_face, comp_0_face, comp_1_face, ..., comp_{n_field-1}_face]
      const unsigned int stride = 2 + n_field;   // doubles per face entry
      std::vector<double> local_packed;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (unsigned int f = 0; f < cell->n_faces(); ++f)
            if (cell->face(f)->at_boundary() &&
                cell->face(f)->boundary_id() == extraction_boundary_id)
              {
                // Point fe_face_values at this cell face so quadrature_point(0)
                // and get_function_values() refer to the face center.
                fe_face_values.reinit(cell, f);

                const double y_face = fe_face_values.quadrature_point(0)[dim-1];
                // Skip faces outside the vertical extraction zone
                if (y_face < extraction_y_min || y_face > extraction_y_max)
                  continue;

                // Pack y-coordinate of face center
                local_packed.push_back(y_face);

                // Pack temperature at face center from the current FE solution
                std::vector<double> T_vals(1);
                fe_face_values[intr.extractors.temperature]
                  .get_function_values(this->get_solution(), T_vals);
                local_packed.push_back(T_vals[0]);

                // Pack each compositional field value at face center
                for (unsigned int c = 0; c < n_field; ++c)
                  {
                    std::vector<double> comp_vals(1);
                    fe_face_values[intr.extractors.compositional_fields[c]]
                      .get_function_values(this->get_solution(), comp_vals);
                    local_packed.push_back(comp_vals[0]);
                  }
              }

      // --- MPI gather: collect local_packed arrays from all processors ---
      // Each processor only owns some boundary faces; we need the full set on every
      // processor so they all integrate the same profile independently.

      const int local_size = static_cast<int>(local_packed.size());
      const unsigned int n_procs =
        Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

      // Step 1: tell every processor how many doubles each other processor has
      std::vector<int> all_sizes(n_procs);
      MPI_Allgather(&local_size, 1, MPI_INT,
                    all_sizes.data(), 1, MPI_INT,
                    this->get_mpi_communicator());

      // Step 2: compute displacements (byte offsets) into the receive buffer
      std::vector<int> displacements(n_procs, 0);
      for (unsigned int i = 1; i < n_procs; ++i)
        displacements[i] = displacements[i-1] + all_sizes[i-1];
      const int total_size = displacements.back() + all_sizes.back();

      // Step 3: all-to-all gather into all_packed (every processor gets the full array)
      std::vector<double> all_packed(total_size);
      MPI_Allgatherv(local_packed.data(), local_size, MPI_DOUBLE,
                     all_packed.data(), all_sizes.data(), displacements.data(),
                     MPI_DOUBLE, this->get_mpi_communicator());

      // --- Sort face centers from shallowest to deepest ---
      // Integration must proceed top-down; faces may arrive in arbitrary order
      // because different processors own different cells.
      const unsigned int n_total = static_cast<unsigned int>(total_size) / stride;
      std::vector<unsigned int> order(n_total);
      for (unsigned int i = 0; i < n_total; ++i) order[i] = i;
      // Sort descending by y (largest y = shallowest = smallest depth)
      std::sort(order.begin(), order.end(),
                [&](unsigned int a, unsigned int b)
                { return all_packed[a*stride] > all_packed[b*stride]; });

      // Unpack into sorted arrays.
      // face_center_depths[i] = y_top - sorted_y[i]  (depth below surface, increasing downward)
      std::vector<double> sorted_y(n_total);
      std::vector<double> sorted_T(n_total);
      std::vector<std::vector<double>> sorted_comp(n_field, std::vector<double>(n_total));
      lith_pressure.resize(n_total);
      face_center_depths.resize(n_total);
      for (unsigned int i = 0; i < n_total; ++i)
        {
          const unsigned int idx = order[i];
          sorted_y[i]           = all_packed[idx*stride + 0];
          sorted_T[i]           = all_packed[idx*stride + 1];
          face_center_depths[i] = y_top - sorted_y[i];   // depth from surface
          for (unsigned int c = 0; c < n_field; ++c)
            sorted_comp[c][i] = all_packed[idx*stride + 2 + c];
        }

      // --- Set up material model inputs/outputs (density only) ---
      typename MaterialModel::Interface<dim>::MaterialModelInputs  in(1, n_field);
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1, n_field);
      // Request only density; skips expensive viscosity/etc. calculations
      in.requested_properties = MaterialModel::MaterialProperties::density;
      in.velocity[0]          = Tensor<1,dim>();

      // x-coordinate of the profile column: on the boundary itself.
      // For a left boundary this is box origin x; for right it is origin x + width.
      const bool is_right = (extraction_boundary_id ==
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("right"));
      Point<dim> p;
      p[0] = is_right ? box->get_origin()[0] + box->get_extents()[0]
                      : box->get_origin()[0];

      // --- Trapezoidal integration from surface down through each face center ---
      // Following the same cumulative scheme as initial_lithostatic_pressure.cc:
      //   lith_pressure[i] = sum + 0.5 * rho[i] * g[i] * dz_i
      //   sum             +=       rho[i] * g[i] * dz_i
      // 'sum' carries the full integral up to the previous face; the 0.5 factor
      // adds a half-step at the current face, giving a midpoint-like correction.
      // Start with surface pressure (typically 0) as the integral at depth = 0.

      // Evaluate density at the shallowest face (index 0) FIRST,
      // because the lines below need out.densities[0] and g0.
      p[dim-1]          = sorted_y[0];
      in.position[0]    = p;
      in.temperature[0] = sorted_T[0];
      for (unsigned int c = 0; c < n_field; ++c)
        in.composition[0][c] = sorted_comp[c][0];
      in.pressure[0] = this->get_surface_pressure();   // initial pressure estimate
      this->get_material_model().evaluate(in, out);    // fills out.densities[0]
      const double g0 = this->get_gravity_model().gravity_vector(p).norm();

      // Initialize running sum to surface pressure, then compute lith_pressure[0]
      // as: P_surface + half the rho*g*dz contribution from surface to face 0.
      double sum = this->get_surface_pressure();
      lith_pressure[0] = sum + face_center_depths[0] * 0.5 * out.densities[0] * g0;
      sum             += face_center_depths[0] * out.densities[0] * g0;

      // Loop over remaining face centers, integrating downward
      for (unsigned int i = 1; i < n_total; ++i)
        {
          // Thickness of the layer between face i-1 and face i
          const double dz = face_center_depths[i] - face_center_depths[i-1];

          // Evaluate solid density at face i using T and composition from FE solution
          p[dim-1]          = sorted_y[i];
          in.position[0]    = p;
          in.temperature[0] = sorted_T[i];
          for (unsigned int c = 0; c < n_field; ++c)
            in.composition[0][c] = sorted_comp[c][i];
          in.pressure[0] = lith_pressure[i-1];   // use previous face pressure as estimate
          this->get_material_model().evaluate(in, out);
          const double g   = this->get_gravity_model().gravity_vector(p).norm();

          // Add half-step at current face, then record full step in sum
          lith_pressure[i] = sum + dz * 0.5 * out.densities[0] * g;
          sum             += dz * out.densities[0] * g;
        }
    }



    template <int dim>
    double
    Density<dim>::interpolate_lith_pressure (const double depth) const
    {
      if (depth <= face_center_depths.front())
        return lith_pressure.front();
      if (depth >= face_center_depths.back())
        return lith_pressure.back();

      // Binary search for the interval containing depth
      const auto it = std::lower_bound(face_center_depths.begin(),
                                        face_center_depths.end(), depth);
      const unsigned int i = static_cast<unsigned int>(
        std::distance(face_center_depths.begin(), it));

      const double d0   = face_center_depths[i-1];
      const double d1   = face_center_depths[i];
      const double frac = (depth - d0) / (d1 - d0);
      return (1.0 - frac) * lith_pressure[i-1] + frac * lith_pressure[i];
    }


    template <int dim>
    void
    Density<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary fluid pressure model");
      {
        prm.enter_subsection("Density");
        {
          prm.declare_entry ("Density formulation", "solid density",
                            Patterns::Selection ("solid density|fluid density|average density|weighted density|side boundary magma extraction"),
                             "The density formulation used to compute the fluid pressure gradient "
                             "at the model boundary."
                             "\n\n"
                             "`solid density' prescribes the gradient of the fluid pressure as "
                             "solid density times gravity (which is the lithostatic "
                             "pressure) and leads to approximately the same pressure in "
                             "the melt as in the solid, so that fluid is only flowing "
                             "in or out due to differences in dynamic pressure."
                             "\n\n"
                             "`fluid density' prescribes the gradient of the fluid pressure as "
                             "fluid density times gravity and causes melt to flow in "
                             "with the same velocity as inflowing solid material, "
                             "or no melt flowing in or out if the solid velocity "
                             "normal to the boundary is zero."
                             "\n\n"
                             "'average density' prescribes the gradient of the fluid pressure as "
                             "the averaged fluid and solid density times gravity "
                             "(which is a better approximation for the lithostatic "
                             "pressure than just the solid density) and leads to approximately the same pressure in "
                             "the melt as in the solid, so that fluid is only flowing "
                                                          "in or out due to differences in dynamic pressure."
                             "\n\n"
                             "'weighted density' prescribes the gradient of the fluid pressure as "
                             "the a weighted average between fluid and solid density times gravity. "
                             "This means that melt can still flow in or out due to differences in "
                             "dynamic pressure, but this flow is impeded depending on the weights "
                             "of the fluid and solid densities as described by the `Fluid density weight'"
                             "parameter."
                             "\n\n"
                             "'side boundary magma extraction' prescribes (grad p_f)·n = "
                             "W*(P_lith - rho_f*g*z)/dx on a user-defined portion of a side boundary, "
                             "where P_lith is the current lithostatic pressure at the boundary, and "
                             "W is the side pressure gradient weight.");

          prm.declare_entry ("Fluid density weight", "1.0",
                             Patterns::Double (0,1),
                             "The value $w$ used to weight the fluid and solid densities in the "
                             "`weighted density' boundary condition. The boundary fluid pressure is "
                             "computed as $w rho\\_{f} + (1-w) rho\\_{s} g$, where $rho\\_{f}$ is the "
                             "fluid density, $rho\\_{s}$ is the solid density and $g$ is the gravity. "
                             "That means that if $w=1$, this formulation is equivalent to the "
                             "`fluid density' formulation, and if $w=0$, this formulation is equivalent "
                             "to the `solid density' formulation. "
                             "This parameter is only used if `weighted density' is selected as "
                             "'Density formulation'.");
          prm.enter_subsection ("Side boundary magma extraction");
          {
            prm.declare_entry ("Extraction boundary", "left",
                               Patterns::Anything(),
                               "Symbolic name of the extraction boundary (e.g. 'left' or 'right').");
            prm.declare_entry ("Extraction y min", "0",
                               Patterns::Double(),
                               "Minimum y-coordinate of the extraction zone. Units: m.");
            prm.declare_entry ("Extraction y max", "1e30",
                               Patterns::Double(),
                               "Maximum y-coordinate of the extraction zone. Units: m.");
            prm.declare_entry ("Extraction horizontal distance scale", "1e3",
                               Patterns::Double(0),
                             "Horizontal distance dx in the formula (P_lith - rho_f*g*z)/dx. , should be ~half the grid size along the boundary Units: m.");
             prm.declare_entry ("Side pressure gradient weight", "0.05",
                               Patterns::Double(0,1),
                               "Fluid pressure gradient on the side is this weight times the difference between"
                               "the lithostatic and magma static pressure/dx");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Density<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary fluid pressure model");
      {
        prm.enter_subsection("Density");
        {
          if (prm.get ("Density formulation") == "solid density")
            density_formulation = DensityFormulation::solid_density;
          else if (prm.get ("Density formulation") == "fluid density")
            density_formulation = DensityFormulation::fluid_density;
          else if (prm.get ("Density formulation") == "average density")
            density_formulation = DensityFormulation::average_density;
          else if (prm.get ("Density formulation") == "weighted density")
          {
            density_formulation = DensityFormulation::weighted_density;
            fluid_density_weight = prm.get_double ("Fluid density weight");
          }
          else if (prm.get ("Density formulation") == "side boundary magma extraction")
            {
              density_formulation = DensityFormulation::side_boundary_magma_extraction;
              prm.enter_subsection ("Side boundary magma extraction");
              {
                extraction_boundary_name      = prm.get        ("Extraction boundary");
                extraction_y_min              = prm.get_double ("Extraction y min");
                extraction_y_max              = prm.get_double ("Extraction y max");
                extraction_dx                 = prm.get_double ("Extraction horizontal distance scale");
                side_pressure_gradient_weight = prm.get_double ("Side pressure gradient weight");
              }
              prm.leave_subsection ();
            }

          else
            AssertThrow (false, ExcNotImplemented());
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryFluidPressure
  {
    ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(Density,
                                                  "density",
                                                  "A plugin that prescribes the fluid pressure gradient at "
                                                  "the boundary based on fluid/solid density from the material "
                                                  "model.")
  }
}
