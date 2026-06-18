/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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

#include <aspect/postprocess/visualization/melting_rate.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      MeltingRate<dim>::MeltingRate()
        :
        DataPostprocessorScalar<dim>("melting_rate", update_values),
        Interface<dim>("1/s")
      {}

      template <int dim>
      void
      MeltingRate<dim>::evaluate_vector_field(
        const DataPostprocessorInputs::Vector<dim> &input_data,
        std::vector<Vector<double>>                &computed_quantities) const
      {
        const unsigned int porosity_component =
          this->introspection().component_indices.compositional_fields
          [this->introspection().compositional_index_for_name("porosity")];

        const double dt = this->get_timestep();

        for (unsigned int q = 0; q < input_data.solution_values.size(); ++q)
          computed_quantities[q](0) = (dt > 0)
                                      ? input_data.solution_values[q][porosity_component] / dt
                                      : 0.0;
      }

      template <int dim>
      bool
      MeltingRate<dim>::needs_reaction_vector() const
      {
        return true;
      }
    }
  }
}

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(
        MeltingRate,
        "melting rate",
        "A visualization output object that outputs the local rate of "
        "change of melting as (d(porosity)/dt in units of 1/s) due to "
        "melting and freezing reactions computed during operator splitting. "
        "Requires operator splitting and melt transport to be enabled. "
        "**Note currently dporosity is erroneously a mass fraction not volume fraction**")
    }
  }
}
