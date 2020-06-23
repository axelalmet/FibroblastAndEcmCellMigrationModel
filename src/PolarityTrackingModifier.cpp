/*

Copyright (c) 2005-2020, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "PolarityTrackingModifier.hpp"
#include "CollagenCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned DIM>
PolarityTrackingModifier<DIM>::PolarityTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mReorientationStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
PolarityTrackingModifier<DIM>::~PolarityTrackingModifier()
{
}

template<unsigned DIM>
double PolarityTrackingModifier<DIM>::GetReorientationStrength()
{
    return mReorientationStrength;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetReorientationStrength(double reorientationStrength)
{
    mReorientationStrength = reorientationStrength;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    // UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Only consider non-collagen cells
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType();

        if (!p_cell_type->template IsType<CollagenCellProliferativeType>())
        {

            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get node index
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index); // Get node

            c_vector<double, DIM> velocity = p_node->rGetAppliedForce(); // The applied force is proportional to the vector

            // Determine the angle the velocity direction makes with the x-axis
            double velocity_angle = atan(velocity[1]/velocity[0]);

            // Adjust for the relevant quadrants
            if (velocity[0] < 0.0) // Quadrant II or III
            {
                velocity_angle += M_PI;
            }
            else
            {
                if (velocity[1] < 0.0) // Quadrant IV
                {
                    velocity_angle += 2.0 * M_PI;
                }
            }

            double phi = cell_iter->GetCellData()->GetItem("direction"); // Get the current cell direction
            double dt = SimulationTime::Instance()->GetTimeStep(); // Get dt
            phi += mReorientationStrength * dt * sin(velocity_angle - phi); // Update the cell polarity
            
            cell_iter->GetCellData()->SetItem("direction", phi);

        }
    }

}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityTrackingModifier<1>;
template class PolarityTrackingModifier<2>;
template class PolarityTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityTrackingModifier)
