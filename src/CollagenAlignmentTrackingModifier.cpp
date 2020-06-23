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

#include "CollagenAlignmentTrackingModifier.hpp"
#include "CollagenCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
CollagenAlignmentTrackingModifier<DIM>::CollagenAlignmentTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mNeighbourhoodRadius(1.5),
      mReorientationStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
CollagenAlignmentTrackingModifier<DIM>::~CollagenAlignmentTrackingModifier()
{
}

template<unsigned DIM>
double CollagenAlignmentTrackingModifier<DIM>::GetNeighbourhoodRadius()
{
    return mNeighbourhoodRadius;
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::SetNeighbourhoodRadius(double neighbourhoodRadius)
{
    mNeighbourhoodRadius = neighbourhoodRadius;
}

template<unsigned DIM>
double CollagenAlignmentTrackingModifier<DIM>::GetReorientationStrength()
{
    return mReorientationStrength;
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::SetReorientationStrength(double reorientationStrength)
{
    mReorientationStrength = reorientationStrength;
}

    /*
     * This method is used to determine whether two fibres may intersect.
     * The way we determine this is by using the collagen orientations, theta_i.
     * 
     * For a cell i's position, given by r_i = (x_i, y_i), we may parametrise a line
     * using the collagen fibre orientation by 
     * 
     * r_i = (x_i, y_i) + t*(cos(theta_i), sin(theta_i)), where t is a real number.
     * 
     * Similarly, for the fibroblast neighbour j, we can write
     * 
     * r_j = (x_j, y_j) + s*(cos(theta_j), sin(theta_j)), where s is a real number.
     * 
     * We can write down the expressions for s and t to determine where they intersect. If
     * both s > 0 and t >0, we say the fibres intersect.
     */ 

template<unsigned DIM>
bool CollagenAlignmentTrackingModifier<DIM>::DoesCollagenFibreIntersectWithFibroblast(AbstractCellPopulation<DIM,DIM>& rCellPopulation, unsigned nodeIndex, unsigned fibroblastIndex)
{
    bool does_fibre_intersect_with_fibroblast = false;

    // Get the considered collagen node's location and orientation
    c_vector<double, 2> current_location = rCellPopulation.GetNode(nodeIndex)->rGetLocation();
    double current_orientation = rCellPopulation.GetCellUsingLocationIndex(nodeIndex)->GetCellData()->GetItem("orientation");

    // Get the neighbouring node's location and orientation
    c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(fibroblastIndex)->rGetLocation();
    double neighbour_orientation = rCellPopulation.GetCellUsingLocationIndex(fibroblastIndex)->GetCellData()->GetItem("direction");

    // Determine the intersecting points, i.e. determine the values of s and t such that the lines intersect.
    double s = (neighbour_location[1] - current_location[1] - (neighbour_location[0] - current_location[0]))/(cos(neighbour_orientation) - sin(neighbour_orientation));
    double t = (cos(neighbour_orientation)*(neighbour_location[1] - current_location[1]) - sin(neighbour_orientation)*(neighbour_location[0] - current_location[0]))/(cos(current_orientation)*(cos(neighbour_orientation) - sin(neighbour_orientation)));

    if ( (s > 0.0)&&(t > 0.0) ) // If s and t are both positive, then the fibres intersect.
    {
        does_fibre_intersect_with_fibroblast = true;
    }

    return does_fibre_intersect_with_fibroblast;

}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // This really only works with node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

    // double neighbourhood_radius = GetNeighbourhoodRadius();

    // We will consider fibres that span from fibroblasts. We define fibres from nodes which
    // express a non-zero amount of collagen (whether or not we then define fibre length based
    // on collagen expression is a future question).
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Only consider collagen cells
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType();

        if (p_cell_type->template IsType<CollagenCellProliferativeType>())
        {

            // Get the collagen amount
            double collagen = cell_iter->GetCellData()->GetItem("collagen");

            if (collagen > 1e-4) // Realise how bad it is to use 0 as a point of comparison
            {
                // Get the node index
                unsigned node_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter);
                    
                // Get the set of neighbouring location indices within a neighbourhood radius.
                // N.B. We may change the neighbourhood radius to reflect the length of the
                // collagen fibre.
                std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringNodeIndices(node_index);

                // We only update the orientation if, well, we can.
                if (!neighbour_indices.empty())
                {
                    // Get the current orientation
                    double current_orientation = cell_iter->GetCellData()->GetItem("orientation");
                    double reorientation_strength = GetReorientationStrength();

                    for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                        iter != neighbour_indices.end();
                        ++iter)
                    {
                        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                        boost::shared_ptr<AbstractCellProperty> p_neighbour_cell_type = p_cell->GetCellProliferativeType();

                        // Only look at fibroblast neighbours.
                        if (p_neighbour_cell_type->template IsType<FibroblastCellProliferativeType>())
                        {

                            // Determine if the collagen fibre intersects with the (migrating) fibroblast
                            bool does_fibre_intersect_with_fibroblast = DoesCollagenFibreIntersectWithFibroblast(rCellPopulation, node_index, *iter);

                            if (does_fibre_intersect_with_fibroblast) // If the fibres intersect, we orient the fibres with respect to the neighbouring orientation
                            {
                                double neighbour_orientation = p_cell->GetCellData()->GetItem("direction");

                                // Adjust the fibroblast orientation so that it lies between -pi/2 and pi/2
                                if ((neighbour_orientation > 0.5*M_PI)&&(neighbour_orientation <= M_PI) ) // 2nd quadrant
                                {
                                    neighbour_orientation -= M_PI;
                                }
                                else if ( (neighbour_orientation > M_PI)&&(neighbour_orientation <= 1.5*M_PI) ) // 3rd quadrant
                                {
                                    neighbour_orientation -= M_PI;
                                }
                                else if ( (neighbour_orientation > 1.5*M_PI)&&(neighbour_orientation <= 2.0*M_PI) ) // 4th quadrant
                                {
                                    neighbour_orientation -= 2.0*M_PI;
                                }

                                double dt = SimulationTime::Instance()->GetTimeStep(); // Get dt
                                current_orientation += reorientation_strength * dt * sin(neighbour_orientation - current_orientation);
                            }

                        }
                    }

                    // Update the fibre orientations                    
                    cell_iter->GetCellData()->SetItem("orientation", current_orientation); // Update scalar data
                }
            }
        }
    }

}

template<unsigned DIM>
void CollagenAlignmentTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CollagenAlignmentTrackingModifier<1>;
template class CollagenAlignmentTrackingModifier<2>;
template class CollagenAlignmentTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CollagenAlignmentTrackingModifier)
