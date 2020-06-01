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

#include "FibreAlignmentBasedMigrationForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "CollagenCellProliferativeType.hpp"

template<unsigned DIM>
FibreAlignmentBasedMigrationForce<DIM>::FibreAlignmentBasedMigrationForce()
    : AbstractForce<DIM>(),
    mReorientationStrength(DOUBLE_UNSET),
    mMigrationForceStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
FibreAlignmentBasedMigrationForce<DIM>::~FibreAlignmentBasedMigrationForce()
{
}

// Method to get the strength of fibre alignment
template<unsigned DIM>
double FibreAlignmentBasedMigrationForce<DIM>::GetReorientationStrength()
{
    return mReorientationStrength;
}

// Method to set the strength of fibre alignment
template<unsigned DIM>
void FibreAlignmentBasedMigrationForce<DIM>::SetReorientationStrength(double reorientationStrength)
{
    mReorientationStrength = reorientationStrength;
}

// Method to get the strength of the migration force
template<unsigned DIM>
double FibreAlignmentBasedMigrationForce<DIM>::GetMigrationForceStrength()
{
    return mMigrationForceStrength;
}

// Method to set the strength of migration force
template<unsigned DIM>
void FibreAlignmentBasedMigrationForce<DIM>::SetMigrationForceStrength(double migrationForceStrength)
{
    mMigrationForceStrength = migrationForceStrength;
}

template<unsigned DIM>
bool FibreAlignmentBasedMigrationForce<DIM>::DoesCollagenFibreIntersectWithFibroblast(AbstractCellPopulation<DIM>& rCellPopulation, unsigned fibroblastIndex, unsigned fibreIndex)
{
    bool does_fibre_intersect_with_fibroblast = false;

    // Get the considered fibroblast node's location and orientation
    c_vector<double, 2> fibroblast_location = rCellPopulation.GetNode(fibroblastIndex)->rGetLocation();
    double fibroblast_direction = rCellPopulation.GetCellUsingLocationIndex(fibroblastIndex)->GetCellData()->GetItem("direction");


    // Get the neighbouring fibre node's location and orientation
    c_vector<double, 2> fibre_location = rCellPopulation.GetNode(fibreIndex)->rGetLocation();
    double fibre_orientation = rCellPopulation.GetCellUsingLocationIndex(fibreIndex)->GetCellData()->GetItem("orientation");

    // Determine the intersecting points, i.e. determine the values of s and t such that the lines intersect.
    double s = (fibre_location[1] - fibroblast_location[1] - (fibre_location[0] - fibroblast_location[0]))/(cos(fibre_orientation) - sin(fibre_orientation));
    double t = (cos(fibre_orientation)*(fibre_location[1] - fibroblast_location[1]) - sin(fibre_orientation)*(fibre_location[0] - fibroblast_location[0]))/(cos(fibroblast_direction)*(cos(fibre_orientation) - sin(fibre_orientation)));

    if ( (s > 1e-4)&&(t > 1e-4) ) // If s and t are both positive, then the fibres intersect.
    {
        does_fibre_intersect_with_fibroblast = true;
    }

    return does_fibre_intersect_with_fibroblast;

}


template<unsigned DIM>
void FibreAlignmentBasedMigrationForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    // For now, this will only apply to node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulation<DIM>* p_tissue = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Get the force parameters
    double reorientation_strength = GetReorientationStrength(); // How strongly fibroblasts align with fibres
    double migration_force_strength = GetMigrationForceStrength(); // Strength of migration force

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        // Only apply this migration force to fibroblasts
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType(); // Get the cell type

        if (p_cell_type->IsType<FibroblastCellProliferativeType>())
        {
            // Get the node index
            unsigned current_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // Get the location of the node
            c_vector<double, DIM> current_location = rCellPopulation.GetNode(current_index)->rGetLocation();

            // Get the neighbouring node indices (fibres span about 0.5um at most, so don't need to set neighbourhood radius)
            std::set<unsigned> neighbouring_indices = p_tissue->GetNeighbouringNodeIndices(current_index);

            // Get the fibroblast migration direction
            double fibroblast_direction = cell_iter->GetCellData()->GetItem("direction");

            // If there are any neighbours
            if (!neighbouring_indices.empty())
            {
                // Iterate over these elements
                for (std::set<unsigned>::iterator elem_iter = neighbouring_indices.begin();
                        elem_iter != neighbouring_indices.end();
                        ++elem_iter)
                {
                    // Get the cell according to the index
                    CellPtr neighbour_cell_iter = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

                    // Only consider fibroblast neighbours
                    boost::shared_ptr<AbstractCellProperty> p_neighbour_cell_type = neighbour_cell_iter->GetCellProliferativeType(); // Get the cell type

                    // Fibres are defined as spanning from collagen neighbours that express 
                    // a positive amount of collagen
                    if (p_neighbour_cell_type->IsType<CollagenCellProliferativeType>())
                    {
                        // Get the neighbouring concentration
                        double collagen = neighbour_cell_iter->GetCellData()->GetItem("collagen");

                        if (collagen > 1e-4)
                        {
                            // Now determine if the fibre intersects with the fibroblast direction
                            bool does_fibre_intersect_with_fibroblast = DoesCollagenFibreIntersectWithFibroblast(rCellPopulation, current_index, *elem_iter);
                            if (does_fibre_intersect_with_fibroblast)
                            {
                                // Get the fibre orientation
                                double fibre_orientation = neighbour_cell_iter->GetCellData()->GetItem("orientation");

                                // We need to adjust the fibre orientation so that the amount of alignment of the 
                                // fibroblast to the collagen fibre is minimised

                                // Note that we define fibre orientations to be only within the 1st or 4th quadrant, as direction does
                                // not matter for collagen fibres. However, as the direction matters for migration, we need to account
                                // for this.

                                if (fibre_orientation >= 0.0 ) // 1st quadrant
                                {
                                    // This is to determine whether the fibroblast better aligns with the fibre in the 1st or 3rd quadrant
                                    if (fabs(fibre_orientation - fibroblast_direction) < fabs(fibre_orientation + M_PI - fibroblast_direction)) // 1st quadrant
                                    {
                                        fibroblast_direction += reorientation_strength * sin(fibre_orientation - fibroblast_direction); // Update fibroblast direction
                                    }
                                    else // 3rd quadrant
                                    {
                                        fibroblast_direction += reorientation_strength * sin(fibre_orientation + M_PI - fibroblast_direction); // Update fibroblast direction
                                    } 
                                }
                                else // 4th quadrant
                                {
                                    // This is to determine whether the fibroblast better aligns with the fibre in the 2nd or 4th quadrant
                                    if (fabs(fibre_orientation - fibroblast_direction) < fabs(fibre_orientation + M_PI - fibroblast_direction)) // 4th quadrant
                                    {
                                        fibroblast_direction += reorientation_strength * sin(fibre_orientation - fibroblast_direction); // Update fibroblast direction
                                    }
                                    else // 
                                    {
                                        fibroblast_direction += reorientation_strength * sin(fibre_orientation + M_PI - fibroblast_direction); // Update fibroblast direction
                                    } 
                                }
                            }
                        }
                    }
                }
            }

        // Update the new fibroblast direction
        cell_iter->GetCellData()->SetItem("direction", fibroblast_direction);

        // Add the force contribution. I'm well aware that this is a 2D vector, despite the templated force direction, but oh well.
        c_vector<double, 2> force_direction; // Define the force direction using the fibroblast direction
        force_direction[0] = cos(fibroblast_direction);
        force_direction[1] = sin(fibroblast_direction);

        rCellPopulation.GetNode(current_index)->AddAppliedForceContribution(migration_force_strength * force_direction);

        }
    }
}

template<unsigned DIM>
void FibreAlignmentBasedMigrationForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
	*rParamsFile <<  "\t\t\t<ReorientationStrength>"<<  mReorientationStrength << "</ReorientationStrength> \n" ;
	*rParamsFile <<  "\t\t\t<MigrationForceStrength>"<<  mMigrationForceStrength << "</MigrationForceStrength> \n" ;

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class FibreAlignmentBasedMigrationForce<1>;
template class FibreAlignmentBasedMigrationForce<2>;
template class FibreAlignmentBasedMigrationForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FibreAlignmentBasedMigrationForce)
