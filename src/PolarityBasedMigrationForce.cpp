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

#include "PolarityBasedMigrationForce.hpp"
#include "FibroblastCellProliferativeType.hpp"

template<unsigned DIM>
PolarityBasedMigrationForce<DIM>::PolarityBasedMigrationForce()
    : AbstractForce<DIM>(),
    mMigrationForceStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
PolarityBasedMigrationForce<DIM>::~PolarityBasedMigrationForce()
{
}

// Method to get the strength of the migration force
template<unsigned DIM>
double PolarityBasedMigrationForce<DIM>::GetMigrationForceStrength()
{
    return mMigrationForceStrength;
}

// Method to set the strength of migration force
template<unsigned DIM>
void PolarityBasedMigrationForce<DIM>::SetMigrationForceStrength(double migrationForceStrength)
{
    mMigrationForceStrength = migrationForceStrength;
}

template<unsigned DIM>
void PolarityBasedMigrationForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

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
            // Get the fibroblast migration direction
            double fibroblast_direction = cell_iter->GetCellData()->GetItem("direction");

            // Add the force contribution. I'm well aware that this is a 2D vector, despite the templated force direction, but oh well.
            c_vector<double, 2> force_direction; // Define the force direction using the fibroblast direction
            force_direction[0] = cos(fibroblast_direction);
            force_direction[1] = sin(fibroblast_direction);

            rCellPopulation.GetNode(current_index)->AddAppliedForceContribution(migration_force_strength * force_direction);

        }
    }
}

template<unsigned DIM>
void PolarityBasedMigrationForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{

	*rParamsFile <<  "\t\t\t<MigrationForceStrength>"<<  mMigrationForceStrength << "</MigrationForceStrength> \n" ;

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityBasedMigrationForce<1>;
template class PolarityBasedMigrationForce<2>;
template class PolarityBasedMigrationForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityBasedMigrationForce)
