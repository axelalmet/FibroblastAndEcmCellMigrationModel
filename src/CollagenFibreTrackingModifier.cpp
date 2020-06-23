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

#include "CollagenFibreTrackingModifier.hpp"
#include "CollagenCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
CollagenFibreTrackingModifier<DIM>::CollagenFibreTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CollagenFibreTrackingModifier<DIM>::~CollagenFibreTrackingModifier()
{
}

template<unsigned DIM>
void CollagenFibreTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    rCellPopulation.Update();
}

template<unsigned DIM>
void CollagenFibreTrackingModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double time = SimulationTime::Instance()->GetTime();
    *mpMarkedSpringsFile << time << " ";

    UpdateCellData(rCellPopulation);

    *mpMarkedSpringsFile << "\n";
}

template<unsigned DIM>
void CollagenFibreTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{

    // Create output file
    OutputFileHandler output_file_handler(outputDirectory + "/", false);
    mpMarkedSpringsFile = output_file_handler.OpenOutputFile("collagenfibrelocations.dat");

    double time = SimulationTime::Instance()->GetTime();

    *mpMarkedSpringsFile << time << " ";

    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);

    *mpMarkedSpringsFile << "\n";
}

template<unsigned DIM>
void CollagenFibreTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // This really only works with node-based cell populations
    assert(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation));

	NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

    std::vector< std::pair<Node<DIM>*, Node<DIM>* > > node_pairs = p_cell_population->rGetNodePairs(); // Get the possible node pairs

    // Iterate over the node pairs
    for (unsigned i = 0; i < node_pairs.size(); i++)
    {
        std::pair<Node<DIM>*, Node<DIM>* > node_pair = node_pairs[i]; // Get the node pair

        Node<DIM>* p_node_A = node_pair.first; // First node
        Node<DIM>* p_node_B = node_pair.second; // Second node

        // Use the node indices to get the cells and hence, the cell types
        unsigned node_index_A = p_node_A->GetIndex(); 
        unsigned node_index_B = p_node_B->GetIndex(); 

        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(node_index_A);
        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(node_index_B);

                // Only record collagen-collagen pairs
        if ( (p_cell_A->GetCellProliferativeType()->IsType<CollagenCellProliferativeType>())
            &&(p_cell_B->GetCellProliferativeType()->IsType<CollagenCellProliferativeType>()) )
        {

            std::pair<CellPtr, CellPtr> cell_pair = p_cell_population->CreateCellPair(p_cell_A, p_cell_B);

            if (p_cell_population->IsMarkedSpring(cell_pair))
            {
                c_vector<double, DIM> cell_location_A = p_node_A->rGetLocation();
                c_vector<double, DIM> cell_location_B = p_node_B->rGetLocation();

                *mpMarkedSpringsFile << node_index_A << " " << node_index_B << " " 
                                        << cell_location_A[0] << " " << cell_location_A[1] << " "
                                        << cell_location_B[0] << " " << cell_location_B[1];
            }

        }

    }

}

template<unsigned DIM>
void CollagenFibreTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CollagenFibreTrackingModifier<1>;
template class CollagenFibreTrackingModifier<2>;
template class CollagenFibreTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CollagenFibreTrackingModifier)
