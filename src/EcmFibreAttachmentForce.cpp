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

#include "EcmFibreAttachmentForce.hpp"
#include "CollagenCellProliferativeType.hpp"
#include "FibroblastCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
EcmFibreAttachmentForce<DIM>::EcmFibreAttachmentForce()
    : AbstractForce<DIM>(),
    mFibreAdhesionStrength(DOUBLE_UNSET),
    mFibreAlignmentStrength(DOUBLE_UNSET),
    mCollagenFibreLength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
EcmFibreAttachmentForce<DIM>::~EcmFibreAttachmentForce()
{
}

// Method to get the strength of the adhesion force
template<unsigned DIM>
double EcmFibreAttachmentForce<DIM>::GetFibreAdhesionStrength()
{
    return mFibreAdhesionStrength;
}

// Method to set the strength of adhesion force
template<unsigned DIM>
void EcmFibreAttachmentForce<DIM>::SetFibreAdhesionStrength(double fibreAdhesionStrength)
{
    mFibreAdhesionStrength = fibreAdhesionStrength;
}

// Method to get the strength of the adhesion force
template<unsigned DIM>
double EcmFibreAttachmentForce<DIM>::GetFibreAlignmentStrength()
{
    return mFibreAlignmentStrength;
}

// Method to set the strength of adhesion force
template<unsigned DIM>
void EcmFibreAttachmentForce<DIM>::SetFibreAlignmentStrength(double fibreAlignmentStrength)
{
    mFibreAlignmentStrength = fibreAlignmentStrength;
}

// Method to get the collagen fibre search radius
template<unsigned DIM>
double EcmFibreAttachmentForce<DIM>::GetCollagenFibreLength()
{
    return mCollagenFibreLength;
}

// Method to set the strength of adhesion force
template<unsigned DIM>
void EcmFibreAttachmentForce<DIM>::SetCollagenFibreLength(double collagenFibreLength)
{
    mCollagenFibreLength = collagenFibreLength;
}

/*
* Get the relevant ECM fibres to consider for the specified node.
*/ 
template<unsigned DIM>
std::vector<std::pair<unsigned, unsigned> > EcmFibreAttachmentForce<DIM>::GetNeighbouringEcmFibres(unsigned nodeIndex, AbstractCellPopulation<DIM>& rCellPopulation)
{
    std::vector<std::pair<unsigned, unsigned> > neighbouring_fibres; // Initialise fibres

    // We need to static-cast the population (assume node-based)
    // to call upon a node-based specific method
    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

        double node_radius = rCellPopulation.GetNode(nodeIndex)->GetRadius(); // Get the node radius to determine fibre interactions
        double phi = rCellPopulation.GetCellUsingLocationIndex(nodeIndex)->GetCellData()->GetItem("direction");

        // Define the cell directions
        c_vector<double, 2> node_direction;
        node_direction[0] = cos(phi);
        node_direction[1] = sin(phi);

        std::vector<unsigned> neighbouring_collagen_indices; // Initialise collagen neighbours

        std::set<unsigned> neighbouring_indices = rCellPopulation.GetNeighbouringNodeIndices(nodeIndex);

        // Iterate through the node indices and only consider collagen nodes.
        for (std::set<unsigned>::iterator elem_iter = neighbouring_indices.begin();
                elem_iter != neighbouring_indices.end();
                ++elem_iter)
        {
            //Get the cell according to the index
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

            //Get the cell type
            boost::shared_ptr<AbstractCellProperty> p_type = p_cell->GetCellProliferativeType();

            if(p_type->IsType<CollagenCellProliferativeType>())
            {
                neighbouring_collagen_indices.push_back(*elem_iter);
            }
        }

        // Now go through each collagen node and determine all the possible fibres to consider.
        for (unsigned i = 0; i < neighbouring_collagen_indices.size(); i++)
        {
            unsigned collagen_index = neighbouring_collagen_indices[i];

            CellPtr collagen_cell = rCellPopulation.GetCellUsingLocationIndex(collagen_index);

            // Now get the neighbouring nodes of this collagen node
            std::set<unsigned> neighbours_of_collagen_node = p_cell_population->GetNodesWithinNeighbourhoodRadius(collagen_index, mCollagenFibreLength);

            for (std::set<unsigned>::iterator elem_iter = neighbours_of_collagen_node.begin();
            elem_iter != neighbours_of_collagen_node.end();
            ++elem_iter)
            {
                CellPtr neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

                //Get the cell type
                boost::shared_ptr<AbstractCellProperty> p_type = neighbour_cell->GetCellProliferativeType();

                if(p_type->IsType<CollagenCellProliferativeType>())
                {
                    std::pair<CellPtr, CellPtr> cell_pair = p_cell_population->CreateCellPair(collagen_cell, neighbour_cell);
                    std::pair<unsigned, unsigned> considered_fibre = std::make_pair(collagen_index, *elem_iter);

                    if (p_cell_population->IsMarkedSpring(cell_pair))
                    {
                        c_vector<double, DIM> orthog_proj = CalculateOrthogonalProjectionOntoFibre(nodeIndex, considered_fibre, rCellPopulation);
                        orthog_proj *= -1.0; // Reverse the direction, as this would be teh actual force direction applied to the cell

                        // If the fibre is close enough, at least orthogonal to the cell, and hasn't been collected, add the fibre to the vector
                        if ( (norm_2(orthog_proj) <= node_radius) // If the fibre is sufficiently close
                                &&(inner_prod(orthog_proj, node_direction) >= 0.0) // If the fibre is at least orthogonal to the cell
                                &&(!HasFibreBeenConsidered(considered_fibre, neighbouring_fibres)) )
                        {
                            neighbouring_fibres.push_back(considered_fibre);
                        }
                    }
                }
            }
        }
    }

    return neighbouring_fibres;
}

/*
* Calculate the orthogonal projection which we will use to calculate the adhesion force strength, as well
* as determine whether or not the cell will interact with this fibre.
*/
template<unsigned DIM>
c_vector<double, DIM> EcmFibreAttachmentForce<DIM>::CalculateOrthogonalProjectionOntoFibre(unsigned nodeIndex, std::pair<unsigned, unsigned> fibre, AbstractCellPopulation<DIM>& rCellPopulation)
{

    c_vector<double, DIM> cell_location = rCellPopulation.GetNode(nodeIndex)->rGetLocation(); // Get the cell location

    // Get the collagen nodes and their locations
    unsigned collagen_node_A = fibre.first;
    unsigned collagen_node_B = fibre.second;

    c_vector<double, DIM> collagen_location_A = rCellPopulation.GetNode(collagen_node_A)->rGetLocation();
    c_vector<double, DIM> collagen_location_B = rCellPopulation.GetNode(collagen_node_B)->rGetLocation();

    // We can define the collagen fibre direction now
    c_vector<double, DIM> fibre_A_to_B = rCellPopulation.rGetMesh().GetVectorFromAtoB(collagen_location_A, collagen_location_B);
    fibre_A_to_B /= norm_2(fibre_A_to_B); // Normalise the fibre vector (just need the direction)

    // Define the vector from node A to the cell
    c_vector<double, DIM> collagen_A_to_cell = rCellPopulation.rGetMesh().GetVectorFromAtoB(collagen_location_A, cell_location);

    return collagen_A_to_cell - inner_prod(collagen_A_to_cell, fibre_A_to_B) * fibre_A_to_B;

}

/*
* Condition to check if an ECM fibre has been considered or not
*/
template<unsigned DIM>
bool EcmFibreAttachmentForce<DIM>::HasFibreBeenConsidered(std::pair<unsigned, unsigned> fibre, std::vector<std::pair<unsigned, unsigned> > neighbouringFibres)
{
    bool has_fibre_been_considered = false;

    for (unsigned i = 0; i < neighbouringFibres.size(); i++)
    {
        std::pair<unsigned, unsigned> neighbouring_fibre = neighbouringFibres[i];

        if ( ( (fibre.first == neighbouring_fibre.first)&&(fibre.second == neighbouring_fibre.second) )
            ||( (fibre.first == neighbouring_fibre.second)&&(fibre.second == neighbouring_fibre.first) ) )
        {
            has_fibre_been_considered = true;
            break;
        }
    }

    return has_fibre_been_considered;
}

template<unsigned DIM>
void EcmFibreAttachmentForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

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
            Node<DIM>* p_node = rCellPopulation.GetNode(current_index);

            std::vector<std::pair<unsigned, unsigned> > neighbouring_collagen_fibres = GetNeighbouringEcmFibres(current_index, rCellPopulation);

            // Get the fibroblast migration direction
            double phi = cell_iter->GetCellData()->GetItem("direction"); // Get the current cell direction
            double dt = SimulationTime::Instance()->GetTimeStep(); // Get dt

            double phi_new = phi;

            for (unsigned i = 0; i < neighbouring_collagen_fibres.size(); i++)
            {
                std::pair<unsigned, unsigned> collagen_fibre = neighbouring_collagen_fibres[i];

                // Get the orthogonal projection
                c_vector<double, DIM> orthog_proj = CalculateOrthogonalProjectionOntoFibre(current_index, collagen_fibre, rCellPopulation);
                double distance_to_fibre = norm_2(orthog_proj);

                c_vector<double, DIM> force_direction = -1.0 * orthog_proj / distance_to_fibre; // The force direction is the reverse of the orthog

                // Calculate the force
                double alpha = 5.0;
                c_vector<double, DIM> force = mFibreAdhesionStrength * distance_to_fibre * exp(-alpha * distance_to_fibre) * force_direction;

                p_node->AddAppliedForceContribution(force);

                // Also align the cell polarity with the fibre.
                unsigned fibre_node_A = collagen_fibre.first;
                unsigned fibre_node_B = collagen_fibre.second;

                c_vector<double, DIM> fibre_location_A = rCellPopulation.GetNode(fibre_node_A)->rGetLocation();
                c_vector<double, DIM> fibre_location_B = rCellPopulation.GetNode(fibre_node_B)->rGetLocation();

                c_vector<double, DIM> fibre_A_to_B = rCellPopulation.rGetMesh().GetVectorFromAtoB(fibre_location_A, fibre_location_B);
                fibre_A_to_B /= norm_2(fibre_A_to_B);

                // Get the angle that the fibre makes with the x-axis
                double theta = atan(fibre_A_to_B[1]/fibre_A_to_B[0]);

                if (fabs(theta - phi) > 0.5*M_PI)
                {
                    if (theta < M_PI)
                    {
                        theta += M_PI;
                    }
                    else
                    {
                        theta -= M_PI;
                    }
                }

                phi_new += dt * mFibreAlignmentStrength * sin(theta - phi);


            }

            cell_iter->GetCellData()->SetItem("direction", phi_new);


        }
    }
}

template<unsigned DIM>
void EcmFibreAttachmentForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{

	*rParamsFile <<  "\t\t\t<FibreAdhesionStrength>"<<  mFibreAdhesionStrength << "</FibreAdhesionStrength> \n" ;
    *rParamsFile <<  "\t\t\t<FibreAlignmentStrength>"<<  mFibreAlignmentStrength << "</FibreAlignmentStrength> \n" ;
	*rParamsFile <<  "\t\t\t<CollagenFibreLength>"<<  mCollagenFibreLength << "</CollagenFibreLength> \n" ;


    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class EcmFibreAttachmentForce<1>;
template class EcmFibreAttachmentForce<2>;
template class EcmFibreAttachmentForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EcmFibreAttachmentForce)
