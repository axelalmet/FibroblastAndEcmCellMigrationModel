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

#include "CellCollagenFibreOrientationWriter.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "UblasVectorInclude.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellCollagenFibreOrientationWriter<ELEMENT_DIM, SPACE_DIM>::CellCollagenFibreOrientationWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellcollagenfibreorientation.dat")
{
    this->mVtkVectorCellDataName = "Cell collagen fibre orientation";
    this->mOutputScalarData = false;
    this->mOutputVectorData = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> CellCollagenFibreOrientationWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    // Only define non-zero-length vectors from cells with non-zero amounts of collagen
    double collagen = pCell->GetCellData()->GetItem("collagen");
    c_vector<double, SPACE_DIM> collagen_fibre_orientation = zero_vector<double>(SPACE_DIM);

    if (dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(pCellPopulation))
    {
        if (collagen > 0.0) // If there is a non-zero amount of collagen
        {
            double orientation = pCell->GetCellData()->GetItem("orientation");
            collagen_fibre_orientation[0] = cos(orientation);
            collagen_fibre_orientation[1] = sin(orientation);
        }
    }

    return collagen_fibre_orientation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellCollagenFibreOrientationWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    double orientation = pCell->GetCellData()->GetItem("orientation");
    c_vector<double, SPACE_DIM> collagen_fibre_orientation = GetVectorCellDataForVtkOutput(pCell, pCellPopulation);

    *this->mpOutStream << location_index << " " << cell_id << " ";
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }
    *this->mpOutStream << orientation << " ";
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << collagen_fibre_orientation[i] << " ";
    }
}

// Explicit instantiation
template class CellCollagenFibreOrientationWriter<1,1>;
template class CellCollagenFibreOrientationWriter<1,2>;
template class CellCollagenFibreOrientationWriter<2,2>;
template class CellCollagenFibreOrientationWriter<1,3>;
template class CellCollagenFibreOrientationWriter<2,3>;
template class CellCollagenFibreOrientationWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellCollagenFibreOrientationWriter)
