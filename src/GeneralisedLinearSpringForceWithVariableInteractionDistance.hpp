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

#ifndef GeneralisedLinearSpringForceWithVariableInteractionDistance_HPP_
#define GeneralisedLinearSpringForceWithVariableInteractionDistance_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * An extension of GeneralisedLinearSpringForce to account for
 * elliptical shapes of fibroblasts and how that changes the assumed
 * rest length of l_{ij} = 2R for overlapping spheres models.
 * The rest length is determined by calculating the range parameter, sigma,
 * from the Gaussian overlap model proposed by Gay and Berne (1972), which is
 * used for liquid crystals. We only use the range parameter from this model,
 * not the Lennard-Jones-type potential that's used to calculate forces and torques.
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class GeneralisedLinearSpringForceWithVariableInteractionDistance : public AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForces;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mCellCellSpringStiffness;
        archive & mMeinekeDivisionRestingSpringLength;
        archive & mMeinekeSpringGrowthDuration;
        archive & mCellEcmStiffnessMultiplicationFactor;
    }

protected:

    /**
     * Spring stiffness.
     *
     * Represented by the parameter mu in the model by Meineke et al (2001) in
     * their off-lattice model of the intestinal crypt
     * (doi:10.1046/j.0960-7722.2001.00216.x).
     */
    double mCellCellSpringStiffness;

    /**
     * Initial resting spring length after cell division.
     * Has units of cell size at equilibrium rest length
     *
     * The value of this parameter should be larger than mDivisionSeparation,
     * because of pressure from neighbouring springs.
     */
    double mMeinekeDivisionRestingSpringLength;

    /**
     * The time it takes for the springs rest length to increase from
     * mMeinekeDivisionRestingSpringLength to its natural length.
     *
     * The value of this parameter is usually the same as the M Phase of the cell cycle and defaults to 1.
     */
    double mMeinekeSpringGrowthDuration;

    /* 
     * Parameter that scales the collagen spring stiffness 
     * as a proxy for buckling.
     */
    double mCellEcmStiffnessMultiplicationFactor;

public:

    /**
     * Constructor.
     */
    GeneralisedLinearSpringForceWithVariableInteractionDistance();

    /**
     * Destructor.
     */
    virtual ~GeneralisedLinearSpringForceWithVariableInteractionDistance();

    /**
     * Return a multiplication factor for the spring constant, which
     * returns a default value of 1.
     *
     * This method may be overridden in subclasses.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    virtual double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                              unsigned nodeBGlobalIndex,
                                                              AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                              bool isCloserThanRestLength);

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes.
     *
     * Note that this assumes they are connected and is called by AddForceContribution()
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /* Calculate the rest length, depending on the cell shapes (as a proxy for cell type).
     * This method is based on the Gay-Berne (1972) model that proposes a Gaussian overlap
     * model for potentials used to study liquid crystals, which accounts for elliptical shapes,
     * which we use to model fibroblasts. 
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * 
     * @return The interaction distance between nodeA and nodeB
     * 
     */
    double CalculateRestLength(unsigned nodeAGlobalIndex,
                            unsigned nodeBGlobalIndex,
                            c_vector<double, SPACE_DIM> unitDifference,
                            AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * @return mMeinekeDivisionRestingSpringLength
     */
    double GetMeinekeDivisionRestingSpringLength();

    /**
     * @return mCellCellSpringStiffness
     */
    double GetCellCellSpringStiffness();

    /**
     * @return mMeinekeSpringGrowthDuration
     */
    double GetMeinekeSpringGrowthDuration();

    /**
     * @return mCellEcmStiffnessMultiplicationFactor
     */
    double GetCellEcmStiffnessMultiplicationFactor();

    /**
     * Set mCellCellSpringStiffness.
     *
     * @param springStiffness the new value of mCellCellSpringStiffness
     */
    void SetCellCellSpringStiffness(double springStiffness);

    /**
     * Set mMeinekeDivisionRestingSpringLength.
     *
     * @param divisionRestingSpringLength the new value of mMeinekeDivisionRestingSpringLength
     */
    void SetMeinekeDivisionRestingSpringLength(double divisionRestingSpringLength);

    /**
     * Set mMeinekeSpringGrowthDuration.
     *
     * @param springGrowthDuration the new value of mMeinekeSpringGrowthDuration
     */
    void SetMeinekeSpringGrowthDuration(double springGrowthDuration);

    /**
     * Set mCellEcmStiffnessMultiplicationFactor
     * 
     * @param CellEcmStiffnessMultiplicationFactor the new value of mCellEcmStiffnessMultiplicationFactor
     */
    void SetCellEcmStiffnessMultiplicationFactor(double CellEcmStiffnessMultiplicationFactor);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(GeneralisedLinearSpringForceWithVariableInteractionDistance)

#endif /*GeneralisedLinearSpringForceWithVariableInteractionDistance_HPP_*/
