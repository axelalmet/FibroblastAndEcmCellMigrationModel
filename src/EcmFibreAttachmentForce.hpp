#ifndef ECMFIBREATTACHMENTFORCE_HPP_
#define ECMFIBREATTACHMENTFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

/**
 * A force class that models the adhesion to and alignment
 * with neighbouring ECM fibres for a migrating cell.
 */
template<unsigned DIM>
class EcmFibreAttachmentForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:


    /*
     * Strength of adhesion to fibres
     */
    double mFibreAdhesionStrength;
    
    /*
     * Strength of alignment to nearby fibres
     */
    double mFibreAlignmentStrength;

    /*
     * Characteristic length of collagen fibres, so that
     * we can determine how far to search
     */
    double mCollagenFibreLength;

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mFibreAdhesionStrength;
        archive & mFibreAlignmentStrength;
        archive & mCollagenFibreLength;
    }

public:

    /**
     * Constructor.
     */
    EcmFibreAttachmentForce();

    /**
     * Destructor.
     */
    ~EcmFibreAttachmentForce();

    /*
     * Get adhesion force strength
     * 
     * @return mAdhesionForceStrength
     */
    double GetFibreAdhesionStrength();

    /*
     * Set the fibre reorientation strength
     * 
     * @param migrationForceStrength
     */
    void SetFibreAdhesionStrength(double fibreAdhesionStrength);


    /*
     * Get fibre reorientation strength
     * 
     * @return mAdhesionForceStrength
     */
    double GetFibreAlignmentStrength();

    /*
     * Set the fibre reorientation strength
     * 
     * @param migrationForceStrength
     */
    void SetFibreAlignmentStrength(double fibreAlignmentStrength);

    /*
     * Set collagen length to determine neighbouring fibres
     */
    double GetCollagenFibreLength();

    void SetCollagenFibreLength(double collagenFibreLength);

    /*
     * Get the relevant ECM fibres to consider for the specified node.
     */ 
    std::vector<std::pair<unsigned, unsigned> > GetNeighbouringEcmFibres(unsigned nodeIndex, AbstractCellPopulation<DIM>& rCellPopulation);

    /*
     * Calculate the orthogonal projection which we will use to calculate the adhesion force strength, as well
     * as determine whether or not the cell will interact with this fibre.
     */
    c_vector<double, DIM> CalculateOrthogonalProjectionOntoFibre(unsigned nodeIndex, std::pair<unsigned, unsigned> fibre, AbstractCellPopulation<DIM>& rCellPopulation);

    /*
     * Condition to check if an ECM fibre has been considered or not
     */
    bool HasFibreBeenConsidered(std::pair<unsigned, unsigned> fibre, std::vector<std::pair<unsigned, unsigned> > neighbouringFibres);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EcmFibreAttachmentForce)

#endif /*EcmFibreAttachmentForce_HPP_*/
