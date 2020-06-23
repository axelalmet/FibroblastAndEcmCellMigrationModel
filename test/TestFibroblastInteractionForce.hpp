#ifndef TESTFIBROBLASTINTERACTIONFORCE_HPP_
#define TESTFIBROBLASTINTERACTIONFORCE_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "NodesOnlyMesh.hpp" // Non-periodic mesh for storing nodes in OS model
#include "GeneralisedLinearSpringForceWithVariableInteractionDistance.hpp" // Standard spring force that implements logarithmic repulsion and exponential attraction for OS models
#include "NoCellCycleModel.hpp" // No cell cycle model
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "FibroblastCellProliferativeType.hpp" // Dermal cell type
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "CellMigrationDirectionWriter.hpp" // Cell writer for migration direction
#include "PolarityTrackingModifier.hpp" // Modifier to update cell polarities
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "FibroblastInteractionForce";
static const double M_DT = 0.005;
static const double M_END_TIME = 5.0;
static const double M_SAMPLING_TIMESTEP = 0.1/M_DT;

/*
* A bunch of tests to make sure that the GeneralisedLinearSpringForceWithVariableInteractionDistance
* does what it's supposed to
*/
class TestFibroblastInteractionForce : public AbstractCellBasedTestSuite
{
public:
    // First, test that the force relaxes two spherical cells to the correct distance.
    void TestRelaxationWithoutAlignment()
    {
        // Set some parameters
        double interaction_radius = 1.5; // Cut-off distance
        double spring_stiffness = 50.0; // Spring stiffness
        double fibroblast_alignment_strength = 5.0; // Remodelling rate of fibroblasts to each other.
        double minor_axis_scale = 0.5; // Scale for minor axis of fibroblast (for now they'll be two circles)

        // Create two nodes
        std::vector<Node<2>*> nodes;

        // Push two nodes into the vector
        nodes.push_back(new Node<2>(0u, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1u, false, 0.25, 0.0));

        // Initialise the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, interaction_radius);
        
        // Create two fibroblasts
        std::vector<CellPtr> cells;

        //Create shared pointers for the cell proliferative type and mutation state
        boost::shared_ptr<AbstractCellProperty> p_fibroblast_type(CellPropertyRegistry::Instance()->Get<FibroblastCellProliferativeType>());        
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        

        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            // Placeholder cell cycle model
            NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); 
            p_cycle_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_fibroblast_type); //Make cell differentiated

            // Random initiate a fibroblast direction
            double fibroblast_direction = 0.5*i*M_PI + 0.5 * M_PI * RandomNumberGenerator::Instance()->ranf();

            // Initialise fibroblast direction and shape scale
            p_cell->GetCellData()->SetItem("direction", fibroblast_direction);

            // Set the shape scale
            p_cell->GetCellData()->SetItem("scale", minor_axis_scale);
            cells.push_back(p_cell); // Push back first cell
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Add cell writer to print the cell migration directions
        cell_population.AddCellWriter<CellMigrationDirectionWriter>();

        OffLatticeSimulation<2> simulator(cell_population);

        //Set output directory and other simulation parameters
        simulator.SetOutputDirectory(M_OUTPUT_DIRECTORY);
        simulator.SetDt(M_DT); // Set dt for Forward Euler method
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); // Pick the frequency at which we sample timesteps
        simulator.SetEndTime(M_END_TIME); // Set the end time.

        // Add linear spring force that accounts for elliptical shapes
        MAKE_PTR(GeneralisedLinearSpringForceWithVariableInteractionDistance<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetFibroblastAlignmentStrength(fibroblast_alignment_strength);
        p_spring_force->SetCutOffLength(interaction_radius);
        simulator.AddForce(p_spring_force);

        // Add simulation modifier to update cell polarities
        MAKE_PTR(PolarityTrackingModifier<2>, p_polarity_tracking_modifier);
        p_polarity_tracking_modifier->SetReorientationStrength(fibroblast_alignment_strength);
        simulator.AddSimulationModifier(p_polarity_tracking_modifier);

        simulator.Solve();
    }
};

#endif /* TESTFIBROBLASTINTERACTIONFORCE_HPP_ */