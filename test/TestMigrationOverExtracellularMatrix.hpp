#ifndef TESTMIGRATIONOVEREXTRACELLULARMATRIX_HPP_
#define TESTMIGRATIONOVEREXTRACELLULARMATRIX_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "GeneralisedLinearSpringForceWithVariableInteractionDistance.hpp" // Version of generalised linear spring force where the spring stiffness depends on cell-cell types
#include "FibreAlignmentBasedMigrationForce.hpp" // Migration force that aligns each fibroblast with surrounding collagen fibres
#include "NoCellCycleModel.hpp" // No cell cycle model for now.
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "NodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "FibroblastCellProliferativeType.hpp" // Dermal cell type
#include "CollagenCellProliferativeType.hpp" // Collagen fibre cell type
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "CollagenAlignmentTrackingModifier.hpp" // Modifier to align fibroblasts with local collagen fibre orientation
#include "CellMigrationDirectionWriter.hpp" // Cell writer for migration direction
#include "CellCollagenFibreOrientationWriter.hpp" // Cell writer for collagen fibre orientations
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "MigrationOverEcm";
static const double M_DT = 0.005;
static const double M_END_TIME = 5.0;
static const double M_SAMPLING_TIMESTEP = 0.1/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into modelling fibroblast migration over an extracellular matrix, which can be generalised
* to both other biological systems, but we will use these components in our wound healing model,
* as we re-think scar formation.s
*/
class TestMigrationOverExtracellularMatrix : public AbstractCellBasedTestSuite
{
private:

    // Generator function to assemble the nodes and cells.
    void GenerateNodesAndCells(std::vector<Node<2>*>& rNodes, std::vector<CellPtr>& rCells,
                                unsigned collagenCellsAcross, unsigned collagenCellsUp,
                                unsigned fibroblastCellsAcross, unsigned fibroblastCellsUp,
                                c_vector<double, 2> fibroblastMeshTranslation, double fibroblastShapeScale)
    {
        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_collagen_type(CellPropertyRegistry::Instance()->Get<CollagenCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_fibroblast_type(CellPropertyRegistry::Instance()->Get<FibroblastCellProliferativeType>());        
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        unsigned node_index = 0; // Initialise node index counter

        // Generate the collagen cells first.
        for (unsigned i = 0; i < collagenCellsUp; i++)
        {
            for (unsigned j = 0; j < collagenCellsAcross; j++)
            {
                // Define the x and y coordinates
                double x = 0.5*(double)j;
                double y = 0.5*(double)i;
                
                bool is_boundary;

                // If we're at the left, right, top, or bottom boundaries the node is a boundary node
                if ( (i == 0)||(i == collagenCellsUp - 1)||(j == 0)||(j == collagenCellsAcross - 1) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false;
                }

                Node<2>* p_node(new Node<2>(node_index, is_boundary, x, y)); // Define the node
                rNodes.push_back(p_node); // Add the node

                // Set up the collagen cell type as well 
                // Set placeholder cell cycle model (i.e. no proliferation)
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                p_cycle_model->SetDimension(2);

                // Create a criss-cross pattern for fibre orientations
                double collagen_fibre_orientation;
                if (i % 2 == 0 ) // Even rows point upwards
                {
                    collagen_fibre_orientation = 0.25*M_PI;
                }
                else // Odd rows point downwards
                {
                    collagen_fibre_orientation = -0.25*M_PI;
                }

                double fibroblast_direction = collagen_fibre_orientation; // For collagen nodes, set the fibroblast migration to be the same as the fibre orientaiton

                double collagen_amount = RandomNumberGenerator::Instance()->ranf(); // Randomise collagen amount (may come in handy later)

                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells

                // Initialise cell data
                p_cell->GetCellData()->SetItem("direction", fibroblast_direction); // Fibroblast migration direction (irrelevant)
                p_cell->GetCellData()->SetItem("orientation", collagen_fibre_orientation); // Collagen fibre orientation
                p_cell->GetCellData()->SetItem("collagen", collagen_amount); // Collagen amount
                p_cell->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

                rCells.push_back(p_cell); // Add the cell

                // Update the node index
                node_index += 1;
            }
        }

        // Now generate the fibroblasts
                // Generate the collagen cells first.
        for (unsigned i = 0; i < fibroblastCellsUp; i++)
        {
            for (unsigned j = 0; j < fibroblastCellsAcross; j++)
            {
                // Define the x and y coordinates
                double x = fibroblastMeshTranslation[0] + 0.5*(double)j;
                double y = fibroblastMeshTranslation[1] + 0.5*(double)i;
                
                bool is_boundary;

                // If we're at the top or bottom boundaries, the node is a boundary node
                if ( (i == 0)||(i == fibroblastCellsUp - 1) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false;
                }

                Node<2>* p_node(new Node<2>(node_index, is_boundary, x, y)); // Define the node
                rNodes.push_back(p_node); // Add the node

                // Set up the collagen cell type as well 
                // Set placeholder cell cycle model (i.e. no proliferation)
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                p_cycle_model->SetDimension(2);

                // Create a criss-cross pattern for fibre orientations
                double fibroblast_direction; 

                if (j == fibroblastCellsAcross - 1) // if we're at the right edge, mark these as the 'leader' cells
                {
                    fibroblast_direction = 0.0;
                }
                else
                {
                    fibroblast_direction = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();
                }
                double collagen_fibre_orientation = fibroblast_direction; // For collagen nodes, set the fibroblast migration to be the same as the fibre orientaiton

                double collagen_amount = 0.0; // Fibroblasts initially express no collagen.

                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_fibroblast_type); // Set cell to be collagen cells

                // Initialise cell data
                p_cell->GetCellData()->SetItem("direction", fibroblast_direction); // Fibroblast migration direction (irrelevant)
                p_cell->GetCellData()->SetItem("orientation", collagen_fibre_orientation); // Collagen fibre orientation
                p_cell->GetCellData()->SetItem("collagen", collagen_amount); // Collagen amount
                p_cell->GetCellData()->SetItem("scale", fibroblastShapeScale); // Shape scale (doesn't matter for collagen nodes)

                rCells.push_back(p_cell); // Add the cell

                // Update the node index
                node_index += 1;
            }
        }
    }

public: 

    void TestMigration()
    {

        //Set the number of cells across and down for the collagen fibre array
        unsigned collagen_cells_across = 12;
        unsigned collagen_cells_up = 10;

        // Do a similar thing for fibroblasts
        unsigned fibroblast_cells_across = 2;
        unsigned fibroblast_cells_up = 10;

        // Set some parameters for node-based cell populations
        double interaction_radius = 1.5; // Radius of interaction to determine neighbourhoods

        // Mechanical parameters
        double spring_stiffness = 30.0; // Spring stiffness
        double fibroblast_alignment_strength = 0.0;
        double migration_force_strength = 10.0*M_DT;
        double reorientation_strength = 2.5*M_DT;

        // Shape parameter for fibroblasts
        double fibroblast_scale = 0.5;

        // We translate the fibroblast mesh across just a smidge
        c_vector<double, 2> fibroblast_mesh_translation = zero_vector<double>(2);
        fibroblast_mesh_translation[0] = 0.25;

        // Generate the nodes and cell vectors
        std::vector<Node<2>*> nodes;
        std::vector<CellPtr> cells;

        // Generate the node positions and cells at the same time.
        GenerateNodesAndCells(nodes, cells, 
                            collagen_cells_across, collagen_cells_up,
                            fibroblast_cells_across, fibroblast_cells_up,
                            fibroblast_mesh_translation, fibroblast_scale);

        // Generate the mesh using the nodes vector
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 2.0*interaction_radius); // Now construct the mesh

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells); // Used for non-periodic

        // Add cell writers
        cell_population.AddCellWriter<CellMigrationDirectionWriter>();
        cell_population.AddCellWriter<CellCollagenFibreOrientationWriter>();

        // Off-lattice simulation class
        OffLatticeSimulation<2> simulator(cell_population);

        //Set output directory and other parameters
        simulator.SetOutputDirectory(M_OUTPUT_DIRECTORY);
        simulator.SetDt(M_DT); // Timestep size for Forward Euler
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        MAKE_PTR(GeneralisedLinearSpringForceWithVariableInteractionDistance<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetFibroblastAlignmentStrength(fibroblast_alignment_strength);
        p_spring_force->SetCutOffLength(interaction_radius);
        simulator.AddForce(p_spring_force);

        // Add the migration force 
        MAKE_PTR(FibreAlignmentBasedMigrationForce<2>, p_migration_force);
        p_migration_force->SetMigrationForceStrength(migration_force_strength);
        p_migration_force->SetReorientationStrength(reorientation_strength);
        simulator.AddForce(p_migration_force);

        // // Create a modifier to realign cell orientations with collagen.
        MAKE_PTR(CollagenAlignmentTrackingModifier<2>, p_collagen_alignment_modifier);
        p_collagen_alignment_modifier->SetNeighbourhoodRadius(interaction_radius);
        p_collagen_alignment_modifier->SetReorientationStrength(reorientation_strength);
		simulator.AddSimulationModifier(p_collagen_alignment_modifier);

        simulator.Solve(); // Run the simulation.

    }
};

#endif /* TESTMIGRATIONOVEREXTRACELLULARMATRIX */