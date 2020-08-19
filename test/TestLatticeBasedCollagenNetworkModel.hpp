#ifndef TESTLATTICEBASEDCOLLAGENNETWORKMODEL_HPP_
#define TESTLATTICEBASEDCOLLAGENNETWORKMODEL_HPP_

#include <algorithm>
#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "CollagenFibreSpringForce.hpp" // Version of generalised linear spring force where the spring stiffness depends on cell-cell types
#include "GeneralisedLinearSpringForceWithVariableInteractionDistance.hpp" // Version of generalised linear spring force where the spring stiffness depends on cell-cell types
#include "NoCellCycleModel.hpp" // No cell cycle model for now.
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "NodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "FibroblastCellProliferativeType.hpp" // Dermal cell type
#include "CollagenCellProliferativeType.hpp" // Collagen fibre cell type
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "CollagenFibreTrackingModifier.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"
#include "MathsCustomFunctions.hpp" // So that I can use SmallPow

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "LatticeBasedCollagenNetworkModel";
static const double M_DT = 0.005;
static const double M_END_TIME = 1.0;
static const double M_SAMPLING_TIMESTEP = 0.1/M_DT;

/*
* A test model to study the various components that we think should be incorporated
* into the dermal extracellular matrix, which is essentially collagen fibres (and other things.)
*/
class TestLatticeBasedEcmRelaxation : public AbstractCellBasedTestSuite
{
private:

    // Generator function to assemble the nodes and cells.
    void GenerateNodesAndCells(std::vector<Node<2>*>& rNodes, std::vector<CellPtr>& rCells,
                                unsigned collagenCellsAcross, unsigned collagenCellsUp, double collagenRadius)
    {
        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_collagen_type(CellPropertyRegistry::Instance()->Get<CollagenCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        unsigned node_index = 0; // Initialise node index counter

        double horizontal_spacing = 2.0 * collagenRadius;
        double vertical_spacing = collagenRadius * sqrt(3.0);

        bool is_boundary = false;

        for (unsigned i = 0; i < collagenCellsUp; i++)
        {
            for (unsigned j = 0; j < collagenCellsAcross; j++)
            {
                // First we create the node

                // This essentially is borrowed from HoneycombMeshGenerator, so that we can create a triangular lattice
                double x = horizontal_spacing * ((double)j + 0.25 * (1.0 + SmallPow(-1.0, i + 1)));
                double y = vertical_spacing * (double)i;

                // Mark the appropriate nodes as boundary nodes
                if ( (i == 0)||(i == collagenCellsUp - 1)||(j == 0)||(j == collagenCellsAcross - 1))
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false;
                }

                // Create the node
                Node<2>* p_node(new Node<2>(node_index, is_boundary, x, y));
                p_node->SetRadius(collagenRadius);

                rNodes.push_back(p_node); // Add to vector of nodes

                node_index += 1; // Update the node index
                
                // Now create the cell
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); // Place-holder cell cycle model
                p_cycle_model->SetDimension(2);

                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                p_cell->GetCellData()->SetItem("direction", 0.0);
                p_cell->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

                rCells.push_back(p_cell); // Add the cell

            }
        }
    }

    // Function to mark collagen fibres
    void MarkCollagenFibres(NodeBasedCellPopulation<2>& rCellPopulation, double collagenRadius, double crosslinkProbability)
    {
        // Get the number of cells, i.e. indices we need to consider. We will randomly sample from the nodes to remove bias.
        unsigned num_cells = rCellPopulation.GetNumNodes();
        std::vector<unsigned> node_indices(num_cells);
        std::iota(node_indices.begin(), node_indices.end(), 0);

        std::random_shuffle(node_indices.begin(), node_indices.end());

        for (unsigned i = 0; i < node_indices.size(); i++)
        {
            unsigned node_index = node_indices[i]; // Get the node index
            // Need to set the node radius
            rCellPopulation.GetNode(node_index)->SetRadius(collagenRadius);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);

            // Get the neighbours
            std::set<unsigned> neighbouring_nodes_set = rCellPopulation.GetNodesWithinNeighbourhoodRadius(node_index, 2.1 * collagenRadius);
            std::vector<unsigned> neighbouring_nodes;

            unsigned fibres_to_mark = 3;
            double probability = RandomNumberGenerator::Instance()->ranf();

            // If p < crossLinkProbability, we choose 3 neighbours. Otherwise, we choose 4.
            if (probability > crosslinkProbability)
            {
                fibres_to_mark += 1;
            }

            unsigned marked_fibres = 0; // Count we initialise for how many fibres we will mark;

            // First check if any fibres have been marked already
            for (std::set<unsigned>::iterator elem_iter = neighbouring_nodes_set.begin();
            elem_iter != neighbouring_nodes_set.end();
            ++elem_iter)
            {                
                CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);
                
                // Chceck if the considered cell pair has been marked or not.
                std::pair<CellPtr, CellPtr> cell_pair = rCellPopulation.CreateCellPair(p_cell, p_neighbour_cell);

                if (rCellPopulation.IsMarkedSpring(cell_pair)) // If unmarked, mark it
                {
                    marked_fibres += 1;
                }
                else
                {
                    neighbouring_nodes.push_back(*elem_iter);
                }

                if (marked_fibres == fibres_to_mark)
                {
                    break;
                }
            }

            if ( (marked_fibres < fibres_to_mark)&&(neighbouring_nodes.size() > 0) ) // If there are still fibres to mark
            {
                // Randomly shuffle the indices
                std::random_shuffle(neighbouring_nodes.begin(), neighbouring_nodes.end());

                for (unsigned j = 0; j < neighbouring_nodes.size(); j++)
                {
                    unsigned neighbour_index = neighbouring_nodes[j];
                    
                    CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);
                    
                    // Chceck if the considered cell pair has been marked or not.
                    std::pair<CellPtr, CellPtr> cell_pair = rCellPopulation.CreateCellPair(p_cell, p_neighbour_cell);

                    if (!rCellPopulation.IsMarkedSpring(cell_pair)) // If unmarked, mark it
                    {
                        rCellPopulation.MarkSpring(cell_pair);
                    }

                    marked_fibres += 1; // Update the counter

                    if (marked_fibres == fibres_to_mark)
                    {
                        break;
                    }
                }
            }
        }
    }

public: 

    void TestRelaxation()
    {

        //Set the number of cells across and down for the collagen fibre array
        unsigned collagen_cells_across = 10;
        unsigned collagen_cells_up = 10;

        double collagen_radius = 0.25; // Set the radius of the collagen node
        double interaction_radius = 0.75; // Radius of interaction to determine neighbourhoods
        double crosslink_probability = 0.2; // Probability of selecting 3 links

        // // // // Mechanical parameters
        double spring_stiffness = 15.0; // Spring stiffness

        std::vector<Node<2>*> nodes; // Vector for collagen ndoes
        std::vector<CellPtr> cells; // Vector for collagen cells    

        // Generate the node positions and cells at the same time.
        GenerateNodesAndCells(nodes, cells, collagen_cells_across, collagen_cells_up, collagen_radius);

        // Generate the mesh using the nodes vector
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, collagen_cells_across * collagen_radius); // Now construct the mesh

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells); // Used for non-periodic

        cell_population.Update(); // Need to call this before searching for neighbours.

        // Mark all the relevant cell pairs
        MarkCollagenFibres(cell_population, collagen_radius, crosslink_probability);

        // Off-lattice simulation class
        OffLatticeSimulation<2> simulator(cell_population);

        // //Set output directory and other parameters
        simulator.SetOutputDirectory(M_OUTPUT_DIRECTORY);
        simulator.SetDt(M_DT); // Timestep size for Forward Euler
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // // Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        // MAKE_PTR(GeneralisedLinearSpringForceWithVariableInteractionDistance<2>, p_spring_force);
        MAKE_PTR(CollagenFibreSpringForce<2>, p_spring_force);
        // p_spring_force->SetCellCellSpringStiffness(spring_stiffness);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetCutOffLength(interaction_radius);
        simulator.AddForce(p_spring_force);

        MAKE_PTR(CollagenFibreTrackingModifier<2>, p_collagen_fibre_tracking_modifier);
        simulator.AddSimulationModifier(p_collagen_fibre_tracking_modifier);

        simulator.Solve(); // Run the simulation.

    }
};

#endif /* TESTLATTICEBASEDCOLLAGENNETWORKMODEL_HPP_ */