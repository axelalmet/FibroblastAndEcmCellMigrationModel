#ifndef TESTMIGRATIONOVEREXTRACELLULARMATRIX_HPP_
#define TESTMIGRATIONOVEREXTRACELLULARMATRIX_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "GeneralisedLinearSpringForceWithVariableInteractionDistance.hpp" // Version of generalised linear spring force where the spring stiffness depends on cell-cell types
#include "PolarityBasedMigrationForce.hpp" // Migration force that aligns each fibroblast with surrounding collagen fibres
#include "NoCellCycleModel.hpp" // No cell cycle model for now.
#include "NodeBasedCellPopulation.hpp" // Overlapping spheres centre-based population
#include "NodesOnlyMesh.hpp" // Mesh with periodic vertical boundaries
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "SmartPointers.hpp" //Enables macros to save typing
#include "FibroblastCellProliferativeType.hpp" // Dermal cell type
#include "CollagenCellProliferativeType.hpp" // Collagen fibre cell type
#include "WildTypeCellMutationState.hpp" // Epidermal mutation state
#include "CellMigrationDirectionWriter.hpp" // Writer to track cell polarities
#include "CollagenFibreTrackingModifier.hpp" // Modifier to align fibroblasts with local collagen fibre orientation
#include "PolarityTrackingModifier.hpp" // Modifier to align fibroblasts with local collagen fibre orientation
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "MigrationOverEcm";
static const double M_DT = 0.005;
static const double M_END_TIME = 10.0;
static const double M_SAMPLING_TIMESTEP = 0.5/M_DT;

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
                            std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> >& rMarkedPairs,
                            std::vector<c_vector<double, 3> >& rPairsAndLowerEndpoint,
                            unsigned collagenCellsAcross, unsigned collagenCellsUp,
                            double collagenRadius, double collagenFibreLength,
                            unsigned fibroblastCellsAcross, unsigned fibroblastCellsUp,
                            c_vector<double, 2> fibroblastMeshTranslation, double fibroblastShapeScale)
    {
        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_collagen_type(CellPropertyRegistry::Instance()->Get<CollagenCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_fibroblast_type(CellPropertyRegistry::Instance()->Get<FibroblastCellProliferativeType>());        
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        unsigned node_index = 0; // Initialise node index counter

        double collagen_mesh_height = 2.0*(double)collagenCellsUp;
        double collagen_mesh_width = 2.0*(double)collagenCellsAcross;

        // Generate the collagen cells first.
        for (unsigned i = 0; i < collagenCellsUp; i++)
        {
            for (unsigned j = 0; j < collagenCellsAcross; j++)
            {
                // Define the x and y coordinates
                unsigned left_start_index, left_end_index, right_start_index, right_end_index;
                double xLeftStart, yLeftStart, xLeftEnd, yLeftEnd, xRightStart, yRightStart, xRightEnd, yRightEnd;
                double theta = 0.25*M_PI;
                bool is_boundary = false;

                // Start point of the left fibre.
                xLeftStart = collagenFibreLength * (double)j;
                yLeftStart = collagenFibreLength * (double)i;

                // If we're at the left, right, top, or bottom boundaries the node is a boundary node
                if ( (xLeftStart == 0.0)||(xLeftStart == collagen_mesh_width)
                        ||(yLeftStart == 0.0)||(yLeftStart == collagen_mesh_height) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false; // Doing this because we're gonna just repeat this later on and it's to be safe
                }

                left_start_index = node_index;

                // Set the node
                Node<2>* p_node_left_start(new Node<2>(left_start_index, is_boundary, xLeftStart, yLeftStart)); // Define the node
                p_node_left_start->SetRadius(collagenRadius);
                rNodes.push_back(p_node_left_start); // Add the node

                node_index += 1;

                // Set the cell
                NoCellCycleModel* p_cycle_model_left_start = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                p_cycle_model_left_start->SetDimension(2);

                CellPtr p_cell_left_start(new Cell(p_wildtype_state, p_cycle_model_left_start));
                p_cell_left_start->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                p_cell_left_start->GetCellData()->SetItem("direction", theta);
                p_cell_left_start->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)
                rCells.push_back(p_cell_left_start); // Add the cell

                // End point of left fibre
                xLeftEnd = xLeftStart + pow(2.0, 0.5) * collagenFibreLength * cos(theta);
                yLeftEnd = yLeftStart + pow(2.0, 0.5) * collagenFibreLength * sin(theta);

                // If we're at the left, right, top, or bottom boundaries the node is a boundary node
                if ( (xLeftEnd == 0.0)||(xLeftEnd == collagen_mesh_width)
                        ||(yLeftEnd == 0.0)||(yLeftEnd == collagen_mesh_height) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false; // Doing this because we're gonna just repeat this later on and it's to be safe
                }

                left_end_index = node_index;

                // Set the node
                Node<2>* p_node_left_end(new Node<2>(left_end_index, is_boundary, xLeftEnd, yLeftEnd)); // Define the node
                p_node_left_end->SetRadius(collagenRadius);
                rNodes.push_back(p_node_left_end); // Add the node

                node_index += 1; // Update the index count

                // Set the cell
                NoCellCycleModel* p_cycle_model_left_end = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                p_cycle_model_left_end->SetDimension(2);

                CellPtr p_cell_left_end(new Cell(p_wildtype_state, p_cycle_model_left_end));
                p_cell_left_end->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                p_cell_left_end->GetCellData()->SetItem("direction", theta);
                p_cell_left_end->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)
                rCells.push_back(p_cell_left_end); // Add the cell
               
                // Mark the fibre
                std::pair<unsigned, unsigned> node_pair = std::make_pair(left_start_index, left_end_index); // Store the pair of indices
                std::vector<unsigned> fibre_indices = {left_start_index, left_end_index}; 

                rMarkedPairs[node_pair] = fibre_indices; 

                c_vector<double, 3> fibre_and_lower_y_coord;
                fibre_and_lower_y_coord[0] = (double)left_start_index;
                fibre_and_lower_y_coord[1] = (double)left_end_index;
                fibre_and_lower_y_coord[2] = yLeftStart;

                rPairsAndLowerEndpoint.push_back(fibre_and_lower_y_coord);

                // Now let's generate the right fibre
                xRightStart = collagenFibreLength * (double)(j + 1);
                yRightStart = collagenFibreLength * (double)i;

                theta = 0.75 * M_PI;

                // If we're at the left, right, top, or bottom boundaries the node is a boundary node
                if ( (xRightStart == 0.0)||(xRightStart == collagen_mesh_width)
                        ||(yRightStart == 0.0)||(yRightStart == collagen_mesh_height) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false; // Doing this because we're gonna just repeat this later on and it's to be safe
                }

                right_start_index = node_index;

                // Set the node
                Node<2>* p_node_right_start(new Node<2>(right_start_index, is_boundary, xRightStart, yRightStart)); // Define the node
                p_node_right_start->SetRadius(collagenRadius);
                rNodes.push_back(p_node_right_start); // Add the node

                node_index += 1;

                // Set the cell
                NoCellCycleModel* p_cycle_model_right_start = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                p_cycle_model_right_start->SetDimension(2);

                CellPtr p_cell_right_start(new Cell(p_wildtype_state, p_cycle_model_right_start));
                p_cell_right_start->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                p_cell_right_start->GetCellData()->SetItem("direction", theta);
                p_cell_right_start->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)
                rCells.push_back(p_cell_right_start); // Add the cell

                // End point of left fibre
                xRightEnd = xRightStart + pow(2.0, 0.5) * collagenFibreLength * cos(theta);
                yRightEnd = yRightStart + pow(2.0, 0.5) * collagenFibreLength * sin(theta);

                // If we're at the left, right, top, or bottom boundaries the node is a boundary node
                if ( (xRightEnd == 0.0)||(xRightEnd == collagen_mesh_width)
                        ||(yRightEnd == 0.0)||(yRightEnd == collagen_mesh_height) )
                {
                    is_boundary = true;
                }
                else
                {
                    is_boundary = false; // Doing this because we're gonna just repeat this later on and it's to be safe
                }

                right_end_index = node_index;

                // Set the node
                Node<2>* p_node_right_end(new Node<2>(right_end_index, is_boundary, xRightEnd, yRightEnd)); // Define the node
                p_node_right_end->SetRadius(collagenRadius);
                rNodes.push_back(p_node_right_end); // Add the node

                node_index += 1; // Update the index count

                // Set the cell
                NoCellCycleModel* p_cycle_model_right_end = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                p_cycle_model_right_end->SetDimension(2);

                CellPtr p_cell_right_end(new Cell(p_wildtype_state, p_cycle_model_right_end));
                p_cell_right_end->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                p_cell_right_end->GetCellData()->SetItem("direction", theta);
                p_cell_right_end->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)
                rCells.push_back(p_cell_right_end); // Add the cell
               
                // Mark the fibre
                node_pair = std::make_pair(right_start_index, right_end_index); // Store the pair of indices
                fibre_indices = {right_start_index, right_end_index}; 

                rMarkedPairs[node_pair] = fibre_indices; 

                fibre_and_lower_y_coord[0] = (double)right_start_index;
                fibre_and_lower_y_coord[1] = (double)right_end_index;
                fibre_and_lower_y_coord[2] = yRightStart;

                rPairsAndLowerEndpoint.push_back(fibre_and_lower_y_coord);

            }
        }

        // Now generate the fibroblasts
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

    // Function to determine all of the fibre intersections
    void DetermineFibreIntersections(std::vector<c_vector<unsigned, 5> >& rFibreIntersections,
                                        std::vector<c_vector<double, 3> >& rPairsAndLowerEndpoint,
                                        std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> >& rMarkedPairs,
                                        NodesOnlyMesh<2>& rMesh,
                                        std::vector<CellPtr>& rCells)
    {

        // Initialise the node index counter to add crosslinks
        unsigned crosslink_node_index = rMesh.GetNumNodes();

        //Create shared pointers for cell and mutation states
        boost::shared_ptr<AbstractCellProperty> p_collagen_type(CellPropertyRegistry::Instance()->Get<CollagenCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_wildtype_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        // Iterate through, based on a rough plane sweeping algorithm
        for (unsigned i = 0; i < rPairsAndLowerEndpoint.size(); i++)
        {
            // Determine the upper y-coordinate, which determines how many segments we may need to consider.
            c_vector<double, 3> node_pair_and_lower_endpoint = rPairsAndLowerEndpoint[i];
            unsigned index_A = (unsigned)node_pair_and_lower_endpoint[0]; // First index of pair
            unsigned index_B = (unsigned)node_pair_and_lower_endpoint[1]; // Second index of pair
            double lower_y_coord = node_pair_and_lower_endpoint[2];

            // Get the locations
            c_vector<double, 2> node_location_A = rMesh.GetNode(index_A)->rGetLocation();
            c_vector<double, 2> node_location_B = rMesh.GetNode(index_B)->rGetLocation();

            double upper_y_coord = std::max(node_location_A[1], node_location_B[1]); // The y-coordinate of the upper end point should tell us when to stop considering neighbours

            std::vector<c_vector<double, 3> > considered_pairs_and_x_coords; // We will sort by the x-coordinate of the intersection with the lower y-coordinate

            for (unsigned j = i; j < rPairsAndLowerEndpoint.size(); j++) // We arguably shouldn't have to consider anything beyond the current segment
            {
                c_vector<double, 3> current_pair_and_lower_endpoint = rPairsAndLowerEndpoint[j];

                unsigned current_index_A = (unsigned)current_pair_and_lower_endpoint[0]; // First index of pair
                unsigned current_index_B = (unsigned)current_pair_and_lower_endpoint[1]; // Second index of pair
                double current_lower_endpoint = current_pair_and_lower_endpoint[2];

                // Get the locations
                c_vector<double, 2> current_node_location_A = rMesh.GetNode(current_index_A)->rGetLocation();
                c_vector<double, 2> current_node_location_B = rMesh.GetNode(current_index_B)->rGetLocation();

                double min_y = std::min(current_node_location_A[1], current_node_location_B[1]);
                double max_y = std::max(current_node_location_A[1], current_node_location_B[1]);

                if ( (lower_y_coord >= min_y)&&(lower_y_coord <= max_y) ) // If the current lower endpoint falls in the interior of a segment, we consider it.
                {
                    c_vector<double, 3> pair_and_x_coord; // Define new pair and x-coord
                    pair_and_x_coord[0] = (double)current_index_A;
                    pair_and_x_coord[1] = (double)current_index_B;

                    // We will sort by the x-coordinate of the intersection with the current lower endpoint
                    double x_intersection = current_node_location_A[0] + (lower_y_coord - current_node_location_A[1])*(current_node_location_B[0] - current_node_location_A[0])/(current_node_location_B[1] - current_node_location_A[1]);
                    pair_and_x_coord[2] = x_intersection;

                    considered_pairs_and_x_coords.push_back(pair_and_x_coord);
                }

                if (current_lower_endpoint > upper_y_coord)
                {
                    break;
                }
            }

            // If non-empty, we can add the current index as well and then determine any intersections
            // on the sorted vector
            if (considered_pairs_and_x_coords.size() > 1) 
            {

                // Sort the vector
                std::sort(considered_pairs_and_x_coords.begin(), considered_pairs_and_x_coords.end(), SortByCoordinate);

                // Now iterate through the vectors and consider only consecutive fibres
                for (unsigned k = 0; k < (considered_pairs_and_x_coords.size() - 1); k++)
                {
                    // Get the start and endpoints of the first fibre
                    c_vector<double, 3> fibre_A = considered_pairs_and_x_coords[k];

                    unsigned start_index_A = (unsigned)fibre_A[0];
                    unsigned end_index_A = (unsigned)fibre_A[1];

                    c_vector<double, 2> fibre_start_A = rMesh.GetNode(start_index_A)->rGetLocation();
                    c_vector<double, 2> fibre_end_A = rMesh.GetNode(end_index_A)->rGetLocation();

                    // Get the start and endpoints of the second fibre
                    c_vector<double, 3> fibre_B = considered_pairs_and_x_coords[k + 1];

                    unsigned start_index_B = (unsigned)fibre_B[0];
                    unsigned end_index_B = (unsigned)fibre_B[1];

                    c_vector<double, 2> fibre_start_B = rMesh.GetNode(start_index_B)->rGetLocation();
                    c_vector<double, 2> fibre_end_B = rMesh.GetNode(end_index_B)->rGetLocation();
                    
                    // Check if they intersect
                    if (DoFibresIntersect(fibre_start_A, fibre_end_A, fibre_start_B, fibre_end_B))
                    {
                        // We've found an intersection!
                        c_vector<unsigned, 5> fibre_intersection;

                        // We always order fibre intersections by the start index, so that it makes it easier to determine
                        // whether or not we've already encountered this intersection
                        if (start_index_A < start_index_B)
                        {
                            fibre_intersection[0] = start_index_A;
                            fibre_intersection[1] = end_index_A;
                            fibre_intersection[2] = start_index_B;
                            fibre_intersection[3] = end_index_B;
                        }
                        else
                        {
                            fibre_intersection[0] = start_index_B;
                            fibre_intersection[1] = end_index_B;
                            fibre_intersection[2] = start_index_A;
                            fibre_intersection[3] = end_index_A;
                        }

                        fibre_intersection[4] = crosslink_node_index;

                        // We only add the intersection if it's new.
                        unsigned count = 0;
                        for (unsigned l = 0; l < rFibreIntersections.size(); l++)
                        {
                            c_vector<unsigned, 5> current_intersection = rFibreIntersections[l];

                            if ( AreFibreIntersectionsEqual(current_intersection, fibre_intersection) )
                            {
                                break;
                            }
                            else
                            {
                                count += 1;
                            }
                        }
                        if (count == rFibreIntersections.size()) // If we're at the end, then this intersection hasn't been accounted for.
                        {
                            rFibreIntersections.push_back(fibre_intersection);

                            // Add the node and cell
                            unsigned start_A = fibre_intersection[0];
                            unsigned end_A = fibre_intersection[1];
                            unsigned start_B = fibre_intersection[2];
                            unsigned end_B = fibre_intersection[3];

                            c_vector<double, 2> start_location_A = rMesh.GetNode(start_A)->rGetLocation();
                            c_vector<double, 2> end_location_A = rMesh.GetNode(end_A)->rGetLocation();
                            c_vector<double, 2> start_location_B = rMesh.GetNode(start_B)->rGetLocation();
                            c_vector<double, 2> end_location_B = rMesh.GetNode(end_B)->rGetLocation(); 

                            c_vector<double, 2> intersection_location = GetFibreIntersection(start_location_A, end_location_A, start_location_B, end_location_B);

                            // Add a new node to the mesh
                            Node<2>* p_node(new Node<2>(crosslink_node_index, false, intersection_location[0], intersection_location[1]));
                            unsigned new_crosslink_index = rMesh.AddNode(p_node);
                            
                            // Add a new cell
                            NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Contact-inhibition-based cycle model yet.
                            p_cycle_model->SetDimension(2);

                            double theta = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

                            CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                            p_cell->SetCellProliferativeType(p_collagen_type); // Set cell to be collagen cells
                            p_cell->GetCellData()->SetItem("direction", theta);
                            p_cell->GetCellData()->SetItem("scale", 1.0); // Shape scale (doesn't matter for collagen nodes)

                            rCells.push_back(p_cell);

                            // Also update the fibre and cross-link maps
                            std::pair<unsigned, unsigned> pair_A = std::make_pair(start_A, end_A);
                            std::pair<unsigned, unsigned> pair_B = std::make_pair(start_B, end_B);

                            std::vector<unsigned> marked_indices_A = rMarkedPairs[pair_A];
                            marked_indices_A.push_back(new_crosslink_index); // Add the crosslink node
                            rMarkedPairs[pair_A] = marked_indices_A; // Update the node pair map

                            std::vector<unsigned> marked_indices_B = rMarkedPairs[pair_B];
                            marked_indices_B.push_back(new_crosslink_index); // Add the crosslink node
                            rMarkedPairs[pair_B] = marked_indices_B; // Update the node pair map

                            crosslink_node_index += 1; // Update the crosslink node index
                        }
                    }
                    
                }
            }

        }
    }

    // Function to compare the coordinates of the fibres. The assumption is that we're comparing
    // x-coordinates to x-coordinates or y-coordinates to y-coordinates
    static bool SortByCoordinate(c_vector<double, 3>& fibreA, c_vector<double, 3>& fibreB)
    {
        bool cond = false;

        if (fibreA[2] != fibreB[2])
        {
            cond = fibreA[2] < fibreB[2];
        }
        else
        {
            cond = fibreA[0] < fibreB[0];
        }

        return cond;
    }

    static bool AreFibreIntersectionsEqual(c_vector<unsigned, 5>& fibreIntersectionA, c_vector<unsigned, 5>& fibreIntersectionB)
    {
        return ( (fibreIntersectionA[0]==fibreIntersectionB[0])&&(fibreIntersectionA[1]==fibreIntersectionB[1])
                    &&(fibreIntersectionA[2]==fibreIntersectionB[2])&&(fibreIntersectionA[3]==fibreIntersectionB[3]) );
    }

    // Function to determine if two fibres intersect, based on a parametrisation from the
    // respective start and endpoints. If both of the resulting parametrisations lie between
    // 0 and 1, they intersect.
    bool DoFibresIntersect(c_vector<double, 2> fibreStartA, c_vector<double, 2> fibreEndA, 
                            c_vector<double, 2> fibreStartB, c_vector<double, 2> fibreEndB)
    {
        // Define the respective x and y coordinates
        double xStartA = fibreStartA[0];
        double xEndA = fibreEndA[0];
        double yStartA = fibreStartA[1];
        double yEndA = fibreEndA[1];

        double xStartB = fibreStartB[0];
        double xEndB = fibreEndB[0];
        double yStartB = fibreStartB[1];
        double yEndB = fibreEndB[1];

        double t_A = ( xEndB*(yStartA - yStartB) + xStartA*(yStartB - yEndB) + xStartB*(yEndB - yStartA) )
                        /( -(xEndB - xStartB)*(yEndA - yStartA) + (xEndA - xStartA)*(yEndB - yStartB) );
        double t_B = ( xStartB*(yStartA - yEndA) + xStartA*(yEndA - yStartB) + xEndA*(yStartB - yStartA) )
                        /( (xEndB - xStartB)*(yEndA - yStartA) - (xEndA - xStartA)*(yEndB - yStartB) );

        return ( (t_A > 0.0)&&(t_A < 1.0)&&(t_B > 0.0)&&(t_B < 1.0) );
    }


    // Function that returns the fibre intersection
    c_vector<double, 2> GetFibreIntersection(c_vector<double, 2> fibreStartA, c_vector<double, 2> fibreEndA, 
                            c_vector<double, 2> fibreStartB, c_vector<double, 2> fibreEndB)
    {
        // Define the respective x and y coordinates
        double xStartA = fibreStartA[0];
        double xEndA = fibreEndA[0];
        double yStartA = fibreStartA[1];
        double yEndA = fibreEndA[1];

        double xStartB = fibreStartB[0];
        double xEndB = fibreEndB[0];
        double yStartB = fibreStartB[1];
        double yEndB = fibreEndB[1];

        double t_A = ( xEndB*(yStartA - yStartB) + xStartA*(yStartB - yEndB) + xStartB*(yEndB - yStartA) )
                        /( -(xEndB - xStartB)*(yEndA - yStartA) + (xEndA - xStartA)*(yEndB - yStartB) );

        c_vector<double, 2> intersection;
        intersection[0] = xStartA + t_A * (xEndA - xStartA);
        intersection[1] = yStartA + t_A * (yEndA - yStartA);

        return intersection;
    }

public: 

    void TestMigration()
    {

        //Set the number of cells across and down for the collagen fibre array
        unsigned collagen_cells_across = 12;
        unsigned collagen_cells_up = 5;

        // Do a similar thing for fibroblasts
        unsigned fibroblast_cells_across = 2;
        unsigned fibroblast_cells_up = 5;

        // Set some parameters for node-based cell populations
        double collagen_radius = 0.5;
        double collagen_fibre_length = 1.0;
        double interaction_radius = 1.5; // Radius of interaction to determine neighbourhoods

        // Mechanical parameters
        double spring_stiffness = 15.0; // Spring stiffness
        double migration_force_strength = 1.0; // Strength of migration force in the direction of polarity
        double reorientation_strength = 2.5; // Rate of remodelling of polarity to velocity
        double cell_ecm_stiffness_factor = 0.1; // Scaling of collagen-collagen spring stiffness when 'buckled'

        // Shape parameter for fibroblasts
        double fibroblast_scale = 0.5;

        // We translate the fibroblast mesh across just a smidge
        c_vector<double, 2> fibroblast_mesh_translation = zero_vector<double>(2);
        fibroblast_mesh_translation[0] = 0.25;
        fibroblast_mesh_translation[1] = 1.5;

        // Generate the nodes and cell vectors
        std::vector<Node<2>*> nodes;
        std::vector<CellPtr> cells;
        std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> > fibres_to_mark;
        std::vector<c_vector<double, 3> > fibres_and_lower_y_coord;

        // Generate the node positions and cells at the same time.
        GenerateNodesAndCells(nodes, cells,
                            fibres_to_mark, fibres_and_lower_y_coord,
                            collagen_cells_across, collagen_cells_up,
                            collagen_radius, collagen_fibre_length, 
                            fibroblast_cells_across, fibroblast_cells_up,
                            fibroblast_mesh_translation, fibroblast_scale);

        // Generate the mesh using the nodes vector
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 2.0*interaction_radius); // Now construct the mesh

        // Let's now determine all the fibre intersections (and their crosslinks)
        std::vector<c_vector<unsigned, 5> > intersecting_fibres;

        DetermineFibreIntersections(intersecting_fibres, fibres_and_lower_y_coord, fibres_to_mark, mesh, cells);

        PRINT_VARIABLE(intersecting_fibres.size());

        //Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells); // Used for non-periodic

        // Mark all the relevant cell pairs
        for (std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> >::iterator map_iter = fibres_to_mark.begin();
        map_iter != fibres_to_mark.end();
        map_iter++)
        {
            std::vector<unsigned> nodes_along_fibre = map_iter->second; // These are all the possible nodes along the fibre

            // We will sort the node indices now by x-coordinate and then mark consecutive pairs
            std::vector<std::pair<double, unsigned> > x_coordinates_and_nodes;

            // Determine the x-coordinate of each node
            for (unsigned j = 0; j < nodes_along_fibre.size(); j++)
            {
                unsigned node_index = nodes_along_fibre[j];
                double x = cell_population.rGetMesh().GetNode(node_index)->rGetLocation()[0];

                std::pair<double, unsigned> x_and_index = std::make_pair(x, node_index);

                x_coordinates_and_nodes.push_back(x_and_index);
            }
            
            std::sort(x_coordinates_and_nodes.begin(), x_coordinates_and_nodes.end()); // Sort the vector by x-coordinate

            // Let's average out the new radii, otherwise putting in these new crosslinks will push everything out.
            double new_node_radius = collagen_radius;

            for (unsigned j = 0; j < (x_coordinates_and_nodes.size() - 1); j++)
            {
                unsigned start_index = x_coordinates_and_nodes[j].second;
                unsigned end_index = x_coordinates_and_nodes[j + 1].second;

                // Update the radii
                cell_population.rGetMesh().GetNode(start_index)->SetRadius(new_node_radius);
                cell_population.rGetMesh().GetNode(end_index)->SetRadius(new_node_radius);

                // Create the cell pair
                CellPtr p_cell_start = cell_population.GetCellUsingLocationIndex(start_index);
                CellPtr p_cell_end = cell_population.GetCellUsingLocationIndex(end_index);

                std::pair<CellPtr, CellPtr> cell_pair = cell_population.CreateCellPair(p_cell_start, p_cell_end);

                cell_population.MarkSpring(cell_pair);
            }

        }

        cell_population.AddCellWriter<CellMigrationDirectionWriter>();

        // Off-lattice simulation class
        OffLatticeSimulation<2> simulator(cell_population);

        //Set output directory and other parameters
        simulator.SetOutputDirectory(M_OUTPUT_DIRECTORY);
        simulator.SetDt(M_DT); // Timestep size for Forward Euler
        simulator.SetSamplingTimestepMultiple(M_SAMPLING_TIMESTEP); //Sample the simulation at every hour
        simulator.SetEndTime(M_END_TIME); //Hopefully this is long enough for a steady state

        // Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
        MAKE_PTR(GeneralisedLinearSpringForceWithVariableInteractionDistance<2>, p_spring_force);
        p_spring_force->SetCellCellSpringStiffness(spring_stiffness);
        p_spring_force->SetCutOffLength(interaction_radius);
        p_spring_force->SetCellEcmStiffnessMultiplicationFactor(cell_ecm_stiffness_factor);
        simulator.AddForce(p_spring_force);

        // Add the migration force 
        MAKE_PTR(PolarityBasedMigrationForce<2>, p_migration_force);
        p_migration_force->SetMigrationForceStrength(migration_force_strength);
        simulator.AddForce(p_migration_force);

        // // Create a modifier to realign cell orientations with collagen.
        MAKE_PTR(CollagenFibreTrackingModifier<2>, p_collagen_fibre_tracking_modifier);
		simulator.AddSimulationModifier(p_collagen_fibre_tracking_modifier);

        // Add polarity tracker
        MAKE_PTR(PolarityTrackingModifier<2>, p_polarity_tracking_modifier);
        p_polarity_tracking_modifier->SetReorientationStrength(reorientation_strength);
        simulator.AddSimulationModifier(p_polarity_tracking_modifier);
        simulator.Solve(); // Run the simulation.

    }
};

#endif /* TESTMIGRATIONOVEREXTRACELLULARMATRIX */