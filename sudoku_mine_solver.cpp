#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
using namespace std;


class SudokuMineSolver{ // start of SudokuMineSolver class definition
    

    public:

        static const int N_ROWS = 9;
        static const int N_COLS = 9;
        static const int N_BLOCKS = 9;
        static const int N_CELLS = 81;


        SudokuMineSolver(){initialize_structure();} // Constructor


        bool load_puzzle(istream &in){ // start of load_puzzle method definition
            
            // reading the puzzle from the input file
            for(int r = 0; r < N_ROWS; ++r){

                for(int c = 0; c < N_COLS; ++c){
                
                    int value;

                    if(!(in >> value)){return false;}

                    int idx = r * N_COLS + c; // calculate the index of the current cell

                    clue[idx] = value; // store the value in the `clue` -- how many allowable mines are in the surrounding cells
                                        // 0 for blank, 1-8 for numbers (how many mines are in the surrounding cells)
                }
            }
            // done reading the puzzle from the input file


            initialize_CSP_state(); // initialize the CSP state --- the domains, constraint flags and adjacency
            return true;

        } // end of load_puzzle method definition


        bool solve(){ // start of solve method definition
            
            if(!check_all_constraints()){return false;} // if the initial board does not satisfy all constraints then there is no solution so return false

            nodes_generated = 0; // initialize the nodes generated counter to 0
            goal_depth = -1; // initialize the depth of the solution node to an invalid

            bool found = backtrack_search(0); // KICK-OFF THE RECURSIVE BACKTRACKING SEARCH AT DEPTH 0
                                                // in Figure5 this is the top part of the pseudo-code
            return found; // return whether a solution was found

        } // end of solve method definition


        void write_output(ostream &out) const{ // start of write_output method definition
            
            out << goal_depth << '\n'; // first line of the output file is d, the depth of the goal node
            out << nodes_generated << '\n'; // second line of the output file is n, the count of nodes generated
            
            // writing lines three to eleven of the output file
            for(int r = 0; r < N_ROWS; ++r){ // loop through each row

                for(int c = 0; c < N_COLS; ++c){ // loop through each column
                    
                    int idx = r * N_COLS + c; // calculate the index of the current cell

                    out << assignment_value[idx]; // OUTPUT 0 OR 1 FOR THIS CELL
                                                    // 0 means no mine, 1 means a mine

                    if(c < N_COLS - 1){ // if this is not the final column of the row ...
                        out << ' '; // print a space
                    }
                    else{ // ... else this is the final column of the row ...
                        out << '\n'; // end the line and move to the next row
                    }
                }
            }
        

        } // end of write_output method definition



    private:

        
        int clue[N_CELLS];   // clue[i] stores the given number at cell i (0-8) which is the number of mines in the surrounding cells
        bool neighbor_sum_constraint_active[N_CELLS]; // neighbor_constraint_active[i] is true if the neighbor-sum constraint is active for cell i


        int row_of[N_CELLS]; // row_of[i] stores the row index (0-8) of cell i
        int col_of[N_CELLS]; // col_of[i] stores the column index (0-8) of cell i
        int block_of[N_CELLS]; // block_of[i] stores the block index (0-8) of cell i

        
        // array of vecs
        vector<int> row_cells[N_ROWS]; // row_cells[r] lists the indices of cells in row r
        vector<int> col_cells[N_COLS]; // col_cells[c] lists the indices of cells in column c
        vector<int> block_cells[N_BLOCKS]; // block_cells[b] lists the indices of cells in block b

        vector<int> neighbor_cells[N_CELLS]; // neighbor_cells[i] lists the indices of the 8-neighbors of cell i
        vector<int> constraint_adjacency[N_CELLS]; // constraint_adjacency[i] lists the indices of cells that are constrained by ANY constraint to cell i
        bool assigned[N_CELLS];  // assigned[i] is true if cell i has an assigned value and false if it does not
        int assignment_value[N_CELLS]; // assignment_value[i] stores the assigned value (0 or 1) of cell i

        bool domain[N_CELLS][2]; // domain[i][v] is true if value v is still allowed for variable (aka cell) i

        long long nodes_generated; // counts how many nodes have been generated during the search
        int goal_depth; // the depth at which the goal was first found


        //------------------------------------------------------------------------------------------------------------------


        void initialize_structure(){ // start of initialize_structure method definition

            
            // compute the row, column and block index for each of the 81 cells
            for(int i = 0; i < N_CELLS; ++i){
                
                int r = i / N_COLS; // calculate the ROW index of cell i
                int c = i % N_COLS; // calculate the COLUMN index of cell i
                row_of[i] = r; // store the row index of cell i
                col_of[i] = c; // store the column index of cell i

                int block_row = r / 3; // is in {0, 1, 2}
                int block_col = c / 3; // is in {0, 1, 2}
                int b = block_row * 3 + block_col; // calculate the BLOCK index of cell i
                block_of[i] = b; // store the block index of cell i


                row_cells[r].push_back(i); // add cell i to the corresponding row's vec
                col_cells[c].push_back(i); // add cell i to the corresponding column's vec
                block_cells[b].push_back(i); // add cell i to the corresponding block's vec
            }


            // compute the eight-neighbors for each cell
            for(int r = 0; r < N_ROWS; ++r){

                for(int c = 0; c < N_COLS; ++c){

                    // inside these two loops we are 'standing at' a specific cell (r,c)

                    int idx = r * N_COLS + c; // calculate the index of the current cell

                    for(int dr = -1; dr <= 1; ++dr){ // we want to visit all delta row and delta column  --- (-1, -1) to (1,1)

                        for(int dc = -1; dc <= 1; ++dc){
                            
                            // inside these two nested loops we are 'standing at' a specific neighbor of (r,c)
                            // so this loop goes over the 3x3 square centered at (r,c)
                            
                            if(dr == 0 && dc == 0){continue;} // skip the current center cell (r,c) itself
                            
                            // grab neighbor's row and column indices
                            int nr = r + dr;
                            int nc = c + dc;

                            if(nr >=0 && nr < N_ROWS && nc >= 0 && nc < N_COLS){ // if neighbor is within bounds...
                                int nidx = nr * N_COLS + nc; // calculate the index of the neighbor
                                neighbor_cells[idx].push_back(nidx); // add the neighbor's index to the list of neighbors for the current cell
                            }

                        }
                    }


                }

            }

        } // end of initialize_structure method definition


        void add_edges_for_constraint_group(const vector<int>& constraint_group, bool adj[N_CELLS][N_CELLS]){ // start of add_edges_for_constraint_group method definition
            
             int n = static_cast<int>(constraint_group.size());  // get the number of cells in this constraint group

             for(int i = 0; i < n; ++i){ // loop through each cell in the constraint group

                int current_cell_idx = constraint_group[i]; // get the index of the current cell in the constraint group

                for(int j = i + 1; j < n; ++j){ // compare the current cell to all subsequent cells in the constraint group

                    int subsequent_cell_idx = constraint_group[j];  // get the index of the subsequent cell in the constraint group

                    // symmetricallyset the entries in the adjacency matrix to true --- indicating that there is a constraint between cells a and b
                    adj[current_cell_idx][subsequent_cell_idx] = true;
                    adj[subsequent_cell_idx][current_cell_idx] = true;
                }

             }

        } // end of add_edges_for_constraint_group method definition


        
        void build_constraint_adjacency(){ // start of build_constraint_adjacency method definition

            bool adj[N_CELLS][N_CELLS]; // instantiate the constraint adjacency matrix for 81 variables
            for(int i = 0; i < N_CELLS; ++i){
                for(int j = 0; j < N_CELLS; ++j){
                    adj[i][j] = false;  // initialize all entries to false --- no edges (aka constraints) between any two cells yet
                }
            }
            

            // tie up row constraints
            for(int r = 0; r < N_ROWS; ++r){
                add_edges_for_constraint_group(row_cells[r], adj); // connect all pairs of cells in the same row
            }
            
            // tie up column constraints
            for(int c = 0; c < N_COLS; ++c){
                add_edges_for_constraint_group(col_cells[c], adj); // connect all pairs of cells in the same column
            }

            // tie up block constraints
            for(int b = 0; b < N_BLOCKS; ++b){
                add_edges_for_constraint_group(block_cells[b], adj); // connect all pairs of cells in the same block
            }

            // tie up neighbor-sum constraints
            for(int idx = 0; idx < N_CELLS; ++idx){
                if(neighbor_sum_constraint_active[idx]){add_edges_for_constraint_group(neighbor_cells[idx], adj);} // connect all pairs of cells that are neighbors and have a neighbor-sum constraint active
            }


            // convert adjacency matrix into ADJACENCY LIST
            for(int i = 0; i < N_CELLS; ++i){ // loop through each cell 

                for(int j = 0; j < N_CELLS; ++j){ // compare i to all other cells j...
                    
                    if(adj[i][j]){constraint_adjacency[i].push_back(j);} // ... if there is an edge (constraint) between i and j (as indicated by the matrix) then add j to the list of constrained neighbors for i
 
                }

            }

        }  // end of build_constraint_adjacency method definition




        void initialize_CSP_state(){ // start of initialize_CSP_state method definition

            /* initializes domains, assignment status, constraint flags and constraint adjacency */
            
            for(int i = 0; i < N_CELLS; ++i){neighbor_sum_constraint_active[i] = (clue[i] > 0);} // initialize the neighbor-sum constraint flag
            
            // initialize assignment status and domain
            for(int i = 0; i < N_CELLS; ++i){

                assigned[i] = false;
                assignment_value[i] = -1;

                domain[i][0] = true; // all variables allow value `0` (no mine)
                if(clue[i] == 0){domain[i][1] = true;}else{domain[i][1] = false;} // if cell has no clue then `1` is allowed otherwise (cell has a number inside) then mines are not allowed
            }

            build_constraint_adjacency(); // build constraint adjacency

        } // end of initialize_CSP_state method definition

        
        int count_domain_size(int var) const{ // start of count_domain_size method definition
            
            /* counts the current domain size (amount of remaining values) for a variable (aka cell)*/

            int size = 0;
            if(domain[var][0]){++size;} // if `0` is still allowed for this cell then increment the domain size
            if(domain[var][1]){++size;} // if `1` is still allowed for this cell then increment the domain size
            return size;
        } // end of count_domain_size method definition


        bool check_all_constraints() const{ // start of check_all_constraints method definition
             
             // check for domain wipeout -- if any cell's domain is empty then this branch of the search tree is dead-end so return false
             for(int i = 0; i < N_CELLS; ++i){
                if(!assigned[i]){
                    if(!domain[i][0] && !domain[i][1]){
                        return false;
                    }
                }
             }


             // check row constraints --- each row must have exactly 3 mines
             for(int r = 0; r < N_ROWS; ++r){

                int min_sum = 0; // given what we know so far, what is the smallest possible amount of mines this row can end up with?
                int max_sum = 0; // given what we know so far, what is the largest possible amount of mines this row can end up with?

                for(int a_cell_idx_in_this_row : row_cells[r]){ // scan through each cell in the current row

                    if(assigned[a_cell_idx_in_this_row]){ // the cell under investigation is already assigned a value
                        int v = assignment_value[a_cell_idx_in_this_row];

                        // fixed contribution to the possible sum trackers
                        min_sum += v;
                        max_sum += v;
                    }
                    else{
                        // for an unassigned cell,
                            // if the domain CAN be 1 then it contribues to max_sum only
                            // if the domain MUST be 1 then it contributes to both min_sum and max_sum
                        if(domain[a_cell_idx_in_this_row][1]){++max_sum;}
                        if(!domain[a_cell_idx_in_this_row][0] && domain[a_cell_idx_in_this_row][1]){++min_sum;}
                    }

                }

                // KEY: we're checking what we can POSSIBLY end up with
                
                if(3 < min_sum || 3 > max_sum){return false;} // if 3 is not with [the smallest possible amount of mines this row can end up with, the largest possible amount of mines this row can end up with] then return false


             }
             // end of checking row constraints


             // check column constraints --- each column must have exactly 3 mines
             for(int c = 0; c < N_COLS; ++c){

                int min_sum = 0; // given what we know so far, what is the smallest possible amount of mines this row can end up with?
                int max_sum = 0; // given what we know so far, what is the largest possible amount of mines this row can end up with?

                for(int a_cell_idx_in_this_col : col_cells[c]){ // scan through each cell in the current row

                    if(assigned[a_cell_idx_in_this_col]){ // the cell under investigation is already assigned a value
                        int v = assignment_value[a_cell_idx_in_this_col];

                        // fixed contribution to the sum
                        min_sum += v;
                        max_sum += v;
                    }
                    else{

                        // for an unassigned cell,
                            // if the domain can be 1 then it contribues to max_sum only
                            // if the domain must be 1 then it contributes to both min_sum and max_sum
                        if(domain[a_cell_idx_in_this_col][1]){++max_sum;}
                        if(!domain[a_cell_idx_in_this_col][0] && domain[a_cell_idx_in_this_col][1]){++min_sum;}
                    }

                }

                // KEY: we're checking what we can POSSIBLY end up with
                
                if(3 < min_sum || 3 > max_sum){return false;}
             }
             // end of checking column constraints


             // check block constraints --- each block must have exactly 3 mines
             for(int b = 0; b < N_BLOCKS; ++b){

                int min_sum = 0; // given what we know so far, what is the smallest possible amount of mines this row can end up with?
                int max_sum = 0; // given what we know so far, what is the largest possible amount of mines this row can end up with?

                for(int a_cell_idx_in_this_block : col_cells[b]){ // scan through each cell in the current row

                    if(assigned[a_cell_idx_in_this_block]){ // the cell under investigation is already assigned a value
                        int v = assignment_value[a_cell_idx_in_this_block];

                        // fixed contribution to the sum
                        min_sum += v;
                        max_sum += v;
                    }
                    else{

                        // for an unassigned cell,
                            // if the domain can be 1 then it contribues to max_sum only
                            // if the domain must be 1 then it contributes to both min_sum and max_sum
                        if(domain[a_cell_idx_in_this_block][1]){++max_sum;}
                        if(!domain[a_cell_idx_in_this_block][0] && domain[a_cell_idx_in_this_block][1]){++min_sum;}
                    }

                }

                // KEY: we're checking what we can POSSIBLY end up with
                
                if(3 < min_sum || 3 > max_sum){return false;}
             }
             // end of checking block constraints

             // check neighbor-sum constraints --- each cell must have the correct number of mines in its neighborhood based on clue
             for(int i = 0; i < N_CELLS; ++i){
                if(!neighbor_sum_constraint_active[i]){continue;} // if the neighbor-sum constraint is not active for this cell then skip it

                int target = clue[i];

                int min_sum = 0;
                int max_sum = 0;

                for(int a_cell_idx_thats_a_neighbor_of_i : neighbor_cells[i]){
                    if(assigned[a_cell_idx_thats_a_neighbor_of_i]){
                        int v = assignment_value[a_cell_idx_thats_a_neighbor_of_i];

                        min_sum += v;
                        max_sum += v;
                    }
                    else{

                        if(domain[a_cell_idx_thats_a_neighbor_of_i][1]){++max_sum;}
                        if(!domain[a_cell_idx_thats_a_neighbor_of_i][0] && domain[a_cell_idx_thats_a_neighbor_of_i][1]){++min_sum;}
                    }

           
                }

                if(target < min_sum || target > max_sum){return false;}



             }
             // end of checking neighbor-sum constraints




             // KEY: all group (e.g. row, col etc) constraint checks are done with the same principle
                   // the distinction is: row, col, block have a fixed target of 3 whereas
                   // neighbor-sum has a variable target based on the clue


             return true; // if we've made it this far then ALL the constraints are satisfied so return true

        } // end of check_all_constraints method definition


        bool is_complete() const{ // start of is_complete method definition
           
            for(int i = 0; i < N_CELLS; ++i){
                if(!assigned[i]){return false;} // if any cell is not assigned then the puzzle is not complete so return false
            }
            return true; // if all cells are assigned then the puzzle is complete so return true

        } // end of is_complete method definition






    // -------------------------------------------------------------------------------------------------------------------


    int select_unassigned_variable() const{ // start of select_unassigned_variable method definition
        
       /* this method implements the MRV and degree heuristic to select the next best variable (aka cell) to assign */


       // initialize some dummies to start with
       int best_var_so_far_idx = -1;
       int best_domain_size_so_far = INT_MAX;
       int best_degree_so_far = -1;

       for(int i = 0; i < N_CELLS; ++i){ // scan through all the variables (aka cells) to try and find the next best variable to assign

          if(assigned[i]){continue;} // if the variable (aka cell) is already assigned then skip it

          int current_candidate_domain_size = count_domain_size(i); // grab the domain size of the current candidate

          if(current_candidate_domain_size < best_domain_size_so_far){ // this is core MRV heuristic logic

               // once we find a better candidate, we update
               best_var_so_far_idx = i;
               best_domain_size_so_far = current_candidate_domain_size;
               

               // we then track the degree of new best candidate
                  // note: as long as a candidate wins MRV, we will automatically grab its degree and call it the best degree so far
                     // which means that `the best degree` is not necessary the largest degree but the one associated with the best MRV candidate
               int current_candidate_degree = 0;
               for(int each_constrained_adjacent_idx : constraint_adjacency[i]){
                   if(!assigned[each_constrained_adjacent_idx]){++current_candidate_degree;}
               }
               best_degree_so_far = current_candidate_degree;

          }
          else if(current_candidate_domain_size == best_domain_size_so_far){ // we have a tie on the MRV... 
               
              // so we need to break it with the degree heuristic


              // we find the degree of the current candidate
              int current_candidate_degree = 0;
              for(int each_constrained_adjacent_idx : constraint_adjacency[i]){
                if(!assigned[each_constrained_adjacent_idx]){++current_candidate_degree;}
              }
              
              // then we compare degrees --- current vs best -- to try to break the tie
              if(current_candidate_degree > best_degree_so_far){
                   best_var_so_far_idx = i;
                   best_degree_so_far = current_candidate_degree;
              }
          }
          // done with comparison logic

       }
       // done scanning through all the variables (aka cells) to try and find the next best variable to assign

       return best_var_so_far_idx;


    } // end of select_unassigned_variable method definition

    bool forward_check(int just_assigned_var_idx, vector<pair<int, int>>& inferences){ // start of forward_check method definition

        /* after we just assigned a value to a variable (aka cell), we need to prune inconsistent values from the domains of the constrained adjacent variables (aka cells) */
        
        for(int constrained_adjacent_idx : constraint_adjacency[just_assigned_var_idx]){ // loop through each constrained adjacent variable (aka cell) of the just-assigned variable

            if(assigned[constrained_adjacent_idx]){continue;} // if the constrained adjacent cell is already assigned then skip it

            for(int v = 0; v < 2; ++v){ // check both values {0,1} for the constrained adjacent cell's domain

                if(!domain[constrained_adjacent_idx][v]){continue;} // if the value is not allowed for the constrained adjacent cell then skip it
                
                // 'HYPOTHETICALLY' assign the value (either 0 or 1) to the constrained adjacent cell for 'testing' purposes
                assigned[constrained_adjacent_idx] = true;
                int old_value = assignment_value[constrained_adjacent_idx];
                assignment_value[constrained_adjacent_idx] = v; // 'hypothetical assignment' of the value (either 0 or 1) to the constrained adjacent cell

                bool ok = check_all_constraints(); // CHECK (ie test) if the 'hypothetical assignment' of the value (either 0 or 1) to the constrained adjacent cell is consistent with all the constraints
                
                // UNDO the 'hypothetical assignment' of the value (either 0 or 1) to the constrained adjacent cell
                assignment_value[constrained_adjacent_idx] = old_value;
                assigned[constrained_adjacent_idx] = false;

                if(!ok){ // if the 'hypothetical assignment' of the value (either 0 or 1) to the constrained adjacent cell is INCONSISTENT (NOT OK) with all the constraints...
                    domain[constrained_adjacent_idx][v] = false; // prune (remove) this value (either 0 or 1) from the constrained adjacent cell's domain
                    inferences.push_back(make_pair(constrained_adjacent_idx, v)); // at this search level, we removed value v from the domain of this constrained adjacent cell
                                                                                     // forward checking modifies domains so if we were to backtrack (go back up a level in the search tree), we would need to restore the domain of this constrained adjacent cell to its original state

                    if(!domain[constrained_adjacent_idx][0] && !domain[constrained_adjacent_idx][1]){return false;} // if the constrained adjacent cell's domain is empty --> domain wipeout (this branch of the search tree is dead-end meaning we have to backtrack) --> return false
                }
            }
            // done checking all the values {0,1} for this particular constrained adjacent cell's domain
        }
        // done checking all the constrained adjacent cells

        return true; // reaching here means that no constrained adjacent cell's domain was wiped out by the forward checking process --> forward checking was successful (we can continue down this branch of the search tree) --> return true

    } // end of forward_check method definition


    // KEY: CORE BACKTRACKING SEARCH
    bool backtrack_search(int depth){ // start of backtrack_search method definition

        /* should match quite close to the figure5 pseudocode */
        
        ++nodes_generated; // increment the number of nodes generated during the search

        if(is_complete()){ // check that all 81 variables (aka cells) have been assigned a value
            goal_depth = depth;
            return true;
        }

        int var_idx_to_assign = select_unassigned_variable(); // use MRV and Degree heuristic to SELECT THE NEXT variable (aka cell) to assign
        

        for(int v = 0; v < 2; ++v){ // for each value in ORDER-DOMAIN-VALUES

            if(!domain[var_idx_to_assign][v]){continue;}
            
            // add {var = value} to the assignment
            assigned[var_idx_to_assign] = true;
            assignment_value[var_idx_to_assign] = v;

            if(check_all_constraints()){ // if v is consistent with assignment

                vector<pair<int, int>> inferences;
                if(forward_check(var_idx_to_assign, inferences)){ // if inference != failure
                    
                    // KEY: recursive call -- go one layer deeper in the search tree
                       // if result != failure then return result
                    if(backtrack_search(depth + 1)){ 
                        return true; 
                    }
                }

                // remove inference from csp
                for(const auto& pruned : inferences){ // undo all the pruned
                    int pruned_var_idx = pruned.first;
                    int pruned_value = pruned.second;
                    domain[pruned_var_idx][pruned_value] = true;
                }

            } // end of checking all constraints if-block


            // remove {var = value} from the assignment
            assigned[var_idx_to_assign] = false;
            assignment_value[var_idx_to_assign] = -1;

        } // end of scanning through all the values in ORDER-DOMAIN-VALUES

        return false; // return failure
                         // no value worked for this variable to assign so we have to backtrack
        
    } // end of backtrack_search method definition


}; // end of SudokuMineSolver class definition




// ======================================================================================================================================= 

int main(int argc, char *argv[]){ // start of main function
                                       // where all the action occurs

    if (argc != 3){ // we expect two command line arguments: the input file and the output file
                       // (note: the 3 includes the name of the program itself)
        cerr << "Error: Wrong number of arguments" << endl;
        return 1; // return 1 to indicate an error
    }

    ifstream in_file(argv[1]); // open the input file specified by the first command line argument
    if(!in_file){
        cerr << "Error: Could not open the input file" << endl;
        return 1; // return 1 to indicate an error
    }

    SudokuMineSolver solver; // instantiate an instance of the SudokuMineSolver class

    if(!solver.load_puzzle(in_file)){ // load the puzzle from the input file
        cerr << "Error: Could not load the puzzle from the input file" << endl;
        return 1; // return 1 to indicate an error
    }

    in_file.close(); 

    // RUN THE BACKTRACKING SEARCH
    if(!solver.solve()){ 
        cerr << "Error: No solution found for this puzzle" << endl;
        return 1; // return 1 to indicate an error
    }

    ofstream out_file(argv[2]);
    if(!out_file){
        cerr << "Error: Could not open the output file" << endl;
        return 1;
    }

    solver.write_output(out_file); // write the solution to the output file

    out_file.close();

    return 0; // return 0 to indicate success

} // end of main function