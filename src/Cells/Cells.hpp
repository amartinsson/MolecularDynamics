#ifndef CELLS_HPP
#define CELLS_HPP

#include <vector>
#include <unordered_map>

#include "LinkedList.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                              CELLS CLASS
// ------------------------------------------------------------------------- //
class Cell
{
public:
    Cell();
    // destructor
    ~Cell();
    // set the neightbour at position i
    void set_neighbour(const unsigned& i, Cell* neighbour);
    // get the neighbour
    Cell* get_neighbour(const unsigned& i);
    // return this cells head for looping over all the particles in this Cell
    ListNode* get_particle_list_head();
    // assign a particle to this cell
    void assign_particle(Particle* particle_pt);
    // add a listnode to this cell
    void add_particle_listnode(ListNode* node_pt);
    // remove a listnode from this cell
    void remove_particle_listnode(ListNode* node_pt);
    // move a node from one cell to another
    void move_list_node(ListNode* node_pt, Cell* new_cell);

private:
    unsigned cell_number;
    // particle list
    LinkedList* particle_list;
    // neightbour list
    // vector<Cell*> neighbour_list;
    unordered_map<int, Cell*> neighbour_list;
};

#endif
