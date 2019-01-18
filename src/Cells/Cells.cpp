#include "Cells.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                              CELLS
// ------------------------------------------------------------------------- //
Cell::Cell()
{
    // initalise the paricle list
    particle_list = new LinkedList;

    // initialise the neightbour list
    neighbour_list.resize(5);
}

// destructor
Cell::~Cell()
{
    // delete the linked list
    delete particle_list;

    // delete the pointers to all the neightbours
    for(unsigned i=0; i<neighbour_list.size(); i++)
    delete neighbour_list[i];
    neighbour_list.clear();
}

// set the neightbour at position i
void Cell::set_neighbour(const unsigned& i, Cell* neighbour)
{
    // +---+---+---+
    // |   |   | 0 |
    // +---+---+---+
    // |   | 1 | 2 |
    // +---+---+---+
    // |   | 3 | 4 |
    // +---+---+---+

    // set the cell
    neighbour_list[i] = neighbour;
}

// get the neighbour
Cell* Cell::get_neighbour(const unsigned& i)
{
    return neighbour_list[i];
}

// return this cells head for looping over all the particles in this Cell
ListNode* Cell::get_particle_list_head()
{
    return particle_list->get_head();
}

// assign a particle to this cell
void Cell::assign_particle(Particle* particle_pt)
{
    particle_list->insert(particle_pt);
}

// add a listnode to this cell
void Cell::add_particle_listnode(ListNode* node_pt)
{
    particle_list->add_node(node_pt);
}

// remove a listnode from this cell
void Cell::remove_particle_listnode(ListNode* node_pt)
{
    particle_list->remove_node(node_pt);
}

// move a node from one cell to another
void Cell::move_list_node(ListNode* node_pt, Cell* new_cell)
{
    // remove the node from this cell
    this->remove_particle_listnode(node_pt);

    // and add the node to the new cell
    new_cell->add_particle_listnode(node_pt);
}
