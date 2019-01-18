#ifndef LINKEDLIST_HPP
#define LINKEDLIST_HPP

#include "Molecules.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                              LINKED LIST NODE
// ------------------------------------------------------------------------- //
// node structure in the linked list
struct ListNode
{
  Particle* particle;
  ListNode* next;
  ListNode* previous;
};

// ------------------------------------------------------------------------- //
//                              LINKED LIST
// ------------------------------------------------------------------------- //
// linked list class
class LinkedList
{
public:
  LinkedList();
  // destructor
  ~LinkedList();
  // node which returns the head of the list
  ListNode* get_head();
  // function for inserting a particle in the list
  void insert(Particle* particle_pt);
  // function which adds a node to a linked list
  void add_node(ListNode* node_pt);
  // function for removing a particle from the list
  void remove_node(ListNode* node_pt);

private:
  ListNode* head;
};

#endif
