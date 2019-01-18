#include "LinkedList.hpp"

using namespace::std;

// ------------------------------------------------------------------------- //
//                              LINKED LIST
// ------------------------------------------------------------------------- //
LinkedList::LinkedList()
{
    // set the head to be NULL
    head = NULL;
}

// destructor
LinkedList::~LinkedList()
{
    // loop through the list
    ListNode* conductor = head;

    while(conductor != NULL)
    {
        ListNode* next = conductor->next;
        delete conductor;
        conductor = next;
    }

    // delete the head
    delete head;
}

// get head of the Linked List
ListNode* LinkedList::get_head() { return head; }

// function for inserting a particle in the list
void LinkedList::insert(Particle* particle_pt)
{
    // construct a new node
    ListNode* new_node = new ListNode;

    // update the previous of the next node
    if(head != NULL)
        head->previous = new_node;

    // update this new node with the correct information
    new_node->particle = particle_pt;
    new_node->next = head;
    new_node->previous = NULL;

    // update the head of the list
    head = new_node;
}

// function which adds a node to a linked list
void LinkedList::add_node(ListNode* node_pt)
{
    // update the next nodes previous pointer
    if(head != NULL)
        head->previous = node_pt;

    // update the next and previous pointer
    node_pt->next = head;
    node_pt->previous = NULL;

    // update the head pointer
    head = node_pt;
}

// function for removing a particle from the list
void LinkedList::remove_node(ListNode* node_pt)
{
    // exchange the nodes
    if(node_pt->next != NULL)
        node_pt->next->previous = node_pt->previous;

    // Check if the previous element is NULL, if it is: then we are removing
    // the first element from the linked list and must therefore reset the
    // head of the linked list
    if(node_pt->previous != NULL)
        node_pt->previous->next = node_pt->next;
    else
        head = node_pt->next;
}
