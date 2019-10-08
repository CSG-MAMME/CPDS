#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct Node {
    Node *prev;
    Node *next;
    int id;
} Node;

typedef struct List {
    Node *first;
    Node *last;
    int size;
} List;

void traverse_list(List l)
{
    Node n;
    for (n = l->first; n; n = n->next)
    {
        #pragma omp task
        printf("Processing element: %i from thread %d\n", n->id, 
                omp_get_thread_num());
    }
}

int main()
{
    List l;
    #pragma omp parallel num_threads(8)
    {
        #pragma single
        {
            traverse_list(l);
        }
    }
}
