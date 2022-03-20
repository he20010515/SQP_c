#include "index_set.h"

Index_set *index_set_alloc(int size)
{
    Index_set *s = (Index_set *)malloc(sizeof(Index_set));
    s->index_range = size;
    s->elem = (int *)malloc(size * sizeof(int));
    for (int i = 0; i < size; i++)
    {
        s->elem[i] = 0;
    }
    return s;
}
void index_set_free(Index_set *set)
{
    free(set->elem);
    free(set);
}

void index_set_append(Index_set *set, const int i)
{
    if (i < 0 OR i >= set->index_range)
    {
        terminate("ERROR: index_set_append index out of range");
    }
    if (set->elem[i] == 1)
    {
        printf("Warning: index_set_append you are trying to append a exitst element");
    }
    set->elem[i] = 1;
}

void index_set_remove(Index_set *set, const int i)
{
    if (i < 0 OR i >= set->index_range)
    {
        terminate("ERROR: index_set_remove index out of range");
    }
    if (set->elem[i] == 0)
    {
        printf("Warning: index_set_remove you are trying to remove a not exitst element");
    }

    set->elem[i] = 0;
    return;
}

int index_set_is_in(const Index_set *set, int i)
{
    if (i < 0 OR i >= set->index_range)
    {
        terminate("ERROR: index_set_insert index out of range");
    }
    return set->elem[i] == 1;
}

void index_set_union(Index_set *A, Index_set *B, Index_set *A_U_B)
{
    if (A->index_range != B->index_range)
    {
        terminate("ERROR: index_set_union set A and B index_range NOT fit");
    }
    for (int i = 0; i < A->index_range; i++)
    {
        if (index_set_is_in(A, i) OR index_set_is_in(B, i))
        {
            index_set_append(A_U_B, i);
        }
    }
    return;
}

void index_set_intersection(const Index_set *A, const Index_set *B, Index_set *A_I_B)
{
    if (A->index_range != B->index_range)
    {
        terminate("ERROR: index_set_union set A and B index_range NOT fit");
    }
    for (int i = 0; i < A->index_range; i++)
    {
        if (index_set_is_in(A, i) AND index_set_is_in(B, i))
        {
            index_set_append(A_I_B, i);
        }
    }
    return;
}

void index_set_print(Index_set *set)
{
    printf("{");
    for (int i = 0; i < set->index_range; i++)
    {
        if (index_set_is_in(set, i))
        {
            printf("%d, ", i);
        }
    }
    printf("}\n");
    return;
}

int index_set_size(const Index_set *set)
{
    int s = 0;
    for (int i = 0; i < set->index_range; i++)
    {
        if (index_set_is_in(set, i))
        {
            s++;
        }
    }
    return s;
}