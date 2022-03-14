#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#pragma once
struct index_set
{
    int index_range;
    int *elem;
};
typedef struct index_set Index_set;
Index_set *index_set_alloc(int size);
void index_set_free(Index_set *set);
void index_set_append(Index_set *set, const int i);
void index_set_remove(Index_set *set, const int i);
int index_set_is_in(const Index_set *set, int i);
void index_set_union(Index_set *A, Index_set *B, Index_set *A_U_B);
void index_set_intersection(Index_set *A, Index_set *B, Index_set *A_I_B);
void index_set_print(Index_set *set);
int index_set_size(Index_set *set);
