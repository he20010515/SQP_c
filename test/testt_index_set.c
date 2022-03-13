#include "index_set.h"

int main(int argc, char const *argv[])
{
    Index_set *set_A = index_set_alloc(4);
    Index_set *set_B = index_set_alloc(4);
    Index_set *A_B = index_set_alloc(set_A->index_range);
    index_set_append(set_A, 0);
    index_set_append(set_A, 1);
    index_set_append(set_B, 2);
    index_set_append(set_B, 3);
    index_set_append(set_B, 2);

    index_set_union(set_A, set_B, A_B);
    Index_set *A_I_B = index_set_alloc(set_A->index_range);
    index_set_intersection(set_A, set_B, A_I_B);
    printf("A Union B\n");
    index_set_print(A_B);
    printf("A intersection B\n");
    index_set_print(A_I_B);

    return 0;
}
