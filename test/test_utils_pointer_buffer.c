#include "pointer_buffer.h"
#include "stdio.h"
int main(int argc, char const *argv[])
{
    double array[64] = {0};
    Pointer_buffer *buffer = Pointer_buffer_alloc();

    for (int i = 0; i < 64; i++)
    {
        array[i] = i * 1.5;
        Pointer_buffer_insert(buffer, &(array[i]));
    }

    for (int i = 0; i < 64; i++)
    {
        double temp = *(double *)Pointer_buffer_get(buffer, i);
        printf("buffer[%d] = %lf\n", i, temp);
    }

    Pointer_buffer_free(buffer);
    return 0;
}
