#include <stdio.h>
#include <math.h>
#include "util.h"

void terminate(char *string)
{
    fprintf(stdout, "\n%s\n", string);
    fprintf(stdout, "The program is exiting now. . . .\n\n");
    exit(-1);
}

int double_equal(double a, double b)
{
    if (abs(a - b) < DOUBLE_ERROR)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}