#include <stdio.h>
void terminate(char *string)
{
    fprintf(stdout, "\n%s\n", string);
    fprintf(stdout, "The program is exiting now. . . .\n\n");
    exit(-1);
}