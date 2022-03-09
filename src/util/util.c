#include <stdio.h>
#define AND &&
#define OR ||

void terminate(char *string)
{
    fprintf(stdout, "\n%s\n", string);
    fprintf(stdout, "The program is exiting now. . . .\n\n");
    exit(-1);
}