#include "random.h"
#include "util.h"
int main(int argc, char const *argv[])
{
    const int start = 10;
    const int end = 20;
    int t;
    for (int i = 0; i < 100; i++)
    {
        t = rand_int(start, end);
        if (t < start OR t >= end)
        {
            terminate("rand_int test faild");
        }
    }

    const double startd = 0.0;
    const double endd = 10.0;
    double td;
    for (int i = 0; i < 100; i++)
    {
        td = rand_double(startd, endd);
        if (td<startd OR td> endd)
        {
            terminate("rand_double test faild");
        }
    }

    return 0;
}
