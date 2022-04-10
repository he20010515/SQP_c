#define LOG_TAG "main"

#include "elog.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>

static void test_elog(void);

int main(void)
{
    /* close printf buffer */
    sqp_init();
    test_elog();

    return EXIT_SUCCESS;
}

/**
 * EasyLogger demo
 */
void test_elog(void)
{
    while (true)
    {
        /* test log output for all level */
        log_a("Hello EasyLogger!");
        log_e("Hello EasyLogger!");
        log_w("Hello EasyLogger!");
        log_i("Hello EasyLogger!");
        log_d("Hello EasyLogger!");
        log_v("Hello EasyLogger!");
        //        elog_raw("Hello EasyLogger!");
        break;
    }
}
