#include "vector.h"
#pragma once
#define AND &&
#define OR ||
#define DOUBLE_ERROR 0.01
void terminate(char *string);
int double_equal(double a, double b);
int vector_any_bigger_equal_than_const(const Vector *v, double a);
void sqp_init(void);

#define MAX(x, y) ({    \
    typeof(x) _x = x;   \
    typeof(y) _y = y;   \
    (void)(&_x == &_y); \
    _x > _y ? _x : _y;  \
})

#define MIN(x, y) ({    \
    typeof(x) _x = x;   \
    typeof(y) _y = y;   \
    (void)(&_x == &_y); \
    _x > _y ? _y : _x;  \
})

#define TRUE 1
#define FALSE 0