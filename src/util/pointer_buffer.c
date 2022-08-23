#include <stdlib.h>
#include <memory.h>
#include "pointer_buffer.h"

Pointer_buffer *Pointer_buffer_alloc()
{
    Pointer_buffer *buffer = (Pointer_buffer *)malloc(sizeof(Pointer_buffer));
    buffer->len = 0;
    buffer->_buffer_size = INIT_BUFFER_SIZE;
    buffer->buffer = malloc(INIT_BUFFER_SIZE * sizeof(void *));
}

int Pointer_buffer_insert(Pointer_buffer *self, void *item)
{
    void *newbuffer = NULL;
    if (self->len >= self->_buffer_size * 0.8)
    {
        newbuffer = malloc(sizeof(void *) * (self->_buffer_size) * 2);
        memcpy(newbuffer, self->buffer, self->_buffer_size * sizeof(void *));
        free(self->buffer);
        self->buffer = newbuffer;
        self->_buffer_size *= 2;
    }

    self->buffer[(self->len)] = item;
    self->len++;
}

void *Pointer_buffer_get(Pointer_buffer *self, int index)
{
    return self->buffer[index];
}

int Pointer_buffer_free(Pointer_buffer *self)
{
    free(self->buffer);
    free(self);
}