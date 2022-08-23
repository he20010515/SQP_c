#define INIT_BUFFER_SIZE 16
typedef struct pointer_buffer
{
    int len;
    int _buffer_size;
    void **buffer;
} Pointer_buffer;

Pointer_buffer *Pointer_buffer_alloc();
int Pointer_buffer_insert(Pointer_buffer *self, void *item);
void *Pointer_buffer_get(Pointer_buffer *self, int index);
int Pointer_buffer_free(Pointer_buffer *self);
