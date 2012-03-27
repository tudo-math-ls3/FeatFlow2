#include <stdio.h>

/**
 * Allocate a block of memory on the heap
 */
int storage_malloc(void **ptr, size_t size)
{
  *ptr = (void*) malloc(size);
  if (*ptr == NULL)
    {
      printf("Unable to allocate memory by malloc()\n");
      return 1;
    } else
    {
      return 0;
    }
  return 0;
}

/**
 * Deallocate a block of memory on the heap.
 */
int storage_free(void **ptr)
{
  free(*ptr);
  *ptr = NULL;
  return 0;
}
