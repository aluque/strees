/* Miscellaneous utilities.
 * Author: Alejandro Luque, 2004
 * License: GPL */

#include <stdlib.h>
#include <stdio.h>

extern char *invok_name;

/*  Allocates memory with check and (eventually) error reporting. */
void*
xmalloc (size_t size)
{
  register void *value = malloc(size);
  if (value == 0) {
    fprintf (stderr, "%s: virtual memory exhausted", invok_name);
    exit (1);
  }
  return value;
}

/* Reallocates memory */
void *
xrealloc (void *ptr, size_t size)
{
  register void *value = realloc (ptr, size);
  if (value == 0) {
    fprintf (stderr, "%s: Virtual memory exhausted", invok_name);
    exit (1);
  }
  return value;
}


/*  Same as before, but now initializes the memory to zero. */
void*
xcalloc (size_t count, size_t size)
{
  register void *value = calloc (count, size);
  if (value == 0){
    fprintf (stderr, "%s: virtual memory exhausted", invok_name);
    exit (1);
  }
  return value;
}
