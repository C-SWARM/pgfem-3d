#include <stdio.h>
#include <stdlib.h>

#include "sort_container.h"
#include "renumber_utils.h"

int compare_integer (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int compare_sort_container (const void * a, const void * b)
{
  return ((*(struct sort_container*)a).value - (*(struct sort_container *)b).value);
}

int compare_sort_container_reverse (const void * a, const void * b)
{
  return ((*(sort_container*)b).value - (*(sort_container *)a).value);
}

int compare_sort_container_revert (const void * a, const void * b)
{
  return ((*(struct sort_container*)a).index - (*(struct sort_container *)b).index);
}
