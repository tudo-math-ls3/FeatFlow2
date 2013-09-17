#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

#include "../feat2macros.h"

/* This function does about the same as "mkdir -p":
 * It assumes *path_name to contain the path to a directory.
 * All the nodes in the path_name-part are created (if they do not yet exist),
 * duplicate slashes are ignored.
 *
 * "/path//to/" , "/path/to/" and "//path//to/" all should give the same result.
 *
 * Function fails on first (unrecoverable) error.
 *
 * Upon return ierr contains:
 * (ierr >= 0) := number of subdirectories created
 * (ierr  < 0) := -errno
 *
 * Based on an implementation discussed in newsgroup comp.unix.programmer in April 2010
 */
void mkdir_recursive(const char *path_name, int* ierr)
{
  char *buff;
  char *slash;
  size_t len;
  int rc = 0;
  int count = 0;
  struct stat statbuff;
#if FEAT2_PP_OS_IS_WIN()
  const char *separator = "/\\";
#else
  const char *separator = "/";
#endif


  len = strlen(path_name);  /* der Wert, der rauskommt, ist zu groß */
  buff = (char *)malloc((len+1) * sizeof(char *));
  if (buff == NULL) {
    *ierr = -EINVAL;
    return;
  }

  /* Create a copy of the given path/file name because we will
   * - add a trailing slash and
   * - blank every separator occurrence later, at least temporarily */
  memcpy(buff, path_name, len);
  /* strncpy(buff, path_name, len); */
  buff[len]   = '/';
  buff[++len] = 0;

  for (slash = buff; (slash = strpbrk(slash+1, separator)); ) {

    /* this is to catch double / in paths */
    if (
#if FEAT2_PP_OS_IS_WIN()
	slash[-1] == '\\' ||
#endif
	slash[-1] == '/') {
      continue;
    }
    *slash = 0;

    rc = stat(buff, &statbuff);
    if (!rc) {
      /* path already exists, skip it */
      *ierr = EEXIST;
    } else {
      rc = mkdir(buff, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      *ierr = (rc) ? errno : 0;
    }

    switch(*ierr) {
      case 0:
	count++;
      case EEXIST:
	break;
      case ENOTDIR:
      case EACCES:
      default:
	goto quit;
    }
    *slash = '/';
  }

quit:
  *ierr = (*ierr) ? -(*ierr) : count;
  free(buff); buff = NULL;
}


void mkdir_recursive_(const char *path_name, int* ierr)
{
  mkdir_recursive(path_name, ierr);
}
