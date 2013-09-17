#include <string.h>

/* Several functions to detect whether a given string identifies a valid
 * file system object.
 * Most Fortran 90 compilers provide an FSTAT intrinsic, but not the macros
 * S_ISLNK(), S_ISDIR(), S_IFREG(). One could use the bits used in the definition
 * of those macros, but (excerpt from 'man 2 stat'):
 *  "POSIX does not describe the S_IFMT, S_IFSOCK, S_IFLNK, S_IFREG, S_IFBLK,
 *   S_IFDIR, S_IFCHR, S_IFIFO, S_ISVTX bits, but instead demands the use of
 *   the macros S_ISDIR(), etc."
 * Confirmed by reading copy of POSIX:2008 specification at
 * http://www.opengroup.org/onlinepubs/9699919799/basedefs/sys_stat.h.html
 */

/* Detect whether given string identifies an existing directory */

#ifdef _MSC_VER /* Windows implementation */

#include <Windows.h>

void isdirectory(char *name, int *isdir, int *ierr)
{
  *ierr = 0;
  *isdir = 0;

  if(strlen(name) <= 0)
    *ierr = 1;
  else
  {
    DWORD attribs = GetFileAttributes(name);
    if(attribs == INVALID_FILE_ATTRIBUTES)
    {
      /* failed to retrieve attibutes */
      *ierr = 2;
    }
    else
    {
      /* if bit 4 is set, the name identifies a directory */
      *isdir = ((attribs & FILE_ATTRIBUTE_DIRECTORY) != 0) ? 1 : 0;
    }
  }
}

#else /* POSIX conformant implementation */

#include <sys/stat.h>

void isdirectory(char *name, int *isdir, int *ierr)
{
  struct stat sb;
  *ierr = 0;
  *isdir = 0;

  if (strlen(name) <= 0)
    *ierr = 1;
  else {
    if (stat(name, &sb) == -1) {
      /* error occurred in stat command */
      *ierr = 2;
    } else {
      /* working code, but does not conform to C99, i.e.
       * gives an error when compiling with 'gcc -std=c99'. */
      /* switch (sb.st_mode & S_IFMT) {
          case S_IFDIR:
            *isdir = 1;
            break;
          case S_IFBLK:
          case S_IFCHR:
          case S_IFIFO:
          case S_IFLNK:
          case S_IFREG:
          case S_IFSOCK:
          default:
            break;
  } */

      /* conforming to C99 */
      *isdir = S_ISDIR(sb.st_mode);
    }
  }
}
#endif

void isdirectory_(char *name, int *isdir, int *ierr)
{
  isdirectory(name, isdir, ierr);
}
