#include "gen_path.h"
#include <sys/stat.h> /* for mkdir etc */
#include <sys/errno.h> /* for errno */
#include <stdlib.h>
#include <string.h>

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

#ifndef ALLOCATION_H
#include "allocation.h"
#endif

int make_dir(const char *path, mode_t mode)
{
  struct stat st;
  int status = 0;

  /* check the status of the file in path */
  if(stat(path,&st) != 0){ /* file does not exist */

    if(mkdir(path,mode) != 0 && errno != EEXIST){
      /* The directory could not be created and the error is something
     other than the file already existing */
      status = -1;
    }

  } else if(!S_ISDIR(st.st_mode)){ /* file exists and is not a directory */
    status = -1;
  }
  return status;
}

int make_path(const char *path, mode_t mode)
{

  struct stat st;
  char *pp;
  char *sp;
  int status, full_path;
  char *copypath = PGFEM_calloc (char, strlen(path)+1);
  strcpy(copypath,path);

  status = 0;
  full_path = 0;
  pp = copypath;

  /* first check to see if the path up to the last dir exists */
  if((sp = strrchr(pp, '/')) == NULL){ /* path in current directory */
    full_path = 1;
  } else if(sp != pp){
    /* not root dir */
    *sp = '\0';
    if(stat(copypath,&st) == 0 && S_ISDIR(st.st_mode)){
      full_path = 1;
    }
    *sp = '/';
  }

  switch (full_path){
  case 0:
    while (status == 0 && (sp = strchr(pp, '/')) != 0){
      if (sp != pp){
    /* Neither root nor double slash in path */

    *sp = '\0'; /* replace '/' with '\0' to terminate the string */

    status = make_dir(copypath, mode);

    *sp = '/'; /* replace the '/' character to regain the full string */
      }
      pp = sp + 1; /* advance the pointer to after the '/' */
    }
    /* Drop through to case 1 */

  case 1:
    if (status == 0)
      status = make_dir(path, mode);
    PGFEM_free(copypath);
    break;
  }

  return (status);
}

/// check directory exist
///
/// \param[in] dir_name directory path
/// \return 1 if input dir_name exists, or 0
int is_directory_exist(const char *dir_name)
{
  struct stat sb;
  if(stat(dir_name, &sb) == 0 && S_ISDIR(sb.st_mode))
    return 1;
  else
    return 0;
}

#ifdef GEN_PATH_TESTING
int main()
{
  int err = 0;
  char out_path[] = "./test_dir";
  char out_dir[500];

  sprintf(out_dir,"%s/VTK",out_path);
  err = make_path(out_dir,0777);
  if(err != 0){
    PGFEM_printerr("Directory not created! (%d: %s)\n",
           errno, strerror(errno));
    abort();
  }

  for(int i=0; i<11; i++){
    sprintf(out_dir,"%s/VTK/step_%0.5d",out_path,i);
    err = make_path(out_dir,0777);
    if(err != 0){
      PGFEM_printerr("Directory not created! (%d: %s)\n",
             errno, strerror(errno));
      abort();
    }
  }
  return 0;
}
#endif
