/* HEADER */
#pragma once 
#ifndef GEN_PATH_H
#define GEN_PATH_H

#include <sys/types.h> /* for stat types */

#define DIR_MODE ((mode_t) 0777)

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

/** Make a single directory. Return 0 on success or -1 on error. */
int make_dir(const char *path, mode_t mode);

/** Make a path. Return 0 on success or -1 on error. */
int make_path(const char *path, mode_t mode);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif
