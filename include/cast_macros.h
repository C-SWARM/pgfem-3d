/* HEADER */
#ifndef CAST_MACROS_H
#define CAST_MACROS_H

/** Define macros to safely cast types. This is mostly to help compile
    without warnings and ensure that const-ness is respected. */

#define CONST_2(type) (const type **)
#define CONST_3(type) (const type ***)
#define CONST_4(type) (const type ****)

#endif /* #ifndef CAST_MACROS_H */
