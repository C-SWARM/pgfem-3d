#ifndef PRINTING_H
#define PRINTING_H

/** Prints Title V1. */
void PrintTitleV1();

/** Prints Title V2. */
void PrintTitleV2();

/** Prints file missing error (with processor info). */
void NoFileProc(char* filename, int proc);

/** Prints file missing error. */
void NoFile(char* filename);

/** Prints output file missing error (with processor info). */
void NoFileOutProc(char* filename, int proc);

/** Prints output file missing error. */
void NoFileOut(char* filename);

#endif
