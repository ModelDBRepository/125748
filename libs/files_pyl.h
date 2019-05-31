//
// files.h
//
// Written by Stefano Fusi and G. La Camera (10/10/2000)
// Used and Revised on 13 Oct 2003 by Michele Giugliano, PhD (info and bug reports to michele@giugliano.info)
//

#include <stdlib.h> 
#include <stdio.h>

// ====================================================
// extract - read the buffer
// returns the number of numbers read from the buffer
// 'buf', otherwise zero if the buffer is empty 
// (by S Fusi)
int extract(double *v,char *buf);


// =============================================
// readline: from 1-line 'file' to vector 'vect' 
// read only 1-line files;
// returns the number of numbers read, 
// otherwise -1 if there's no such file to open;
// uses extract();
int readline(char *file, double vect[]);

// readmtrix: read a matrix-formatted file
// read from the file 'file', and put it into the matrix 'matr'
// which has 'rownum' rows and 'colnum' columns
int readmatrix(char *file, double **matr, int *rownum, int *colnum);
