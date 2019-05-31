//
// files.c
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

int extract(double *v,char *buf)
{
  int n=0;
  char *s;

  s=buf;

  // reach the first number (if any) of the current buffer
  while(*s && (*s==' ' || *s=='\t')) s++; // skip blanks
  if(*s==0 || *s=='\n') return 0; // no numbers

  // read the number
  while(*s)
    {
      v[n]=atof(s);
      n++;
      while(*s && !(*s==' ' || *s=='\t')) s++; // skip non-blanks
      while(*s && (*s==' ' || *s=='\t' || *s=='\n')) s++; // skip blanks
      if(*s==0)  break;
    }
  return n;  // Number of columns of the datafile
}


// =============================================
// readline: from 1-line 'file' to vector 'vect' 
// read only 1-line files;
// returns the number of numbers read, 
// otherwise -1 if there's no such file to open;
// uses extract();

int readline(char *file, double vect[])
{
   FILE *dev;
   int num;
   char buf[1000];  // tmp buffer

   dev = fopen(file,"r");
   if (dev == NULL)
      {
         printf("\n\n>>> ERROR: can't open '%s'\n\n",file);
         return -1;
      }
   
   fgets(buf,1000,dev);
   num = extract(vect,buf); //
   
   fclose(dev);

   return num; // normal exit
}

// readmtrix: read a matrix-formatted file
// read from the file 'file', and put it into the matrix 'matr'
// which has 'rownum' rows and 'colnum' columns
int readmatrix(char *file, double **matr, int *rownum, int *colnum)
{
   FILE *dev;
   int i=0;
   char buf[1000]; // tmp buffer
   
   dev = fopen(file,"r");
   if (dev == NULL)
      {
         printf("\n\n>>> ERROR: can't open '%s'\n\n",file);
         return -1;
      }

   while(fgets(buf,1000,dev))
      {
         *colnum = extract(matr[i],buf);
         i++;
      }
   *rownum = i;

   fclose(dev); 
   return 0; // normal exit  
}
