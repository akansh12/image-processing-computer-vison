#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>

#define PI 3.1415927 
#define EPSILON 10E-7

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                 CANNY EDGE DETECTOR FOR GREY VALUE IMAGES                */
/*                                                                          */
/* (Copyright by Markus Mainberger, Pascal Peter, Joachim Weickert 5/2022)  */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*
  Canny edge detector for grey value images
*/

/*--------------------------------------------------------------------------*/

void alloc_double_vector

     (double **vector,   /* vector */
      long   n1)         /* size */

/*
  allocates memory for a double format vector of size n1
*/

{
*vector = (double *) malloc (n1 * sizeof(double));

if (*vector == NULL)
   {
   printf("alloc_double_vector: not enough memory available\n");
   exit(1);
   }

return;

}  /* alloc_double_vector */

/*--------------------------------------------------------------------------*/

void alloc_double_matrix

     (double ***matrix,  /* matrix */
      long   n1,         /* size in direction 1 */
      long   n2)         /* size in direction 2 */

/*
  allocates memory for a double format matrix of size n1 * n2
*/

{
long i;    /* loop variable */

*matrix = (double **) malloc (n1 * sizeof(double *));

if (*matrix == NULL)
   {
   printf("alloc_double_matrix: not enough memory available\n");
   exit(1);
   }

for (i=0; i<n1; i++)
    {
    (*matrix)[i] = (double *) malloc (n2 * sizeof(double));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_double_matrix: not enough memory available\n");
       exit(1);
       }
    }

return;

}  /* alloc_double_matrix */

/*--------------------------------------------------------------------------*/

void free_double_vector

     (double  *vector,    /* vector */
      long    n1)         /* size */

/*
  frees memory for a double format vector of size n1
*/

{

free(vector);
return;

}  /* free_double_vector */

/*--------------------------------------------------------------------------*/

void free_double_matrix

     (double  **matrix,   /* matrix */
      long    n1,         /* size in direction 1 */
      long    n2)         /* size in direction 2 */

/*
  frees memory for a double format matrix of size n1 * n2
*/

{
long i;   /* loop variable */

for (i=0; i<n1; i++)
free(matrix[i]);

free(matrix);

return;

}  /* free_double_matrix */

/*--------------------------------------------------------------------------*/

void read_string

     (char *v)         /* string to be read */

/*
  reads a string v
*/

{
if (fgets (v, 80, stdin) == NULL)
{
   printf("could not read string, aborting\n");
   exit(1);
}

if (v[strlen(v)-1] == '\n')
   v[strlen(v)-1] = 0;

return;

}  /* read_string */

/*--------------------------------------------------------------------------*/

void read_double

     (double *v)         /* value to be read */

/*
  reads a double value v
*/

{
char   row[80];    /* string for reading data */

if (fgets (row, 80, stdin) == NULL)
{
   printf("could not read string, aborting\n");
   exit(1);
}

if (row[strlen(row)-1] == '\n')
   row[strlen(row)-1] = 0;
sscanf(row, "%lf", &*v);

return;

}  /* read_double */

/*--------------------------------------------------------------------------*/

void skip_white_space_and_comments

     (FILE *inimage)  /* input file */

/*
  skips over white space and comments while reading the file
*/

{

int   ch = 0;   /* holds a character */
char  row[80];  /* for reading data */

/* skip spaces */
while (((ch = fgetc(inimage)) != EOF) && isspace(ch));

/* skip comments */
if (ch == '#')
   {
   if (fgets(row, sizeof(row), inimage))
      skip_white_space_and_comments (inimage);
   else
      {
      printf("skip_white_space_and_comments: cannot read file\n");
      exit(1);
      }
   }
else
   fseek (inimage, -1, SEEK_CUR);

return;

} /* skip_white_space_and_comments */

/*--------------------------------------------------------------------------*/

void read_pgm_to_double

     (const char  *file_name,    /* name of pgm file */
      long        *nx,           /* image size in x direction, output */
      long        *ny,           /* image size in y direction, output */
      double      ***u)          /* image, output */

/*
  reads a greyscale image that has been encoded in pgm format P5 to
  an image u in double format;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
char  row[80];      /* for reading data */
long  i, j;         /* image indices */
long  max_value;    /* maximum color value */
FILE  *inimage;     /* input file */

/* open file */
inimage = fopen (file_name, "rb");
if (inimage == NULL)
   {
   printf ("read_pgm_to_double: cannot open file '%s'\n", file_name);
   exit(1);
   }

/* read header */
if (fgets(row, 80, inimage) == NULL)
   {
   printf ("read_pgm_to_double: cannot read file\n");
   exit(1);
   }

/* image type: P5 */
if ((row[0] == 'P') && (row[1] == '5'))
   {
   /* P5: grey scale image */
   }
else
   {
   printf ("read_pgm_to_double: unknown image format\n");
   exit(1);
   }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", nx))
   {
   printf ("read_pgm_to_double: cannot read image size nx\n");
   exit(1);
   }

/* read image size in x direction */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", ny))
   {
   printf ("read_pgm_to_double: cannot read image size ny\n");
   exit(1);
   }

/* read maximum grey value */
skip_white_space_and_comments (inimage);
if (!fscanf (inimage, "%ld", &max_value))
   {
   printf ("read_pgm_to_double: cannot read maximal value\n");
   exit(1);
   }
fgetc(inimage);

/* allocate memory */
alloc_double_matrix (u, (*nx)+2, (*ny)+2);

/* read image data row by row */
for (j=1; j<=(*ny); j++)
 for (i=1; i<=(*nx); i++)
     (*u)[i][j] = (double) getc(inimage);

/* close file */
fclose (inimage);

return;

}  /* read_pgm_to_double */

/*--------------------------------------------------------------------------*/

void comment_line

     (char* comment,       /* comment string (output) */
      char* lineformat,    /* format string for comment line */
      ...)                 /* optional arguments */

/* 
  Adds a line to the comment string comment. The string line can contain 
  plain text and format characters that are compatible with sprintf.
  Example call: 
  print_comment_line(comment, "Text %lf %ld", double_var, long_var).
  If no line break is supplied at the end of the input string, it is 
  added automatically.
*/

{
char     line[80];
va_list  arguments;

/* get list of optional function arguments */
va_start (arguments, lineformat);

/* convert format string and arguments to plain text line string */
vsprintf (line, lineformat, arguments);

/* add line to total commentary string */
strncat (comment, line, 80);

/* add line break if input string does not end with one */
if (line[strlen(line)-1] != '\n')
   strncat (comment, "\n", 1); 

/* close argument list */
va_end (arguments);

return;

}  /* comment_line */

/*--------------------------------------------------------------------------*/

void write_double_to_pgm

     (double  **u,          /* image, unchanged */
      long    nx,           /* image size in x direction */
      long    ny,           /* image size in y direction */
      char    *file_name,   /* name of pgm file */
      char    *comments)    /* comment string (set 0 for no comments) */

/*
  writes a greyscale image in double format into a pgm P5 file
*/

{
FILE           *outimage;  /* output file */
long           i, j;       /* loop variables */
double         aux;        /* auxiliary variable */
unsigned char  byte;       /* for data conversion */

/* open file */
outimage = fopen (file_name, "wb");
if (NULL == outimage)
   {
   printf("could not open file '%s' for writing, aborting\n", file_name);
   exit(1);
   }

/* write header */
fprintf (outimage, "P5\n");                  /* format */
if (comments != 0)
   fputs (comments, outimage);               /* comments */
fprintf (outimage, "%ld %ld\n", nx, ny);     /* image size */
fprintf (outimage, "255\n");                 /* maximal value */

/* write image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     aux = u[i][j] + 0.499999;    /* for correct rounding */
     if (aux < 0.0)
        byte = (unsigned char)(0.0);
     else if (aux > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(aux);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }

/* close file */
fclose (outimage);

return;

}  /* write_double_to_pgm */

/*--------------------------------------------------------------------------*/

void dummies_double

     (double **u,        /* image */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

/*
  creates dummy boundaries for a double format image u by mirroring
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    u[i][0]    = u[i][1];
    u[i][ny+1] = u[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    u[0][j]    = u[1][j];
    u[nx+1][j] = u[nx][j];
    }

return;

}  /* dummies_double */

/*--------------------------------------------------------------------------*/

void zero_boundaries

     (double **u,        /* image matrix */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

/*
  sets boundary pixels to zero
*/

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    u[i][0]    = 0.0;
    u[i][ny+1] = 0.0;
    }

for (j=0; j<=ny+1; j++)
    {
    u[0][j]    = 0.0;
    u[nx+1][j] = 0.0;
    }

return;

}  /* zero_boundaries */

/*--------------------------------------------------------------------------*/

void gauss_conv

    (double   sigma,     /* standard deviation of the Gaussian */
     long     btype,     /* type of boundary condition */
     double   prec,      /* cutoff at precision * sigma */
     long     nx,        /* image dimension in x direction */
     long     ny,        /* image dimension in y direction */
     double   hx,        /* pixel size in x direction */
     double   hy,        /* pixel size in y direction */
     double   **u)       /* input: original image ;  output: smoothed */


/*
  Gaussian convolution with a truncated and resampled Gaussian
*/


{
long    i, j, k, l, p;        /* loop variables */
long    length;               /* convolution vector: 0..length */
long    pmax;                 /* upper bound for p */
double  aux1, aux2;           /* time savers */
double  sum;                  /* for summing up */
double  *conv;                /* convolution vector */
double  *help;                /* row or column with dummy boundaries */


/* ----------------------- convolution in x direction -------------------- */

/* compute length of convolution vector */
length = (long)(prec * sigma / hx) + 1;

/* allocate memory for convolution vector */
alloc_double_vector (&conv, length+1);

/* compute entries of convolution vector */
aux1 = 1.0 / (sigma * sqrt(2.0 * 3.1415927));
aux2 = (hx * hx) / (2.0 * sigma * sigma);
for (i=0; i<=length; i++)
    conv[i] = aux1 * exp (- i * i * aux2);

/* normalisation */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate memory for a row */
alloc_double_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy u in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = u[i][j];

    /* extend signal according to the boundary conditions */
    k = length;
    l = length + nx - 1;
    while (k > 0)
          {
          /* pmax = min (k, nx) */
          if (k < nx)
             pmax = k;
          else
             pmax = nx;

          /* extension on both sides */
          if (btype == 0)
             /* reflecting b.c.: symmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = help[k+p-1];
                 help[l+p] = help[l-p+1];
                 }
          else
             /* Dirichlet b.c.: antisymmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = - help[k+p-1];
                 help[l+p] = - help[l-p+1];
                 }

          /* update k and l */
          k = k - nx;
          l = l + nx;
          }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* compute convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        u[i-length+1][j] = sum;
        }
    } /* for j */

/* free memory */
free_double_vector (help, nx+length+length);
free_double_vector (conv, length + 1);


/* ----------------------- convolution in y direction -------------------- */

/* compute length of convolution vector */
length = (long)(prec * sigma / hy) + 1;

/* allocate memory for convolution vector */
alloc_double_vector (&conv, length + 1);

/* compute entries of convolution vector */
aux1 = 1.0 / (sigma * sqrt(2.0 * 3.1415927));
aux2 = (hy * hy) / (2.0 * sigma * sigma);
for (j=0; j<=length; j++)
    conv[j] = aux1 * exp (- j * j * aux2);

/* normalisation */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate memory for a row */
alloc_double_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy u in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = u[i][j];

    /* extend signal according to the boundary conditions */
    k = length;
    l = length + ny - 1;
    while (k > 0)
          {
          /* pmax = min (k, ny) */
          if (k < ny)
             pmax = k;
          else
             pmax = ny;

          /* extension on both sides */
          if (btype == 0)
             /* reflecting b.c.: symmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = help[k+p-1];
                 help[l+p] = help[l-p+1];
                 }
          else
             /* Dirichlet b.c.: antisymmetric extension */
             for (p=1; p<=pmax; p++)
                 {
                 help[k-p] = - help[k+p-1];
                 help[l+p] = - help[l-p+1];
                 }

          /* update k and l */
          k = k - ny;
          l = l + ny;
          }

    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* compute convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        u[i][j-length+1] = sum;
        }
    } /* for i */

/* free memory */
free_double_vector (help, ny+length+length);
free_double_vector (conv, length+1);

return;

} /* gauss_conv */

/*--------------------------------------------------------------------------*/

void sobel

     (double  **u,        /* image, unchanged */
      long    nx,         /* x-dimension of u */
      long    ny,         /* y-dimension of u */
      double  hx,         /* pixel size in x-direction */
      double  hy,         /* pixel size in y-direction */
      double  **ux,       /* derivatives in x-direction */
      double  **uy)       /* derivatives in y-direction */

/*
  computes the derivatives in x and y direction with Sobel operators
*/

{
long    i, j;         /* loop variables */
double  aux1, aux2;   /* time savers */
double  aux3, aux4;   /* time savers */

/* mirror image boundaries */
dummies_double (u, nx, ny);

/* compute time savers */
aux1 = 0.25 / (2.0 * hx);
aux2 = 0.50 / (2.0 * hx);
aux3 = 0.25 / (2.0 * hy);
aux4 = 0.50 / (2.0 * hy);

/* compute derivatives with Sobel operators */
/*
  INSERT CODE
*/

return;

}  /* sobel */

/*--------------------------------------------------------------------------*/

double get_direction

     (double  x,            /* first component of vector */
      double  y)            /* second component of vector */

/*
  computes the direction of a vector (x,y) in radian;
  output in range [-PI/2, PI/2]
*/

{
double  direction;     /* direction */

if (abs(x) < EPSILON)
   if (abs(y) < EPSILON)
      direction = 0.0;
   else
      direction = 0.5 * PI;
else
   direction = atan (y / x);

return direction;

}  /* get_direction */

/*--------------------------------------------------------------------------*/

void heapsort

     (long    n,
      double  ra[])

/*
  Sorts an array ra[1..n] into ascending numerical order using the
  Heapsort algorithm.
  n is input; ra is replaced on output by its sorted rearrangement.
  Ref.:  Press et al: Numerical recipes in C. Second Edition,
         Section 8.3, Cambridge University Press, 1992.
*/

{
long    i, ir, j, l;
double  rra;

if (n < 2)
   return;
l = (n >> 1) + 1;
ir = n;

/* The index l will be decremented from its initial value down to 1  */
/* during the “hiring” (heap creation) phase. Once it reaches 1, the */
/* index ir will be decremented from its initial value down to 1     */
/* during the “retirement-and-promotion” (heap selection) phase.     */

for (;;)
    {
    if (l > 1)
       {
       /* Still in hiring phase. */
       rra = ra[--l];
       }
    else

       {
       /* ---- In retirement-and-promotion phase. ---- */

       /* Clear a space at end of array. */
       rra = ra[ir];

       /* Retire the top of the heap into it. */
       ra[ir] = ra[1];

       if (--ir == 1)
          /* Done with the last promotion.      */
          /* The least competent worker of all! */
          {
          ra[1] = rra;
          break;
          }
       } /* else */

    /* Whether in the hiring phase or promotion phase, we here */
    /* set up to sift down element rra to its proper level.    */
    i = l;
    j = l + l;

    while (j <= ir)
          {
          /* Compare to the better underling */
          if ((j < ir) && (ra[j] < ra[j+1]))
             j++;

          /* Demote rra. */
          if (rra < ra[j])
             {
             ra[i] = ra[j];
             i = j;
             j <<= 1;
             }
          else
             /* Found rra’s level. Terminate the sift-down */
             break;
          } /* while */

    /* Put rra into its slot. */
    ra[i] = rra;
    } /* for */

}  /* heapsort */

/*--------------------------------------------------------------------------*/

void nms

     (double  fn1,       /* value of neighbour pixel 1 */
      double  fc,        /* value of central pixel */
      double  fn2,       /* value of neighbour pixel 2 */ 
      double  *uc)       /* new value of central pixel, output */

/*
  applies nonmaxima suppression to the three values fn1, fc, fn2:
  sets result uc to 0 if fn1 or fn2 is larger than fc
*/

{
if ((fn1 > fc) || (fn2 > fc))
   *uc = 0.0;

return;

}  /* nms */

/*--------------------------------------------------------------------------*/

void nonmaxima_suppression

     (double   **u,        /* image */
      double   **ux,       /* x-derivatives */
      double   **uy,       /* y-derivatives */
      long     nx,         /* x-dimension of u */
      long     ny)         /* y-dimension of u */

/*
  performs nonmaxima suppression on all pixels of the image u
*/

{
long    i, j;       /* loop variables */
double  **f;        /* work copy of image u */
double  direction;  /* directional information */

/* allocate memory */
alloc_double_matrix (&f, nx+2, ny+2);

/* copy image u to f */
for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     f[i][j] = u[i][j];

/* nonmaxima suppression */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (f[i][j] != 0.0)
        {
        /* get edge direction in pixel (i,j) with values in [-PI/2, PI/2] */
        direction = get_direction (ux[i][j], uy[i][j]);

        /* apply nonmaxima suppression perpendicular to the edge */
        if ((direction >= - 0.125 * PI) && (direction <= 0.125 * PI))
           /* nonmaxima suppression in x direction */
           nms (f[i-1][j], f[i][j], f[i+1][j], &u[i][j]);

        else if ((direction > 0.125 * PI) && (direction <= 0.375 * PI))
           /* nonmaxima suppression in first diagonal direction */
           nms (f[i-1][j-1], f[i][j], f[i+1][j+1], &u[i][j]);

        else if ((direction < -0.125 * PI) && (direction >= -0.375 * PI))
           /* nonmaxima suppression in second diagonal direction */
           nms (f[i+1][j-1], f[i][j], f[i-1][j+1], &u[i][j]);

        else 
           /* nonmaxima suppression in y direction */
           nms (f[i][j-1], f[i][j], f[i][j+1], &u[i][j]);
        }

/* free memory */
free_double_matrix (f, nx+2, ny+2);

return;

}  /* nonmaxima_suppression */

/*--------------------------------------------------------------------------*/

void trace_edge

     (long    i,          /* x component of current pixel */
      long    j,          /* y component of current pixel */
      double  **u,        /* image */
      double  T1,         /* lower threshold */
      double  T2)         /* higher threshold */

/*
  Consider all 8 neighbours of a pixel (i,j) that is assumed to be >=T2. 
  If a neighbour is below the higher threshold T2 and above the lower 
  threshold T1, add this pixel to the edge pixels and repeat all steps 
  for that pixel.
*/

{
/* pixel u(i,j) is assumed to be >= T2 and thus an edge pixel */
u[i][j] = 255;

/* check neighbour (i+1,j) */
if ((u[i+1][j] <= T2) && (u[i+1][j] > T1))
   trace_edge (i+1, j, u, T1, T2);

/* check neighbour (i+1,j+1) */
if ((u[i+1][j+1] <= T2) && (u[i+1][j+1] > T1))
   trace_edge (i+1, j+1, u, T1, T2);

/* check neighbour (i,j+1) */
if ((u[i][j+1] <= T2) && (u[i][j+1] > T1))
   trace_edge (i, j+1, u, T1, T2);

/* check neighbour (i-1,j+1) */
if ((u[i-1][j+1] <= T2) && (u[i-1][j+1] > T1))
   trace_edge (i-1, j+1, u, T1, T2);

/* check neighbour (i-1,j) */
if ((u[i-1][j] <= T2) && (u[i-1][j] > T1))
   trace_edge (i-1, j, u, T1, T2);

/* check neighbour (i-1,j-1) */
if ((u[i-1][j-1] <= T2) && (u[i-1][j-1] > T1))
   trace_edge (i-1, j-1, u, T1, T2);

/* check neighbour (i,j-1) */
if ((u[i][j-1] <= T2) && (u[i][j-1] > T1))
   trace_edge (i, j-1, u, T1, T2);

/* check neighbour (i+1,j-1) */
if ((u[i+1][j-1] <= T2) && (u[i+1][j-1] > T1))
   trace_edge (i+1, j-1, u, T1, T2);

return;

} /* trace_edge */

/*--------------------------------------------------------------------------*/

void hysteresis_thresholding

     (double   **u,        /* image */
      double   T1,         /* lower threshold */
      double   T2,         /* higher threshold */
      long     nx,         /* x-dimension of u */
      long     ny)         /* y-dimension of u */

/*
  performs hysteresis thresholding and creates the binary edge image
*/

{
long  i, j;    /* loop variables */

/* use pixels that are larger than T2 as seed points to find
   additional edge pixels that are larger than T1 */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (u[i][j] >= T2)
        trace_edge (i, j, u, T1, T2);

/* create a binary image */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (u[i][j] < 255)
        u[i][j] = 0.0;

return;

}  /* hysteresis_thresholding */

/*--------------------------------------------------------------------------*/

void canny

     (double   **f,        /* input image, unchanged */
      long     nx,         /* x-dimension of u */
      long     ny,         /* y-dimension of u */
      double   hx,         /* pixel size in x-direction */
      double   hy,         /* pixel size in y-direction */
      double   sigma,      /* Gaussian standard deviation */
      double   q1,         /* lower quantile */
      double   q2,         /* higher quantile */
      double   *T1,        /* lower threshold, output */
      double   *T2,        /* higher threshold, output */
      double   **u)        /* output image with Canny edges */

/*
  Canny edge detection 
*/

{
long    i, j;    /* loop variables */
double  *v;      /* vector with all gradient magnitudes */
double  **ux;    /* derivative in x direction */
double  **uy;    /* derivative in y direction */

/* allocate memory */
alloc_double_vector (&v, nx*ny+1); 
alloc_double_matrix (&ux, nx+2, ny+2);
alloc_double_matrix (&uy, nx+2, ny+2);

/* copy image f into u */
for (i=0; i<=nx+1; i++)
 for (j=0; j<=ny+1; j++)
     u[i][j] = f[i][j];

/* perform Gaussian smoothing of u */
if (sigma > 0.0)
   gauss_conv (sigma, 0, 5.0, nx, ny, hx, hy, u);

/* compute Sobel derivatives of u */
sobel (u, nx, ny, hx, hx, ux, uy);

/* set boundaries of u to zero */
zero_boundaries (u, nx, ny);

/* compute the gradient magnitude */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = sqrt (ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j]);

/* write all gradient magnitudes in a vector v, starting with index 1 */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     v[(i-1)*ny+j] = u[i][j];

/* order it with heapsort */
heapsort (nx * ny, v);

/* find thresholds T1, T2 that belong to quantiles q1, q2 */
*T1 = v[(long) (ceil (q1 * (nx * ny)))];
*T2 = v[(long) (ceil (q2 * (nx * ny)))];

/* apply nonmaxima suppression */
nonmaxima_suppression (u, ux, uy, nx, ny);

/* apply hysteresis thresholding */
hysteresis_thresholding (u, *T1, *T2, nx, ny);

/* free memory */
free_double_vector (v, nx*ny+1); 
free_double_matrix (ux, nx+2, ny+2);
free_double_matrix (uy, nx+2, ny+2);

return;

}  /* canny */

/*--------------------------------------------------------------------------*/

int main ()

{
char    in[80];               /* for reading data */
char    out[80];              /* for reading data */
double  **f;                  /* original image */
double  **u;                  /* output image */
long    nx, ny;               /* image size in x, y direction */
double  sigma;                /* Gaussian standard deviation */
double  q1;                   /* lower quantile */
double  q2;                   /* higher quantile */
double  T1;                   /* lower threshold */
double  T2;                   /* higher threshold */
char    comments[1600];       /* string for comments */

printf ("\n");
printf ("CANNY EDGE DETECTOR FOR GREYSCALE IMAGES\n\n");
printf ("*********************************************************\n\n");
printf ("    Copyright 2022 by                                    \n");
printf ("    Markus Mainberger, Pascal Peter, Joachim Weickert    \n");
printf ("    Dept. of Mathematics and Computer Science            \n");
printf ("    Saarland University, Saarbruecken, Germany           \n\n");
printf ("    All rights reserved. Unauthorized usage,             \n");
printf ("    copying, hiring, and selling prohibited.             \n\n");
printf ("    Send bug reports to                                  \n");
printf ("    weickert@mia.uni-saarland.de                         \n\n");
printf ("*********************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

printf ("input image (pgm):                      ");
read_string (in);
read_pgm_to_double (in, &nx, &ny, &f);   /* also allocates memory for f */


/* ---- read parameters ---- */

printf ("Gaussian standard deviation sigma:      ");
read_double (&sigma);

printf ("higher quantile (<=1.00):               ");
read_double (&q2);

printf ("lower  quantile (<=%4.2lf):               ", q2);
read_double (&q1);

printf ("output image (pgm):                     ");
read_string (out);
printf ("\n");


/* ---- allocate memory ---- */

alloc_double_matrix (&u,  nx+2, ny+2);


/* ---- apply Canny edge detection ---- */

canny (f, nx, ny, 1.0, 1.0, sigma, q1, q2, &T1, &T2, u);
printf ("threshold T1: %5.2lf\n",   T1);
printf ("threshold T2: %5.2lf\n\n", T2);


/* ---- write output image (pgm format P5) ---- */

/* generate comment string */
comments[0]='\0';
comment_line (comments, "# Canny edge detector for greyscale images\n");
comment_line (comments, "# input image:   %s\n", in);
comment_line (comments, "# sigma: %7.2lf\n", sigma);
comment_line (comments, "# q1:    %7.2lf\n", q1);
comment_line (comments, "# q2:    %7.2lf\n", q2);
comment_line (comments, "# T1:    %7.2lf\n", T1);
comment_line (comments, "# T2:    %7.2lf\n", T2);

/* write image */
write_double_to_pgm (u, nx, ny, out, comments);
printf ("output image %s successfully written\n\n", out);


/* ---- free memory  ---- */

free_double_matrix (u, nx+2, ny+2);
free_double_matrix (f, nx+2, ny+2);

return(0);

}  /* main */

