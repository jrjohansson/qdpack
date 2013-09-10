//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

#define _GNU_SOURCE

#include <stdio.h>

#include <qdpack/qdpack.h>

/*
 * Print complex matrix with element format "(real, imag)"
 *
 */
void
qdpack_operator_print(qdpack_operator_t *m)
{
    double val_re, val_im;
    int i, j;
    
    for (i = 0; i < m->m; i++) 
    {
           for (j = 0; j < m->n; j++)
        {    
            val_re = QDPACK_REAL(qdpack_operator_get(m, i, j));
            val_im = QDPACK_IMAG(qdpack_operator_get(m, i, j));
            printf (" (%2.5f,%2.5f)", val_re, val_im);
        }
        printf("\n");
    }
}
/*
 *
 *
 */
void
qdpack_operator_print_real(qdpack_operator_t *m)
{
    double val;
    int i, j;
    
    for (i = 0; i < m->m; i++) 
    {
           for (j = 0; j < m->n; j++)
        {    
            val = QDPACK_REAL(qdpack_operator_get(m, i, j));
            if (val >= 0.0)
                printf (" %2.5f  ", val);
            else
                printf ("%2.5f  ", val);
        }
        printf("\n");
    }
}

/*
 *
 *
 */
void
qdpack_operator_print_imag(qdpack_operator_t *m)
{
    int i, j;
    
    for (i = 0; i < m->m; i++) 
    {
           for (j = 0; j < m->n; j++)
            printf ("%2.5f  ", QDPACK_IMAG(qdpack_operator_get(m, i, j)));
        printf("\n");
    }
}

/*
 *
 *
 */
void
qdpack_operator_print_abs2(qdpack_operator_t *m)
{
    qdpack_complex z;
    int i, j;
    
    for (i = 0; i < m->m; i++) 
    {
           for (j = 0; j < m->n; j++)
        {
            z = qdpack_operator_get(m, i, j);
            printf("%2.5f  ", QDPACK_REAL(z)*QDPACK_REAL(z) + QDPACK_IMAG(z)*QDPACK_IMAG(z));
        }
        printf("\n");
    }
}


/*
 *
 *
 *
 */
int
qdpack_file_row_add(const char *filename, const char *row)
{
    FILE *file;

//    while ((file = fopen(filename, "ax")) == NULL)
//    {
//        if (errno == EEXIST)
//            continue;
//        
//        fprintf(stderr, "qdpack_file_row_add: Failed to open file %s\n", filename);
//        return -1;
//    }


    if ((file = fopen(filename, "a")) == NULL)
    {
        fprintf(stderr, "qdpack_file_row_add: Failed to open file %s\n", filename);
        return -1;
    }

    fprintf(file, "%s\n", row);
    
    fclose(file);


/*
    {
        int fd;

        if ((fd = open(file, O_APPEND)) == -1)
        {
            fprintf(stderr, "qdpack_file_row_add: Failed to open file %s\n", filename);
                    return -1;
        }

        while (flock(fd, LOCK_EX) == -1)
        {
            ;
        }

        write(fd, 

        flock(fd, LOCK_UN);
        close(fd);
    }
*/    

    return 0;
}

/*
 * Write/read matrix to/from file:
 *
 * Format
 *
 * a11 a12 a13 ... a1N 
 * a21 a22 a23 ... a2N 
 * ... 
 * aN1 aN2 aN3 ... aNN
 *
 * Arguments:
 *
 *  m         = matrix to write to
 *  filename  = name of file to read from
 *  n         = dimension of (square) matrix
 *
 */
int
qdpack_operator_write(char *filename, qdpack_operator_t *m)
{
    int i, j;
    FILE *file;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_write: Failed to open file %s\n", filename);
        return -1;
    }

    for (i = 0; i < m->m; i++)
    {
        for (j = 0; j < m->n; j++)
        {
            qdpack_complex z;
            z = qdpack_operator_get(m, i, j);
            fprintf(file, "(%f, %f) ", QDPACK_REAL(z), QDPACK_IMAG(z));
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return 0;
}

int
qdpack_operator_write_real(char *filename, qdpack_operator_t *m)
{
    int i, j;
    FILE *file;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_write: Failed to open file %s\n", filename);
        return -1;
    }

    for (i = 0; i < m->m; i++)
    {
        for (j = 0; j < m->m; j++)
        {
            qdpack_complex z;
            z = qdpack_operator_get(m, i, j);
            fprintf(file, "%f ", QDPACK_REAL(z));
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return 0;
}

int
qdpack_operator_write_imag(char *filename, qdpack_operator_t *m)
{
    int i, j;
    FILE *file;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_write: Failed to open file %s\n", filename);
        return -1;
    }

    for (i = 0; i < m->m; i++)
    {
        for (j = 0; j < m->n; j++)
        {
            qdpack_complex z;
            z = qdpack_operator_get(m, i, j);
            fprintf(file, "%f ", QDPACK_IMAG(z));
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return 0;
}

/*
 * Write/read matrix to/from file (either real or imaginary part):
 *
 * Format
 *
 * a11 a12 a13 ... a1N 
 * a21 a22 a23 ... a2N 
 * ... 
 * aN1 aN2 aN3 ... aNN
 *
 * Arguments:
 *
 *  m         = matrix to write to
 *  filename  = name of file to read from
 *  n         = dimension of (square) matrix
 *
 */
int
qdpack_operator_write_double(char *filename, qdpack_operator_t *m, int real)
{
    int i, j;
    FILE *file;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_write: Failed to open file %s\n", filename);
        return -1;
    }

    for (i = 0; i < m->m; i++)
    {
        for (j = 0; j < m->m; j++)
        {
            qdpack_complex z;
            z = qdpack_operator_get(m, i, j);
            if (real)
            {
                fprintf(file, "%f\t", QDPACK_REAL(z));
            }
            else
            {
                fprintf(file, "%f\t", QDPACK_IMAG(z));
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return 0;
}

int
qdpack_operator_read(char *filename, qdpack_operator_t *m)
{
    int i = 0, j;
    char *linebuff = NULL; //[100000];
    FILE *file;
    size_t size;

    if ((file = fopen(filename, "r")) == NULL)
    {
            fprintf(stderr, "qdpack_operator_read: Failed to open file %s\n", filename);
            return -1;
    }

    while (getline(&linebuff, &size, file) > 0)
    {
        qdpack_complex z;
        double zre, zim;
        char *endptr, *ptr;
        j = 0;
        

        ptr = &linebuff[1];

        while (j < m->n)
        {
            zre = strtod(ptr, &endptr);
            ptr = endptr + 1;
            zim = strtod(ptr, &endptr);
            ptr = endptr + 3;

            QDPACK_SET_COMPLEX(&z, zre, zim);
            
            qdpack_operator_set(m, i, j, z);

            j++;
        }
        i++;

        //free(linebuff);
    }

    return 0;    
}



/*
 * Write/read diagonal matrix to/from file:
 *
 * Format
 *
 * a11 a22 a33 ... aNN 
 *
 * Arguments:
 *
 *  m         = matrix to write to
 *  filename  = name of file to read from
 *  n         = dimension of (square) matrix
 *
 */
int
qdpack_operator_diag_write(char *filename, qdpack_operator_t *m, char *msg)
{
    int i;
    FILE *file;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_read_write: Failed to open file %s\n", filename);
        return -1;
    }

    fprintf(file, "%s", msg);
    for (i = 0; i < m->n; i++)
    {
        qdpack_complex z;
        z = qdpack_operator_get(m, i, i);
        fprintf(file, "(%f, %f) ", QDPACK_REAL(z), QDPACK_IMAG(z));
    }
    fprintf(file, "\n");

    fclose(file);
    return 0;
}

int
qdpack_operator_diag_read(char *filename, qdpack_operator_t *m)
{
    int j;
    char *linebuff = NULL;
    FILE *file;
    size_t size;

    if ((file = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_diag_read: Failed to open file %s\n", filename);
        return -1;
    }

    qdpack_operator_set_zero(m);

    if (getline(&linebuff, &size, file) > 0)
    {
        qdpack_complex z;
        double zre, zim;
        char *endptr, *ptr;
        j = 0;

        ptr = &linebuff[1];

        while (j < m->n)
        {
            zre = strtod(ptr, &endptr);
            ptr = endptr + 1;
            zim = strtod(ptr, &endptr);
            ptr = endptr + 3;

            QDPACK_SET_COMPLEX(&z, zre, zim);
            
            qdpack_operator_set(m, j, j, z);

            j++;
        }
        //free(linebuff);
    }

    return 0;    
}

int
qdpack_operator_diag_real_write(char *filename, qdpack_operator_t *m, char *msg)
{
    int i;
    FILE *file;

    if ((file = fopen(filename, "w")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_diag_real_write: Failed to open file %s\n", filename);
        return -1;
    }

    fprintf(file, "%s", msg);
    for (i = 0; i < m->n; i++)
    {
        qdpack_complex z;
        z = qdpack_operator_get(m, i, i);
        fprintf(file, "%f ", QDPACK_REAL(z));
    }
    fprintf(file, "\n");

    fclose(file);
    return 0;
}


int
qdpack_operator_diag_real_read(char *filename, qdpack_operator_t *m)
{
    int j;
    char *linebuff = NULL;
    FILE *file;
    size_t size;

    if ((file = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "qdpack_operator_diag_real_read: Failed to open file %s\n", filename);
        return -1;
    }

    qdpack_operator_set_zero(m);

    if (getline(&linebuff, &size, file) > 0)
    {
        qdpack_complex z;
        double zre;
        char *endptr, *ptr;
        j = 0;

        ptr = linebuff;

        while (j < m->n)
        {
            zre = strtod(ptr, &endptr);
            ptr = endptr;

            QDPACK_SET_COMPLEX(&z, zre, 0.0);
            
            qdpack_operator_set(m, j, j, z);

            j++;
        }
        //free(linebuff);
    }

    return 0;    
}



// -----------------------------------------------------------------------------
// Vector functions
//
//
void
qdpack_state_print_real(qdpack_state_t *v)
{
    int i;
    
    for (i = 0; i < v->n; i++) 
    {
        printf("%2.5f ", QDPACK_REAL(qdpack_state_get(v, i)));
    }
    printf("\n");
}

void
qdpack_state_print_imag(qdpack_state_t *v)
{
    int i;
    
    for (i = 0; i < v->n; i++) 
    {
        printf("%2.5f ", QDPACK_IMAG(qdpack_state_get(v, i)));
    }
    printf("\n");
}






