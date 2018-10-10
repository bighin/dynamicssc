#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>

#include "auxx.h"
#include "spline.h"

/*
	Spline interpolation, given a grid of points, essentially a nice wrapper
	over the code contained in Numerical Recipes (NR).

	Note that the code in NR is __not__ thread safe, the 'static' keyword has to
	be removed in the splint() function.
*/

struct interpolation_t *init_interpolation(double *x,double *y,int n)
{
	struct interpolation_t *ret;
	
	if(!(ret=malloc(sizeof(struct interpolation_t))))
		return NULL;

	if(!(ret->x=malloc(sizeof(double)*n)))
	{
		if(ret)
			free(ret);
	
		return NULL;
	}

	if(!(ret->y=malloc(sizeof(double)*n)))
	{
		if(ret)
		{
			if(ret->x)
				free(ret->x);

			free(ret);
		}
		
		return NULL;
	}

	if(!(ret->y2=malloc(sizeof(double)*n)))
	{
		if(ret)
		{
			if(ret->x)
				free(ret->x);

			if(ret->y)
				free(ret->y);
		
			free(ret);
		}

		return NULL;
	}

	memcpy(ret->x,x,sizeof(double)*n);
	memcpy(ret->y,y,sizeof(double)*n);
	ret->n=n;
	
	spline(ret->x,ret->y,ret->n,0.0f,0.0f,ret->y2);

	return ret;
}

void fini_interpolation(struct interpolation_t *it)
{
	if(it)
	{
		if(it->x)
			free(it->x);

		if(it->y)
			free(it->y);

		if(it->y2)
			free(it->y2);
	
		free(it);
	}
}

void copy_interpolation(struct interpolation_t *dst,struct interpolation_t *src)
{
	assert(dst);
	assert(src);

	dst->n=src->n;

	memcpy(dst->x,src->x,sizeof(double)*dst->n);
	memcpy(dst->y,src->y,sizeof(double)*dst->n);
	memcpy(dst->y2,src->y2,sizeof(double)*dst->n);
}

double get_point(struct interpolation_t *it,double x)
{
	double y;
	
	splint(it->x,it->y,it->y2,it->n,x,&y);
	
	return y;
}

void fcopy(FILE *f1,FILE *f2)
{
	char buffer[1024];
	size_t n;

	while((n=fread(buffer,sizeof(char),1024,f1))>0)
	{
		if (fwrite(buffer,sizeof(char),n,f2)!=n)
			perror("Write failed!\n");
	}
}

FILE *fopen_mkdir(const char *name,const char *mode)
{
	char *mname=strdup(name);

	for(int i=0;mname[i]!='\0';i++)
	{
		if((i>0)&&((mname[i]=='\\')||(mname[i]=='/')))
		{
			char slash=mname[i];

			mname[i]='\0';

			if(!mkdir(mname,S_IRWXU|S_IRWXG|S_IRWXO))
			{
				free(mname);
				return NULL;
			}
			
			mname[i]=slash;
		}
	}

	free(mname);
	return fopen(name,mode);
}

double arccot(double x)
{
	return (M_PI/2.0f)-atan(x);
}

/*
 * From: https://stackoverflow.com/questions/4833347/removing-substring-from-a-string
 * 
 * Description:
 *   Find and replace text within a string.
 *
 * Parameters:
 *   src  (in) - pointer to source string
 *   from (in) - pointer to search text
 *   to   (in) - pointer to replacement text
 *
 * Returns:
 *   Returns a pointer to dynamically-allocated memory containing string
 *   with occurences of the text pointed to by 'from' replaced by with the
 *   text pointed to by 'to'.
 */

char *find_and_replace(const char *src, const char *from, const char *to)
{
   /*
    * Find out the lengths of the source string, text to replace, and
    * the replacement text.
    */
   size_t size    = strlen(src) + 1;
   size_t fromlen = strlen(from);
   size_t tolen   = strlen(to);
   /*
    * Allocate the first chunk with enough for the original string.
    */
   char *value = malloc(size);
   /*
    * We need to return 'value', so let's make a copy to mess around with.
    */
   char *dst = value;
   /*
    * Before we begin, let's see if malloc was successful.
    */
   if ( value != NULL )
   {
      /*
       * Loop until no matches are found.
       */
      for ( ;; )
      {
         /*
          * Try to find the search text.
          */
         const char *match = strstr(src, from);
         if ( match != NULL )
         {
            /*
             * Found search text at location 'match'. :)
             * Find out how many characters to copy up to the 'match'.
             */
            size_t count = match - src;
            /*
             * We are going to realloc, and for that we will need a
             * temporary pointer for safe usage.
             */
            char *temp;
            /*
             * Calculate the total size the string will be after the
             * replacement is performed.
             */
            size += tolen - fromlen;
            /*
             * Attempt to realloc memory for the new size.
             */
            temp = realloc(value, size);
            if ( temp == NULL )
            {
               /*
                * Attempt to realloc failed. Free the previously malloc'd
                * memory and return with our tail between our legs. :(
                */
               free(value);
               return NULL;
            }
            /*
             * The call to realloc was successful. :) But we'll want to
             * return 'value' eventually, so let's point it to the memory
             * that we are now working with. And let's not forget to point
             * to the right location in the destination as well.
             */
            dst = temp + (dst - value);
            value = temp;
            /*
             * Copy from the source to the point where we matched. Then
             * move the source pointer ahead by the amount we copied. And
             * move the destination pointer ahead by the same amount.
             */
            memmove(dst, src, count);
            src += count;
            dst += count;
            /*
             * Now copy in the replacement text 'to' at the position of
             * the match. Adjust the source pointer by the text we replaced.
             * Adjust the destination pointer by the amount of replacement
             * text.
             */
            memmove(dst, to, tolen);
            src += fromlen;
            dst += tolen;
         }
         else /* No match found. */
         {
            /*
             * Copy any remaining part of the string. This includes the null
             * termination character.
             */
            strcpy(dst, src);
            break;
         }
      }
   }
   return value;
}

