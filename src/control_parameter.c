
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"


ControlParameterOps control_parameter_int = { control_parameter_set_int, control_parameter_list_int };
ControlParameterOps control_parameter_float = { control_parameter_set_float, control_parameter_list_float };
ControlParameterOps control_parameter_double = { control_parameter_set_double, control_parameter_list_double };
ControlParameterOps control_parameter_string = { control_parameter_set_string, control_parameter_list_string };


#define MAX_NAMES   10


typedef struct
{
  int set;
  void *ptr;
  const char *help;
  const char *names[MAX_NAMES];
  ControlParameterSetter setter;
  ControlParameterLister lister;
} ControlParameter;

ControlParameter *control_parameters = 0;
int control_parameters_size = 0;
int num_control_parameters = 0;


void control_parameter_add_worker(ControlParameterSetter setter, ControlParameterLister lister, void *ptr, int num_names, const char **names, const char *help)
{
  int i, j;
  ControlParameter *tmp;

  cart_assert(num_names>0 && names!=NULL && names[0]!=NULL && ptr!=NULL && setter!=NULL && lister!=NULL);

  if(num_control_parameters == control_parameters_size)
    {
      control_parameters_size += 100;
      tmp = (ControlParameter *)malloc(control_parameters_size*sizeof(ControlParameter));
      if(control_parameters != 0)
	{
	  memcpy(tmp,control_parameters,num_control_parameters*sizeof(ControlParameter));
	  free(control_parameters);
	}
      control_parameters = tmp;
    }

  /*
  //  Check that pointers are unique.
  */
  for(i=0; i<num_control_parameters; i++)
    {
      if(control_parameters[i].ptr == ptr)
	{
	  cart_error("Pointer to parameter %s is not unique (%d,%d)",names[0],i,num_control_parameters);
	}
    }

  /*
  //  Check that primary names are unique.
  */
  for(i=0; i<num_control_parameters; i++)
    {
      if(strcmp(control_parameters[i].names[0],names[0]) == 0)
	{
	  cart_error("Primary parameter name %s is not unique (%d,%d)",names[0],i,num_control_parameters);
	}
    }

  control_parameters[num_control_parameters].set = 0;
  control_parameters[num_control_parameters].ptr = ptr;
  control_parameters[num_control_parameters].help = help;
  control_parameters[num_control_parameters].setter = setter;
  control_parameters[num_control_parameters].lister = lister;

  for(j=0; j<num_names; j++)
    {
      control_parameters[num_control_parameters].names[j] = names[j];
    }
  for(j=num_names; j<MAX_NAMES; j++)
    {
      control_parameters[num_control_parameters].names[j] = NULL;
    }

  num_control_parameters++;
}


int control_parameter_is_set(const char *name)
{
  int i, j;

  for(i=0; i<num_control_parameters; i++)
    {
      for(j=0; j<MAX_NAMES && control_parameters[i].names[j]!=NULL; j++)
	{
	  if(strcmp(name,control_parameters[i].names[j]) == 0)
	    {
	      return control_parameters[i].set;
	    }
	}
    }

  cart_error("String '%s' is not a valid name for a control parameter.",name);
  return 0;
}


void control_parameter_add(ControlParameterOps ops, void *ptr, const char *name, const char *help)
{
  control_parameter_add_worker(ops.setter,ops.lister,ptr,1,&name,help);
}


void control_parameter_add2(ControlParameterOps ops, void *ptr, const char *name1, const char *name2, const char *help)
{
  const char *names[2] = { name1, name2 };
  control_parameter_add_worker(ops.setter,ops.lister,ptr,2,names,help);
}


void control_parameter_add3(ControlParameterOps ops, void *ptr, const char *name1, const char *name2, const char *name3, const char *help)
{
  const char *names[3] = { name1, name2, name3 };
  control_parameter_add_worker(ops.setter,ops.lister,ptr,3,names,help);
}


void control_parameter_print_name(FILE *f, const char *name)
{
  const int tab = 40;
  int j, n, offset;

  offset = 0;
  if(name[0] == '@') offset++;
  if(name[0] == '%') offset++;

  fprintf(f,"%s ",name+offset);
  n = tab - strlen(name);
  for(j=0; j<n; j++) fprintf(f," ");
}


void control_parameters_print_worker(FILE *f, int with_help, int hidden)
{
  const int wrap = 70;
  int i, j, k, lastws;
  char str[wrap+1];
  
  for(i=0; i<num_control_parameters; i++) if(control_parameters[i].names[0][0]!='%' && ((hidden && control_parameters[i].names[0][0]=='@') || (!hidden && control_parameters[i].names[0][0]!='@')))
    {
      control_parameter_print_name(f,control_parameters[i].names[0]);
      control_parameters[i].lister(f,control_parameters[i].ptr);
      fprintf(f,"\n");

      if(with_help && control_parameters[i].help!=NULL)
	{
	  j = 0;
	  while(j < strlen(control_parameters[i].help))
	    {
	      strncpy(str,control_parameters[i].help+j,wrap);
	      /*
	      //  Find the break point
	      */
	      lastws = -1;
	      for(k=0; k<wrap; k++)
		{
		  if(str[k]=='\n' || str[k]==0) break;
		  if(str[k] == ' ') lastws = k;
		}

	      /* Break at whitespace */
	      if(k==wrap && lastws>-1) k = lastws;

	      str[k] = 0;
	      fprintf(f,"\t# %s\n",str);
	      j += ++k;
	    }
	}
    }
}


void control_parameter_print(FILE *f, int with_help)
{
  control_parameters_print_worker(f,with_help,0);
}


void control_parameter_print_hidden(FILE *f, int with_help)
{
  control_parameters_print_worker(f,with_help,1);
}


void control_parameter_set_int(const char *value, void *ptr, int ind)
{
  if(sscanf(value,"%d",(int *)ptr) != 1) cart_error("Unable to read INT parameter from string '%s'",value);
}


void control_parameter_set_float(const char *value, void *ptr, int ind)
{
  if(sscanf(value,"%g",(float *)ptr) != 1) cart_error("Unable to read FLOAT parameter from string '%s'",value);
}


void control_parameter_set_double(const char *value, void *ptr, int ind)
{
  if(sscanf(value,"%lg",(double *)ptr) != 1) cart_error("Unable to read DOUBLE parameter from string '%s'",value);
}


void control_parameter_set_string(const char *value, void *ptr, int ind)
{
  strncpy((char *)ptr,value,CONTROL_PARAMETER_STRING_LENGTH);
  ((char *)ptr)[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
}


void control_parameter_list_int(FILE *stream, const void *ptr)
{
  fprintf(stream,"%-d",*(int *)ptr);
}


void control_parameter_list_float(FILE *stream, const void *ptr)
{
  fprintf(stream,"%-g",*(float *)ptr);
}


void control_parameter_list_double(FILE *stream, const void *ptr)
{
  fprintf(stream,"%-lg",*(double *)ptr);
}


void control_parameter_list_string(FILE *stream, const void *ptr)
{
  fprintf(stream,"%s",(char *)ptr);
}


int control_parameter_read(const char *tag, const char *value)
{
  int i, j, offset;
  /*
  //  unknown tag 
  */
  int ret = 1;
  
  for(i=0; i<num_control_parameters; i++)
    {
      for(j=0; j<MAX_NAMES && control_parameters[i].names[j]!=NULL; j++)
	{
	  /*
	  //  Skip control char
	  */
	  offset = 0;
	  if(j == 0)
	    {
	      if(control_parameters[i].names[j][0] == '@') offset++;
	      if(control_parameters[i].names[j][0] == '%') offset++;
	    }
	  if(strcmp(tag,control_parameters[i].names[j]+offset) == 0)
	    {
	      /*
	      //  Issue warning if setting hidden or deprecated parameter
	      */
	      if(offset > 0)
		{
		  if(control_parameters[i].names[j][0] == '@')
		    {
		      cart_debug("Parameter <%s> is for internal use only. Changing this parameter value can break the code and is not recommended.",control_parameters[i].names[j]+offset);
		    }
		  if(control_parameters[i].names[j][0] == '%')
		    {
		      cart_debug("Parameter <%s> is deprecated. Please avoid using this parameter in the config file in the future.",control_parameters[i].names[j]+offset);
		    }
		}
	      control_parameters[i].set = 1;
	      control_parameters[i].setter(value,control_parameters[i].ptr,j);
	      /*
	      //  Permit setting more than one parameter with common 
	      //  secondary names (for backward compatibility).
	      */
	      if(j == 0) return 0; else ret = 0;
	    }
	}
    }

  return ret;
}

