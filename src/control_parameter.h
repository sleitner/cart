#ifndef __CONTROL_PARAMETER_H__
#define __CONTROL_PARAMETER_H__


#include <stdio.h>


#define CONTROL_PARAMETER_STRING_LENGTH    256

typedef void (*ControlParameterSetter)(const char *value, void *ptr, int ind);
typedef void (*ControlParameterLister)(FILE *stream, const void *ptr);

typedef struct 
{
  ControlParameterSetter setter;
  ControlParameterLister lister;
}
ControlParameterOps;


void control_parameter_set_int(const char *value, void *ptr, int ind);
void control_parameter_set_float(const char *value, void *ptr, int ind);
void control_parameter_set_double(const char *value, void *ptr, int ind);
void control_parameter_set_string(const char *value, void *ptr, int ind);

void control_parameter_list_int(FILE *stream, const void *ptr);
void control_parameter_list_float(FILE *stream, const void *ptr);
void control_parameter_list_double(FILE *stream, const void *ptr);
void control_parameter_list_string(FILE *stream, const void *ptr);


extern ControlParameterOps control_parameter_int;
extern ControlParameterOps control_parameter_float;
extern ControlParameterOps control_parameter_double;
extern ControlParameterOps control_parameter_string;


/*
//  The primary (first) name must be unique. Secondary names do not 
//  have to be unique, and all parameters with the same secondary name
//  will be set by a single read_control_parameter(...) call.
//  Parameters with primary names starting with @ are considered 
//  hidden ("for internal use only") and are not listed by default.
//  Separate function must be used to list these parameters.
//  Parameters with primary names starting with % are considered 
//  deprecated and are not listed. They are retained for backward 
//  compatibility only and their use is discouraged.
*/
void control_parameter_add(ControlParameterOps ops, void *ptr, const char *name, const char *help);
void control_parameter_add2(ControlParameterOps ops, void *ptr, const char *name1, const char *name2, const char *help);
void control_parameter_add3(ControlParameterOps ops, void *ptr, const char *name1, const char *name2, const char *name3, const char *help);
void control_parameter_add4(ControlParameterOps ops, void *ptr, const char *name1, const char *name2, const char *name3, const char *name4, const char *help);

int control_parameter_is_set(const char *name);

int control_parameter_read(const char *tag, const char *value);

void control_parameter_print(FILE *f, int with_help);
void control_parameter_print_hidden(FILE *f, int with_help);

#endif
