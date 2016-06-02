//
// YOCTO_CMDLINE: utilities for writing command line applications, mostly
// a simple to use command line parser.
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. create a parser object
//     - an option for printing help is automatically added
//     parser = yc_init_parser(argc, argv, program description)
// 2. for each option, parse it calling the functions yc_parse_optXXX
//     - options are parsed on the flh and a comprehensive help is
//       automatically generated
//     - supports bool (flags), int, float, double, string, enums
//     - options names are "--longname" for longname and "-s" for short
//     - command line format is "--longname value", "-s v" for all but flags
//     opt = yc_parse_opti(parser, longname, flag, description, default value)
// 3. for each unnamed argument, parse it calling the functions yc_parse_argXXX
//     - names are only used for help
//     - supports above type and all a special one that comsumes all args
//     saving them in a string array
//     arg = yc_parse_args(parser, longname, description, deafult value,
//               required)
// 4. end cmdline parsing with yc_done_parser
//    - check for unsued values, missing arguments and print help if needed
//    yc_done_parser(parser)
// 5. since arguments are parsed immediately, one can easily implement
//    subcommands by just branching the command line code based on a read
//    argument without any need for complex syntax
// 6. contains also utilities to load whole files in memory and split paths
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Notes: the end of this file contains a test function that also
// illustreates the library usage.
//

//
// COMPILATION:
//
// All functions in this library are inlined by default for ease of use in
// C/C++.
// To use the library as a .h/.cpp pair do the following:
// - to use as a .h, just #define YGL_DECLARATION before including this file
// - to build as a .cpp, just #define YGL_IMPLEMENTATION before including this
// file into only one file that you can either link directly or pack as a lib.
//
// This file depends on yocto_math.h.
//

//
// HISTORY:
// - v 0.1: C++ implementation
// - v 0.0: initial release in C99
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#ifndef _YC_H_
#define _YC_H_

// compilation options
#ifdef __cplusplus
#ifndef YGL_DECLARATION
#define YGL_API inline
#define YGLC_API inline
#else
#define YGL_API
#define YGLC_API extern "C"
#endif
#endif

#ifndef __cplusplus
#define YGLC_API extern
#include <stdbool.h>
#endif

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

// COMMAND LINE PARSING --------------------------------------------------------

//
// Forward declaration of an opaque parser object.
//
typedef struct yc_parser yc_parser;

//
// Inits the parser with command line arguments and a usage string.
//
YGLC_API yc_parser* yc_init_parser(int argc, const char** argv,
                                   const char* help);

//
// Finish parsing, checking for errors and rinting help if needed.
//
YGLC_API bool yc_done_parser(yc_parser* parser);

//
// Parses an optional flag as in the intro.
//
YGLC_API bool yc_parse_optb(yc_parser* parser, const char* longname,
                            const char* shortname, const char* help, bool def);

//
// Parses an optional int as described in the intro.
//
YGLC_API int yc_parse_opti(yc_parser* parser, const char* longname,
                           const char* shortnamee, const char* help, int def);

//
// Parses an optional float as described in the intro.
//
YGLC_API float yc_parse_optf(yc_parser* parser, const char* longname,
                             const char* shortname, const char* help,
                             float def);

//
// Parses an optional double as described in the intro.
//
YGLC_API double yc_parse_optd(yc_parser* parser, const char* longname,
                              const char* shortname, const char* help,
                              double def);

//
// Parses an optional string with the argument described in the intro.
// No additional memory is allocated.
//
YGLC_API const char* yc_parse_opts(yc_parser* parser, const char* longname,
                                   const char* shortname, const char* help,
                                   const char* def);

//
// Parses an optional enum defined by the NULL-terminated array of strings
// choice and returns the index of the string in that array.
//
YGLC_API int yc_parse_opte(yc_parser* parser, const char* longname,
                           const char* shortnamee, const char* help, int def,
                           const char** choices);

//
// Parses an unnamed int argument as described in the intro.
//
YGLC_API int yc_parse_argi(yc_parser* parser, const char* longname,
                           const char* help, int def, bool required);

//
// Parses an unnamed string argument as described in the intro.
//
YGLC_API const char* yc_parse_args(yc_parser* parser, const char* longname,
                                   const char* help, const char* def,
                                   bool required);

//
// Parses an unnamed array of string argument as described in the intro.
//
YGLC_API const char** yc_parse_argsa(yc_parser* parser, const char* longname,
                                     const char* help, const char** def,
                                     bool required, const char* values[],
                                     int max_values);

// FILE LOADING ----------------------------------------------------------------

//
// Loads the contents of a binary file in an in-memory array.
//
YGLC_API void* yc_load_binfile(const char* filename, int* size);

//
// Loads the contents of a text file into a NULL-terminated string.
//
YGLC_API char* yc_load_txtfile(const char* filename);

// PATH MANIPULATION -----------------------------------------------------------

//
// Splits paths into directory name (including '/'), basename, and extension
// (including '.'). dirname, basename and ext have to preallocated and have
// enough memory to contain the data (e.g. strlen(path)+1 would do).
// If any the pointers is missing, that component will be skipped.
//
YGLC_API void yc_split_path(const char* filename, char* dirname, char* basename,
                            char* ext);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if !defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION)

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//
// limits for parser data
//
#define yc__maxargs 4096  // max cmdline arguments
#define yc__maxopts 64    // max options to parse
#define yc__maxline 256   // max help line length

//
// type of parser data
//
enum {
    yc__bool = 1,
    yc__int = 2,
    yc__float = 3,
    yc__double = 4,
    yc__string = 5,
    yc__string_array = 6
};

//
// Command line parser
//
struct yc_parser {
    // saved arguments ----------------------------------------
    int argc;                       // saved numer of arguments
    const char* argv[yc__maxargs];  // saved arguments

    // help strings -------------------------------------------
    char help_prog[yc__maxline];                // program description
    char help_usage[yc__maxline];               // sage lines
    char help_opts[yc__maxopts * yc__maxline];  // options lines
    char help_args[yc__maxopts * yc__maxline];  // unnamed arguments lines

    // parse data ---------------------------------------------
    bool error;            // whether a parsing error occurred
    char error_msg[4096];  // error message
    bool exit_on_error;    // exit on error
    bool print_help;       // print help

    // helper values ------------------------------------------
    char* stras_values[yc__maxargs];  // helper values
};

//
// Inits the parser.
//
YGLC_API yc_parser* yc_init_parser(int argc, const char** argv,
                                   const char* help) {
    // Allocates parser and copy argument data
    yc_parser* parser = new yc_parser();
    assert(argc < yc__maxargs);
    parser->argc = argc - 1;
    for (int i = 0; i < parser->argc; i++) parser->argv[i] = argv[i + 1];
    parser->error = false;
    parser->exit_on_error = true;
    strcpy(parser->help_prog, argv[0]);
    strcpy(parser->help_usage, help);
    parser->help_opts[0] = 0;
    parser->help_args[0] = 0;
    parser->print_help =
        yc_parse_optb(parser, "--help", "-h", "print help", false);
    return parser;
}

//
// Print help based on the help lines collected during parsing.
//
static inline void yc__print_help(yc_parser* parser) {
    printf("usage: %s [options] %s\n", parser->help_prog,
           strlen(parser->help_args) ? "<arguments>" : "");
    if (strlen(parser->help_usage)) printf("  %s\n", parser->help_usage);
    printf("\n");
    if (strlen(parser->help_opts)) printf("options:\n%s\n", parser->help_opts);
    if (strlen(parser->help_args))
        printf("arguments:\n%s\n", parser->help_args);
}

//
// Ends parsing checking for error for unused options or arguments.
// Exit if needed.
//
YGLC_API bool yc_done_parser(yc_parser* parser) {
    // check for error
    if (!parser->error && parser->argc > 0) {
        parser->error = true;
        if (parser->argv[0][0] == '-')
            sprintf(parser->error_msg, "unknown option %s", parser->argv[0]);
        else
            strcpy(parser->error_msg, "unsued values");
    };
    // check whether we need to print help and exit
    if (parser->error) printf("error: %s\n", parser->error_msg);
    if (parser->print_help || parser->error) {
        yc__print_help(parser);
        if (parser->exit_on_error) exit(EXIT_FAILURE);
    }
    return !parser->error;
}

//
// Check if an option name is valid
//
static inline void yc__check_name(yc_parser* parser, const char* longname,
                                  const char* shortname, int type) {
    // check name
    assert(longname);
    if (shortname)
        assert(strlen(shortname) > 1 && shortname[0] == '-' &&
               shortname[1] != '-');
    assert(strlen(longname) > 1);
    bool flag = longname[0] == '-';
    if (flag)
        assert(strlen(longname) > 2 && longname[0] == '-' &&
               longname[1] == '-' && longname[2] != '-');
    else
        assert(type != yc__bool);
}

//
// Add a formatted help line
//
static inline void yc__add_help(yc_parser* parser, const char* longname,
                                const char* shortname, const char* help,
                                int type, const void* def,
                                const void* choices) {
    // get whether it is an opt
    bool opt = longname[0] == '-';

    // initialize help str
    char val_str[256] = {0};
    val_str[0] = 0;
    const char* type_str = 0;
    char choice_str[256 * 10] = {0};

    // set strings
    switch (type) {
        case yc__bool:
            if (def) sprintf(val_str, "[%s]", (*(bool*)def) ? "true" : "false");
            type_str = "";
            break;
        case yc__int:
            if (def) sprintf(val_str, "[%d]", *(int*)def);
            type_str = "<int>";
            break;
        case yc__float:
            if (def) sprintf(val_str, "[%f]", *(float*)def);
            type_str = "<num>";
            break;
        case yc__double:
            if (def) sprintf(val_str, "[%lf]", *(double*)def);
            type_str = "<num>";
            break;
        case yc__string:
            if (def) sprintf(val_str, "[%s]", *(const char**)def);
            type_str = "<str>";
            if (choices) {
                const char** choice = (const char**)choices;
                strcat(choice_str, "(");
                while (*choice) {
                    if (choice != choices) strcat(choice_str, ", ");
                    strcat(choice_str, *choice);
                    choice++;
                }
                strcat(choice_str, ")");
            }
            break;
        case yc__string_array:
            if (def) {
                if (!((char**)def)[0])
                    sprintf(val_str, "[<empty>]");
                else if (!((char**)def)[1])
                    sprintf(val_str, "[%s]", ((char**)def)[0]);
                else
                    sprintf(val_str, "[%s,...]", ((char**)def)[0]);
                type_str = "<str>*";
            } else {
                type_str = "<str>+";
            }
            break;
        default: assert(false); break;
    }

    // print full name
    char help_line[yc__maxline], help_fullname[yc__maxline];
    if (shortname)
        sprintf(help_fullname, "%s/%s %s", longname, shortname, type_str);
    else
        sprintf(help_fullname, "%s %s", longname, type_str);

    // print help line
    sprintf(help_line, "  %-24s  %s %s\n", help_fullname, help, val_str);
    if (strlen(choice_str)) {
        char buf[4096];
        sprintf(buf, "  %-24s  %s\n", "", choice_str);
        strcat(help_line, buf);
    }

    // add line to proper help
    strcat((opt) ? parser->help_opts : parser->help_args, help_line);
}

//
// Set default values for each argument type
//
static inline void yc__set_default(yc_parser* parser, int type, const void* def,
                                   void* value) {
    // if not present, set default and exit
    switch (type) {
        case yc__bool: *(bool*)value = *(bool*)def; break;
        case yc__int: *(int*)value = *(int*)def; break;
        case yc__float: *(float*)value = *(float*)def; break;
        case yc__double: *(double*)value = *(double*)def; break;
        case yc__string: *(const char**)value = *(const char**)def; break;
        case yc__string_array:
            for (int i = 0;; i++) {
                ((char**)value)[i] = ((char**)def)[i];
                if (((char**)def)[i]) break;
            }
            break;
        default: assert(false);
    }
}

//
// Main parsing routine working in generic data.
//
static inline void yc__parse_val(yc_parser* parser, const char* longname,
                                 const char* shortname, const char* help,
                                 int type, const void* def, void* value,
                                 void* choices) {
    // check whether the name is good
    yc__check_name(parser, longname, shortname, type);

    // add help line
    yc__add_help(parser, longname, shortname, help, type, def, choices);

    // skip if alreasy in error
    if (parser->error) return;

    // set value to default if present to handle errors
    if (def) yc__set_default(parser, type, def, value);

    // check whether it is an option
    bool opt = longname[0] == '-';

    // grab opt or arg value
    int val_pos = -1;
    if (opt) {
        // find the option
        for (int i = 0; i < parser->argc && val_pos == -1; i++) {
            const char* sn = (shortname) ? shortname : "";
            if (!strcmp(sn, parser->argv[i]) ||
                !strcmp(longname, parser->argv[i]))
                val_pos = i;
        }

        // check if argument is present
        if (val_pos > parser->argc - 1) {
            parser->error = true;
            sprintf(parser->error_msg, "missing value for option %s", longname);
            return;
        }

        // remove opt name from the list if found
        if (val_pos >= 0) {
            for (int i = val_pos; i < parser->argc - 1; i++)
                parser->argv[i] = parser->argv[i + 1];
            parser->argc -= 1;
        }
    } else {
        // set arg position to first param if still there
        if (parser->argc) {
            // check tht this is not an option
            if (parser->argv[0][0] == '-') {
                parser->error = true;
                sprintf(parser->error_msg, "unknown option %s",
                        parser->argv[0]);
                return;
            } else
                val_pos = 0;
        }
    }

    // if not present, error if not default, otherwise exit
    if (val_pos < 0) {
        if (!def) {
            parser->error = true;
            sprintf(parser->error_msg, "missing value for %s", longname);
        }
        return;
    }

    // parse boolean case
    if (type == yc__bool) {
        *(bool*)value = (def) ? !*(bool*)def : true;
    }

    // parse array of strings
    if (type == yc__string_array) {
        int nstr = parser->argc - val_pos;
        for (int i = 0; i < nstr; i++)
            ((const char**)value)[i] = parser->argv[i + val_pos];
        ((char**)value)[nstr] = 0;
        parser->argc -= nstr;
    }

    // parse all standard values
    if (type != yc__bool && type != yc__string_array) {
        // parse value
        const char* val_str = parser->argv[val_pos];
        switch (type) {
            case yc__int:
                if (!sscanf(val_str, "%i", (int*)value)) parser->error = true;
                break;
            case yc__float:
                if (!sscanf(val_str, "%f", (float*)value)) parser->error = true;
                break;
            case yc__double:
                if (!sscanf(val_str, "%lf", (double*)value))
                    parser->error = true;
                break;
            case yc__string: *(const char**)value = val_str; break;
            default: assert(false); break;
        }

        // if parsing error, exit with a message
        if (parser->error) {
            sprintf(parser->error_msg, "bad value for opt %s", longname);
            return;
        }

        // remove argument
        for (int i = val_pos; i < parser->argc - 1; i++)
            parser->argv[i] = parser->argv[i + 1];
        parser->argc -= 1;
    }

    // done
    return;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API bool yc_parse_optb(yc_parser* parser, const char* longname,
                            const char* shortname, const char* help, bool def) {
    bool value = false;
    yc__parse_val(parser, longname, shortname, help, yc__bool, &def, &value, 0);
    return value;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API int yc_parse_opti(yc_parser* parser, const char* longname,
                           const char* shortname, const char* help, int def) {
    int value = 0;
    yc__parse_val(parser, longname, shortname, help, yc__int, &def, &value, 0);
    return value;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API float yc_parse_optf(yc_parser* parser, const char* longname,
                             const char* shortname, const char* help,
                             float def) {
    float value = 0;
    yc__parse_val(parser, longname, shortname, help, yc__float, &def, &value,
                  0);
    return value;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API double yc_parse_optd(yc_parser* parser, const char* longname,
                              const char* shortname, const char* help,
                              double def) {
    double value = 0;
    yc__parse_val(parser, longname, shortname, help, yc__double, &def, &value,
                  0);
    return value;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API const char* yc_parse_opts(yc_parser* parser, const char* longname,
                                   const char* shortname, const char* help,
                                   const char* def) {
    char* value = 0;
    yc__parse_val(parser, longname, shortname, help, yc__string, &def, &value,
                  0);
    return value;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API int yc_parse_opte(yc_parser* parser, const char* longname,
                           const char* shortname, const char* help, int def,
                           const char** choices) {
    char* values = 0;
    yc__parse_val(parser, longname, shortname, help, yc__string,
                  &(choices[def]), &values, choices);
    if (parser->error) return -1;
    for (int count = 0; *choices; count++, choices++) {
        if (!strcmp(*choices, values)) return count;
    }
    assert(false);
    return -1;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API int yc_parse_argi(yc_parser* parser, const char* longname,
                           const char* help, int def, bool required) {
    int value = 0;
    yc__parse_val(parser, longname, 0, help, yc__int, (required) ? 0 : &def,
                  &value, 0);
    return value;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API const char* yc_parse_args(yc_parser* parser, const char* longname,
                                   const char* help, const char* def,
                                   bool required) {
    char* value = 0;
    yc__parse_val(parser, longname, 0, help, yc__string, (required) ? 0 : &def,
                  &value, 0);
    return value;
}

//
// Just a convenience wrapper for yc__parse_val. Public interface.
//
YGLC_API const char** yc_parse_argsa(yc_parser* parser, const char* longname,
                                     const char* help, const char** def,
                                     bool required, const char* values[],
                                     int max_values) {
    yc__parse_val(parser, longname, 0, help, yc__string_array,
                  (required) ? 0 : def, parser->stras_values, 0);
    if (parser->error) return 0;
    int nvalues = 0;
    for (char** c = parser->stras_values; *c; c++) nvalues++;
    if (!values)
        values = (const char**)calloc((nvalues + 1), sizeof(char*));
    else
        assert(nvalues <= max_values);
    memcpy(values, parser->stras_values, sizeof(char*) * (nvalues + 1));
    return values;
}

//
// Parse a binary file a chuck at a time and loads its content full in memory.
// Does not attempt to determine file size upfront to handle all cases.
//
YGLC_API void* yc_load_binfile(const char* filename, int* size) {
    FILE* file = fopen(filename, "rb");
    if (!file) return 0;
    int len = 0;
    char* ret = 0;
    char buf[4096];
    int bufn;
    while ((bufn = (int)fread(buf, 1, sizeof(buf), file))) {
        ret = (char*)realloc(ret, len + bufn);
        memcpy(ret + len, buf, bufn);
        len += bufn;
    }
    if (size) *size = len;
    return ret;
}

//
// Parse a text file a chuck at a time and loads its content full in memory.
// Does not attempt to determine file size upfront to handle all cases.
//
YGLC_API char* yc_load_txtfile(const char* filename) {
    FILE* file = fopen(filename, "rt");
    if (!file) return 0;
    int len = 0;
    char* ret = 0;
    char buf[4096];
    int bufn;
    while ((bufn = (int)fread(buf, 1, sizeof(buf), file))) {
        ret = (char*)realloc(ret, len + bufn + 1);
        memcpy(ret + len, buf, bufn);
        ret[len + bufn] = 0;
        len += bufn;
    }
    return ret;
}

//
// Splits a path. Public interface.
//
YGLC_API void yc_split_path(const char* filename, char* dirname, char* basename,
                            char* ext) {
    // walk till end keeping the position of '/', '\\' and '.'
    const char *path_sep = 0, *ext_sep = 0;
    for (const char* p = filename; *p; p++) {
        if (*p == '/' || *p == '\\') path_sep = p;
        if (*p == '.') ext_sep = p;
    }

    // copy strings
    if (dirname) {
        if (path_sep) {
            strncpy(dirname, filename, 1 + path_sep - filename);
            if (path_sep) dirname[1 + path_sep - filename] = 0;
        } else
            strcpy(dirname, "");
    }
    if (basename) {
        const char* start = (path_sep) ? path_sep + 1 : filename;
        if (ext_sep) {
            strncpy(basename, start, ext_sep - filename);
            if (ext_sep) basename[ext_sep - start] = 0;
        } else
            strcpy(basename, start);
    }
    if (ext) {
        if (ext_sep)
            strcpy(ext, ext_sep);
        else
            strcpy(ext, "");
    }
}

#endif

// -----------------------------------------------------------------------------

#ifdef SGLC_TEST
#include <math.h>

int main(int argc, char** argv) {
    // test empty
    int test0_argc = 1;
    const char* test0_argv[] = {"test0"};
    yc_parser* parser0 = parser0 =
        yc_init_parser(test0_argc, test0_argv, "test0");
    parser0->exit_on_error = false;
    assert(yc_done_parser(parser0) == true);

    // test exit on help
    int test1_argc = 2;
    const char* test1_argv[] = {"test1", "-h"};
    yc_parser* parser1 = parser1 =
        yc_init_parser(test1_argc, test1_argv, "test1");
    parser1->exit_on_error = false;
    assert(yc_done_parser(parser1) == true);

    // test opts
    int test2_argc = 10;
    const char* test2_argv[] = {"test2", "--int",    "10",   "--float",
                                "3.14",  "--double", "6.28", "--str",
                                "bob",   "--flag"};
    yc_parser* parser2 = parser2 =
        yc_init_parser(test2_argc, test2_argv, "test2");
    parser2->exit_on_error = false;
    assert(yc_parse_optb(parser2, "--flag", 0, "", false) == true);
    assert(yc_parse_opti(parser2, "--int", 0, "", 0) == 10);
    assert(fabsf(yc_parse_optf(parser2, "--float", 0, "", 0) - 3.14f) < 0.01);
    assert(fabs(yc_parse_optd(parser2, "--double", 0, "", 0) - 6.28) < 0.01);
    assert(!strcmp(yc_parse_opts(parser2, "--str", 0, "", "mike"), "bob"));
    assert(yc_parse_optb(parser2, "--flag_def", 0, "", false) == false);
    assert(yc_parse_opti(parser2, "--int_def", 0, "", 5) == 5);
    assert(yc_parse_optf(parser2, "--float_def", 0, "", 2.67f) == 2.67f);
    assert(yc_parse_optd(parser2, "--double_def", 0, "", 9.54) == 9.54);
    assert(!strcmp(yc_parse_opts(parser2, "--str_def", 0, "", "alex"), "alex"));
    assert(yc_done_parser(parser2) == true);

    // test args
    int test3_argc = 3;
    const char* test3_argv[] = {"test3", "10", "bob"};
    yc_parser* parser3 = yc_init_parser(test3_argc, test3_argv, "test3");
    parser3->exit_on_error = false;
    assert(yc_parse_argi(parser3, "int", "", 0, true) == 10);
    assert(!strcmp(yc_parse_args(parser3, "str", "", "mike", true), "bob"));
    assert(yc_done_parser(parser3) == true);

    // test bad opts
    int test4_argc = 3;
    const char* test4_argv[] = {"test4", "--int", "bob"};
    yc_parser* parser4 = yc_init_parser(test4_argc, test4_argv, "test4");
    parser4->exit_on_error = false;
    assert(yc_parse_opti(parser4, "int", 0, "", 0) == 0);
    assert(yc_done_parser(parser4) == false);
}

#endif

#endif
