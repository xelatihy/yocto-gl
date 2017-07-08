# Yocto/Utils

Utilities for writing command line applications, mostly
a simple to use command line parser, a logger, a thread-pool,
and string, path and file functions. All functions are defined
in separate namespaces.

## Usage for Command Line Parsing

1. namespace logging
2. create a parser object
    - an option for printing help is automatically added
    `parser = make_parser(argc, argv, program description)`
3. for each option, parse it calling the functions `parse_opt()`
    - options are parsed on the fly and a comprehensive help is
      automatically generated
    - supports bool (flags), int, float, double, std::string, enums
    - options names are "--longname" for longname and "-s" for short
    - command line format is "--longname value", "-s v" for all but flags
    - values are parsed with `iostream <<` operators
    - for general use `opt = parse_opt<type>()`
    - for boolean flags is `parse_flag()`
    - for enums use `parse_opte()`
4. for each unnamed argument, parse it calling the functions parse_arg()
    - names are only used for help
    - supports types as above
    - for general use `arg = parse_arg<type>()`
    - to parse all remaining values use args = parse_arga<type>(...)
5. end cmdline parsing with `check_parser()` to check for unsued values,
   missing arguments and print help if needed
6. since arguments are parsed immediately, one can easily implement
   subcommands by just branching the command line code based on a read
   argument without any need for complex syntax

Notes: the end of this file contains a test function that also
illustreates the library usage.

## Usage for Logging

1. namespace logging
2. create loggers with `make_file_logger()`, `make_stderr_logger()`,
   `make_stdout_logger()`
3. you can set default loggers with `get_default_loggers()`; note that none
   are set by default
4. write log messages with `log_msg()` and its variants.

## Usage for Concurrent Execution

1. namespace concurrent
2. either create a thread pool `make_thread_pool` or use the global one
3. run tasks in parallel `thread_pool_for()`
4. run tasks asynchronously `thread_pool_async()

## Utilities

1. filename splitting functions in namespace path
2. loading and save entire files in namespace file
3. Python-line string manipulation in namespace string
5. Python-like operators for standard containers in namespace stl_operators


## History

- v 0.17: renamed to yocto utils
- v 0.16: split into namespaces
- v 0.15: remove inline compilation
- v 0.14: Python-like operator for std::vector
- v 0.13: more file and string utilities
- v 0.12: better thread pool implementation
- v 0.11: added a few more path utilities
- v 0.10: changed default name for help option; better help printing
- v 0.9: C-like string formatting
- v 0.8: switch to .h/.cpp pair (templated functions are specialized)
- v 0.7: logging
- v 0.6: doxygen comments
- v 0.5: added a few python-like std::string manipulation functions.
- v 0.4: [major API change] move to modern C++ interface
- v 0.3: adding a few functions for path splitting.
- v 0.2: [API change] C++ API
- v 0.1: C++ implementation
- v 0.0: initial release in C99

## Namespace yu

General utilities to write command line applications

## Namespace cmdline

Command line parsing

### Struct parser

~~~ .cpp
struct parser;
~~~

Command line parser.

### Function make_parser()

~~~ .cpp
parser* make_parser(
    const std::vector<std::string>& args, const std::string& help = "");
~~~

Inits a command line parser.

### Function make_parser()

~~~ .cpp
parser* make_parser(int argc, char* argv[], const char* help);
~~~

Inits a command line parser.

### Function check_parser()

~~~ .cpp
bool check_parser(parser* prs);
~~~

Ends parsing checking for error for unused options or arguments.
Exit if needed.

### Function parse_flag()

~~~ .cpp
bool parse_flag(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def = false);
~~~

Parses an optional flag as described in the intro.

### Function parse_opt()

~~~ .cpp
template <typename T>
T parse_opt(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    bool required = false, const std::vector<T>& choices =;
~~~

Parses an option as described in the intro.

### Function parse_opts()

~~~ .cpp
std::string parse_opts(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help,
    const std::string& def, bool required = false,
    const std::vector<std::string>& choices =;
~~~

Specialization of parse_opt()

### Function parse_opti()

~~~ .cpp
int parse_opti(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    bool required = false, const std::vector<int>& choices =;
~~~

Specialization of parse_opt()

### Function parse_optf()

~~~ .cpp
float parse_optf(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, float def,
    bool required = false, const std::vector<float>& choices =;
~~~

Specialization of parse_opt()

### Function parse_optd()

~~~ .cpp
double parse_optd(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, double def,
    bool required = false, const std::vector<double>& choices =;
~~~

Specialization of parse_opt()

### Function parse_optb()

~~~ .cpp
bool parse_optb(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def,
    bool required = false);
~~~

Specialization of parse_opt()

### Function parse_opte()

~~~ .cpp
int parse_opte(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    const std::vector<std::pair<std::string, int>>& vals,
    bool required = false);
~~~

Parses an option enum as described in the intro.

### Function parse_opta()

~~~ .cpp
template <typename T>
T parse_opta(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    int nargs, bool required = false, const std::vector<T>& choices =;
~~~

Parses an option array as described in the intro.

### Function parse_arg()

~~~ .cpp
template <typename T>
T parse_arg(parser* par, const std::string& longname, const std::string& help,
    const T& def, bool required = false, const std::vector<T>& choices =;
~~~

Parses an argument as described in the intro.

### Function parse_args()

~~~ .cpp
std::string parse_args(parser* par, const std::string& longname,
    const std::string& help, const std::string& def, bool required = false,
    const std::vector<std::string>& choices =;
~~~

Specialization of parse_arg()

### Function parse_argi()

~~~ .cpp
int parse_argi(parser* par, const std::string& longname,
    const std::string& help, int def, bool required = false,
    const std::vector<int>& choices =;
~~~

Specialization of parse_arg()

### Function parse_argf()

~~~ .cpp
float parse_argf(parser* par, const std::string& longname,
    const std::string& help, float def, bool required = false,
    const std::vector<float>& choices =;
~~~

Specialization of parse_arg()

### Function parse_argd()

~~~ .cpp
double parse_argd(parser* par, const std::string& longname,
    const std::string& help, double def, bool required = false,
    const std::vector<double>& choices =;
~~~

Specialization of parse_arg()

### Function parse_arga()

~~~ .cpp
template <typename T>
std::vector<T> parse_arga(parser* par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs = -1,
    bool required = false, const std::vector<T>& choices =;
~~~

Parses an argument array as described in the intro.

### Function parse_argas()

~~~ .cpp
std::vector<std::string> parse_argas(parser* par, const std::string& longname,
    const std::string& help, const std::vector<std::string>& def,
    int nargs = -1, bool required = false,
    const std::vector<std::string>& choices =;
~~~

Specialization of parse_arga()

### Function parse_argai()

~~~ .cpp
std::vector<int> parse_argai(parser* par, const std::string& longname,
    const std::string& help, const std::vector<int>& def, int nargs = -1,
    bool required = false, const std::vector<int>& choices =;
~~~

Specialization of parse_arga()

### Function parse_argaf()

~~~ .cpp
std::vector<float> parse_argaf(parser* par, const std::string& longname,
    const std::string& help, const std::vector<float>& def, int nargs = -1,
    bool required = false, const std::vector<float>& choices =;
~~~

Specialization of parse_arga()

### Function parse_argad()

~~~ .cpp
std::vector<double> parse_argad(parser* par, const std::string& longname,
    const std::string& help, const std::vector<double>& def, int nargs = -1,
    bool required = false, const std::vector<double>& choices =;
~~~

Specialization of parse_arga()

## Namespace file

File loading and saving

### Function load_binfile()

~~~ .cpp
std::vector<unsigned char> load_binfile(const std::string& filename);
~~~

Loads the contents of a binary file in an in-memory array.

### Function load_txtfile()

~~~ .cpp
std::string load_txtfile(const std::string& filename);
~~~

Loads the contents of a text file into a std::string.

### Function save_binfile()

~~~ .cpp
void save_binfile(
    const std::string& filename, const std::vector<unsigned char>& data);
~~~

Saves binary data to a file.

### Function save_txtfile()

~~~ .cpp
void save_txtfile(const std::string& filename, const std::string& str);
~~~

Saves a string to a text file.

## Namespace path

Path manipulation

### Function get_dirname()

~~~ .cpp
std::string get_dirname(const std::string& filename);
~~~

Get directory name (including '/').

### Function get_basename()

~~~ .cpp
std::string get_basename(const std::string& filename);
~~~

Get file basename.

### Function get_extension()

~~~ .cpp
std::string get_extension(const std::string& filename);
~~~

Get extension (including '.').

### Function get_filename()

~~~ .cpp
std::string get_filename(const std::string& filename);
~~~

Get filename without directory (equiv to get_basename() +
get_extension()).

### Function replace_extension()

~~~ .cpp
std::string replace_extension(
    const std::string& filename, const std::string& ext);
~~~

Replace extension.

### Function prepend_extension()

~~~ .cpp
std::string prepend_extension(
    const std::string& filename, const std::string& prep);
~~~

Prepend a string to the extension.

### Function split_path()

~~~ .cpp
void split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext);
~~~

Splits a path calling the above functions.

## Namespace string

String manipulation

### Function starts_with()

~~~ .cpp
bool starts_with(const std::string& str, const std::string& substr);
~~~

Checks if a std::string starts with a prefix.

### Function ends_with()

~~~ .cpp
bool ends_with(const std::string& str, const std::string& substr);
~~~

Checks if a std::string ends with a prefix.

### Function contains()

~~~ .cpp
bool contains(const std::string& str, const std::string& substr);
~~~

Check is a string contains a substring.

### Function split_lines()

~~~ .cpp
std::vector<std::string> split_lines(
    const std::string& str, bool keep_newline = false);
~~~

Splits a std::string into lines at the '\n' character. The line
terminator is kept if keep_newline. This function does not work on
Window if keep_newline is true.

### Function partition_str()

~~~ .cpp
std::vector<std::string> partition_str(
    const std::string& str, const std::string& split);
~~~

Partition the string.

### Function split_str()

~~~ .cpp
std::vector<std::string> split_str(const std::string& str);
~~~

Splits the string.

### Function strip_str()

~~~ .cpp
std::string strip_str(const std::string& str);
~~~

Strip the string.

### Function rstrip_str()

~~~ .cpp
std::string rstrip_str(const std::string& str);
~~~

Strip the string.

### Function lstrip_str()

~~~ .cpp
std::string lstrip_str(const std::string& str);
~~~

Strip the string.

### Function join_strings()

~~~ .cpp
std::string join_strings(
    const std::vector<std::string>& strs, const std::string& sep);
~~~

Joins a list of std::string with a std::string as separator.

### Function to_lower()

~~~ .cpp
std::string to_lower(const std::string& str);
~~~

Converts an ASCII string to lowercase.

### Function to_upper()

~~~ .cpp
std::string to_upper(const std::string& str);
~~~

Converts an ASCII string to uppercase.

### Function is_space()

~~~ .cpp
bool is_space(const std::string& str);
~~~

Strung is space.

### Function replace_str()

~~~ .cpp
std::string replace_str(
    const std::string& str, const std::string& s1, const std::string& s2);
~~~

Replace s1 with s2 in str.

### Function format_str()

~~~ .cpp
std::string format_str(const char* fmt, va_list args);
~~~

C-like string formatting. This is only meant for short strings with max
length 10000 chars. Memory corruption will happen for longer strings.

### Function format_str()

~~~ .cpp
std::string format_str(const char* fmt, ...);
~~~

C-like string formatting. See format_str(fmt,args);

## Namespace operators

Python-like STL operators

### Function operator+()

~~~ .cpp
template <typename T>
inline std::vector<T> operator+(const std::vector<T>& v, const T& vv);
~~~

Append an element to a vector

### Function operator+=()

~~~ .cpp
template <typename T>
inline std::vector<T>& operator+=(std::vector<T>& v, const T& vv);
~~~

Append an element to a vector

### Function operator+()

~~~ .cpp
template <typename T, typename ET>
inline std::vector<T> operator+(const std::vector<T>& v, const ET& vv);
~~~

Append an element to a vector

### Function operator+=()

~~~ .cpp
template <typename T, typename ET>
inline std::vector<T>& operator+=(std::vector<T>& v, const ET& vv);
~~~

Append an element to a vector

### Function operator+()

~~~ .cpp
template <typename T>
inline std::vector<T> operator+(
    const std::vector<T>& v, const std::vector<T>& vv);
~~~

Append a vector to a vector

### Function operator+=()

~~~ .cpp
template <typename T>
inline std::vector<T>& operator+=(std::vector<T>& v, const std::vector<T>& vv);
~~~

Append a vector to a vector

## Namespace logging

Simple Logging

### Enum log_level : int

~~~ .cpp
enum struct log_level : int {
    verbose = -3,
    trace = -2,
    debug = -1,
    info = 0,
    warning = 1,
    error = 2,
    fatal = 3,
}
~~~

Logging level

- Values:
    - verbose:      verbose
    - trace:      trace
    - debug:      debug
    - info:      info
    - warning:      warning
    - error:      error
    - fatal:      fatal


### Struct logger

~~~ .cpp
struct logger;
~~~

Logger object

### Function make_file_logger()

~~~ .cpp
logger* make_file_logger(const std::string& filename, bool append,
    log_level output_level = log_level::info,
    log_level flush_level = log_level::info);
~~~

Make a file logger. See set_logger() for Parameters.

- Parameters:
    - filename: logger filename or stderr if empty
    - apend: append or write open mode for file logger

### Function make_stderr_logger()

~~~ .cpp
logger* make_stderr_logger(log_level output_level = log_level::info,
    log_level flush_level = log_level::info);
~~~

Make a stderr logger. See set_logger() for Parameters.

### Function make_stdout_logger()

~~~ .cpp
logger* make_stdout_logger(log_level output_level = log_level::info,
    log_level flush_level = log_level::info);
~~~

Make a stderr logger. See set_logger() for Parameters.

### Function clear_logger()

~~~ .cpp
void clear_logger(logger* lgr);
~~~

Clear a logger

### Function get_default_loggers()

~~~ .cpp
std::vector<logger*>* get_default_loggers();
~~~

Get default loggers. This is a modifiable reference.

### Function set_logger()

~~~ .cpp
void set_logger(logger* lgr, log_level output_level,
    log_level flush_level = log_level::error);
~~~

Set logger level

- Parameters:
    - lgr: logger
    - name : logger name
    - output_level: output level
    - flush_level: output level

### Function log_msg()

~~~ .cpp
void log_msg(logger* lgr, log_level level, const std::string& name,
    const std::string& msg);
~~~

Log a message

- Parameters:
    - lgr: logger
    - level: message level
    - code: message code (5 chars)
    - msg: message

### Function log_msg()

~~~ .cpp
void log_msg(logger* lgr, log_level level, const char* name, const char* msg);
~~~

Log a message

- Parameters:
    - lgr: logger
    - level: message level
    - code: message code (5 chars)
    - msg: message

### Function log_msgf()

~~~ .cpp
void log_msgf(
    logger* lgr, log_level level, const char* name, const char* msg, ...);
~~~

Log a message formatted ala printf.

- Parameters:
    - lgr: logger
    - level: message level
    - code: message code (5 chars)
    - msg: message

### Function log_msgfv()

~~~ .cpp
void log_msgfv(logger* lgr, log_level level, const char* name, const char* msg,
    va_list args);
~~~

Logs a message to the default loggers

### Function log_msg()

~~~ .cpp
void log_msg(log_level level, const std::string& name, const std::string& msg);
~~~

Logs a message to the default loggers

### Function log_msgf()

~~~ .cpp
void log_msgf(log_level level, const char* name, const char* msg, ...);
~~~

Logs a message to the default loggers

### Function log_msgfv()

~~~ .cpp
void log_msgfv(
    log_level level, const char* name, const char* msg, va_list args);
~~~

Logs a message to the default loggers

## Namespace concurrent

Concurrent  and async execution based on a thread pool

### Struct thread_pool

~~~ .cpp
struct thread_pool;
~~~

Forward declaration of thread pool.

### Function make_pool()

~~~ .cpp
thread_pool* make_pool(int nthread = 0);
~~~

Initialize a thread pool with a certain number of threads (0 for
defatul).

### Function clear_pool()

~~~ .cpp
void clear_pool(thread_pool* pool);
~~~

Clear thread pool

### Function wait_pool()

~~~ .cpp
void wait_pool(thread_pool* pool);
~~~

Wait for all jobs to finish

### Function parallel_for()

~~~ .cpp
void parallel_for(
    thread_pool* pool, int count, const std::function<void(int idx)>& task);
~~~

Parallel for implementation

### Function run_async()

~~~ .cpp
std::shared_future<void> run_async(
    thread_pool* pool, const std::function<void()>& task);
~~~

Runs a task asynchronously onto a thread pool

### Function wait_pool()

~~~ .cpp
void wait_pool();
~~~

Wait for all jobs to finish on a global thread pool

### Function run_async()

~~~ .cpp
std::shared_future<void> run_async(const std::function<void()>& task);
~~~

Runs a task asynchronously onto a global thread pool

### Function parallel_for()

~~~ .cpp
void parallel_for(int count, const std::function<void(int idx)>& task);
~~~

Parallel for implementation on a global thread pool

