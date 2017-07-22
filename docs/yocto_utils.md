# Yocto/Utils

Utilities for writing command line applications, mostly
a simple to use command line parser, a logger, a thread-pool,
and string, path and file functions, and a timer. All functions are defined
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
    - to parse all remaining values use `args = parse_arga<type>(...)`
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
4. run tasks asynchronously `thread_pool_async()`

## Utilities

1. filename splitting functions in namespace `path`
2. loading and save entire files in namespace `file`
3. Python-line string manipulation in namespace `string`
4. Python-like operators for standard containers in namespace `operators`
5. simple timer in namespace `timer`
6. hashing functions in the `hashing` namespace (SHA1 and xxHash)


## History

- v 0.22: seimpler logging
- v 0.21: move to header-only mode
- v 0.20: simpler logging
- v 0.19: some containers ops
- v 0.18: timer
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
inline parser* make_parser(const std::vector<std::string>& args,
    const std::string& name = "", const std::string& help = "");
~~~

Inits a command line parser.

### Function make_parser()

~~~ .cpp
inline parser* make_parser(
    int argc, char* argv[], const char* name, const char* help);
~~~

Inits a command line parser.

### Function check_parser()

~~~ .cpp
inline bool check_parser(parser* prs);
~~~

Ends parsing checking for error for unused options or arguments.
Exit if needed.

### Function parse_flag()

~~~ .cpp
inline bool parse_flag(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, bool def = false);
~~~

Parses an optional flag as described in the intro.

### Function parse_opt()

~~~ .cpp
template <typename T>
inline T parse_opt(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    bool required = false, const std::vector<T>& choices =;
~~~

Parses an option as described in the intro.

### Function parse_opts()

~~~ .cpp
inline std::string parse_opts(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help,
    const std::string& def, bool required = false,
    const std::vector<std::string>& choices =;
~~~

Specialization of parse_opt()

### Function parse_opti()

~~~ .cpp
inline int parse_opti(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    bool required = false, const std::vector<int>& choices =;
~~~

Specialization of parse_opt()

### Function parse_optf()

~~~ .cpp
inline float parse_optf(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, float def,
    bool required = false, const std::vector<float>& choices =;
~~~

Specialization of parse_opt()

### Function parse_optd()

~~~ .cpp
inline double parse_optd(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, double def,
    bool required = false, const std::vector<double>& choices =;
~~~

Specialization of parse_opt()

### Function parse_opte()

~~~ .cpp
inline int parse_opte(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, int def,
    const std::vector<std::pair<std::string, int>>& vals,
    bool required = false);
~~~

Parses an option enum as described in the intro.

### Function parse_opte()

~~~ .cpp
template <typename T>
inline T parse_opte(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, T def,
    const std::vector<std::pair<std::string, T>>& vals, bool required = false);
~~~

Parses an option enum as described in the intro.

### Function parse_opta()

~~~ .cpp
template <typename T>
inline T parse_opta(parser* par, const std::string& longname,
    const std::string& shortname, const std::string& help, const T& def,
    int nargs, bool required = false, const std::vector<T>& choices =;
~~~

Parses an option array as described in the intro.

### Function parse_arg()

~~~ .cpp
template <typename T>
inline T parse_arg(parser* par, const std::string& longname,
    const std::string& help, const T& def, bool required = false,
    const std::vector<T>& choices =;
~~~

Parses an argument as described in the intro.

### Function parse_args()

~~~ .cpp
inline std::string parse_args(parser* par, const std::string& longname,
    const std::string& help, const std::string& def, bool required = false,
    const std::vector<std::string>& choices =;
~~~

Specialization of parse_arg()

### Function parse_argi()

~~~ .cpp
inline int parse_argi(parser* par, const std::string& longname,
    const std::string& help, int def, bool required = false,
    const std::vector<int>& choices =;
~~~

Specialization of parse_arg()

### Function parse_argf()

~~~ .cpp
inline float parse_argf(parser* par, const std::string& longname,
    const std::string& help, float def, bool required = false,
    const std::vector<float>& choices =;
~~~

Specialization of parse_arg()

### Function parse_argd()

~~~ .cpp
inline double parse_argd(parser* par, const std::string& longname,
    const std::string& help, double def, bool required = false,
    const std::vector<double>& choices =;
~~~

Specialization of parse_arg()

### Function parse_arga()

~~~ .cpp
template <typename T>
inline std::vector<T> parse_arga(parser* par, const std::string& longname,
    const std::string& help, const std::vector<T>& def, int nargs = -1,
    bool required = false, const std::vector<T>& choices =;
~~~

Parses an argument array as described in the intro.

### Function parse_argas()

~~~ .cpp
inline std::vector<std::string> parse_argas(parser* par,
    const std::string& longname, const std::string& help,
    const std::vector<std::string>& def, int nargs = -1, bool required = false,
    const std::vector<std::string>& choices =;
~~~

Specialization of parse_arga()

## Namespace file

File loading and saving

### Function load_binfile()

~~~ .cpp
inline std::vector<unsigned char> load_binfile(const std::string& filename);
~~~

Loads the contents of a binary file in an in-memory array.

### Function load_txtfile()

~~~ .cpp
inline std::string load_txtfile(const std::string& filename);
~~~

Loads the contents of a text file into a std::string.

### Function save_binfile()

~~~ .cpp
inline bool save_binfile(
    const std::string& filename, const std::vector<unsigned char>& data);
~~~

Saves binary data to a file.

### Function save_txtfile()

~~~ .cpp
inline bool save_txtfile(const std::string& filename, const std::string& str);
~~~

Saves a string to a text file.

## Namespace path

Path manipulation

### Function get_dirname()

~~~ .cpp
inline std::string get_dirname(const std::string& filename);
~~~

Get directory name (including '/').

### Function get_extension()

~~~ .cpp
inline std::string get_extension(const std::string& filename);
~~~

Get extension (including '.').

### Function get_basename()

~~~ .cpp
inline std::string get_basename(const std::string& filename);
~~~

Get file basename.

### Function get_filename()

~~~ .cpp
inline std::string get_filename(const std::string& filename);
~~~

Get filename without directory (equiv to get_basename() +
get_extension()).

### Function replace_extension()

~~~ .cpp
inline std::string replace_extension(
    const std::string& filename, const std::string& ext);
~~~

Replace extension.

### Function prepend_extension()

~~~ .cpp
inline std::string prepend_extension(
    const std::string& filename, const std::string& prep);
~~~

Prepend a string to the extension.

### Function split_path()

~~~ .cpp
inline void split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext);
~~~

Splits a path calling the above functions.

## Namespace string

String manipulation

### Function startswith()

~~~ .cpp
inline bool startswith(const std::string& str, const std::string& substr);
~~~

Checks if a std::string starts with a prefix.

### Function endswith()

~~~ .cpp
inline bool endswith(const std::string& str, const std::string& substr);
~~~

Checks if a std::string ends with a prefix.

### Function contains()

~~~ .cpp
inline bool contains(const std::string& str, const std::string& substr);
~~~

Check is a string contains a substring.

### Function splitlines()

~~~ .cpp
inline std::vector<std::string> splitlines(
    const std::string& str, bool keep_newline = false);
~~~

Splits a std::string into lines at the '\n' character. The line
terminator is kept if keep_newline. This function does not work on
Window if keep_newline is true.

### Function partition()

~~~ .cpp
inline std::vector<std::string> partition(
    const std::string& str, const std::string& split);
~~~

Partition the string.

### Function split()

~~~ .cpp
inline std::vector<std::string> split(const std::string& str);
~~~

Splits the string.

### Function rstrip()

~~~ .cpp
inline std::string rstrip(const std::string& str);
~~~

Strip the string.

### Function lstrip()

~~~ .cpp
inline std::string lstrip(const std::string& str);
~~~

Strip the string.

### Function strip()

~~~ .cpp
inline std::string strip(const std::string& str);
~~~

Strip the string.

### Function join()

~~~ .cpp
inline std::string join(
    const std::vector<std::string>& strs, const std::string& sep);
~~~

Joins a list of std::string with a std::string as separator.

### Function lower()

~~~ .cpp
inline std::string lower(const std::string& str);
~~~

Converts an ASCII string to lowercase.

### Function upper()

~~~ .cpp
inline std::string upper(const std::string& str);
~~~

Converts an ASCII string to uppercase.

### Function isspace()

~~~ .cpp
inline bool isspace(const std::string& str);
~~~

Check if a string is space.

### Function replace()

~~~ .cpp
inline std::string replace(
    const std::string& str, const std::string& s1, const std::string& s2);
~~~

Replace s1 with s2 in str.

### Function format()

~~~ .cpp
template <typename... Args>
inline std::string format(const char* fmt, Args&&... args);
~~~

C-like string formatting. This is only meant for short strings with max
length 10000 chars. Memory corruption will happen for longer strings.

## Namespace containers

Python-like STL functions

### Function contains()

~~~ .cpp
template <typename T>
inline bool contains(const std::vector<T>& v, const T& vv);
~~~

Checks if a containers contains a value

### Function contains()

~~~ .cpp
template <typename K, typename V>
inline bool contains(const std::map<K, V>& v, const K& vv);
~~~

Checks if a containers contains a value

### Function contains()

~~~ .cpp
template <typename K, typename V>
inline bool contains(const std::unordered_map<K, V>& v, const K& vv);
~~~

Checks if a containers contains a value

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

Logger object. A logger can output messages to multiple streams.
Use add streams commands for it.

### Function make_logger()

~~~ .cpp
inline logger* make_logger(const char* name, bool add_console_stream = false);
~~~

Make a logger with an optional console stream.

### Function set_logger_name()

~~~ .cpp
inline void set_logger_name(logger* lgr, const char* name);
~~~

Set logger default name

### Function free_logger()

~~~ .cpp
inline void free_logger(logger*& lgr);
~~~

Free logger

### Function add_file_stream()

~~~ .cpp
inline bool add_file_stream(logger* lgr, const std::string& filename,
    bool append, bool short_message = false,
    log_level output_level = log_level::info,
    log_level flush_level = log_level::info);
~~~

Add a file stream to a logger.

- Parameters:
    - lgr: logger
    - filename: filename
    - append: append or write open mode for file logger
    - short_message: whether to use a short message version
    - output_level: output level
    - flush_level: output level
- Returns:
    - true if ok

### Function add_console_stream()

~~~ .cpp
inline bool add_console_stream(logger* lgr, bool use_std_error = false,
    bool short_message = true, log_level output_level = log_level::info,
    log_level flush_level = log_level::info);
~~~

Add a console stream to a logger.

- Parameters:
    - lgr: logger
    - filename: logger filename or stderr if empty
    - use_std_error: use standard error instead of standard out
    - short_message: whether to use a short message version
    - output_level: output level
    - flush_level: output level
- Returns:
    - true if ok

### Function get_default_logger()

~~~ .cpp
inline logger* get_default_logger();
~~~

Get default logger.
By default a non-verbose stdout logger is creater.

### Function set_logger_name()

~~~ .cpp
inline void set_logger_name(const char* name);
~~~

Set default logger name

### Function add_file_stream()

~~~ .cpp
inline void add_file_stream(const std::string& filename, bool append,
    bool short_message = false, log_level output_level = log_level::info,
    log_level flush_level = log_level::info);
~~~

Add a file logger to the default loggers.

### Function log_msg()

~~~ .cpp
inline void log_msg(
    logger* lgr, log_level level, const char* name, const char* msg);
~~~

Log a message

- Parameters:
    - lgr: logger
    - level: message level
    - code: message code (5 chars)
    - msg: message

### Function log_msg()

~~~ .cpp
template <typename... Args>
inline void log_msg(logger* lgr, log_level level, const char* name,
    const char* msg, const Args&... args);
~~~

Log a message formatted ala printf.

- Parameters:
    - lgr: logger
    - level: message level
    - code: message code (5 chars)
    - msg: message

### Function log_msg()

~~~ .cpp
template <typename... Args>
inline void log_msg(
    log_level level, const char* name, const char* msg, const Args&... args);
~~~

Logs a message to the default loggers

### Function log_info()

~~~ .cpp
template <typename... Args>
inline void log_info(const char* msg, const Args&... args);
~~~

Logs a message to the default loggers

### Function log_error()

~~~ .cpp
template <typename... Args>
inline void log_error(const char* msg, const Args&... args);
~~~

Logs a message to the default loggers

### Function log_fatal()

~~~ .cpp
template <typename... Args>
inline void log_fatal(const char* msg, const Args&... args);
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
inline thread_pool* make_pool(int nthread = 0);
~~~

Initialize a thread pool with a certain number of threads (0 for
defatul).

### Function free_pool()

~~~ .cpp
inline void free_pool(thread_pool*& pool);
~~~

Free the thread pool

### Function wait_pool()

~~~ .cpp
inline void wait_pool(thread_pool* pool);
~~~

Wait for all jobs to finish

### Function clear_pool()

~~~ .cpp
inline void clear_pool(thread_pool* pool);
~~~

Clear all jobs

### Function parallel_for()

~~~ .cpp
inline void parallel_for(
    thread_pool* pool, int count, const std::function<void(int idx)>& task);
~~~

Parallel for implementation

### Function run_async()

~~~ .cpp
inline std::shared_future<void> run_async(
    thread_pool* pool, const std::function<void()>& task);
~~~

Runs a task asynchronously onto a thread pool

### Function wait_pool()

~~~ .cpp
inline void wait_pool();
~~~

Wait for all jobs to finish on a global thread pool

### Function clear_pool()

~~~ .cpp
inline void clear_pool();
~~~

Clear all jobs on a global thread pool

### Function run_async()

~~~ .cpp
inline std::shared_future<void> run_async(const std::function<void()>& task);
~~~

Runs a task asynchronously onto a global thread pool

### Function parallel_for()

~~~ .cpp
inline void parallel_for(int count, const std::function<void(int idx)>& task);
~~~

Parallel for implementation on a global thread pool

## Namespace timer

Simple timer for performance measumrents.

### Struct timer

~~~ .cpp
struct timer {
    timer(bool autostart = true); 
    void start(); 
    void stop(); 
    double elapsed(); 
}
~~~

A simple wrapper for std::chrono.

- Members:
    - timer():      initialize a timer and start it if necessary
    - start():      start a timer
    - stop():      stops a timer
    - elapsed():      elapsed time


