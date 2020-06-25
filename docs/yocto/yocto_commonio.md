# Yocto/CommonIO: Utilities for writing command-line apps

Yocto/CommonIO is a collection of utilities used in writing command-line
applications, including parsing command line arguments, simple path
manipulation, file lading and saving, and printing values, timers and
progress bars.

## Printing values

Use `print_info(message)` to print a message, and `print_fatal(message)`
to print and exit. To time a block of code use `print_timed(message)`
to use an RIIA timer or call `print_elapsed(timer)` to print the elapsed
time as needed. Use `print_progress(message, current, total)` to print
a progress bar with a message for a given current and total number of tasks.

```cpp
print_info("Message");                      // print message
print_fatal("Error and exit");              // print error and exit
{ auto timer = print_timed("Timer"); ... }  // time a block of code
auto timer = print_timed("Timer");          // start timer
print_elapsed(timer);                       // Print elapsed time
for(auto task : range(10))                  // iterate over tasks
  print_progress("Progress bar", task, 10); // print progress bar
```

## Command-Line Parsing

Yocto/CommonIO includes a simple command-line parser that supports optional
and positional arguments, automatic help generation, and error checking.

The command-line parser is initialized with `make_cli(name, help)` that takes
a program name and a help string. Then add command-line options with
`add_option(cli, name, value, help, req)` where `name` is the option name,
`value` is a reference to the variable that we want to set, `help` is a usage
message and `req` is a flag that determines whether or not the argument is
required. After adding all arguments, use `parse_cli(cli)` to parse the
command-line. If an error occurs, the parser exists. A help flag is also
handled in this function.

The argument name determines the argument type. If the arguments name starts
with '--' or '-' then it is an option, otherwise it is a positional argument.
An argument may have multiple names separated by commas. Options and arguments
may be intermixed. The type of each option is determined by the passed reference.
The parser supports integers, floating point numbers, strings, bool options
and arrays of strings. Boolean flags are indicated with a pair of names
"--name/--no-name", so that we have both options available.

```cpp
auto samples = 10; auto flag = false;                // state
auto out = ""s, name = ""s;
auto names = vector<string>{};
auto cli = make_cli("app", "testing cli");           // initialize cli
add_option(cli, "--out", out, "out", true);          // required argument
add_option(cli, "--samples,-s", samples, "samples"); // optional argument
add_option(cli, "--flag/--no-flag", flag, "flag");   // optional flag
add_option(cli, "name", name, "name", true);         // positional argument
add_option(cli, "names", names, "names");            // positional arguments
parse_cli(cli);
```

## Text and binary serialization

Text files are loaded with `load_text(filename, text, error)` and saved with
`save_text(filename, text, error)`.
Binary files are loaded with `load_binary(filename, binary, error)` and saved
with `save_binary(filename, binary, error)`.
Both loading and saving take a filename, a text or binary buffer, and return
whether or not the file was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.

```cpp
auto error = string{};                  // error buffer
auto text  = string{};                  // text buffer
if(!load_text(filename, text, error))   // load a text file
  print_error(error);                   // check and print error
if(!save_text(filename, text, error))   // save a text file
  print_error(error);                   // check and print error
auto data  = vector<byte>{};            // data buffer
if(!load_binary(filename, data, error)) // load a binary file
  print_error(error);                   // check and print error
if(!save_binary(filename, data, error)) // save a binary file
  print_error(error);                   // check and print error
```
