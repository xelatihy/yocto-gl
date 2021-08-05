# Yocto/Cli: Utilities for writing command-line apps

Yocto/Cli is a collection of utilities used in writing command-line
applications, including parsing command line arguments, printing values,
timers and progress bars.
Yocto/Cli is implemented in `yocto_cli.h` and `yocto_cli.cpp`.

## Printing values

Use `print_info(message)` to print a message, and `print_fatal(message)`
to print and exit. To time a block of code use `print_timed(message)`
to use an RIIA timer or call `print_elapsed(timer)` to print the elapsed
time as needed. Use `print_progress(message, current, total)` to print
a progress bar with a message for a given current and total number of tasks,
or `print_progress_start(message, total)`, `print_progress_next()` and
`print_progress_end()` to use a progress bar without managing manually the
progress counter.

```cpp
print_info("Message");                      // print message
print_fatal("Error and exit");              // print error and exit
{ auto timer = print_timed("Timer"); ... }  // time a block of code
auto timer = print_timed("Timer");          // start timer
print_elapsed(timer);                       // Print elapsed time
for(auto task : range(10))                  // iterate over tasks
  print_progress("Progress bar", task, 10); // print progress bar
print_progress_start("Progress bar", 10);   // start progress bar
for(auto task : range(10))                  // iterate over tasks
  print_progress_next();                    // print progress bar
```

## Command-Line Parsing

Yocto/Cli includes a simple command-line parser that supports optional
and positional arguments, automatic help generation, and error checking.

The command-line parser is initialized with `make_cli(name, help)` that takes
a program name and a help string. Then add command-line options with
`add_option(cli, name, value, help, check, alt, req)` for optional arguments and
`add_argument(cli, name, value, help, check, req)` for positional arguments.
In these commands, `name` is the option name, `value` is a reference to the
variable that we want to set, `help` is a usage message, `alt` is the short
option flag, and `req` determines whether the option is required, and `check`
is either a range of valid values er a list of valid choices.
A help flag is automatically added.

After adding all arguments, use `parse_cli(cli, args)` to parse the
command-line. If an error occurs, the parser throws a `cli_error` exception,
that includes a message that can be printed to inform the user.
If you cannot use exceptions, use the function `parse_cli(cli, args, error)`
that returns whether an error has occurred as a boolean and the error string.

For positional arguments, the argument name is only used for printing help.
For optional arguments, the command line is parsed for options named `--` 
followed by the option `name` or `-` followed by the `alt` short name.
The type of each option is determined by the passed reference.
The parser supports integers, floating point numbers, strings, boolean flags,
enums and arrays of strings.

An help command is added automatically and checked for printing help and exiting
without errors.

```cpp
auto samples = 10; auto flag = false;                // state
auto out = ""s, scene = ""s;
auto scenes = vector<string>{};
auto cli = make_cli("app", "testing cli");      // initialize cli
add_option(cli, "out", out, "out");             // optional argument
add_option(cli, "samples", samples, "samples"); // optional argument
add_option(cli, "flag", flag, "flag");          // optional flag
add_argument(cli, "scene", scene, "scene");     // positional argument
add_option(cli, "scenes", scenes, "scenes");    // positional arguments
auto args = make_cli_args(argc, argv);          // make a vector of args
try {
  parse_cli(cli, args);                         // parse args
} catch(const cli_error& error) {
  print_info(error.what());                     // exit with error
  exit(1);
} catch(const cli_help& help) {
  print_info(help.what());                      // exit with help
  exit(0);
} catch (const std::exception& error) {               
  print_info(error.what());                     // exit with other error
  exit(1);
}
```

You can avoid using writing the exception chain above by using the 
`handle_errors(func)` or `handle_errors(func, args)` function that executes a 
function and handle all errors internally.

```cpp
void app(const vector<string>& args) { ... }    // parse cli and run
handle_errors(app, make_cli_args(argc, argv));  // run and handle errors
```

## Command-Line Commands

The command line parser also support commands, that are created by
`add_command(cli, name, help)`, where `name` is the command name and
`help` is the command usage message. You can add options to commands
using the functions above. Commands can have any optional arguments,
but support either sub-commands or positional arguments. To get the selected
command, set a string variable with `set_command_var(cli, var)` that will
receive the command name.

```cpp
auto samples = 10; auto scene = ""s;                // state
auto command = string{};                            // selected command 
auto cli = make_cli("app", "testing cli");          // initialize cli
set_command_var(cli, command);                      // command variable
auto render = add_command(cli, "render", "render"); // command
add_option(render, "samples", samples, "samples");  // optional argument
add_argument(render, "scene", scene, "scene");      // positional argument
auto convert = add_command(cli, "convert", "convert"); // command
add_argument(convert, "scene", scene, "scene");     // positional argument
parse_cli(cli, make_cli_args(argc, argv));          // parse cli
if (command == "render") { ... }                    // execute commands
else if (command == "samples") { ... }
```

## Command-Line JSON Serialization

Command-line arguments can be loaded from JSON config files by specifying the 
`--config <filename>` option in the command line. The semantic of the 
config file is that configs are loaded first, and then overridden by 
the value specified in the terminal.

The JSON file should contain a value for each specified option, of the same
type as the option itself. For commands, use a JSON object to specify their
values. The command variable is called `command`.

Alternatively, you can specify an option that, if specified, indicates a
filename in a directory where the config might be found, using 
`add_option_with_config(cli, name, value, usage, config, alt)`.
Here `config` is the filename of the optional config file.

## Command-Line Positional Arguments And Flags

While for now positional arguments are supported, their use is deprecated.
The same is true for the use of short flags, like `-h`.
Th main reasons for this is to make the input to CLI tools more readable by 
using long names for options and avoiding unnamed positional arguments.
This is also a better match for the JSON config files and since internally
the CLI uses JSON as a data model to process values.
In future releases, positional arguments will be deprecated and eventually
removed. 