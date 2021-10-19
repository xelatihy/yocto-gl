# Yocto/Cli: Utilities for writing command-line apps

Yocto/Cli is a collection of utilities used in writing command-line
applications, including parsing command line arguments, printing values,
and timers.
Yocto/Cli is implemented in `yocto_cli.h`.

## Printing values

Use `print_info(message)` to print a message, and `print_error(message)`
to print an error. To time a block of code create a `simple_timer` and
get the elapsed time with `elapsed_format(timer)`.

```cpp
print_info("Message");          // print message
print_fatal("Error and exit");  // print error and exit
auto timer = simple_timer{};    // start timer
print_info("elapsed: {}",       // print elapsed time
  elapsed_formatted(timer));
```

## Command-Line Parsing

Yocto/Cli includes a simple command-line parser that supports optional
arguments with automatic help generation and error checking. This parser
is meant to be simple and force a command-line style where an app has several
commands, all of which are options with default values.

The command-line parser is initialized with `make_cli(name, help)` that takes
a program name and a help string. Then add command-line options with
`add_option(cli, name, value, help, check)` for optional arguments.
In these commands, `name` is the option name, `value` is a reference to the
variable that we want to set, `help` is a usage message, and `check`
is either a range of valid values er a list of valid choices.
A help flag is automatically added.

The command line is parsed for options named `--` followed by the option `name`.
The type of each option is determined by the passed reference.
The parser supports integers, floating point numbers, strings, boolean flags,
enums, arrays and vectors.

After adding all arguments, use `ok = parse_cli(cli, args, error)` to parse the
command-line. If an error occurs, the parser throws `false` and stores in
`error` a message that can be printed to inform the user. You can also use 
`parse_cli(cli, args)` in which case a `cli_error` exception is thrown if
errors occur.

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
auto error = string{};                          // error string
if (!parse_cli(cli, args, error))               // parse args
  print_error(error);                           // print error
try {
  parse_cli(cli, args);                         // parse args
} catch(const cli_error& error) {
  print_error(error);                           // print error
}
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
