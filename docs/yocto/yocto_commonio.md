# Yocto/CommonIO: Common IO functionality

Yocto/CommonIO is a collection of utilities used in writing IO functionality,
including file IO, Json IO, and path manipulation.
Yocto/CommonIO is implemented in `yocto_commonio.h` and `yocto_commonio.cpp`, 
and depends on `json.hpp` for Json serialization and number printing, and
`fast_float.h` for number parsing.

## Errors handling in IO functions

IO functions in Yocto/GL have a dual interface to support their use with and 
without exceptions. IO functions with exception are written as 
`load_<type>(filename, data, <options>)` and 
`save_<type>(filename, data, <options>)` where `<type>` is the data type 
read or written, and `<options>` is an optional list of IO options. 
A convenience shortcut is provided for all loading functions that returns the 
data directly, written as `data = load_<type>(filename, <options>)`.
Upon errors, an `io_error` is thrown from all IO functions.
This makes the error-handling code more uniform across the library. 
`io_error` has fields for `filename` and `message` that can retrieved directly, 
and a message for users that can be retrieved with the `.what()` method.

```cpp
auto text = string{};
try {
  load_text("input_file.txt",  text);   // load text
  text = load_text("input_file.txt");   // alternative load
  save_text("output_file.txt", text);   // save text
} catch(const io_error& error) {
  print_info(error.what()); exit(1);    // handle error
}
```

IO functions without exceptions are written as
`load_<type>(filename, data, error, <options>)` and 
`save_<type>(filename, data, error, <options>)` where `error` is a string 
that contains a message if an error occurred. These functions return a boolean 
flag that indicates whether the operation succeeded.  

```cpp
auto text = string{};
auto error = string{};
if(!load_text("input_file.txt",  text))      // load text
  print_info(error); exit(1);                // handle error
if(!save_text("output_file.txt", text))      // save text
  print_info(error); exit(1);                // handle error
```

## Text and binary serialization

Yocto/CommonIO supports reading and writing whole files of either binary
data or text. Use `load_text(filename, text)` or `text = load_text(filename)` 
to load text, and `save_text(filename, text)` to save it.
Use `load_binary(filename, data)` or `data = load_binary(filename)` 
to load text, and `save_binary(filename, data)` to save it. 
Text is stored as a string and binary data is stored as an array of bytes.
Use without exception is supported as described above.

```cpp
auto text = string{};
load_text("input_file.txt",  text);          // load text
save_text("output_file.txt", text);          // save text
auto text1 = load_text("input_file.txt");    // alternative load
auto data = vector<byte>{};
load_binary("input_file.bin",  data);        // load data
save_binary("output_file.bin", data);        // save data
auto data1 = load_binary("input_file.bin");  // alternative load
```

## Path manipulation utilities

Yocto/CommonIO contains several helper function to manipulate paths. Paths
are encoded in UTF8 across the library and these functions make it easier to
handle UTF8-encoded paths across operating systems, by wrapping 
`std::filesystem` with a string interface.

Use `path_dirname(filename)`, `path_extension(filename)`,  
`path_filename(filename)`, `path_basename(fillename)`
to extract the directory, extension, filename and basename from a path.
Use `path_join(patha,pathb)` to joins paths and
`replace_extension(filename,ext)` to replace a path extension.
Use `path_exists(filename)` to check if a path exists and
`path_isdir(filename)` and `path_isfile(filename)` to check whether
it is a directory ot file respectively.
Use `list_directory(dirname)` to list directory contents, 
`make_directory(dirname)` to create a directory, and
`path_current()` to get the current directory.

## Json data type

Throughout the library, JSON data is used extensively as a simple and flexible
IO format. Most JSON libraries we reviewed use header-only formulation that
slow down compile time significantly. For this reason, we provide a simple
JSON data type modeled after `nlohmann::json`. Here we review quickly its
functionality, for those familiar with the JSON format.

We model JSON values with the `json_value` type. Most functionality in this
type mimic the C++ `vector` and `map` functions. In fact, if you know how 
to use those type and you know how to handle JSON data, then it should be 
easy to use `json_value`. Json

Json values are unions of either null, strings, booleans, numbers, arrays 
of Json values and dictionaries of string to Json values. To avoid loosing
precision, we store numbers with 64 bits per value and as either signed
integers, unsigned integers, or floating point. Json arrays are stored as
`vector<json_value>`, while Json dictionaries are stored as ordered maps,
which are roughly equivalent to `vector<pair<string, json_value>>`. 
The choice off dictionary storage handle well small dictionaries, which is 
mostly the case in Json.

The type of a `json_value` is queried by the `json.type()` method that returns 
a `json_type` value. More conveniently, you can query thee type with 
`json.is_null()`, `json.is_integer()`, `json.is_number()`, `json.is_boolean()`, 
`json.is_string()`,`json.is_array()`, and `json.is_object()`. 
`json.is_integer()` checks if a number is signed or unsigned integer, 
while `json.is_number()` returns true for any number type.

You can convert values to a `json_value` by either constructing it,
assigning to it and calling the `json.set(value)` method with values of all 
builtin types, arrays of builtin types and enums, if they define the trait 
`json_enum_trait`.

```cpp
auto json = json_value{5};        // construct a value of type integer
json = "text";                    // converts to a value of type string
json.set(vector<float>{...});     // converts to a value of type array
```

You can convert a `json_value` to values by either casting to the type of the
value or calling the `json.get(value)` and `value = json.get<type>()` methods. 
A `json_error` is thrown if the conversion cannot be done.

```cpp
auto json = json_value{"hello"};  // construct a json_value
auto text = (string)json;         // convert to string
json.get(text);                   // convert to string
```

Besides conversions, Json arrays and dictionaries can be manipulated directly.
Json arrays are constructed by assigning a `json_array` to a `json_value`, 
either constructed directly or via the `json_value::array()` and 
`json_value::array(size)`. Array sizes are queries via the `json.empty()` and 
`json.size()` methods. Arrays can be resized with `json.resize(size)`, or by
appending values at the array end with `json.push_back(value)` and 
`json.emplace_back(args...)`. Array items are accessed by index with the 
`[index]` operator and the `json.at(index)` method. Arrays can be iterated
by using the `json.begin()` and `json.end()` methods.

```cpp
auto json = json_value::array();  // construct a value containing an array
json.push_back(json_value{...});  // append a value
json.push_back(5);                // append a value converted from an int
auto size = json.size();          // get size
auto converted = (int)json.at(5); // access value by index
for(auto& item : json) { ... }    // iterate over the array
```

Json dictionaries are constructed by assigning a `json_object` to a `json_value`, 
either constructed directly or via the `json_value::object()`. Object sizes 
are queries via the `json.empty()` and `json.size()` methods. Values are
inserted in dictionaries by using the `[key]` operator, or using the 
`insert_back(key, value)` method for a faster insertion that does not check
for duplicate keys. Object items are accessed by key with the `[key]` operator 
and the `json.at(key)` method. The existence of a key is checked with the 
`json.contains(key)` method. Json dictionaries are iterated with the 
`json.items().begin()` and `json.items().end()`.
Missing keys are so common in Json that Json dictionaries support the 
`json.try_get(key, value)` that assigns to the value only if the key is present,
and `json.get_or(key, default)` that returns the value if the key is present or
the default value specified.

```cpp
auto json = json_value::object();     // construct a value containing a dict.
json["key"] = json_value{...};        // insert a value
json["key"] = 5;                      // insert a value converted from an int
auto size = json.size();              // get size
auto converted = (int)json.at("key"); // access value by key
for(auto& [key, item] : json.items) { ... }  // iterate over the dictionary
auto value = json.get_or("key", 5);   // get value or default
json.try_get("key", value);           // get value if present
```

## Json IO

Use `load_json(filename, json)` or `json = load_json(filename)` to load json 
data from a file, and `save_json(filename, json)` to save json data to a file. 
Use without exception is supported as described above.

```cpp
auto json = json_value{};
load_json("input_file.json",  json);       // load json
save_json("output_file.json", json);       // save json
auto json1 = load_json("input_file.txt");  // alternative load
```

Use `parse_json(text, json)`  or `json = parse_json(text)` to parse json data 
from a string, and `format_json(text, json)` or `text = format_json(json)` to 
write json data to a string.
Use without exception is supported as described above.

```cpp
auto input_string  = string{...};
auto json = json_value{};
parse_json(input_string, json);            // parse json
auto output_string = string{};
format_json(output_string, json);          // format json
auto json1 = parse_json(input_string);     // alternative parse
auto output_string1 = format_json(json1);  // alternative format
```

## Ordered map

Json dictionaries are stored in an `ordered_map` structure. Ordered maps are
unsorted just arrays of key-value pairs, with methods added to match the 
map data structures in the standard library, that we do not report here since
it would be a repetition.

We use unsorted arrays for dictionaries since dictionary implementations
tradeoff memory consumption, number of allocation, and insertion and
query speeds. For th specific use of Json values in Yocto/GL, a good trade off
is to reduce memory pressure, while optimizing the lookups only for small 
dictionaries, as is generally the case for the Json values in Yocto/GL.
One exception are large dictionaries that need to be only iterated. 
In that case, we can skip the key checks and insert values at the end directly.
