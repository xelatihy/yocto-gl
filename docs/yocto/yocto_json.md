# Yocto/Json: Json functionality (deprecated)

Yocto/Json is an implementation of a Json type and related IO functionality.
Yocto/Json is implemented in `yocto_json.h` and `yocto_json.cpp`, 
and depends on `json.hpp` for Json serialization.

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

## Json serialization

Use `load_json(filename, json)` or `json = load_json(filename)` to load json 
data from a file, and `save_json(filename, json)` to save json data to a file. 
Upon errors, an `json_error` is thrown from all IO functions.

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
