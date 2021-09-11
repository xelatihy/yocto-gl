#! /usr/bin/env python3 -B

import jsonschema, json, typer

def main(input: str, schema: str):
  with open(input) as f:
    json_input = json.load(f)
  with open(schema) as f:
    json_schema = json.load(f)
  validator = jsonschema.Draft7Validator(json_schema)
  for error in validator.iter_errors(json_input):
     print(error.message)

if __name__ == '__main__':
  typer.run(main)
