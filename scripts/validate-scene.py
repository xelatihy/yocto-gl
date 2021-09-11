#! /usr/bin/env python3 -B

import jsonschema, json, typer
from typing import List

def main(schema: str, inputs: List[str]):
  with open(schema) as f:
    json_schema = json.load(f)
  for input in inputs:
    print(f'--- {input}')
    with open(input) as f:
      json_input = json.load(f)
    validator = jsonschema.Draft7Validator(json_schema)
    for error in validator.iter_errors(json_input):
      print(error.message)

if __name__ == '__main__':
  typer.run(main)
