# my-proyect-add-one-chain

Script to generate and assemble random polymers (peptides) onto a protein at the N-ter or C-ter end.

---

## 1. Dependencia de `random/`

The script generates polymers following a RIS model that generates 100,000 random configurations.
The subroutine is at `src/random/polymer.x`.
The default conditions are already set in `src/random/polymer.x`:

- `lseg = 0.38`  (nm) , distance between segments
- `nrot = 100`        , rotations included
- `cuantas = 100000`  , Total number of generated configurations. 100000/100 = 1000 generated configurations

These parameters can be modified in the aux-main.f95 script.

---

## 2. Scripts principales

| file                            | principal function                                            |
|---------------------------------|---------------------------------------------------------------|
| `get_polymers.py`               | Generates all polymers that fit to a protein from a sequence  |
| `get_protein-polymer.py`        | Assembles the all protein-polymer complex                     |
| `my_function/`                  | Auxiliary utilities                                           |

---

## 5. Cómo correr

get all polymer sequences:

`python3 get_polymers.py.py -input input_filename.xyz -output output_filename.xyz -xter Nter sequence`

It combines all the polymers with the protein, and generates a number of random conformations of the protein-polymer conjugates (plus the total)

`python3 get_protein-polymers.py.py -input input_filename.xyz -output output_filename.xyz -xter Nter nconfs sequence`

