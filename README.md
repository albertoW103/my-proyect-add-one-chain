# my-proyect-add-one-chain

Script to generate and assemble random polymers (peptides) onto a protein at the N-ter or C-ter end.

---

## 1. Dependencia de `random/` to generate configuration of polymers

The script generates polymers following a Rotational Isomeric State (RIS) model that generates 100,000 random configurations.
The subroutine is at `src/random/polymer.x`.
The default conditions are already set in `src/random/polymer.x`:

- `lseg = 0.38`  (nm) , distance between segments for alpha carbons in proteins
- `nrot = 100`        , rotations included
- `cuantas = 100000`  , total number of generated configurations. 100000/100 = 1000 generated configurations (the rest are rotations)

These parameters can be modified in the `aux-main.f95` script.

---

## 3. Dependencia de `random/` to generate configuration of polymers

The script employed protein in the format xyz for molecular theory.

---

## 3. Scripts principales

| file                            | principal function                                            |
|---------------------------------|---------------------------------------------------------------|
| `get_polymers.py`               | Generates all polymers that fit to a protein                  |
| `get_protein-polymer.py`        | Assembles the all protein-polymer complexes                   |
| `my_functions/my_functions.py`  | Auxiliary utilities                                           |

---

## 4. How to run

get all polymer sequences:

`python3 get_polymers.py.py -input input_filename.xyz -output output_filename.xyz -xter Nter sequence`

It combines all the polymers with the protein, and generates a number of random conformations of the protein-polymer conjugates (plus the total)

`python3 get_protein-polymers.py.py -input input_filename.xyz -output output_filename.xyz -xter Nter nconfs sequence`












