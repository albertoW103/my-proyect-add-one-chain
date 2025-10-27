# my-proyect-add-one-chain

The script generates a random set of protein-polymer configurations.
The polymer is attached to the N-ter or C-ter end of the protein.

---

## 1. How polymer configurations are generated

The script generates polymers configurations following a Rotational Isomeric State (RIS) model.
This model generate 100000 random configurations.

By default, polymer configurations are generated following the parameters:
- `lseg = 0.38`  (nm) , distance between segments for alpha carbons in proteins.
- `nrot = 100`        , numer of rotations that are included (following Euler angles).
- `cuantas = 100000`  , total number of configurations to generated. 100000/100 = 1000 generated configurations (the rest are rotations).


This parameters can be modified in the file `aux-main.f95` at `/random`.
The executable `polymer.x` is compiled each time the code is employed.


---

## 2. How the protein file must be

The script employed protein in the format xyz for molecular theory.

---

## 3. How the protein and polymer are merged

The polymer configuration is first superposed onto the protein by placing its first bead at the selected terminal (N-ter or C-ter).
Any configuration where a polymer bead is closer than 0.38 nm (this cutoff can be modified in`get_polymers.py`) to the protein is discarded to avoid steric clash.
All remaining configuration form the valid protein–polymer set (this set is store for comparison).
From this filtered set, a random subset (see `get_protein-polymers.py` to change the random seed) is selected to be used in simulations.

---

## 3. Scripts

These are the mains scripts:

| file                               | principal function                                            |
|------------------------------------|---------------------------------------------------------------|
| `get_polymers.py`                  | Generates all polymers that fit to a protein                  |
| `get_protein-polymer.py`           | Assembles the all protein-polymer complexes                   |
| `my_functions/my_functions.py`     | Auxiliary utilities                                           |

---

## 4. How to run

First, we generate all polymer sequences that fit to the protein:

`python3 get_polymers.py --input 1J05.xyz --xter Nter --seq GGSGGSGGS`

Then, we combine those polymers with the protein and generates a number of random conformations of the protein-polymer conjugates (plus the total)

`python3 get_protein-polymer.py --input 1J05.xyz --output 1J05 --xter Nter --nconfs '5 10 20' --seq GGSGGSGGS`












