# my-proyect-add-one-chain

Pipeline para generar y ensamblar polímeros aleatorios (péptidos) sobre una proteína en el N-ter o C-ter, con filtrado geométrico por distancia mínima.

---

## 1. Dependencia de `random/`

El script necesita la subrutina `src/random/` para generar las configuraciones del polímero.  
En `src/random/` ya están seteadas las condiciones por defecto:

- `lseg = 0.38`  (nm)
- `nrot = 1`
- `cuantas = 1000000`  (configuraciones generadas)

Estas se pueden modificar directamente dentro de esa carpeta.

---

## 2. Scripts principales

| Archivo                         | Función principal                                  |
|----------------------------------|----------------------------------------------------|
| `get_polymers.py`               | Genera polímeros a partir de una secuencia         |
| `get_protein-polymer.py`        | Lee proteína, filtra polímeros y arma el complejo  |
| `my_function/`                 | Utilidades auxiliares                              |
| `run.sh`                        | Script de ejemplo para ejecución completa          |

---

## 3. Scripts de gráficas (opcional)

En `plots/` se incluyen herramientas para analizar/visualizar resultados:

- `plot_histo-end2end_confs.py`
- `temp_origin_for-peptides.mpstyle`

---

## 4. Inputs (definidos en `run.sh`)

Todos los parámetros de entrada (proteína, secuencia, terminal, número de conformaciones) se ajustan dentro de `src/run.sh`.

---

## 5. Cómo correr

### Opción A — usando `run.sh`

```bash
cd src
bash run.sh






