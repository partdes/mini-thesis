# structure_builder — Usage Guide

## Activate environment

```powershell
& C:\Users\Yasemin\environments\minitez\Scripts\Activate.ps1
```

## 1. Create a spec file

Copy the template below to `data/specs/<your-name>.yaml` and fill it in.
Use `/article-structure-extractor` to generate a spec automatically from article text.

```yaml
name: cubic_TiAlN_2x2x2_random
base_structure: data/structures/cubic_TiN.cif   # or mp-492
type: compositional                              # or: sqs
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.5
    Al: 0.5
  # nonmetal:                                   # omit if pure nitride
  #   N: 0.75
  #   C: 0.25
output:
  formats: [vasp, cif]
  directory: data/structures/
source:
  article: "Author et al., Journal, Year"
  description: "brief quote from article"
compositional_options:
  arrangement: random                            # random | ordered | explicit
  # ordered_type: layered                        # layered | l10 | l12 | b2 | checkerboard
  # layer_axis: 2                                # 0=a, 1=b, 2=c
```

For SQS, replace `compositional_options` with:

```yaml
sqs_options:
  iterations: 50000
  backend: auto     # auto | pymatgen | atat
```

## 2. Run from the terminal

```powershell
cd C:\Users\Yasemin\Desktop\mini-thesis
python -m structure_builder data/specs/cubic_TiAlN_2x2x2_random.yaml
```

Multiple specs at once:

```powershell
python -m structure_builder data/specs/*.yaml
```

## 3. Run from a Jupyter notebook

```python
from structure_builder.builder import build_from_spec
written = build_from_spec("data/specs/cubic_TiAlN_2x2x2_random.yaml")
for f in written: print(f)
```

## 4. Composition rules

- `metal` and `nonmetal` fractions must each sum to **1.0**
- Omit the `nonmetal` block entirely if the non-metal sublattice is unchanged from the base structure
- Species names must match pymatgen element symbols exactly (e.g. `Ti`, `Al`, `N`, `C`)

## 5. Ordering options

| `ordered_type` | Suitable for | Notes |
|---|---|---|
| `layered` | any N species | Alternates species layer by layer along `layer_axis` |
| `l10` | binary | Alternating (001) planes — same as layered for 2 species |
| `l12` | binary | Minority at face-center positions |
| `b2` | binary | Two interpenetrating BCC-like sublattices |
| `checkerboard` | binary | Alternates in the ab-plane |
| `explicit` | any | Specify every site by index in `site_occupancies` |

## 6. Using a Materials Project ID as base structure

Set `base_structure: mp-492` (for example). Requires `MP_API_KEY` in `0.settings/.env`.
