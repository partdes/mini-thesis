import sys
from pathlib import Path

from .builder import build_from_spec


def main():
    if len(sys.argv) < 2:
        print("Usage: python -m structure_builder <spec.yaml> [<spec2.yaml> ...]")
        sys.exit(1)

    for spec_path in sys.argv[1:]:
        path = Path(spec_path)
        print(f"Building from {path.name}...")
        written = build_from_spec(path)
        for f in written:
            print(f"  -> {f}")


if __name__ == "__main__":
    main()
