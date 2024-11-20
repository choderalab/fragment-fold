from plip_analysis import PLIntReport
from pathlib import Path
import argparse
import yaml

def parse_args():
    parser = argparse.ArgumentParser(description="Get PLIP interactions")
    parser.add_argument("--yaml_input", type=Path, help="Path to input yaml file containing name: path pairs")
    parser.add_argument("--output-dir", type=Path, help="Path to output directory")
    return parser.parse_args()


def main():
    args = parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True)

    with open(args.yaml_input, "r") as f:
        input_dict = yaml.safe_load(f)

    for name, structure_dir in input_dict.items():
        structure_dir = Path(structure_dir)
        if not structure_dir.exists():
            raise FileNotFoundError(f"{structure_dir} does not exist")

        print(f"Loading all pdb structures in {structure_dir}")
        structures = [structure for structure in structure_dir.glob("*.pdb")]

        print(f"Analyzing {len(structures)} structures")

        for structure in structures:

            interactions = PLIntReport.from_complex_path(complex_path=structure)

            interactions.to_csv(output_dir / f"{name}_{structure.stem}_interactions.csv")


if __name__ == "__main__":
    main()