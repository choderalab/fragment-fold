from plip_analysis import PLIntReport
from pathlib import Path
import argparse
import yaml
from pebble import ProcessPool
from functools import partial

def parse_args():
    parser = argparse.ArgumentParser(description="Get PLIP interactions")
    parser.add_argument("--yaml_input", type=Path, help="Path to input yaml file containing name: path pairs")
    parser.add_argument("--output-dir", type=Path, help="Path to output directory")
    parser.add_argument("--ncpus", type=int, default=1, help="Number of cpus to use for parallel processing")
    return parser.parse_args()


def analyze_structure(structure: Path, name: str, output_dir: Path):
    outpath = output_dir / f"{name}_{structure.stem}_interactions.csv"
    interactions = PLIntReport.from_complex_path(complex_path=structure, create_pymol_session=True, pymol_session_path=outpath.with_suffix(".pse"))
    interactions.to_csv(outpath)
    return outpath

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

        analyze_structure_partial = partial(analyze_structure, name=name, output_dir=output_dir)

        # parallelize with pebble
        with ProcessPool(max_workers=args.ncpus) as pool:
            result = pool.map(analyze_structure_partial, structures)
        print(result.result())

if __name__ == "__main__":
    main()