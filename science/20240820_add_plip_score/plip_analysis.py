from pydantic import BaseModel
from enum import Enum, auto
from plip.structure.preparation import PDBComplex, PLInteraction
import pandas as pd
from pathlib import Path
import json

class InteractionType(Enum):
    """
    Protein-Ligand Interaction Types. Names refer to the role of the Ligand
    """
    HydrogenBondDonor = 'HydrogenBondDonor'
    HydrogenBondAcceptor = 'HydrogenBondAcceptor'
    HydrophobicInteraction = 'HydrophobicInteraction'
    PiStacking = 'PiStacking'
    HalogenBond = 'HalogenBond'
    SaltBridge = 'SaltBridge'

    def __str__(self):
        return self.value

class AtomType(str, Enum):
    NONE = 'none'
    H = 'H'
    C = 'C'
    N = 'N'
    O = 'O'
    S = 'S'
    P = 'P'
    F = 'F'
    Cl = 'Cl'
    Br = 'Br'

    def __str__(self):
        return self.value

class ProteinResidueType(str, Enum):
    ALA = 'ALA'
    ARG = 'ARG'
    ASN = 'ASN'
    ASP = 'ASP'
    CYS = 'CYS'
    GLN = 'GLN'
    GLU = 'GLU'
    GLY = 'GLY'
    HIS = 'HIS'
    ILE = 'ILE'
    LEU = 'LEU'
    LYS = 'LYS'
    MET = 'MET'
    PHE = 'PHE'
    PRO = 'PRO'
    SER = 'SER'
    THR = 'THR'
    TRP = 'TRP'
    TYR = 'TYR'
    VAL = 'VAL'

    def __str__(self):
        return self.value

class FormalCharge(str, Enum):
    Positive = 'positive'
    Negative = 'negative'
    Neutral = 'neutral'

    def __str__(self):
        return self.name

class ProteinLigandInteraction(BaseModel):
    interaction_type: InteractionType
    protein_residue_number: int
    protein_residue_type: ProteinResidueType
    to_sidechain: bool
    count: int = 1 # Default to 1
    ligand_atom_type: AtomType = 'none'
    protein_atom_type: AtomType = 'none'
    ligand_charge: FormalCharge = 'neutral'
    protein_charge: FormalCharge = 'neutral'

    def to_json_dict(self):
        return {k: v.value if k in [self.interaction_type, self.protein_residue_type] else v for k, v in self.dict().items()}



def plip_constructor(plip: PLInteraction) -> list[ProteinLigandInteraction]:
    """
    Function to get a tidy dataframe of interactions from a plip fingerprint
    :param plip:
    :return:
    """

    interactions = []
    for interaction in plip.all_itypes:
        interaction_type = type(interaction).__name__

        protein_atom_type = AtomType.NONE
        ligand_atom_type = AtomType.NONE
        ligand_charge = FormalCharge.Neutral
        protein_charge = FormalCharge.Neutral

        if interaction_type == 'hbond':
            if interaction.protisdon:
                output_interaction_type = InteractionType.HydrogenBondAcceptor
                protein_atom_type = AtomType(interaction.dtype[0])
                ligand_atom_type = AtomType(interaction.atype[0])

            else:
                output_interaction_type = InteractionType.HydrogenBondDonor
                protein_atom_type = AtomType(interaction.atype[0])
                ligand_atom_type = AtomType(interaction.dtype[0])

        elif interaction_type == 'halogenbond':
            output_interaction_type = InteractionType.HalogenBond
            protein_atom_type = AtomType(interaction.acctype[0])
            ligand_atom_type = AtomType(interaction.donortype[0])

        elif interaction_type == 'saltbridge':
            output_interaction_type = InteractionType.SaltBridge
            if interaction.protispos:
                ligand_charge = FormalCharge.Negative
                protein_charge = FormalCharge.Positive
            else:
                ligand_charge = FormalCharge.Positive
                protein_charge = FormalCharge.Negative

        elif interaction_type == 'hydroph_interaction':
            output_interaction_type = InteractionType.HydrophobicInteraction

        elif interaction_type == 'pication':
            output_interaction_type = InteractionType.PiStacking
            protein_charge = FormalCharge.Positive if interaction.protcharged else FormalCharge.Negative

        elif interaction_type == 'pistack':
            output_interaction_type = InteractionType.PiStacking

        else:
            raise NotImplementedError(f"Interaction type not implemented: {interaction_type} for interaction: {interaction}")

        protein_residue_number = interaction.resnr
        protein_residue_type = interaction.restype

        # These interactions must be to the sidechain, and therefore PLIP does not provide the sidechain boolean
        if interaction_type in ['hydroph_interaction', 'pistack', 'saltbridge', 'pication']:
            to_sidechain = True
        else:
            to_sidechain = interaction.sidechain

        interactions.append(ProteinLigandInteraction(interaction_type=output_interaction_type,
                                                     protein_residue_number=protein_residue_number,
                                                     protein_residue_type=protein_residue_type,
                                                     to_sidechain=to_sidechain,
                                                     ligand_atom_type=ligand_atom_type,
                                                     protein_atom_type=protein_atom_type,
                                                     ligand_charge=ligand_charge,
                                                     protein_charge=protein_charge
                                                     ))
    # check if there are any interactions that are the same

    return interactions


def collect_duplicates(plints: list[ProteinLigandInteraction]):
    """
    This function takes a list of ProteinLigandInteractions and
    returns a list of unique interactions with the count updated if any duplicates are found.
    :param plints:
    :return:
    """
    obj_dict = {}

    for obj in plints:
        if obj.json() in obj_dict:
            obj_dict[obj.json()].count += 1
        else:
            obj_dict[obj.json()] = obj

    return list(obj_dict.values())


class PLIntReport(BaseModel):
    """
    Class to store a report of PLIP interactions
    """
    structure: str
    interactions: list[ProteinLigandInteraction]

    @classmethod
    def from_complex_path(cls, complex_path: str | Path, ligand_id="UNK") -> "PLIntReport":
        my_mol = PDBComplex()
        my_mol.load_pdb(str(complex_path))
        my_mol.analyze()

        binding_site = None

        for k, v in my_mol.interaction_sets.items():
            if ligand_id in k:
                binding_site = k
                break

        # Convert the PLIP interaction analysis of the binding site to a list of ProteinLigandInteractions and then
        # Collect any duplicate interactions
        interactions = collect_duplicates(plip_constructor(my_mol.interaction_sets[binding_site]))
        return cls(structure=str(complex_path), interactions=interactions)

    def to_csv(self, path: str | Path):
        df = pd.DataFrame.from_records([json.loads(interaction.json()) for interaction in self.interactions])
        df.to_csv(path, index=False)

    @classmethod
    def from_csv(cls, path: str | Path):
        df = pd.read_csv(path)
        interactions = [ProteinLigandInteraction(**row) for _, row in df.iterrows()]
        return PLIntReport(structure=path.stem, interactions=interactions)


class FingerprintLevel(Enum):
    """
    Enum to specify the level of detail for the fingerprint
    """
    ByInteractionType = 'ByInteractionType'
    ByInteractionTypeAndResidueType = 'ByInteractionTypeAndResidueType'
    ByInteractionTypeAndAtomTypes = 'ByInteractionTypeAndAtomTypes'
    ByInteractionTypeAndResidueTypeAndBBorSC = 'ByInteractionTypeAndResidueTypeAndBBorSC'
    ByInteractionTypeAndResidueTypeAndNumber = 'ByInteractionTypeAndResidueTypeAndNumber'
    ByEverything = 'ByEverything'

    def __str__(self):
        return self.value


def calculate_fingerprint(plint_report: PLIntReport, level: FingerprintLevel) -> dict:
    """
    Calculate a fingerprint of the interactions in a PLIntReport
    :param plint_report:
    :param level:
    :return: dict
    """
    fingerprint_dict = {}
    for interaction in plint_report.interactions:
        if level == FingerprintLevel.ByInteractionType:
            key = interaction.interaction_type.value
        elif level == FingerprintLevel.ByInteractionTypeAndResidueType:
            key = f"{interaction.interaction_type.value}_{interaction.protein_residue_type}"
        elif level == FingerprintLevel.ByInteractionTypeAndAtomTypes:
            key = f"{interaction.interaction_type.value}_Protein_{interaction.protein_atom_type}_Ligand_{interaction.ligand_atom_type}"
        elif level == FingerprintLevel.ByInteractionTypeAndResidueTypeAndBBorSC:
            key = f"{interaction.interaction_type.value}_{interaction.protein_residue_type}_{'SC' if interaction.to_sidechain else 'BB'}"
        elif level == FingerprintLevel.ByInteractionTypeAndResidueTypeAndNumber:
            key = f"{interaction.interaction_type.value}_{interaction.protein_residue_type}{interaction.protein_residue_number}"
        elif level == FingerprintLevel.ByEverything:
            key = "_".join([f"{k}_{str(v)}" for k, v in interaction.dict().items()])
        else:
            raise ValueError("Invalid Fingerprint Level")
        original_count = fingerprint_dict.get(key, 0)
        fingerprint_dict[key] = original_count + interaction.count
    return fingerprint_dict


class SimilarityScore(BaseModel):
    """
    Class to store a similarity score between two fingerprints with some record of how it was calculated
    """
    provenance: str
    score: float
    number_of_interactions_in_reference: int
    number_of_interactions_in_query: int
    number_of_interactions_in_intersection: int
    number_of_interactions_in_union: int


def calculate_tversky(fingerprint1: dict, fingerprint2: dict, alpha: float = 1, beta: float = 0) -> SimilarityScore:
    """
    Calculate the Tversky Index between two fingerprints.
    To calculate the Tanimoto Coefficient, set alpha=beta=1
    To calculate the Dice Coefficient, set alpha=beta=0.5

    :param fingerprint1:
    :param fingerprint2:
    :param alpha:
    :param beta:
    :return:
    """

    # First get the union of the keys
    fp_types = set(fingerprint1.keys()).union(set(fingerprint2.keys()))

    # Calculate the number of interactions that are in both fingerprints
    matched = sum([min(fingerprint1.get(a, 0), fingerprint2.get(a, 0)) for a in fp_types])

    # Calculate the number of interactions in the union
    union = sum(fingerprint1.values()) + sum(fingerprint2.values()) - 2 * matched

    # Calculate the Tversky Index
    score = matched / (matched + alpha * (sum(fingerprint1.values()) - matched) + beta * (
                sum(fingerprint2.values()) - matched))
    return SimilarityScore(provenance=f"Tversky_alpha{alpha:.2f}_beta{beta:.2f}", score=score,
                           number_of_interactions_in_reference=sum(fingerprint1.values()),
                           number_of_interactions_in_query=sum(fingerprint2.values()),
                           number_of_interactions_in_intersection=matched, number_of_interactions_in_union=union, )


class InteractionScore(BaseModel):
    """
    Class to store the results of comparing two PLIntReports.
    Automatically calculates the Tanimoto Coefficient (alpha=beta=1) and the Tversky Index (alpha=1, beta=0)
    """
    provenance: str
    number_of_interactions_in_query: int
    number_of_interactions_in_reference: int
    number_of_interactions_in_intersection: int
    number_of_interactions_in_union: int
    tanimoto_coefficient: float
    tversky_index: float
    reference_fingerprint: dict
    query_fingerprint: dict

    @classmethod
    def from_fingerprints(cls, reference: PLIntReport, query: PLIntReport,
                          level: FingerprintLevel) -> "InteractionScore":

        # Calculate the fingerprints according to the specified level of detail
        reference = calculate_fingerprint(reference, level)
        query = calculate_fingerprint(query, level)

        # Calculate the Tanimoto Coefficient and Tversky Index
        tanimoto = calculate_tversky(reference, query, alpha=1, beta=1)
        tversky_index = calculate_tversky(reference, query)

        return cls(provenance=level.value, number_of_interactions_in_query=tanimoto.number_of_interactions_in_query,
                   number_of_interactions_in_reference=tanimoto.number_of_interactions_in_reference,
                   number_of_interactions_in_intersection=tanimoto.number_of_interactions_in_intersection,
                   number_of_interactions_in_union=tanimoto.number_of_interactions_in_union,
                   tanimoto_coefficient=tanimoto.score,
                   tversky_index=tversky_index.score,
                   reference_fingerprint=reference,
                   query_fingerprint=query)