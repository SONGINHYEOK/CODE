""" Module containing all types and type imports
"""
# pylint: disable=unused-import
from typing import (
    Any,
    Callable,  # noqa
    Dict,
    List,  # noqa
    Iterable,  # noqa
    Optional,
    Sequence,  # noqa
    Set,  # noqa
    Tuple,
    TypeVar,  # noqa
    Union,
)
import abc
from PIL.Image import Image
from rdkit import Chem
from rdkit import Chem, DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.Chem import AllChem, Descriptors
import io
import sys
from PIL import Image, ImageDraw
from rdkit.Chem import Draw
import hashlib
import pandas as pd
import numpy as np

import operator
import json
import abc
import hashlib
from typing import TYPE_CHECKING

import networkx as nx
from networkx.algorithms.traversal.depth_first_search import dfs_tree

StrDict = Dict[str, Any]
RdMol = Chem.rdchem.Mol
RdReaction = Chem.rdChemReactions.ChemicalReaction
BitVector = ExplicitBitVect
PilImage = Image
PilColor = Union[str, Tuple[int, int, int]]
FrameColors = Optional[Dict[bool, PilColor]]

class MoleculeException(Exception):
    """An exception that is raised by molecule class"""
    
class Molecule:
    """
    A base class for molecules. Encapsulate an RDKit mol object and
    functions that can be applied to such a molecule.

    The objects of this class is hashable by the inchi key and hence
    comparable with the equality operator.

    :ivar rd_mol: the RDkit mol object that is encapsulated
    :ivar smiles: the SMILES representation of the molecule

    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    def __init__(
        self, rd_mol: RdMol = None, smiles: str = None, sanitize: bool = False
    ) -> None:
        if not rd_mol and not smiles:
            raise MoleculeException(
                "Need to provide either a rdkit Mol object or smiles to create a molecule"
            )

        if rd_mol:
            self.rd_mol = rd_mol
            self.smiles = Chem.MolToSmiles(rd_mol)
        else:
            self.smiles = smiles
            self.rd_mol = Chem.MolFromSmiles(smiles, sanitize=False)

        self._inchi_key: Optional[str] = None
        self._inchi: Optional[str] = None
        self._fingerprints: Dict[Union[Tuple[int, int], Tuple[int]], np.ndarray] = {}
        self._is_sanitized: bool = False

        # Atom mapping -> atom index dictionary
        self._atom_mappings: Dict[int, int] = {}
        # Atom index -> atom mapping dictionary
        self._reverse_atom_mappings: Dict[int, int] = {}

        if sanitize:
            self.sanitize()

    def __hash__(self) -> int:
        return hash(self.inchi_key)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Molecule):
            return False
        return self.inchi_key == other.inchi_key

    def __len__(self) -> int:
        return self.rd_mol.GetNumAtoms()

    def __str__(self) -> str:
        return self.smiles

    @property
    def inchi(self) -> str:
        """
        The inchi representation of the molecule
        Created by lazy evaluation. Will cause the molecule to be sanitized.

        :return: the inchi
        """
        if not self._inchi:
            self.sanitize(raise_exception=False)
            self._inchi = Chem.MolToInchi(self.rd_mol)
            if self._inchi is None:
                raise MoleculeException("Could not make InChI")
        return self._inchi

    @property
    def inchi_key(self) -> str:
        """
        The inchi key representation of the molecule
        Created by lazy evaluation. Will cause the molecule to be sanitized.

        :return: the inchi key
        """
        if not self._inchi_key:
            self.sanitize(raise_exception=False)
            self._inchi_key = Chem.MolToInchiKey(self.rd_mol)
            if self._inchi_key is None:
                raise MoleculeException("Could not make InChI key")
        return self._inchi_key

    @property
    def index_to_mapping(self) -> Dict[int, int]:
        """Return a dictionary mapping to atom indices to atom mappings"""
        if not self._reverse_atom_mappings:
            self._reverse_atom_mappings = {
                index: mapping for mapping, index in self.mapping_to_index.items()
            }
        return self._reverse_atom_mappings

    @property
    def mapping_to_index(self) -> Dict[int, int]:
        """Return a dictionary mapping to atom mappings to atom indices"""
        if not self._atom_mappings:
            self._atom_mappings = {
                atom.GetAtomMapNum(): atom.GetIdx()
                for atom in self.rd_mol.GetAtoms()
                if atom.GetAtomMapNum()
            }
        return self._atom_mappings

    @property
    def weight(self) -> float:
        """Return the exact molecular weight of the molecule"""
        self.sanitize(raise_exception=False)
        return Descriptors.ExactMolWt(self.rd_mol)

    def basic_compare(self, other: "Molecule") -> bool:
        """
        Compare this molecule to another but only to
        the basic part of the inchi key, thereby ignoring stereochemistry etc

        :param other: the molecule to compare to
        :return: True if chemical formula and connectivity is the same
        """
        return self.inchi_key[:14] == other.inchi_key[:14]

    def fingerprint(self, radius: int, nbits: int = 2048) -> np.ndarray:
        """
        Returns the Morgan fingerprint of the molecule

        :param radius: the radius of the fingerprint
        :param nbits: the length of the fingerprint
        :return: the fingerprint
        """
        key = radius, nbits

        if key not in self._fingerprints:
            self.sanitize()
            bitvect = AllChem.GetMorganFingerprintAsBitVect(self.rd_mol, *key)
            array = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(bitvect, array)
            self._fingerprints[key] = array

        return self._fingerprints[key]

    def has_atom_mapping(self) -> bool:
        """
        Determines if a the molecule has atom mappings

        :return: True if at least one atom has a mapping
        """
        for atom in self.rd_mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                return True
        return False

    def make_unique(self) -> "UniqueMolecule":
        """
        Returns an instance of the UniqueMolecule class that
        is representing the same molecule but is not hashable or comparable.

        :return: the unique molecule
        """
        return UniqueMolecule(rd_mol=self.rd_mol)

    def remove_atom_mapping(self, exceptions: Sequence[int] = None) -> None:
        """
        Remove all mappings of the atoms and update the smiles

        :param exceptions: keep the listed atom mappings
        """
        exceptions = exceptions or []
        for atom in self.rd_mol.GetAtoms():
            if exceptions and atom.GetAtomMapNum() in exceptions:
                continue
            atom.SetAtomMapNum(0)
        self.smiles = Chem.MolToSmiles(self.rd_mol)
        self._clear_cache()

    def sanitize(self, raise_exception: bool = True) -> None:
        """
        Sanitizes the molecule if it has not been done before.

        :param raise_exception: if True will raise exception on failed sanitation
        :raises MoleculeException: if the molecule could not be sanitized
        """
        if self._is_sanitized:
            return

        try:
            AllChem.SanitizeMol(self.rd_mol)
        # pylint: disable=bare-except
        except:  # noqa, there could be many reasons why the molecule cannot be sanitized
            if raise_exception:
                raise MoleculeException(f"Unable to sanitize molecule ({self.smiles})")
            self.rd_mol = Chem.MolFromSmiles(self.smiles, sanitize=False)
            return

        self.smiles = Chem.MolToSmiles(self.rd_mol)
        self._clear_cache()
        self._is_sanitized = True

    def _clear_cache(self):
        self._inchi = None
        self._inchi_key = None
        self._fingerprints = {}
        self._atom_mappings = {}
        self._reverse_atom_mappings = {}

class UniqueMolecule(Molecule):
    """
    A special molecule with the hash set to the id of the object.
    Therefore no two instances of this class will be comparable.

    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    def __init__(
        self, rd_mol: RdMol = None, smiles: str = None, sanitize: bool = False
    ) -> None:
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)

    def __hash__(self) -> int:
        return id(self)

    def __eq__(self, _) -> bool:
        return False

def none_molecule() -> UniqueMolecule:
    """Return an empty molecule"""
    return UniqueMolecule(rd_mol=Chem.MolFromSmiles(""))

class _ReactionInterfaceMixin:
    """
    Mixin class to define a common interface for all reaction class

    The methods `_products_getter` and `_reactants_getter` needs to be implemented by subclasses
    """

    def fingerprint(self, radius: int, nbits: int = None) -> np.ndarray:
        """
        Returns a difference fingerprint

        :param radius: the radius of the fingerprint
        :param nbits: the length of the fingerprint. If not given it will use RDKit default, defaults to None
        :return: the fingerprint
        """
        product_fp = sum(
            mol.fingerprint(radius, nbits) for mol in self._products_getter()  # type: ignore
        )
        reactants_fp = sum(
            mol.fingerprint(radius, nbits) for mol in self._reactants_getter()  # type: ignore
        )
        return reactants_fp - product_fp  # type: ignore

    def hash_list(self) -> List[str]:
        """
        Return all the products and reactants as hashed SMILES

        :return: the hashes of the SMILES string
        """
        mols = self.reaction_smiles().replace(".", ">>").split(">>")
        return [hashlib.sha224(mol.encode("utf8")).hexdigest() for mol in mols]

    def hash_key(self) -> str:
        """
        Return a code that can be use to identify the reaction

        :return: the hash code
        """
        reactants = sorted([mol.inchi_key for mol in self._reactants_getter()])  # type: ignore
        products = sorted([mol.inchi_key for mol in self._products_getter()])  # type: ignore
        hash_ = hashlib.sha224()
        for item in reactants + [">>"] + products:
            hash_.update(item.encode())
        return hash_.hexdigest()

    def rd_reaction_from_smiles(self) -> RdReaction:
        """
        The reaction as a RDkit reaction object but created from the reaction smiles
        instead of the SMARTS of the template.

        :return: the reaction object
        """
        return AllChem.ReactionFromSmarts(self.reaction_smiles(), useSmiles=True)

    def reaction_smiles(self) -> str:
        """
        Get the reaction SMILES, i.e. the SMILES of the reactants and products joined together

        :return: the SMILES
        """
        reactants = ".".join(mol.smiles for mol in self._reactants_getter())  # type: ignore
        products = ".".join(mol.smiles for mol in self._products_getter())  # type: ignore
        return "%s>>%s" % (reactants, products)

class FixedRetroReaction(_ReactionInterfaceMixin):
    """
    A retrosynthesis reaction that has the same interface as `RetroReaction`
    but it is fixed so it does not support SMARTS application or any creation of reactants.

    The reactants are set by using the `reactants` property.

    :ivar mol: the UniqueMolecule object that this reaction is applied to
    :ivar smiles: the SMILES representation of the RDKit reaction
    :ivar metadata: meta data associated with the reaction
    :ivar reactants: the reactants of this reaction

    :param mol: the molecule
    :param smiles: the SMILES of the reaction
    :param metadata: some meta data
    """

    def __init__(
        self, mol: UniqueMolecule, smiles: str = "", metadata: StrDict = None
    ) -> None:
        self.mol = mol
        self.smiles = smiles
        self.metadata = metadata or {}
        self.reactants: Tuple[Tuple[UniqueMolecule, ...], ...] = ()

    def copy(self) -> "FixedRetroReaction":
        """
        Shallow copy of this instance.

        :return: the copy
        """
        new_reaction = FixedRetroReaction(self.mol, self.smiles, self.metadata)
        new_reaction.reactants = tuple(mol_list for mol_list in self.reactants)
        return new_reaction

    def _products_getter(self) -> Tuple[UniqueMolecule, ...]:
        return self.reactants[0]

    def _reactants_getter(self) -> List[UniqueMolecule]:
        return [self.mol]
    
class TreeMolecule(Molecule):
    """
    A special molecule that keeps a reference to a parent molecule.
    If the class is instantiated without specifying the `transform` argument,
    it is computed by increasing the value of the `parent.transform` variable.
    If no parent is provided the atoms with atom mapping number are tracked
    and inherited to children.
    :ivar mapped_mol: the tracked molecule with atom mappings
    :ivar mapped_smiles: the SMILES of the tracked molecule with atom mappings
    :ivar original_smiles: the SMILES as passed when instantiating the class
    :ivar parent: parent molecule
    :ivar transform: a numerical number corresponding to the depth in the tree
    :ivar tracked_atom_indices: tracked atom indices and what indices they correspond to in this molecule
    :param parent: a TreeMolecule object that is the parent
    :param transform: the transform value, defaults to None
    :param rd_mol: a RDKit mol object to encapsulate, defaults to None
    :param smiles: a SMILES to convert to a molecule object, defaults to None
    :param sanitize: if True, the molecule will be immediately sanitized, defaults to False
    :param mapping_update_callback: if given will call this method before setting up the `mapped_smiles`
    :raises MoleculeException: if neither rd_mol or smiles is given, or if the molecule could not be sanitized
    """

    # pylint: disable=too-many-arguments
    def __init__(
        self,
        parent: Optional["TreeMolecule"],
        transform: int = None,
        rd_mol: RdMol = None,
        smiles: str = None,
        sanitize: bool = False,
        mapping_update_callback: Callable[["TreeMolecule"], None] = None,
    ) -> None:
        super().__init__(rd_mol=rd_mol, smiles=smiles, sanitize=sanitize)
        self.parent = parent
        if transform is None and parent and parent.transform is not None:
            self.transform: int = parent.transform + 1
        else:
            self.transform = transform or 0

        self.original_smiles = smiles
        self.tracked_atom_indices: Dict[int, Optional[int]] = {}
        self.mapped_mol = Chem.Mol(self.rd_mol)
        if not self.parent:
            self._init_tracking()
        elif mapping_update_callback is not None:
            mapping_update_callback(self)

        AllChem.SanitizeMol(self.mapped_mol)
        self.mapped_smiles = Chem.MolToSmiles(self.mapped_mol)

        if self.parent:
            self.remove_atom_mapping()
            self._update_tracked_atoms()

    @property
    def mapping_to_index(self) -> Dict[int, int]:
        """Return a dictionary mapping to atom mappings to atom indices"""
        if not self._atom_mappings:
            self._atom_mappings = {
                atom.GetAtomMapNum(): atom.GetIdx()
                for atom in self.mapped_mol.GetAtoms()
                if atom.GetAtomMapNum()
            }
        return self._atom_mappings

    def _init_tracking(self):
        self.tracked_atom_indices = dict(self.mapping_to_index)
        for idx, atom in enumerate(self.mapped_mol.GetAtoms()):
            atom.SetAtomMapNum(idx + 1)
        self._atom_mappings = {}

    def _update_tracked_atoms(self) -> None:
        if self.parent is None:
            return

        if not self.parent.tracked_atom_indices:
            return

        parent2child_map = {
            atom_index: self.mapping_to_index.get(mapping_index)
            for mapping_index, atom_index in self.parent.mapping_to_index.items()
        }

        self.tracked_atom_indices = {
            tracked_index: parent2child_map[parent_index]  # type: ignore
            for tracked_index, parent_index in self.parent.tracked_atom_indices.items()
        }

class RetroReaction(abc.ABC, _ReactionInterfaceMixin):
    """
    A retrosynthesis reaction. Only a single molecule is the reactant.
    This is an abstract class and child classes needs to implement the `_apply` and `_make_smiles` functions
    that should create the reactants molecule objects and the reaction SMILES representation, respectively.
    :ivar mol: the TreeMolecule object that this reaction is applied to
    :ivar index: a unique index of this reaction,
                 to count for the fact that a reaction can produce more than one outcome
    :ivar metadata: meta data associated with the reaction
    :param mol: the molecule
    :param index: the index, defaults to 0
    :param metadata: some meta data
    :params kwargs: any extra parameters for child classes
    """

    _required_kwargs: List[str] = []

    def __init__(
        self, mol: TreeMolecule, index: int = 0, metadata: StrDict = None, **kwargs: Any
    ) -> None:
        if any(name not in kwargs for name in self._required_kwargs):
            raise KeyError(
                f"A {self.__class__.__name__} class needs to be initiated "
                f"with keyword arguments: {', '.join(self._required_kwargs)}"
            )
        self.mol = mol
        self.index = index
        self.metadata: StrDict = metadata or {}
        self._reactants: Optional[Tuple[Tuple[TreeMolecule, ...], ...]] = None
        self._smiles: Optional[str] = None
        self._kwargs: StrDict = kwargs

    @classmethod
    def from_serialization(
        cls, init_args: StrDict, reactants: List[List[TreeMolecule]]
    ) -> "RetroReaction":
        """
        Create an object from a serialization. It does
        1) instantiate an object using the `init_args` and
        2) set the reactants to a tuple-form of `reactants
        :param init_args: the arguments passed to the `__init__` method
        :param reactants: the reactants
        :return: the deserialized object
        """
        obj = cls(**init_args)
        obj._reactants = tuple(tuple(mol for mol in lst_) for lst_ in reactants)
        return obj

    def __str__(self) -> str:
        return f"reaction on molecule {self.mol.smiles}"

    @property
    def reactants(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        """
        Returns the reactant molecules.
        Apply the reaction if necessary.
        :return: the products of the reaction
        """
        if not self._reactants:
            self._reactants = self._apply()
        return self._reactants

    @property
    def smiles(self) -> str:
        """
        The reaction as a SMILES
        :return: the SMILES
        """
        if self._smiles is None:
            try:
                self._smiles = self._make_smiles()
            except ValueError:
                self._smiles = ""  # noqa
        return self._smiles

    @property
    def unqueried(self) -> bool:
        """
        Return True if the reactants has never been retrieved
        """
        return self._reactants is None

    def copy(self, index: int = None) -> "RetroReaction":
        """
        Shallow copy of this instance.
        :param index: new index, defaults to None
        :return: the copy
        """
        # pylint: disable=protected-access
        index = index if index is not None else self.index
        new_reaction = self.__class__(self.mol, index, self.metadata, **self._kwargs)
        new_reaction._reactants = tuple(mol_list for mol_list in self._reactants or [])
        new_reaction._smiles = self._smiles
        return new_reaction

    def mapped_reaction_smiles(self) -> str:
        """
        Get the mapped reaction SMILES if it exists
        :return: the SMILES
        """
        reactants = self.mol.mapped_smiles
        products = ".".join(mol.mapped_smiles for mol in self._products_getter())
        return reactants + ">>" + products

    def to_dict(self) -> StrDict:
        """
        Return the retro reaction as dictionary
        This dictionary is not suitable for serialization, but is used by other serialization routines
        The elements of the dictionary can be used to instantiate a new reaction object
        """
        return {
            "mol": self.mol,
            "index": self.index,
            "metadata": dict(self.metadata),
        }

    @abc.abstractmethod
    def _apply(self) -> Tuple[Tuple[TreeMolecule, ...], ...]:
        pass

    @abc.abstractmethod
    def _make_smiles(self) -> str:
        pass

    def _products_getter(self) -> Tuple[TreeMolecule, ...]:
        return self.reactants[self.index]

    def _reactants_getter(self) -> List[TreeMolecule]:
        return [self.mol]

    @staticmethod
    def _update_unmapped_atom_num(mol: TreeMolecule, exclude_nums: Set[int]) -> None:
        mapped_nums = {num for num in mol.mapping_to_index.keys() if 0 < num < 900}
        offset = max(mapped_nums) + 1 if mapped_nums else 1
        for atom in mol.mapped_mol.GetAtoms():
            if 0 < atom.GetAtomMapNum() < 900:
                continue
            while offset in exclude_nums:
                offset += 1
            atom.SetAtomMapNum(offset)
            exclude_nums.add(offset)


    """
    Encapsulation of a bipartite reaction tree of a single route.
    The nodes consists of either FixedRetroReaction or UniqueMolecule objects.

    The reaction tree is initialized at instantiation and is not supposed to
    be updated.

    :ivar graph: the bipartite graph
    :ivar is_solved: if all of the leaf nodes are in stock
    :ivar root: the root of the tree
    :ivar created_at_iteration: iteration the reaction tree was created
    """

    def __init__(self) -> None:
        self.graph = nx.DiGraph()
        self.root = none_molecule()
        self.is_solved: bool = False
        self.created_at_iteration: Optional[int] = None

    @classmethod
    def from_dict(cls, tree_dict: StrDict) -> "ReactionTree":
        """
        Create a new ReactionTree by parsing a dictionary.

        This is supposed to be the opposite of ``to_dict``,
        but because that format loses information, the returned
        object is not a full copy as the stock will only contain
        the list of molecules marked as ``in_stock`` in the dictionary.

        The returned object should be sufficient to e.g. generate an image of the route.

        :param tree_dict: the dictionary representation
        :returns: the reaction tree
        """
        return ReactionTreeFromDict(tree_dict).tree

    @property
    def metadata(self) -> StrDict:
        """Return a dicitionary with route metadata"""
        return {
            "created_at_iteration": self.created_at_iteration,
            "is_solved": self.is_solved,
        }

    def depth(self, node: Union[UniqueMolecule, FixedRetroReaction]) -> int:
        """
        Return the depth of a node in the route

        :param node: the query node
        :return: the depth
        """
        return self.graph.nodes[node].get("depth", -1)

    def distance_to(self, other: "ReactionTree", content: str = "both") -> float:
        """
        Calculate the distance to another reaction tree

        This is a tree edit distance, with unit cost to
        insert and deleted nodes, and the Jaccard distance for substituting nodes

        :param other: the reaction tree to compare to
        :param content: determine what part of the tree to include in the calculation
        :return: the distance between the routes
        """
        if not SUPPORT_DISTANCES:
            raise ValueError(
                "Distance calculations are not supported by this installation."
                " Please install aizynthfinder with extras dependencies."
            )
        calculator = route_distances_calculator("ted", content=content)
        distances = calculator([self.to_dict(), other.to_dict()])
        return distances[0, 1]

    def hash_key(self) -> str:
        """
        Calculates a hash code for the tree using the sha224 hash function recursively

        :return: the hash key
        """
        return self._hash_func(self.root)

    def in_stock(self, node: Union[UniqueMolecule, FixedRetroReaction]) -> bool:
        """
        Return if a node in the route is in stock

        Note that is a property set on creation and as such is not updated.

        :param node: the query node
        :return: if the molecule is in stock
        """
        return self.graph.nodes[node].get("in_stock", False)

    def is_branched(self) -> bool:
        """
        Returns if the route is branched

        i.e. checks if the maximum depth is not equal to the number
        of reactions.
        """
        nsteps = len(list(self.reactions()))
        max_depth = max(self.depth(leaf) for leaf in self.leafs())
        return nsteps != max_depth // 2

    def leafs(self) -> Iterable[UniqueMolecule]:
        """
        Generates the molecules nodes of the reaction tree that has no predecessors,
        i.e. molecules that has not been broken down

        :yield: the next leaf molecule in the tree
        """
        for node in self.graph:
            if isinstance(node, UniqueMolecule) and not self.graph[node]:
                yield node

    def molecules(self) -> Iterable[UniqueMolecule]:
        """
        Generates the molecule nodes of the reaction tree

        :yield: the next molecule in the tree
        """
        for node in self.graph:
            if isinstance(node, UniqueMolecule):
                yield node

    def reactions(self) -> Iterable[FixedRetroReaction]:
        """
        Generates the reaction nodes of the reaction tree

        :yield: the next reaction in the tree
        """
        for node in self.graph:
            if not isinstance(node, Molecule):
                yield node

    def subtrees(self) -> Iterable[ReactionTree]:
        """
        Generates the subtrees of this reaction tree a
        subtree is a reaction treee starting at a molecule node that has children.

        :yield: the next subtree
        """

        def create_subtree(root_node):
            subtree = ReactionTree()
            subtree.root = root_node
            subtree.graph = dfs_tree(self.graph, root_node)
            for node in subtree.graph:
                prop = dict(self.graph.nodes[node])
                prop["depth"] -= self.graph.nodes[root_node].get("depth", 0)
                if "transform" in prop:
                    prop["transform"] -= self.graph.nodes[root_node].get("transform", 0)
                subtree.graph.nodes[node].update(prop)
            subtree.is_solved = all(subtree.in_stock(node) for node in subtree.leafs())
            return subtree

        for node in self.molecules():
            if node is not self.root and self.graph[node]:
                yield create_subtree(node)

    def to_dict(self, include_metadata=False) -> StrDict:
        """
        Returns the reaction tree as a dictionary in a pre-defined format.
        :param include_metadata: if True include metadata
        :return: the reaction tree
        """
        return self._build_dict(self.root, include_metadata=include_metadata)

    def to_image(
        self,
        in_stock_colors: FrameColors = None,
        show_all: bool = True,
    ) -> PilImage:
        """
        Return a pictorial representation of the route

        :param in_stock_colors: the colors around molecules, defaults to {True: "green", False: "orange"}
        :param show_all: if True, also show nodes that are marked as hidden
        :return: the image of the route
        """
        factory = RouteImageFactory(
            self.to_dict(), in_stock_colors=in_stock_colors, show_all=show_all
        )
        return factory.image

    def to_json(self, include_metadata=False) -> str:
        """
        Returns the reaction tree as a JSON string in a pre-defined format.

        :return: the reaction tree
        """
        return json.dumps(
            self.to_dict(include_metadata=include_metadata), sort_keys=False, indent=2
        )

    def _build_dict(
        self,
        node: Union[UniqueMolecule, FixedRetroReaction],
        dict_: StrDict = None,
        include_metadata=False,
    ) -> StrDict:
        if dict_ is None:
            dict_ = {}

        if node is self.root and include_metadata:
            dict_["route_metadata"] = self.metadata

        dict_["type"] = "mol" if isinstance(node, Molecule) else "reaction"
        dict_["hide"] = self.graph.nodes[node].get("hide", False)
        dict_["smiles"] = node.smiles
        if isinstance(node, UniqueMolecule):
            dict_["is_chemical"] = True
            dict_["in_stock"] = self.in_stock(node)
        elif isinstance(node, FixedRetroReaction):
            dict_["is_reaction"] = True
            dict_["metadata"] = dict(node.metadata)
        else:
            raise ValueError(
                f"This is an invalid reaction tree. Unknown node type {type(node)}"
            )

        dict_["children"] = []

        children = list(self.graph.successors(node))
        if isinstance(node, FixedRetroReaction):
            children.sort(key=operator.attrgetter("weight"))
        for child in children:
            child_dict = self._build_dict(child)
            dict_["children"].append(child_dict)

        if not dict_["children"]:
            del dict_["children"]
        return dict_

    def _hash_func(self, node: Union[FixedRetroReaction, UniqueMolecule]) -> str:
        if isinstance(node, UniqueMolecule):
            hash_ = hashlib.sha224(node.inchi_key.encode())
        else:
            hash_ = hashlib.sha224(node.hash_key().encode())
        child_hashes = sorted(
            self._hash_func(child) for child in self.graph.successors(node)
        )
        for child_hash in child_hashes:
            hash_.update(child_hash.encode())
        return hash_.hexdigest()

class ReactionTreeLoader(abc.ABC):
    """
    Base class for classes that creates a reaction tree object

    This class makes sure that node attributes are set after the
    graph is generated, and provides utility methods.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._unique_mols: Dict[int, UniqueMolecule] = {}
        self._unique_reactions: Dict[int, FixedRetroReaction] = {}
        self.tree = ReactionTree()
        self._load(*args, **kwargs)

        self.tree.is_solved = all(
            self.tree.in_stock(node) for node in self.tree.leafs()
        )

    def _add_node(
        self,
        node: Union[UniqueMolecule, FixedRetroReaction],
        depth: int = 0,
        transform: int = 0,
        in_stock: bool = False,
        hide: bool = False,
    ) -> None:
        attributes = {
            "hide": hide,
            "depth": depth,
        }
        if isinstance(node, Molecule):
            attributes.update({"transform": transform, "in_stock": in_stock})
            if not self.tree.root:
                self.tree.root = node
        self.tree.graph.add_node(node, **attributes)

    @abc.abstractmethod
    def _load(self, *args: Any, **kwargs: Any) -> None:
        pass

    def _unique_mol(self, molecule: Molecule) -> UniqueMolecule:
        id_ = id(molecule)
        if id_ not in self._unique_mols:
            self._unique_mols[id_] = molecule.make_unique()
        return self._unique_mols[id_]

    def _unique_reaction(self, reaction: RetroReaction) -> FixedRetroReaction:
        id_ = id(reaction)
        if id_ not in self._unique_reactions:
            metadata = dict(reaction.metadata)
            if ":" in reaction.mapped_reaction_smiles():
                metadata["mapped_reaction_smiles"] = reaction.mapped_reaction_smiles()
            self._unique_reactions[id_] = FixedRetroReaction(
                self._unique_mol(reaction.mol),
                smiles=reaction.smiles,
                metadata=metadata,
            )
        return self._unique_reactions[id_]


class ReactionTreeFromDict(ReactionTreeLoader):
    """Creates a reaction tree object from a dictionary"""

    def _load(self, tree_dict: StrDict) -> None:  # type: ignore
        if tree_dict.get("route_metadata"):
            self.tree.created_at_iteration = tree_dict["route_metadata"].get(
                "created_at_iteration"
            )
        self._parse_tree_dict(tree_dict)

    def _parse_tree_dict(self, tree_dict: StrDict, ncalls: int = 0) -> UniqueMolecule:
        product_node = UniqueMolecule(smiles=tree_dict["smiles"])
        self._add_node(
            product_node,
            depth=2 * ncalls,
            transform=ncalls,
            hide=tree_dict.get("hide", False),
            in_stock=tree_dict["in_stock"],
        )

        rxn_tree_dict = tree_dict.get("children", [])
        if not rxn_tree_dict:
            return product_node

        rxn_tree_dict = rxn_tree_dict[0]
        reaction_node = FixedRetroReaction(
            product_node,
            smiles=rxn_tree_dict["smiles"],
            metadata=rxn_tree_dict.get("metadata", {}),
        )
        self._add_node(
            reaction_node, depth=2 * ncalls + 1, hide=rxn_tree_dict.get("hide", False)
        )
        self.tree.graph.add_edge(product_node, reaction_node)

        reactant_nodes = []
        for reactant_tree in rxn_tree_dict.get("children", []):
            reactant_node = self._parse_tree_dict(reactant_tree, ncalls + 1)
            self.tree.graph.add_edge(reaction_node, reactant_node)
            reactant_nodes.append(reactant_node)
        reaction_node.reactants = (tuple(reactant_nodes),)

        return product_node


class ReactionTreeFromExpansion(ReactionTreeLoader):
    """
    Create a ReactionTree from a single reaction

    This is mainly intended as a convenience function for the expander interface
    """

    def _load(self, reaction: RetroReaction) -> None:  # type: ignore
        root = self._unique_mol(reaction.mol)
        self._add_node(root)

        rxn = self._unique_reaction(reaction)
        if hasattr(reaction, "smarts"):
            rxn.metadata["smarts"] = reaction.smarts  # type: ignore
        self._add_node(rxn)
        self.tree.graph.add_edge(root, rxn)

        reactant_nodes = []
        for reactant in reaction.reactants[0]:
            reactant_node = self._unique_mol(reactant)
            reactant_nodes.append(reactant_node)
            self._add_node(reactant_node)
            self.tree.graph.add_edge(rxn, reactant_node)
        rxn.reactants = (tuple(reactant_nodes),)
         
class ReactionTreeLoader(abc.ABC):
    """
    Base class for classes that creates a reaction tree object

    This class makes sure that node attributes are set after the
    graph is generated, and provides utility methods.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self._unique_mols: Dict[int, UniqueMolecule] = {}
        self._unique_reactions: Dict[int, FixedRetroReaction] = {}
        self.tree = ReactionTree()
        self._load(*args, **kwargs)

        self.tree.is_solved = all(
            self.tree.in_stock(node) for node in self.tree.leafs()
        )

    def _add_node(
        self,
        node: Union[UniqueMolecule, FixedRetroReaction],
        depth: int = 0,
        transform: int = 0,
        in_stock: bool = False,
        hide: bool = False,
    ) -> None:
        attributes = {
            "hide": hide,
            "depth": depth,
        }
        if isinstance(node, Molecule):
            attributes.update({"transform": transform, "in_stock": in_stock})
            if not self.tree.root:
                self.tree.root = node
        self.tree.graph.add_node(node, **attributes)

    @abc.abstractmethod
    def _load(self, *args: Any, **kwargs: Any) -> None:
        pass

    def _unique_mol(self, molecule: Molecule) -> UniqueMolecule:
        id_ = id(molecule)
        if id_ not in self._unique_mols:
            self._unique_mols[id_] = molecule.make_unique()
        return self._unique_mols[id_]

    def _unique_reaction(self, reaction: RetroReaction) -> FixedRetroReaction:
        id_ = id(reaction)
        if id_ not in self._unique_reactions:
            metadata = dict(reaction.metadata)
            if ":" in reaction.mapped_reaction_smiles():
                metadata["mapped_reaction_smiles"] = reaction.mapped_reaction_smiles()
            self._unique_reactions[id_] = FixedRetroReaction(
                self._unique_mol(reaction.mol),
                smiles=reaction.smiles,
                metadata=metadata,
            )
        return self._unique_reactions[id_]

class ReactionTreeFromDict(ReactionTreeLoader):
    """Creates a reaction tree object from a dictionary"""

    def _load(self, tree_dict: StrDict) -> None:  # type: ignore
        if tree_dict.get("route_metadata"):
            self.tree.created_at_iteration = tree_dict["route_metadata"].get(
                "created_at_iteration"
            )
        self._parse_tree_dict(tree_dict)

    def _parse_tree_dict(self, tree_dict: StrDict, ncalls: int = 0) -> UniqueMolecule:
        product_node = UniqueMolecule(smiles=tree_dict["smiles"])
        self._add_node(
            product_node,
            depth=2 * ncalls,
            transform=ncalls,
            hide=tree_dict.get("hide", False),
            in_stock=tree_dict["in_stock"],
        )

        rxn_tree_dict = tree_dict.get("children", [])
        if not rxn_tree_dict:
            return product_node

        rxn_tree_dict = rxn_tree_dict[0]
        reaction_node = FixedRetroReaction(
            product_node,
            smiles=rxn_tree_dict["smiles"],
            metadata=rxn_tree_dict.get("metadata", {}),
        )
        self._add_node(
            reaction_node, depth=2 * ncalls + 1, hide=rxn_tree_dict.get("hide", False)
        )
        self.tree.graph.add_edge(product_node, reaction_node)

        reactant_nodes = []
        for reactant_tree in rxn_tree_dict.get("children", []):
            reactant_node = self._parse_tree_dict(reactant_tree, ncalls + 1)
            self.tree.graph.add_edge(reaction_node, reactant_node)
            reactant_nodes.append(reactant_node)
        reaction_node.reactants = (tuple(reactant_nodes),)
        
        print(product_node)

        return product_node
        


data = pd.read_json("/Users/song-inhyeok/CODING/short_test/output.json", orient="table")
all_trees = data.trees.values  # This contains a list of all the trees for all the compounds
trees_for_first_target = all_trees[0]

for itree, tree in enumerate(trees_for_first_target):
    ReactionTreeFromDict(tree)