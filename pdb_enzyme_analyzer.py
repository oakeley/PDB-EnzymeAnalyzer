#!/usr/bin/env python3
"""
PDB Enzymatic Function and Structural Analysis Pipeline

Analyzes PDB files to predict:
- Enzymatic catalytic functions
- Binding sites and active sites
- Structural properties and differences
- Physicochemical characteristics
- Functional annotations with confidence scores

Requirements:
- biopython
- numpy, scipy, matplotlib
- sklearn
- requests (for database queries)
"""

import argparse
import sys
import time
import logging
import json
import base64
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO
import warnings
warnings.filterwarnings('ignore')

# Biopython imports
from Bio.PDB import PDBParser, DSSP, NeighborSearch, Selection
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import PDBIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
import Bio.PDB.vectors as vectors

# Scientific computing
from scipy.spatial.distance import cdist, pdist
from scipy.stats import entropy
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

# For enzyme prediction
import re
import subprocess
from dataclasses import dataclass, asdict
from collections import defaultdict, Counter

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress matplotlib font manager debug messages
logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

def three_to_one(resname):
    """Convert a three-letter amino acid code to one-letter."""
    return protein_letters_3to1.get(resname.capitalize(), 'X')  # 'X' for unknowns

@dataclass
class CatalyticSite:
    """Container for catalytic site information"""
    residues: List[Tuple[str, int, str]]  # (residue_name, residue_number, chain_id)
    center_coords: np.ndarray
    site_type: str
    confidence: float
    description: str
    supporting_evidence: List[str]

@dataclass
class BindingSite:
    """Container for binding site information"""
    residues: List[Tuple[str, int, str]]
    center_coords: np.ndarray
    volume: float
    surface_area: float
    hydrophobicity: float
    charge: float
    binding_likelihood: float
    site_classification: str

@dataclass
class EnzymeFunction:
    """Container for predicted enzyme function"""
    ec_number: Optional[str]
    enzyme_class: str
    function_description: str
    confidence: float
    supporting_features: List[str]
    catalytic_mechanism: str
    substrate_prediction: List[str]
    is_biosynthetic: bool = False
    pathway_type: Optional[str] = None  # e.g., "Polyketide", "NRPS", "Terpene"
    microbial_origin: Optional[str] = None  # "Bacterial", "Fungal", "Unknown"

@dataclass
class StructuralFeatures:
    """Container for structural characteristics"""
    secondary_structure: Dict[str, float]
    geometric_properties: Dict[str, float]
    surface_properties: Dict[str, float]
    domain_architecture: List[Dict]
    unusual_features: List[str]

@dataclass
class ProteinAnalysisResult:
    """Complete analysis result for a protein"""
    pdb_id: str
    pdb_file: str
    sequence: str
    length: int
    catalytic_sites: List[CatalyticSite]
    binding_sites: List[BindingSite]
    enzyme_functions: List[EnzymeFunction]
    structural_features: StructuralFeatures
    physicochemical_properties: Dict[str, float]
    quality_metrics: Dict[str, float]
    analysis_confidence: float

class CatalyticSitePredictor:
    """Predicts catalytic sites using structural patterns and motifs"""
    
    def __init__(self):
        # Common catalytic residues
        self.catalytic_residues = {
            'SER': 'Serine protease/esterase',
            'CYS': 'Cysteine protease/thiol chemistry',
            'HIS': 'Histidine catalysis/proton transfer',
            'ASP': 'Aspartic acid catalysis/nucleophile',
            'GLU': 'Glutamic acid catalysis/nucleophile',
            'LYS': 'Lysine catalysis/base',
            'ARG': 'Arginine catalysis/base',
            'TYR': 'Tyrosine radical chemistry',
            'TRP': 'Tryptophan radical chemistry'
        }
        
        # Metal coordination patterns
        self.metal_binding_patterns = {
            'zinc_binding': ['HIS', 'CYS', 'ASP', 'GLU'],
            'iron_binding': ['HIS', 'CYS', 'MET'],
            'calcium_binding': ['ASP', 'GLU', 'ASN', 'GLN'],
            'magnesium_binding': ['ASP', 'GLU', 'ASN']
        }
        
        # Known catalytic motifs
        self.catalytic_motifs = {
            'serine_protease_triad': [('SER', 'HIS', 'ASP'), 6.0],
            'cysteine_protease_triad': [('CYS', 'HIS', 'ASN'), 6.0],
            'zinc_peptidase': [('HIS', 'HIS', 'GLU'), 8.0],
            'phosphatase_motif': [('ASP', 'ASP', 'GLY'), 5.0]
        }
    
    def find_catalytic_sites(self, structure, chain_id: str = None) -> List[CatalyticSite]:
        """Find potential catalytic sites in the structure"""
        sites = []
        
        # Get all chains or specific chain
        chains = [structure[0][chain_id]] if chain_id else structure[0]
        
        for chain in chains:
            # Find catalytic residue clusters
            catalytic_residues = self._get_catalytic_residues(chain)
            if len(catalytic_residues) < 2:
                continue
            
            # Find spatial clusters of catalytic residues
            clusters = self._find_spatial_clusters(catalytic_residues)
            
            for cluster in clusters:
                if len(cluster) >= 2:  # Minimum 2 residues for a catalytic site
                    site = self._analyze_catalytic_cluster(cluster, chain.id)
                    if site:
                        sites.append(site)
        
        # Look for metal binding sites
        metal_sites = self._find_metal_binding_sites(structure)
        sites.extend(metal_sites)
        
        return sites
    
    def _get_catalytic_residues(self, chain) -> List[Tuple]:
        """Get all potentially catalytic residues with coordinates"""
        catalytic_residues = []
        
        for residue in chain:
            if residue.get_resname() in self.catalytic_residues:
                try:
                    if residue.get_resname() in ['SER', 'TYR', 'THR', 'CYS']:
                        coord = residue['OG'].coord if 'OG' in residue else residue['SG'].coord
                    elif residue.get_resname() in ['HIS']:
                        coord = residue['ND1'].coord if 'ND1' in residue else residue['NE2'].coord
                    elif residue.get_resname() in ['ASP', 'GLU']:
                        coord = residue['OD1'].coord if 'OD1' in residue else residue['OE1'].coord
                    elif residue.get_resname() in ['LYS']:
                        coord = residue['NZ'].coord
                    elif residue.get_resname() in ['ARG']:
                        coord = residue['NH1'].coord if 'NH1' in residue else residue['NH2'].coord
                    else:
                        coord = residue['CA'].coord
                    
                    catalytic_residues.append((
                        residue.get_resname(),
                        residue.get_id()[1],
                        coord,
                        residue
                    ))
                except KeyError:
                    # Missing atoms, skip
                    continue
        
        return catalytic_residues
    
    def _find_spatial_clusters(self, residues: List[Tuple], max_distance: float = 8.0) -> List[List]:
        """Find spatial clusters of catalytic residues"""
        if len(residues) < 2:
            return []
        
        # Extract coordinates
        coords = np.array([res[2] for res in residues])
        
        # Use DBSCAN clustering
        clustering = DBSCAN(eps=max_distance, min_samples=2).fit(coords)
        
        clusters = []
        for cluster_id in set(clustering.labels_):
            if cluster_id != -1:  # -1 is noise
                cluster_residues = [residues[i] for i in range(len(residues)) 
                                   if clustering.labels_[i] == cluster_id]
                clusters.append(cluster_residues)
        
        return clusters
    
    def _analyze_catalytic_cluster(self, cluster: List[Tuple], chain_id: str) -> Optional[CatalyticSite]:
        """Analyze a cluster of catalytic residues"""
        if len(cluster) < 2:
            return None
        
        residue_names = [res[0] for res in cluster]
        residue_numbers = [res[1] for res in cluster]
        coords = np.array([res[2] for res in cluster])
        center = np.mean(coords, axis=0)
        
        # Check for known catalytic motifs
        motif_type = self._identify_catalytic_motif(residue_names)
        
        # Calculate confidence based on known patterns
        confidence = self._calculate_catalytic_confidence(residue_names, coords)
        
        # Create description
        description = f"Potential catalytic site with {len(cluster)} residues: {', '.join(residue_names)}"
        
        supporting_evidence = []
        if motif_type != "Unknown":
            supporting_evidence.append(f"Matches {motif_type} pattern")
        
        # Check spatial arrangement
        if len(cluster) >= 3:
            distances = pdist(coords)
            if np.max(distances) < 10.0:  # Tight cluster
                supporting_evidence.append("Tight spatial clustering (<10Å)")
        
        return CatalyticSite(
            residues=[(res[0], res[1], chain_id) for res in cluster],
            center_coords=center,
            site_type=motif_type,
            confidence=confidence,
            description=description,
            supporting_evidence=supporting_evidence
        )
    
    def _identify_catalytic_motif(self, residue_names: List[str]) -> str:
        """Identify known catalytic motifs"""
        residue_set = set(residue_names)
        
        # Check for serine protease triad
        if {'SER', 'HIS', 'ASP'}.issubset(residue_set):
            return "Serine protease triad"
        elif {'CYS', 'HIS', 'ASN'}.issubset(residue_set):
            return "Cysteine protease triad"
        elif residue_names.count('HIS') >= 2 and 'GLU' in residue_names:
            return "Zinc binding motif"
        elif {'ASP', 'GLU'}.intersection(residue_set) and len(residue_set) >= 2:
            return "Acid catalysis motif"
        elif {'LYS', 'ARG'}.intersection(residue_set):
            return "Base catalysis motif"
        else:
            return "Unknown"
    
    def _calculate_catalytic_confidence(self, residue_names: List[str], coords: np.ndarray) -> float:
        """Calculate confidence score for catalytic site prediction"""
        confidence = 0.0
        
        # Bonus for known catalytic residues
        known_catalytic = sum(1 for res in residue_names if res in self.catalytic_residues)
        confidence += known_catalytic * 0.2
        
        # Bonus for spatial arrangement
        if len(coords) >= 3:
            distances = pdist(coords)
            avg_distance = np.mean(distances)
            if 3.0 <= avg_distance <= 8.0:  # Optimal range
                confidence += 0.3
        
        # Bonus for known motifs
        motif = self._identify_catalytic_motif(residue_names)
        if motif != "Unknown":
            confidence += 0.4
        
        return min(1.0, confidence)
    
    def _find_metal_binding_sites(self, structure) -> List[CatalyticSite]:
        """Find potential metal binding sites"""
        sites = []
        
        for chain in structure[0]:
            # Look for metal coordination patterns
            metal_coordinating = []
            
            for residue in chain:
                resname = residue.get_resname()
                if resname in ['HIS', 'CYS', 'ASP', 'GLU', 'MET', 'ASN', 'GLN']:
                    try:
                        if resname == 'HIS':
                            coord = residue['ND1'].coord if 'ND1' in residue else residue['NE2'].coord
                        elif resname == 'CYS':
                            coord = residue['SG'].coord
                        elif resname in ['ASP', 'GLU']:
                            coord = residue['OD1'].coord if 'OD1' in residue else residue['OE1'].coord
                        elif resname == 'MET':
                            coord = residue['SD'].coord
                        else:
                            coord = residue['CA'].coord
                        
                        metal_coordinating.append((resname, residue.get_id()[1], coord, residue))
                    except KeyError:
                        continue
            
            # Find clusters that could coordinate metals
            if len(metal_coordinating) >= 3:
                coords = np.array([res[2] for res in metal_coordinating])
                clustering = DBSCAN(eps=6.0, min_samples=3).fit(coords)
                
                for cluster_id in set(clustering.labels_):
                    if cluster_id != -1:
                        cluster_residues = [metal_coordinating[i] for i in range(len(metal_coordinating)) 
                                           if clustering.labels_[i] == cluster_id]
                        
                        if len(cluster_residues) >= 3:
                            residue_names = [res[0] for res in cluster_residues]
                            center = np.mean([res[2] for res in cluster_residues], axis=0)
                            
                            # Determine metal type
                            metal_type = self._predict_metal_type(residue_names)
                            
                            site = CatalyticSite(
                                residues=[(res[0], res[1], chain.id) for res in cluster_residues],
                                center_coords=center,
                                site_type=f"Metal binding ({metal_type})",
                                confidence=0.7,
                                description=f"Potential {metal_type} binding site",
                                supporting_evidence=[f"Metal coordinating residues: {', '.join(residue_names)}"]
                            )
                            sites.append(site)
        
        return sites
    
    def _predict_metal_type(self, residue_names: List[str]) -> str:
        """Predict the type of metal based on coordinating residues"""
        residue_set = set(residue_names)
        
        if residue_names.count('HIS') >= 2 and ('GLU' in residue_set or 'ASP' in residue_set):
            return "Zinc"
        elif residue_names.count('CYS') >= 2:
            return "Iron-sulfur"
        elif residue_names.count('HIS') >= 2 and 'MET' in residue_set:
            return "Copper"
        elif ('ASP' in residue_set or 'GLU' in residue_set) and residue_names.count('ASP') + residue_names.count('GLU') >= 2:
            return "Calcium/Magnesium"
        else:
            return "Unknown metal"

class BindingSiteAnalyzer:
    """Analyzes potential binding sites and cavities"""
    
    def __init__(self):
        self.probe_radius = 1.4  # Water probe radius
        
    def find_binding_sites(self, structure) -> List[BindingSite]:
        """Find potential binding sites using geometric analysis"""
        binding_sites = []
        
        for chain in structure[0]:
            # Get all atoms
            atoms = []
            for residue in chain:
                for atom in residue:
                    if atom.element != 'H':  # Skip hydrogens
                        atoms.append(atom)
            
            if len(atoms) < 10:
                continue
            
            # Find cavities using alpha shape approach
            cavities = self._find_cavities(atoms)
            
            for cavity in cavities:
                binding_site = self._analyze_cavity(cavity, chain)
                if binding_site:
                    binding_sites.append(binding_site)
        
        return binding_sites
    
    def _find_cavities(self, atoms: List, grid_spacing: float = 1.0) -> List[Dict]:
        """Find cavities using a simplified grid-based approach"""
        # Get bounding box
        coords = np.array([atom.coord for atom in atoms])
        min_coords = np.min(coords, axis=0) - 5.0
        max_coords = np.max(coords, axis=0) + 5.0
        
        # Create grid points
        x_points = np.arange(min_coords[0], max_coords[0], grid_spacing)
        y_points = np.arange(min_coords[1], max_coords[1], grid_spacing)
        z_points = np.arange(min_coords[2], max_coords[2], grid_spacing)
        
        cavity_points = []
        
        # Check each grid point
        for x in x_points[::2]:  # Subsample for speed
            for y in y_points[::2]:
                for z in z_points[::2]:
                    point = np.array([x, y, z])
                    
                    # Check if point is in a cavity
                    if self._is_cavity_point(point, coords):
                        cavity_points.append(point)
        
        if len(cavity_points) < 5:
            return []
        
        # Cluster cavity points
        cavity_points = np.array(cavity_points)
        clustering = DBSCAN(eps=2.0, min_samples=5).fit(cavity_points)
        
        cavities = []
        for cluster_id in set(clustering.labels_):
            if cluster_id != -1:
                cluster_points = cavity_points[clustering.labels_ == cluster_id]
                if len(cluster_points) >= 5:
                    cavities.append({
                        'points': cluster_points,
                        'center': np.mean(cluster_points, axis=0),
                        'volume': len(cluster_points) * (grid_spacing ** 3)
                    })
        
        return cavities
    
    def _is_cavity_point(self, point: np.ndarray, atom_coords: np.ndarray, 
                        min_distance: float = 1.5, max_distance: float = 8.0) -> bool:
        """Check if a point is in a cavity"""
        distances = np.linalg.norm(atom_coords - point, axis=1)
        min_dist = np.min(distances)
        
        # Point should be away from atoms but not too far (interior cavity)
        if min_dist < min_distance or min_dist > max_distance:
            return False
        
        # Check if point is surrounded by atoms (cavity-like)
        nearby_atoms = np.sum(distances < max_distance)
        return nearby_atoms >= 8  # Minimum surrounding atoms
    
    def _analyze_cavity(self, cavity: Dict, chain) -> Optional[BindingSite]:
        """Analyze a cavity to determine if it's a binding site"""
        center = cavity['center']
        volume = cavity['volume']
        
        # Get nearby residues
        nearby_residues = []
        for residue in chain:
            try:
                ca_coord = residue['CA'].coord
                distance = np.linalg.norm(ca_coord - center)
                if distance < 8.0:  # Within 8Å
                    nearby_residues.append(residue)
            except KeyError:
                continue
        
        if len(nearby_residues) < 3:
            return None
        
        # Calculate properties
        hydrophobicity = self._calculate_hydrophobicity(nearby_residues)
        charge = self._calculate_charge(nearby_residues)
        surface_area = self._estimate_surface_area(cavity)
        
        # Calculate binding likelihood
        binding_likelihood = self._calculate_binding_likelihood(
            volume, hydrophobicity, charge, len(nearby_residues)
        )
        
        if binding_likelihood < 0.3:
            return None
        
        # Classify binding site
        site_classification = self._classify_binding_site(hydrophobicity, charge, volume)
        
        return BindingSite(
            residues=[(res.get_resname(), res.get_id()[1], chain.id) for res in nearby_residues],
            center_coords=center,
            volume=volume,
            surface_area=surface_area,
            hydrophobicity=hydrophobicity,
            charge=charge,
            binding_likelihood=binding_likelihood,
            site_classification=site_classification
        )
    
    def _calculate_hydrophobicity(self, residues: List) -> float:
        """Calculate hydrophobicity score for residues"""
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET', 'PRO'}
        hydrophobic_count = sum(1 for res in residues if res.get_resname() in hydrophobic_residues)
        return hydrophobic_count / len(residues) if residues else 0.0
    
    def _calculate_charge(self, residues: List) -> float:
        """Calculate net charge of residues"""
        positive_residues = {'LYS', 'ARG', 'HIS'}
        negative_residues = {'ASP', 'GLU'}
        
        positive_count = sum(1 for res in residues if res.get_resname() in positive_residues)
        negative_count = sum(1 for res in residues if res.get_resname() in negative_residues)
        
        return (positive_count - negative_count) / len(residues) if residues else 0.0
    
    def _estimate_surface_area(self, cavity: Dict) -> float:
        """Estimate surface area of cavity"""
        # Simplified estimation based on volume
        volume = cavity['volume']
        radius = (3 * volume / (4 * np.pi)) ** (1/3)
        return 4 * np.pi * radius ** 2
    
    def _calculate_binding_likelihood(self, volume: float, hydrophobicity: float, 
                                    charge: float, residue_count: int) -> float:
        """Calculate likelihood that this is a functional binding site"""
        score = 0.0
        
        # Volume score (prefer moderate volumes)
        if 100 < volume < 2000:
            score += 0.3
        elif volume < 100:
            score += 0.1
        
        # Hydrophobicity score (moderate hydrophobicity is good)
        if 0.2 < hydrophobicity < 0.8:
            score += 0.2
        
        # Charge score (some charge is good for binding)
        if abs(charge) > 0.1:
            score += 0.2
        
        # Residue count score
        if residue_count >= 5:
            score += 0.3
        
        return min(1.0, score)
    
    def _classify_binding_site(self, hydrophobicity: float, charge: float, volume: float) -> str:
        """Classify the type of binding site"""
        if hydrophobicity > 0.6:
            return "Hydrophobic binding site"
        elif charge > 0.3:
            return "Positively charged binding site"
        elif charge < -0.3:
            return "Negatively charged binding site"
        elif volume > 1000:
            return "Large binding cavity"
        elif volume < 200:
            return "Small binding pocket"
        else:
            return "Mixed binding site"

class BiosyntheticEnzymeAnalyzer:
    """Specialized analyzer for biosynthetic enzymes in natural product pathways"""
    
    def __init__(self):
        # Biosynthetic enzyme families with specific motifs and characteristics
        self.biosynthetic_families = {
            'polyketide_synthase_typeI': {
                'domains': ['KS', 'AT', 'ACP', 'TE', 'KR', 'DH', 'ER'],
                'motifs': {
                    'KS_active_site': ['CYS', 'HIS', 'HIS'],  # Ketosynthase catalytic triad
                    'AT_active_site': ['SER', 'HIS'],  # Acyltransferase
                    'KR_active_site': ['SER', 'TYR', 'LYS']  # Ketoreductase
                },
                'cofactors': ['CoA', 'NADPH', 'malonyl-CoA'],
                'description': 'Type I Polyketide Synthase',
                'pathway_role': 'Polyketide backbone assembly',
                'products': ['Polyketides', 'Macrolides', 'Polyenes']
            },
            'polyketide_synthase_typeII': {
                'domains': ['KS_alpha', 'KS_beta', 'ACP'],
                'motifs': {
                    'KS_active_site': ['CYS', 'HIS', 'ASN']
                },
                'description': 'Type II Polyketide Synthase',
                'pathway_role': 'Aromatic polyketide synthesis',
                'products': ['Anthracyclines', 'Tetracyclines', 'Aureolic acids']
            },
            'polyketide_synthase_typeIII': {
                'motifs': {
                    'chalcone_synthase_active': ['CYS', 'HIS', 'ASN']
                },
                'description': 'Type III Polyketide Synthase',
                'pathway_role': 'Plant-like polyketide synthesis',
                'products': ['Chalcones', 'Pyrones', 'Resorcinols']
            },
            'nonribosomal_peptide_synthetase': {
                'domains': ['A', 'PCP', 'C', 'TE', 'E'],
                'motifs': {
                    'adenylation_domain': ['TRP', 'LYS', 'ASP'],  # A domain signature
                    'condensation_domain': ['HIS', 'ASP'],  # C domain catalytic dyad
                    'thioesterase_domain': ['SER', 'HIS', 'ASP']  # TE domain triad
                },
                'cofactors': ['ATP', 'CoA', 'Mg2+'],
                'description': 'Nonribosomal Peptide Synthetase',
                'pathway_role': 'Peptide antibiotic synthesis',
                'products': ['Cyclic peptides', 'Lipopeptides', 'Siderophores']
            },
            'terpene_synthase': {
                'motifs': {
                    'class_I_active_site': ['ASP', 'ASP', 'GLU'],  # DDXXE motif
                    'class_II_active_site': ['ASP', 'CYS', 'TYR']  # DXDD motif region
                },
                'cofactors': ['Mg2+', 'Mn2+', 'FPP', 'GGPP'],
                'description': 'Terpene Synthase',
                'pathway_role': 'Terpenoid biosynthesis',
                'products': ['Monoterpenes', 'Sesquiterpenes', 'Diterpenes']
            },
            'cytochrome_p450': {
                'motifs': {
                    'heme_binding': ['CYS'],  # Cysteine heme ligand
                    'oxygen_binding': ['THR', 'ASP']  # I-helix region
                },
                'cofactors': ['Heme', 'NADPH', 'O2'],
                'description': 'Cytochrome P450 Monooxygenase',
                'pathway_role': 'Oxidative tailoring',
                'products': ['Hydroxylated metabolites', 'Epoxides', 'Dealkylated products']
            },
            'methyltransferase': {
                'motifs': {
                    'sam_binding': ['GLY', 'GLY', 'GLY'],  # Rossmann fold
                    'catalytic_base': ['ASP', 'GLU', 'HIS']
                },
                'cofactors': ['SAM', 'SAH'],
                'description': 'S-Adenosyl Methionine Methyltransferase',
                'pathway_role': 'Methylation modification',
                'products': ['O-methylated compounds', 'N-methylated compounds', 'C-methylated compounds']
            },
            'halogenase': {
                'motifs': {
                    'flavin_binding': ['GLY', 'TYR', 'SER'],
                    'halide_binding': ['TRP', 'PHE', 'LYS']
                },
                'cofactors': ['FAD', 'NADH', 'Cl-', 'Br-'],
                'description': 'Flavin-dependent Halogenase',
                'pathway_role': 'Halogenation',
                'products': ['Chlorinated metabolites', 'Brominated metabolites']
            },
            'glycosyltransferase': {
                'motifs': {
                    'udp_binding': ['ASP', 'ASP', 'ASP'],  # DXD motif
                    'acceptor_binding': ['HIS', 'ASP']
                },
                'cofactors': ['UDP-sugar', 'Mn2+', 'Mg2+'],
                'description': 'Glycosyltransferase',
                'pathway_role': 'Glycosylation modification',
                'products': ['Glycosides', 'Nucleoside antibiotics']
            },
            'aminotransferase': {
                'motifs': {
                    'plp_binding': ['LYS'],  # PLP Schiff base
                    'substrate_binding': ['TYR', 'PHE', 'ARG']
                },
                'cofactors': ['PLP', 'PMP'],
                'description': 'Pyridoxal Phosphate Aminotransferase',
                'pathway_role': 'Amino group transfer',
                'products': ['Amino acids', 'Aminoglycosides']
            },
            'flavin_monooxygenase': {
                'motifs': {
                    'flavin_binding': ['GLY', 'GLY', 'GLY'],  # Rossmann fold
                    'nadph_binding': ['GLY', 'ALA', 'GLY']
                },
                'cofactors': ['FAD', 'FMN', 'NADPH'],
                'description': 'Flavin Monooxygenase',
                'pathway_role': 'Oxidative modifications',
                'products': ['Hydroxylated products', 'N-oxides', 'Sulfoxides']
            }
        }
        
        # Bacterial and fungal specific features
        self.microbial_features = {
            'bacterial_specific': {
                'iron_sulfur_clusters': ['CYS', 'CYS', 'CYS', 'CYS'],
                'molybdenum_cofactor': ['CYS', 'SER', 'GLY'],
                'cobalamin_binding': ['HIS', 'ASP', 'MET'],
                'pqq_binding': ['TRP', 'TYR', 'ASP']
            },
            'fungal_specific': {
                'nonheme_iron': ['HIS', 'ASP', 'GLU'],
                'copper_binding': ['HIS', 'MET', 'CYS'],
                'flavin_binding_fungal': ['TYR', 'PHE', 'ARG']
            }
        }
        
        # Product classification patterns
        self.product_patterns = {
            'polyketides': ['fatty_acid_like', 'aromatic_rings', 'lactone_rings'],
            'peptides': ['amino_acid_chains', 'cyclic_structures', 'modified_residues'],
            'terpenoids': ['isoprene_units', 'cyclic_structures', 'hydroxyl_groups'],
            'alkaloids': ['nitrogen_heterocycles', 'complex_rings'],
            'phenylpropanoids': ['aromatic_rings', 'propyl_chains'],
            'fatty_acids': ['long_carbon_chains', 'carboxyl_groups']
        }

class EnzymeFunctionPredictor:
    """Predicts enzyme function based on structural features"""
    
    def __init__(self):
        # EC class patterns (simplified)
        self.ec_patterns = {
            'EC 1': {  # Oxidoreductases
                'catalytic_residues': ['CYS', 'HIS', 'TYR'],
                'cofactors': ['FAD', 'NAD', 'NADP', 'FMN'],
                'metal_binding': True,
                'description': 'Oxidation-reduction reactions'
            },
            'EC 2': {  # Transferases
                'catalytic_residues': ['SER', 'THR', 'TYR', 'LYS'],
                'cofactors': ['ATP', 'SAM', 'CoA'],
                'description': 'Transfer of functional groups'
            },
            'EC 3': {  # Hydrolases
                'catalytic_residues': ['SER', 'HIS', 'ASP', 'GLU', 'CYS'],
                'water_binding': True,
                'description': 'Hydrolysis reactions'
            },
            'EC 4': {  # Lyases
                'catalytic_residues': ['LYS', 'ARG', 'HIS', 'CYS'],
                'description': 'Addition or removal of groups to form double bonds'
            },
            'EC 5': {  # Isomerases
                'catalytic_residues': ['SER', 'CYS', 'HIS'],
                'description': 'Intramolecular rearrangements'
            },
            'EC 6': {  # Ligases
                'catalytic_residues': ['LYS', 'ARG', 'ASP', 'GLU'],
                'cofactors': ['ATP', 'GTP'],
                'metal_binding': True,
                'description': 'Formation of bonds with ATP cleavage'
            }
        }
        
        # Initialize biosynthetic analyzer
        self.biosynthetic_analyzer = BiosyntheticEnzymeAnalyzer()
        
        # Specific enzyme families
        self.enzyme_families = {
            'serine_protease': {
                'motif': ['SER', 'HIS', 'ASP'],
                'ec_class': 'EC 3.4',
                'description': 'Serine endopeptidase',
                'mechanism': 'Nucleophilic attack by serine hydroxyl'
            },
            'cysteine_protease': {
                'motif': ['CYS', 'HIS'],
                'ec_class': 'EC 3.4',
                'description': 'Cysteine endopeptidase',
                'mechanism': 'Nucleophilic attack by cysteine thiol'
            },
            'metallopeptidase': {
                'motif': ['HIS', 'GLU'],
                'metal_required': True,
                'ec_class': 'EC 3.4',
                'description': 'Metal-dependent peptidase'
            },
            'kinase': {
                'motif': ['LYS', 'ASP'],
                'cofactors': ['ATP', 'Mg2+'],
                'ec_class': 'EC 2.7',
                'description': 'Protein kinase'
            },
            'phosphatase': {
                'motif': ['ASP', 'ASP'],
                'ec_class': 'EC 3.1',
                'description': 'Phosphatase'
            }
        }
    
    def predict_enzyme_function(self, catalytic_sites: List[CatalyticSite], 
                              binding_sites: List[BindingSite],
                              structure) -> List[EnzymeFunction]:
        """Predict enzyme functions based on structural analysis"""
        functions = []
        
        # Analyze for biosynthetic enzymes first (more specific)
        biosynthetic_functions = self._predict_biosynthetic_functions(structure, catalytic_sites, binding_sites)
        functions.extend(biosynthetic_functions)
        
        # Analyze catalytic sites for general enzyme functions
        for site in catalytic_sites:
            function = self._analyze_catalytic_site_function(site)
            if function and not self._is_duplicate_function(function, functions):
                functions.append(function)
        
        # General functional prediction based on overall features
        general_function = self._predict_general_function(catalytic_sites, binding_sites, structure)
        if general_function and not self._is_duplicate_function(general_function, functions):
            functions.append(general_function)
        
        return functions
    
    def _analyze_catalytic_site_function(self, site: CatalyticSite) -> Optional[EnzymeFunction]:
        """Analyze a specific catalytic site for function"""
        residue_names = [res[0] for res in site.residues]
        
        # Check against known enzyme families
        for family_name, family_data in self.enzyme_families.items():
            motif = family_data['motif']
            if all(res in residue_names for res in motif):
                return EnzymeFunction(
                    ec_number=family_data.get('ec_class'),
                    enzyme_class=family_name.replace('_', ' ').title(),
                    function_description=family_data['description'],
                    confidence=0.8,
                    supporting_features=[f"Contains {family_name} motif: {', '.join(motif)}"],
                    catalytic_mechanism=family_data.get('mechanism', 'Unknown mechanism'),
                    substrate_prediction=self._predict_substrates(family_name)
                )
        
        return None
    
    def _predict_general_function(self, catalytic_sites: List[CatalyticSite],
                                binding_sites: List[BindingSite], structure) -> Optional[EnzymeFunction]:
        """Predict general enzyme function based on all features"""
        if not catalytic_sites:
            return None
        
        # Collect all catalytic residues
        all_catalytic_residues = []
        for site in catalytic_sites:
            all_catalytic_residues.extend([res[0] for res in site.residues])
        
        catalytic_counter = Counter(all_catalytic_residues)
        
        # Predict EC class
        ec_scores = {}
        for ec_class, patterns in self.ec_patterns.items():
            score = 0.0
            required_residues = patterns['catalytic_residues']
            
            for residue in required_residues:
                if residue in catalytic_counter:
                    score += catalytic_counter[residue]
            
            ec_scores[ec_class] = score
        
        if max(ec_scores.values()) > 0:
            predicted_ec = max(ec_scores, key=ec_scores.get)
            confidence = min(0.9, ec_scores[predicted_ec] / 10.0)
            
            return EnzymeFunction(
                ec_number=predicted_ec,
                enzyme_class=predicted_ec.replace('EC ', 'Enzyme Class '),
                function_description=self.ec_patterns[predicted_ec]['description'],
                confidence=confidence,
                supporting_features=[
                    f"Catalytic residues: {', '.join(catalytic_counter.keys())}",
                    f"Number of catalytic sites: {len(catalytic_sites)}",
                    f"Number of binding sites: {len(binding_sites)}"
                ],
                catalytic_mechanism="General enzymatic catalysis",
                substrate_prediction=["Unknown substrates"]
            )
        
        return None
    
    def _predict_substrates(self, enzyme_family: str) -> List[str]:
        """Predict likely substrates based on enzyme family"""
        substrate_predictions = {
            'serine_protease': ['Peptides', 'Proteins', 'Peptide bonds'],
            'cysteine_protease': ['Peptides', 'Proteins'],
            'metallopeptidase': ['Peptides', 'Proteins', 'Metalloprotein substrates'],
            'kinase': ['Proteins', 'Nucleotides', 'ATP', 'Serine/threonine/tyrosine residues'],
            'phosphatase': ['Phosphoproteins', 'Phosphate esters']
        }
        
        return substrate_predictions.get(enzyme_family, ['Unknown substrates'])
    
    def _predict_biosynthetic_functions(self, structure, catalytic_sites: List[CatalyticSite], 
                                      binding_sites: List[BindingSite]) -> List[EnzymeFunction]:
        """Predict biosynthetic enzyme functions specific to natural product pathways"""
        biosynthetic_functions = []
        
        # Extract all residues from the structure
        all_residues = []
        for chain in structure[0]:
            for residue in chain:
                all_residues.append(residue)
        
        # Analyze for each biosynthetic enzyme family
        for family_name, family_data in self.biosynthetic_analyzer.biosynthetic_families.items():
            confidence, supporting_features = self._analyze_biosynthetic_family(
                all_residues, family_data, catalytic_sites, binding_sites
            )
            
            if confidence > 0.3:  # Threshold for biosynthetic enzyme prediction
                # Determine EC number based on family
                ec_number = self._get_biosynthetic_ec_number(family_name)
                
                # Predict pathway role and products
                pathway_role = family_data.get('pathway_role', 'Secondary metabolite biosynthesis')
                products = family_data.get('products', ['Secondary metabolites'])
                
                # Determine pathway type
                pathway_type = self._determine_pathway_type(family_name)
                
                # Create enzyme function
                function = EnzymeFunction(
                    ec_number=ec_number,
                    enzyme_class=f"Biosynthetic: {family_data['description']}",
                    function_description=f"{pathway_role} - {family_data['description']}",
                    confidence=confidence,
                    supporting_features=supporting_features,
                    catalytic_mechanism=self._get_biosynthetic_mechanism(family_name),
                    substrate_prediction=self._predict_biosynthetic_substrates(family_name, family_data),
                    is_biosynthetic=True,
                    pathway_type=pathway_type,
                    microbial_origin="Unknown"
                )
                
                biosynthetic_functions.append(function)
        
        # Check for microbial-specific features
        microbial_functions = self._analyze_microbial_features(all_residues, structure)
        biosynthetic_functions.extend(microbial_functions)
        
        return biosynthetic_functions
    
    def _analyze_biosynthetic_family(self, residues: List, family_data: Dict, 
                                   catalytic_sites: List[CatalyticSite], 
                                   binding_sites: List[BindingSite]) -> Tuple[float, List[str]]:
        """Analyze structure for specific biosynthetic enzyme family features"""
        confidence = 0.0
        supporting_features = []
        
        # Get residue names for analysis
        residue_names = [res.get_resname() for res in residues]
        residue_counter = Counter(residue_names)
        
        # Check for family-specific motifs
        if 'motifs' in family_data:
            for motif_name, motif_residues in family_data['motifs'].items():
                motif_score = self._check_motif_presence(residues, motif_residues)
                if motif_score > 0.5:
                    confidence += motif_score * 0.3
                    supporting_features.append(f"Contains {motif_name}: {', '.join(motif_residues)}")
        
        # Check for domain architecture (for modular enzymes like PKS/NRPS)
        if 'domains' in family_data:
            domain_score = self._analyze_domain_architecture(residues, family_data['domains'])
            if domain_score > 0.3:
                confidence += domain_score * 0.4
                supporting_features.append(f"Domain architecture consistent with {family_data['description']}")
        
        # Check for cofactor binding sites
        if 'cofactors' in family_data:
            cofactor_score = self._analyze_cofactor_binding(residues, family_data['cofactors'])
            if cofactor_score > 0.2:
                confidence += cofactor_score * 0.2
                supporting_features.append(f"Cofactor binding sites: {', '.join(family_data['cofactors'])}")
        
        # Check for characteristic amino acid composition
        composition_score = self._analyze_aa_composition_biosynthetic(residue_counter, family_data)
        confidence += composition_score * 0.1
        
        return min(1.0, confidence), supporting_features
    
    def _check_motif_presence(self, residues: List, motif_residues: List[str]) -> float:
        """Check for presence of specific motifs in the structure"""
        # Get coordinates of motif residues
        motif_coords = []
        
        for res in residues:
            if res.get_resname() in motif_residues:
                try:
                    if 'CA' in res:
                        motif_coords.append((res.get_resname(), res['CA'].coord, res.get_id()[1]))
                except KeyError:
                    continue
        
        if len(motif_coords) < len(motif_residues):
            return 0.0
        
        # Look for spatial clusters of motif residues
        best_cluster_score = 0.0
        
        for i in range(len(motif_coords)):
            for j in range(i+1, len(motif_coords)):
                for k in range(j+1, len(motif_coords)):
                    if len(motif_residues) <= 3:
                        # Check if these three residues could form the motif
                        coords = [motif_coords[i][1], motif_coords[j][1], motif_coords[k][1]]
                        residue_names = [motif_coords[i][0], motif_coords[j][0], motif_coords[k][0]]
                        
                        # Check spatial proximity
                        max_distance = max(np.linalg.norm(coords[0] - coords[1]),
                                         np.linalg.norm(coords[1] - coords[2]),
                                         np.linalg.norm(coords[0] - coords[2]))
                        
                        if max_distance < 12.0:  # Within reasonable distance
                            # Check if residue types match motif
                            motif_match = sum(1 for req_res in motif_residues if req_res in residue_names)
                            cluster_score = motif_match / len(motif_residues)
                            
                            if cluster_score > best_cluster_score:
                                best_cluster_score = cluster_score
        
        return best_cluster_score
    
    def _analyze_domain_architecture(self, residues: List, expected_domains: List[str]) -> float:
        """Analyze domain architecture for modular enzymes"""
        # This is a simplified domain analysis
        # In reality, you'd use domain prediction tools like InterPro
        
        protein_length = len(residues)
        
        # Estimate domains based on protein length and composition
        domain_score = 0.0
        
        # Large proteins are more likely to be modular
        if protein_length > 1000:  # Typical for Type I PKS/NRPS
            domain_score += 0.5
        elif protein_length > 500:
            domain_score += 0.3
        
        # Check for signature sequences
        residue_names = [res.get_resname() for res in residues]
        
        # Look for conserved cysteine (phosphopantetheine attachment)
        if 'CYS' in residue_names:
            cys_count = residue_names.count('CYS')
            if cys_count >= 2:  # Multiple cysteines suggest ACP domains
                domain_score += 0.3
        
        # Look for GGXS motif (simplified)
        if 'GLY' in residue_names and 'SER' in residue_names:
            domain_score += 0.2
        
        return min(1.0, domain_score)
    
    def _analyze_cofactor_binding(self, residues: List, cofactors: List[str]) -> float:
        """Analyze potential cofactor binding based on residue patterns"""
        cofactor_score = 0.0
        residue_names = [res.get_resname() for res in residues]
        
        for cofactor in cofactors:
            if cofactor in ['CoA', 'Acetyl-CoA', 'malonyl-CoA']:
                # Look for CoA binding patterns
                if 'GLY' in residue_names and 'SER' in residue_names and 'THR' in residue_names:
                    cofactor_score += 0.3
            elif cofactor in ['NADPH', 'NADH']:
                # Look for Rossmann fold patterns
                gly_count = residue_names.count('GLY')
                if gly_count >= 6:  # Rossmann folds are glycine-rich
                    cofactor_score += 0.3
            elif cofactor in ['SAM', 'SAH']:
                # S-adenosyl methionine binding
                if 'GLY' in residue_names and 'ASP' in residue_names:
                    cofactor_score += 0.3
            elif cofactor in ['FAD', 'FMN']:
                # Flavin binding
                if 'TYR' in residue_names and 'PHE' in residue_names:
                    cofactor_score += 0.3
        
        return min(1.0, cofactor_score)
    
    def _analyze_aa_composition_biosynthetic(self, residue_counter: Counter, family_data: Dict) -> float:
        """Analyze amino acid composition specific to biosynthetic enzymes"""
        composition_score = 0.0
        total_residues = sum(residue_counter.values())
        
        if total_residues == 0:
            return 0.0
        
        # Biosynthetic enzymes often have specific amino acid preferences
        biosynthetic_preferences = {
            'GLY': 0.08,  # Flexibility for conformational changes
            'ALA': 0.08,  # Small, non-polar
            'VAL': 0.06,  # Hydrophobic cores
            'LEU': 0.09,  # Hydrophobic
            'ILE': 0.05,  # Hydrophobic
            'SER': 0.07,  # Hydroxyl groups for binding
            'THR': 0.05,  # Hydroxyl groups
            'CYS': 0.03,  # Metal coordination, disulfides
            'MET': 0.02,  # Hydrophobic, sulfur
            'PHE': 0.04,  # Aromatic
            'TYR': 0.03,  # Aromatic, hydroxyl
            'TRP': 0.01,  # Aromatic, bulky
            'HIS': 0.02,  # Metal coordination, catalysis
            'LYS': 0.06,  # Positive charge
            'ARG': 0.05,  # Positive charge
            'ASP': 0.05,  # Negative charge
            'GLU': 0.06,  # Negative charge
            'ASN': 0.04,  # Polar
            'GLN': 0.04,  # Polar
            'PRO': 0.05   # Structural
        }
        
        # Calculate deviation from expected composition
        for aa, expected_freq in biosynthetic_preferences.items():
            actual_freq = residue_counter.get(aa, 0) / total_residues
            deviation = abs(actual_freq - expected_freq)
            if deviation < 0.02:  # Close to expected
                composition_score += 0.1
        
        return min(1.0, composition_score)
    
    def _analyze_microbial_features(self, residues: List, structure) -> List[EnzymeFunction]:
        """Analyze for bacterial and fungal specific enzyme features"""
        microbial_functions = []
        residue_names = [res.get_resname() for res in residues]
        
        # Check for bacterial-specific features
        bacterial_score = 0.0
        bacterial_features = []
        
        # Iron-sulfur clusters (common in bacterial enzymes)
        cys_count = residue_names.count('CYS')
        if cys_count >= 4:
            bacterial_score += 0.3
            bacterial_features.append("Multiple cysteine residues suggest iron-sulfur clusters")
        
        # Molybdenum cofactor binding
        if 'CYS' in residue_names and 'SER' in residue_names and 'GLY' in residue_names:
            mo_pattern_score = self._check_molybdenum_pattern(residues)
            if mo_pattern_score > 0.3:
                bacterial_score += 0.4
                bacterial_features.append("Molybdenum cofactor binding pattern")
        
        # Cobalamin (B12) binding
        if 'HIS' in residue_names and 'ASP' in residue_names and 'MET' in residue_names:
            bacterial_score += 0.2
            bacterial_features.append("Potential cobalamin binding")
        
        if bacterial_score > 0.4:
            function = EnzymeFunction(
                ec_number="EC 1.-.-",
                enzyme_class="Bacterial Oxidoreductase",
                function_description="Bacterial-specific redox enzyme with metal cofactors",
                confidence=bacterial_score,
                supporting_features=bacterial_features,
                catalytic_mechanism="Metal-mediated electron transfer",
                substrate_prediction=["Organic substrates", "Metal ions", "Electron donors/acceptors"],
                is_biosynthetic=True,
                pathway_type="Bacterial metabolism",
                microbial_origin="Bacterial"
            )
            microbial_functions.append(function)
        
        # Check for fungal-specific features
        fungal_score = 0.0
        fungal_features = []
        
        # Nonheme iron enzymes (common in fungi)
        if 'HIS' in residue_names and 'ASP' in residue_names and 'GLU' in residue_names:
            fungal_score += 0.3
            fungal_features.append("Nonheme iron binding pattern")
        
        # Copper-containing enzymes
        if 'HIS' in residue_names and 'MET' in residue_names:
            copper_score = self._check_copper_binding_pattern(residues)
            if copper_score > 0.3:
                fungal_score += 0.4
                fungal_features.append("Copper binding site")
        
        if fungal_score > 0.4:
            function = EnzymeFunction(
                ec_number="EC 1.-.-",
                enzyme_class="Fungal Oxidoreductase",
                function_description="Fungal-specific enzyme with copper or iron cofactors",
                confidence=fungal_score,
                supporting_features=fungal_features,
                catalytic_mechanism="Copper/iron-mediated oxidation",
                substrate_prediction=["Aromatic compounds", "Phenolic substrates", "Oxygen"],
                is_biosynthetic=True,
                pathway_type="Fungal metabolism",
                microbial_origin="Fungal"
            )
            microbial_functions.append(function)
        
        return microbial_functions
    
    def _determine_pathway_type(self, family_name: str) -> str:
        """Determine the biosynthetic pathway type"""
        pathway_mapping = {
            'polyketide_synthase_typeI': 'Type I Polyketide',
            'polyketide_synthase_typeII': 'Type II Polyketide', 
            'polyketide_synthase_typeIII': 'Type III Polyketide',
            'nonribosomal_peptide_synthetase': 'NRPS Peptide',
            'terpene_synthase': 'Terpenoid',
            'cytochrome_p450': 'P450 Oxidation',
            'methyltransferase': 'Methylation',
            'halogenase': 'Halogenation',
            'glycosyltransferase': 'Glycosylation',
            'aminotransferase': 'Amino Transfer',
            'flavin_monooxygenase': 'Flavin Oxidation'
        }
        return pathway_mapping.get(family_name, 'Unknown pathway')
    
    def _check_molybdenum_pattern(self, residues: List) -> float:
        """Check for molybdenum cofactor binding patterns"""
        # Look for specific spatial arrangements of Cys, Ser, Gly
        mo_residues = []
        for res in residues:
            if res.get_resname() in ['CYS', 'SER', 'GLY']:
                try:
                    mo_residues.append((res.get_resname(), res['CA'].coord))
                except KeyError:
                    continue
        
        if len(mo_residues) < 3:
            return 0.0
        
        # Check for spatial clustering
        max_score = 0.0
        for i in range(len(mo_residues)):
            for j in range(i+1, len(mo_residues)):
                for k in range(j+1, len(mo_residues)):
                    coords = [mo_residues[i][1], mo_residues[j][1], mo_residues[k][1]]
                    max_dist = max(np.linalg.norm(coords[0] - coords[1]),
                                  np.linalg.norm(coords[1] - coords[2]),
                                  np.linalg.norm(coords[0] - coords[2]))
                    
                    if max_dist < 10.0:
                        score = 1.0 - (max_dist / 10.0)
                        max_score = max(max_score, score)
        
        return max_score
    
    def _check_copper_binding_pattern(self, residues: List) -> float:
        """Check for copper binding patterns typical in fungal enzymes"""
        copper_residues = []
        for res in residues:
            if res.get_resname() in ['HIS', 'MET', 'CYS']:
                try:
                    copper_residues.append((res.get_resname(), res['CA'].coord))
                except KeyError:
                    continue
        
        if len(copper_residues) < 2:
            return 0.0
        
        # Check for close histidine pairs (typical in copper enzymes)
        his_pairs = []
        for i, (res1_name, coord1) in enumerate(copper_residues):
            for j, (res2_name, coord2) in enumerate(copper_residues[i+1:], i+1):
                if res1_name == 'HIS' and res2_name == 'HIS':
                    distance = np.linalg.norm(coord1 - coord2)
                    if distance < 8.0:
                        his_pairs.append(distance)
        
        if his_pairs:
            return 0.8  # Strong evidence for copper binding
        
        return 0.3  # Weaker evidence
    
    def _get_biosynthetic_ec_number(self, family_name: str) -> Optional[str]:
        """Get EC number for biosynthetic enzyme families"""
        ec_mapping = {
            'polyketide_synthase_typeI': 'EC 2.3.1.-',
            'polyketide_synthase_typeII': 'EC 2.3.1.-',
            'polyketide_synthase_typeIII': 'EC 2.3.1.-',
            'nonribosomal_peptide_synthetase': 'EC 6.3.2.-',
            'terpene_synthase': 'EC 4.2.3.-',
            'cytochrome_p450': 'EC 1.14.13.-',
            'methyltransferase': 'EC 2.1.1.-',
            'halogenase': 'EC 1.14.19.-',
            'glycosyltransferase': 'EC 2.4.1.-',
            'aminotransferase': 'EC 2.6.1.-',
            'flavin_monooxygenase': 'EC 1.14.13.-'
        }
        return ec_mapping.get(family_name)
    
    def _get_biosynthetic_mechanism(self, family_name: str) -> str:
        """Get catalytic mechanism for biosynthetic enzymes"""
        mechanisms = {
            'polyketide_synthase_typeI': 'Sequential decarboxylative condensation of malonyl-CoA units',
            'polyketide_synthase_typeII': 'Iterative condensation with separate enzymes for each step',
            'polyketide_synthase_typeIII': 'Chalcone synthase-like condensation mechanism',
            'nonribosomal_peptide_synthetase': 'Amino acid activation, thioesterification, and condensation',
            'terpene_synthase': 'Isoprenoid cyclization via carbocation intermediates',
            'cytochrome_p450': 'Heme-mediated hydroxylation using oxygen and electrons',
            'methyltransferase': 'SAM-dependent methyl transfer to acceptor substrates',
            'halogenase': 'Flavin-mediated halogenation with molecular halides',
            'glycosyltransferase': 'Sugar transfer from nucleotide sugar donors',
            'aminotransferase': 'PLP-mediated amino group transfer',
            'flavin_monooxygenase': 'Flavin-mediated oxygenation reactions'
        }
        return mechanisms.get(family_name, 'Unknown biosynthetic mechanism')
    
    def _predict_biosynthetic_substrates(self, family_name: str, family_data: Dict) -> List[str]:
        """Predict substrates for biosynthetic enzymes"""
        # Use the products from family_data as potential substrates/products
        products = family_data.get('products', [])
        cofactors = family_data.get('cofactors', [])
        
        # Combine and add common biosynthetic substrates
        substrates = []
        substrates.extend(products)
        substrates.extend(cofactors)
        
        # Add family-specific substrates
        family_substrates = {
            'polyketide_synthase_typeI': ['Malonyl-CoA', 'Acetyl-CoA', 'Propionyl-CoA'],
            'nonribosomal_peptide_synthetase': ['Amino acids', 'ATP', 'CoA'],
            'terpene_synthase': ['Farnesyl pyrophosphate', 'Geranylgeranyl pyrophosphate'],
            'cytochrome_p450': ['Organic substrates', 'NADPH', 'Oxygen'],
            'methyltransferase': ['S-adenosyl methionine', 'Acceptor molecules'],
            'halogenase': ['Tryptophan', 'Aromatic compounds', 'Halide ions']
        }
        
        specific_substrates = family_substrates.get(family_name, [])
        substrates.extend(specific_substrates)
        
        return list(set(substrates))  # Remove duplicates
    
    def _is_duplicate_function(self, new_function: EnzymeFunction, existing_functions: List[EnzymeFunction]) -> bool:
        """Check if a function is already in the list to avoid duplicates"""
        for existing in existing_functions:
            if (existing.enzyme_class == new_function.enzyme_class or 
                existing.ec_number == new_function.ec_number):
                return True
        return False

class StructuralAnalyzer:
    """Analyzes structural properties and features"""

    def __init__(self, output_dir: str = "pdb_analysis_results"):
        self.output_dir = Path(output_dir)
        self.dssp_pdb_dir = self.output_dir / "dsspPDB"
        self.reports_dir = self.output_dir / "individual_reports"
        self.dssp_available = self._check_dssp()

        self.dssp_pdb_dir.mkdir(parents=True, exist_ok=True)
        self.reports_dir.mkdir(parents=True, exist_ok=True)

    def _check_dssp(self) -> bool:
        try:
            result = subprocess.run(['mkdssp', '--version'], capture_output=True, timeout=5)
            if result.returncode == 0:
                logger.info("Found mkdssp (DSSP v4.x)")
                return True
        except Exception:
            pass
        logger.warning("DSSP not found or not working. Secondary structure analysis will use fallback phi/psi method.")
        return False

    def _write_clean_pdb(self, structure, pdb_filename):
        output_path = self.dssp_pdb_dir / Path(pdb_filename).name

        # Write temporary raw PDB
        tmp = output_path.with_suffix(".tmp")
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(tmp))

        # Read contents, add HEADER line, write to final output
        with open(tmp, "r") as f:
            lines = f.readlines()

        header = f"HEADER    {Path(pdb_filename).stem[:66]:<66}\n"
        with open(output_path, "w") as f:
            f.write(header)
            f.writelines(lines)

        tmp.unlink()  # Clean up
        return output_path

    def _write_dssp_output(self, cleaned_path: Path):
        try:
            dssp_output = subprocess.run([
                'mkdssp', str(cleaned_path), '-'
            ], capture_output=True, text=True, timeout=10)

            if dssp_output.returncode == 0 and dssp_output.stdout:
                dssp_report_path = self.reports_dir / (cleaned_path.stem + '_dssp.txt')
                with open(dssp_report_path, 'w') as f:
                    f.write(dssp_output.stdout)
        except Exception as e:
            logger.warning(f"Failed to write DSSP output for {cleaned_path.name}: {e}")

    def analyze_structure(self, structure, pdb_file: str) -> StructuralFeatures:
        secondary_structure = self._analyze_secondary_structure(structure, pdb_file)
        geometric_properties = self._calculate_geometric_properties(structure)
        surface_properties = self._calculate_surface_properties(structure)
        domain_architecture = self._analyze_domain_architecture(structure)
        unusual_features = self._identify_unusual_features(structure)
        return StructuralFeatures(
            secondary_structure=secondary_structure,
            geometric_properties=geometric_properties,
            surface_properties=surface_properties,
            domain_architecture=domain_architecture,
            unusual_features=unusual_features
        )

    def _analyze_secondary_structure(self, structure, pdb_file: str) -> Dict[str, float]:
        ss_content = {'helix': 0.0, 'sheet': 0.0, 'loop': 0.0}
        if self.dssp_available:
            try:
                cleaned_path = self._write_clean_pdb(structure, pdb_file)
                self._write_dssp_output(cleaned_path)
                dssp = DSSP(structure[0], str(cleaned_path), dssp='mkdssp')
                ss_codes = {'H': 'helix', 'B': 'sheet', 'E': 'sheet', 'G': 'helix',
                            'I': 'helix', 'T': 'loop', 'S': 'loop', '-': 'loop'}
                ss_count = Counter()
                for residue in dssp:
                    ss_code = residue[2]
                    ss_type = ss_codes.get(ss_code, 'loop')
                    ss_count[ss_type] += 1
                total = sum(ss_count.values())
                if total > 0:
                    for ss_type in ss_content:
                        ss_content[ss_type] = ss_count[ss_type] / total
            except Exception as e:
                logger.warning(f"DSSP analysis failed on {pdb_file}: {e}")
                ss_content = self._simple_secondary_structure_analysis(structure)
        else:
            ss_content = self._simple_secondary_structure_analysis(structure)
        return ss_content
    
    def _simple_secondary_structure_analysis(self, structure) -> Dict[str, float]:
        """Simple secondary structure analysis without DSSP"""
        ppb = PPBuilder()
        helix_count = 0
        sheet_count = 0
        total_residues = 0
        
        for chain in structure[0]:
            for polypeptide in ppb.build_peptides(chain):
                phi_psi = polypeptide.get_phi_psi_list()
                
                for phi, psi in phi_psi:
                    if phi is not None and psi is not None:
                        total_residues += 1
                        
                        # Simple phi/psi classification
                        phi_deg = np.degrees(phi)
                        psi_deg = np.degrees(psi)
                        
                        # Alpha helix region
                        if -90 <= phi_deg <= -30 and -70 <= psi_deg <= 50:
                            helix_count += 1
                        # Beta sheet region
                        elif -180 <= phi_deg <= -90 and 90 <= psi_deg <= 180:
                            sheet_count += 1
        
        if total_residues > 0:
            helix_frac = helix_count / total_residues
            sheet_frac = sheet_count / total_residues
            loop_frac = 1.0 - helix_frac - sheet_frac
            return {'helix': helix_frac, 'sheet': sheet_frac, 'loop': loop_frac}
        
        return {'helix': 0.0, 'sheet': 0.0, 'loop': 1.0}
    
    def _calculate_geometric_properties(self, structure) -> Dict[str, float]:
        """Calculate geometric properties"""
        properties = {}
        
        # Get all CA atoms
        ca_atoms = []
        for chain in structure[0]:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
        
        if len(ca_atoms) < 3:
            return properties
        
        ca_coords = np.array([atom.coord for atom in ca_atoms])
        
        # Radius of gyration
        center = np.mean(ca_coords, axis=0)
        rg = np.sqrt(np.mean(np.sum((ca_coords - center)**2, axis=1)))
        properties['radius_of_gyration'] = float(rg)
        
        # Maximum distance
        distances = pdist(ca_coords)
        properties['max_distance'] = float(np.max(distances))
        properties['mean_distance'] = float(np.mean(distances))
        
        # Asphericity (shape measure)
        # Calculate gyration tensor
        coords_centered = ca_coords - center
        gyration_tensor = np.dot(coords_centered.T, coords_centered) / len(ca_coords)
        eigenvals = np.linalg.eigvals(gyration_tensor)
        eigenvals = np.sort(eigenvals)[::-1]  # Sort descending
        
        if eigenvals[0] > 0:
            asphericity = eigenvals[0] - 0.5 * (eigenvals[1] + eigenvals[2])
            properties['asphericity'] = float(asphericity / eigenvals[0])
        
        return properties
    
    def _calculate_surface_properties(self, structure) -> Dict[str, float]:
        """Calculate surface-related properties"""
        properties = {}
        
        # Simplified surface area calculation
        all_atoms = []
        for chain in structure[0]:
            for residue in chain:
                for atom in residue:
                    if atom.element != 'H':
                        all_atoms.append(atom)
        
        if len(all_atoms) < 10:
            return properties
        
        # Estimate accessible surface area
        coords = np.array([atom.coord for atom in all_atoms])
        # Very simplified - just use convex hull volume as proxy
        try:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(coords)
            properties['approximate_surface_area'] = float(hull.area)
            properties['approximate_volume'] = float(hull.volume)
        except:
            # Fallback to bounding box
            ranges = np.max(coords, axis=0) - np.min(coords, axis=0)
            properties['bounding_box_volume'] = float(np.prod(ranges))
        
        return properties
    
    def _analyze_domain_architecture(self, structure) -> List[Dict]:
        """Analyze domain architecture (simplified)"""
        domains = []
        
        # Simple domain prediction based on chain breaks and compactness
        for chain in structure[0]:
            residues = list(chain.get_residues())
            if len(residues) < 50:
                continue
            
            # Look for potential domain boundaries based on CA distance breaks
            ca_coords = []
            residue_numbers = []
            
            for residue in residues:
                if 'CA' in residue:
                    ca_coords.append(residue['CA'].coord)
                    residue_numbers.append(residue.get_id()[1])
            
            if len(ca_coords) < 50:
                continue
            
            # Find large gaps in CA-CA distances
            domain_boundaries = [0]
            for i in range(1, len(ca_coords)):
                distance = np.linalg.norm(ca_coords[i] - ca_coords[i-1])
                if distance > 10.0:  # Large gap
                    domain_boundaries.append(i)
            domain_boundaries.append(len(ca_coords))
            
            # Create domain entries
            for i in range(len(domain_boundaries) - 1):
                start_idx = domain_boundaries[i]
                end_idx = domain_boundaries[i + 1]
                
                if end_idx - start_idx >= 30:  # Minimum domain size
                    domain_coords = ca_coords[start_idx:end_idx]
                    center = np.mean(domain_coords, axis=0)
                    
                    domain = {
                        'domain_id': f"Domain_{i+1}",
                        'start_residue': residue_numbers[start_idx],
                        'end_residue': residue_numbers[end_idx-1],
                        'length': end_idx - start_idx,
                        'center_coords': center.tolist(),
                        'chain_id': chain.id
                    }
                    domains.append(domain)
        
        return domains
    
    def _identify_unusual_features(self, structure) -> List[str]:
        """Identify unusual structural features"""
        unusual_features = []
        
        for chain in structure[0]:
            residues = list(chain.get_residues())
            
            # Check for unusual amino acids
            standard_aa = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                          'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                          'THR', 'TRP', 'TYR', 'VAL'}
            
            for residue in residues:
                resname = residue.get_resname()
                if resname not in standard_aa and len(resname) == 3:
                    unusual_features.append(f"Non-standard amino acid: {resname}")
            
            # Check for potential disulfide bonds
            cys_residues = [res for res in residues if res.get_resname() == 'CYS']
            if len(cys_residues) >= 2:
                for i, cys1 in enumerate(cys_residues):
                    for cys2 in cys_residues[i+1:]:
                        if 'SG' in cys1 and 'SG' in cys2:
                            distance = cys1['SG'] - cys2['SG']
                            if distance < 3.0:  # Potential disulfide bond
                                unusual_features.append(
                                    f"Potential disulfide bond: CYS{cys1.get_id()[1]}-CYS{cys2.get_id()[1]}"
                                )
            
            # Check for very long loops
            ca_atoms = [res['CA'] for res in residues if 'CA' in res]
            if len(ca_atoms) > 10:
                for i in range(1, len(ca_atoms)):
                    distance = ca_atoms[i] - ca_atoms[i-1]
                    if distance > 10.0:
                        unusual_features.append(f"Large structural gap at position {i}")
        
        return unusual_features

class PDBAnalysisPipeline:
    """Main pipeline for analyzing PDB files"""
    
    def __init__(self, output_dir: str = "pdb_analysis_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize analyzers
        self.catalytic_predictor = CatalyticSitePredictor()
        self.binding_analyzer = BindingSiteAnalyzer()
        self.enzyme_predictor = EnzymeFunctionPredictor()
        self.structural_analyzer = StructuralAnalyzer(output_dir)
        
        # PDB parser
        self.pdb_parser = PDBParser(QUIET=True)
    
    def analyze_pdb_folder(self, pdb_folder: str) -> Dict:
        """Analyze all PDB files in a folder"""
        pdb_folder = Path(pdb_folder)
        pdb_files = list(pdb_folder.glob("*.pdb"))
        
        if not pdb_files:
            raise ValueError(f"No PDB files found in {pdb_folder}")
        
        logger.info(f"Found {len(pdb_files)} PDB files to analyze")
        
        results = []
        failed_files = []
        
        for pdb_file in pdb_files:
            try:
                logger.info(f"Analyzing {pdb_file.name}...")
                result = self.analyze_single_pdb(str(pdb_file))
                results.append(result)
            except Exception as e:
                logger.error(f"Failed to analyze {pdb_file.name}: {e}")
                failed_files.append(str(pdb_file))
        
        # Compile final results
        analysis_results = {
            'analysis_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_files': len(pdb_files),
            'successful_analyses': len(results),
            'failed_analyses': len(failed_files),
            'failed_files': failed_files,
            'protein_analyses': results,
            'comparative_analysis': self._perform_comparative_analysis(results)
        }
        
        # Save results
        self._save_results(analysis_results)
        
        return analysis_results
    
    def analyze_single_pdb(self, pdb_file: str) -> ProteinAnalysisResult:
        """Analyze a single PDB file"""
        pdb_path = Path(pdb_file)
        
        # Parse structure
        structure = self.pdb_parser.get_structure(pdb_path.stem, pdb_file)
        
        # Extract sequence
        sequence = self._extract_sequence(structure)
        
        # Find catalytic sites
        catalytic_sites = self.catalytic_predictor.find_catalytic_sites(structure)
        
        # Find binding sites
        binding_sites = self.binding_analyzer.find_binding_sites(structure)
        
        # Predict enzyme functions
        enzyme_functions = self.enzyme_predictor.predict_enzyme_function(
            catalytic_sites, binding_sites, structure
        )
        
        # Structural analysis
        structural_features = self.structural_analyzer.analyze_structure(structure, pdb_file)
        
        # Physicochemical properties
        physicochemical_props = self._calculate_physicochemical_properties(sequence)
        
        # Quality metrics
        quality_metrics = self._calculate_quality_metrics(structure)
        
        # Overall confidence
        analysis_confidence = self._calculate_analysis_confidence(
            catalytic_sites, binding_sites, enzyme_functions, quality_metrics
        )
        
        return ProteinAnalysisResult(
            pdb_id=pdb_path.stem,
            pdb_file=str(pdb_file),
            sequence=sequence,
            length=len(sequence),
            catalytic_sites=catalytic_sites,
            binding_sites=binding_sites,
            enzyme_functions=enzyme_functions,
            structural_features=structural_features,
            physicochemical_properties=physicochemical_props,
            quality_metrics=quality_metrics,
            analysis_confidence=analysis_confidence
        )
    
    def _extract_sequence(self, structure) -> str:
        """Extract protein sequence from structure"""
        ppb = PPBuilder()
        sequences = []
        
        for chain in structure[0]:
            for polypeptide in ppb.build_peptides(chain):
                seq_str = str(polypeptide.get_sequence())
                sequences.append(seq_str)
        
        return ''.join(sequences)
    
    def _calculate_physicochemical_properties(self, sequence: str) -> Dict[str, float]:
        """Calculate physicochemical properties of the sequence"""
        if not sequence:
            return {}
        
        try:
            analyzer = ProteinAnalysis(sequence)
            
            properties = {
                'molecular_weight': analyzer.molecular_weight(),
                'isoelectric_point': analyzer.isoelectric_point(),
                'instability_index': analyzer.instability_index(),
                'aromaticity': analyzer.aromaticity(),
                'gravy': analyzer.gravy(),  # Grand average of hydropathy
            }
            
            # Amino acid composition
            aa_percent = analyzer.get_amino_acids_percent()
            for aa, percent in aa_percent.items():
                properties[f'percent_{aa}'] = percent
            
            # Secondary structure fraction
            ss_fraction = analyzer.secondary_structure_fraction()
            properties.update({
                'helix_fraction': ss_fraction[0],
                'turn_fraction': ss_fraction[1], 
                'sheet_fraction': ss_fraction[2]
            })
            
            return properties
            
        except Exception as e:
            logger.warning(f"Failed to calculate physicochemical properties: {e}")
            return {}
    
    def _calculate_quality_metrics(self, structure) -> Dict[str, float]:
        """Calculate structure quality metrics"""
        metrics = {}
        
        # Count atoms and residues
        total_atoms = 0
        total_residues = 0
        missing_atoms = 0
        
        for chain in structure[0]:
            for residue in chain:
                total_residues += 1
                residue_atoms = len(list(residue.get_atoms()))
                total_atoms += residue_atoms
                
                # Check for missing backbone atoms
                if residue.get_resname() in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                    required_atoms = ['N', 'CA', 'C', 'O']
                    for atom_name in required_atoms:
                        if atom_name not in residue:
                            missing_atoms += 1
        
        metrics['total_atoms'] = total_atoms
        metrics['total_residues'] = total_residues
        metrics['missing_atoms'] = missing_atoms
        metrics['completeness'] = 1.0 - (missing_atoms / (total_residues * 4)) if total_residues > 0 else 0.0
        
        # Resolution (if available in header)
        try:
            resolution = structure.header.get('resolution')
            if resolution:
                metrics['resolution'] = float(resolution)
        except:
            pass
        
        return metrics
    
    def _calculate_analysis_confidence(self, catalytic_sites: List[CatalyticSite],
                                     binding_sites: List[BindingSite],
                                     enzyme_functions: List[EnzymeFunction],
                                     quality_metrics: Dict[str, float]) -> float:
        """Calculate overall confidence in the analysis"""
        confidence = 0.0
        
        # Structure quality contribution
        completeness = quality_metrics.get('completeness', 0.0)
        confidence += completeness * 0.3
        
        # Catalytic sites contribution
        if catalytic_sites:
            avg_catalytic_confidence = np.mean([site.confidence for site in catalytic_sites])
            confidence += avg_catalytic_confidence * 0.3
        
        # Binding sites contribution
        if binding_sites:
            avg_binding_confidence = np.mean([site.binding_likelihood for site in binding_sites])
            confidence += avg_binding_confidence * 0.2
        
        # Function prediction contribution
        if enzyme_functions:
            avg_function_confidence = np.mean([func.confidence for func in enzyme_functions])
            confidence += avg_function_confidence * 0.2
        
        return min(1.0, confidence)
    
    def _perform_comparative_analysis(self, results: List[ProteinAnalysisResult]) -> Dict:
        """Perform comparative analysis across all proteins"""
        if not results:
            return {}
        
        comparative = {
            'total_proteins': len(results),
            'proteins_with_catalytic_sites': sum(1 for r in results if r.catalytic_sites),
            'proteins_with_binding_sites': sum(1 for r in results if r.binding_sites),
            'proteins_with_predicted_functions': sum(1 for r in results if r.enzyme_functions),
            'average_confidence': np.mean([r.analysis_confidence for r in results]),
            'length_distribution': {
                'min': min(r.length for r in results),
                'max': max(r.length for r in results),
                'mean': np.mean([r.length for r in results]),
                'median': np.median([r.length for r in results])
            }
        }
        
        # Enzyme class distribution
        ec_classes = []
        pathway_types = []
        microbial_origins = []
        
        for result in results:
            for func in result.enzyme_functions:
                if func.ec_number:
                    ec_classes.append(func.ec_number.split('.')[0])
                if func.is_biosynthetic and func.pathway_type:
                    pathway_types.append(func.pathway_type)
                if func.microbial_origin and func.microbial_origin != "Unknown":
                    microbial_origins.append(func.microbial_origin)
        
        ec_distribution = Counter(ec_classes)
        comparative['enzyme_class_distribution'] = dict(ec_distribution)
        
        # Biosynthetic pathway distribution
        pathway_distribution = Counter(pathway_types)
        comparative['biosynthetic_pathway_distribution'] = dict(pathway_distribution)
        
        # Microbial origin distribution
        origin_distribution = Counter(microbial_origins)
        comparative['microbial_origin_distribution'] = dict(origin_distribution)
        
        # Count biosynthetic enzymes
        biosynthetic_count = sum(1 for result in results 
                               if any(func.is_biosynthetic for func in result.enzyme_functions))
        comparative['proteins_with_biosynthetic_functions'] = biosynthetic_count
        
        # Most common catalytic residues
        all_catalytic_residues = []
        for result in results:
            for site in result.catalytic_sites:
                all_catalytic_residues.extend([res[0] for res in site.residues])
        
        catalytic_distribution = Counter(all_catalytic_residues)
        comparative['catalytic_residue_distribution'] = dict(catalytic_distribution.most_common(10))
        
        return comparative
    
    def _save_results(self, results: Dict):
        """Save analysis results in multiple formats"""
        # Save JSON results
        json_file = self.output_dir / "pdb_analysis_results.json"
        with open(json_file, 'w') as f:
            # Convert dataclasses to dicts for JSON serialization
            json_results = self._prepare_for_json(results)
            json.dump(json_results, f, indent=2, default=str)
        
        logger.info(f"Saved JSON results to {json_file}")
        
        # Save individual protein reports
        self._save_individual_reports(results['protein_analyses'])
        
        # Generate HTML report
        self._generate_html_report(results)
        
        # Generate summary report
        self._generate_summary_report(results)
    
    def _prepare_for_json(self, obj):
        """Prepare object for JSON serialization"""
        if isinstance(obj, list):
            return [self._prepare_for_json(item) for item in obj]
        elif isinstance(obj, dict):
            return {key: self._prepare_for_json(value) for key, value in obj.items()}
        elif hasattr(obj, '__dict__'):  # Dataclass or similar
            return self._prepare_for_json(asdict(obj) if hasattr(obj, '__dataclass_fields__') else obj.__dict__)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        else:
            return obj
    
    def _save_individual_reports(self, protein_analyses: List[ProteinAnalysisResult]):
        """Save individual protein analysis reports"""
        reports_dir = self.output_dir / "individual_reports"
        reports_dir.mkdir(exist_ok=True)
        
        for analysis in protein_analyses:
            report_file = reports_dir / f"{analysis.pdb_id}_analysis.txt"
            
            report_lines = []
            report_lines.append(f"PROTEIN ANALYSIS REPORT: {analysis.pdb_id}")
            report_lines.append("=" * 60)
            report_lines.append(f"PDB File: {analysis.pdb_file}")
            report_lines.append(f"Sequence Length: {analysis.length} amino acids")
            report_lines.append(f"Analysis Confidence: {analysis.analysis_confidence:.3f}")
            report_lines.append("")
            
            # Catalytic sites
            if analysis.catalytic_sites:
                report_lines.append("CATALYTIC SITES:")
                report_lines.append("-" * 20)
                for i, site in enumerate(analysis.catalytic_sites, 1):
                    report_lines.append(f"Site {i}: {site.site_type}")
                    report_lines.append(f"  Confidence: {site.confidence:.3f}")
                    report_lines.append(f"  Residues: {', '.join([f'{res[0]}{res[1]}' for res in site.residues])}")
                    report_lines.append(f"  Description: {site.description}")
                    if site.supporting_evidence:
                        report_lines.append(f"  Evidence: {'; '.join(site.supporting_evidence)}")
                    report_lines.append("")
            
            # Binding sites
            if analysis.binding_sites:
                report_lines.append("BINDING SITES:")
                report_lines.append("-" * 20)
                for i, site in enumerate(analysis.binding_sites, 1):
                    report_lines.append(f"Site {i}: {site.site_classification}")
                    report_lines.append(f"  Binding likelihood: {site.binding_likelihood:.3f}")
                    report_lines.append(f"  Volume: {site.volume:.1f} Ų")
                    report_lines.append(f"  Hydrophobicity: {site.hydrophobicity:.3f}")
                    report_lines.append(f"  Charge: {site.charge:.3f}")
                    report_lines.append("")
            
            # Enzyme functions
            if analysis.enzyme_functions:
                report_lines.append("PREDICTED ENZYME FUNCTIONS:")
                report_lines.append("-" * 30)
                for i, func in enumerate(analysis.enzyme_functions, 1):
                    report_lines.append(f"Function {i}: {func.enzyme_class}")
                    if func.ec_number:
                        report_lines.append(f"  EC Number: {func.ec_number}")
                    report_lines.append(f"  Confidence: {func.confidence:.3f}")
                    report_lines.append(f"  Description: {func.function_description}")
                    report_lines.append(f"  Mechanism: {func.catalytic_mechanism}")
                    if func.substrate_prediction:
                        report_lines.append(f"  Predicted substrates: {', '.join(func.substrate_prediction)}")
                    report_lines.append("")
            
            with open(report_file, 'w') as f:
                f.write('\n'.join(report_lines))
        
        logger.info(f"Saved individual reports to {reports_dir}")
    
    def _generate_html_report(self, results: Dict):
        """Generate comprehensive HTML report"""
        html_file = self.output_dir / "pdb_analysis_report.html"
        
        # Create visualizations
        charts = self._create_analysis_charts(results)
        
        html_content = self._build_html_content(results, charts)
        
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Saved HTML report to {html_file}")
    
    def _create_analysis_charts(self, results: Dict) -> Dict[str, str]:
        """Create analysis charts for the HTML report"""
        charts = {}
        
        protein_analyses = results['protein_analyses']
        comparative = results['comparative_analysis']
        
        try:
            # Chart 1: Analysis confidence distribution
            confidences = [analysis.analysis_confidence for analysis in protein_analyses]
            if confidences:
                plt.figure(figsize=(10, 6))
                plt.hist(confidences, bins=20, color='skyblue', alpha=0.7, edgecolor='black')
                plt.xlabel('Analysis Confidence')
                plt.ylabel('Number of Proteins')
                plt.title('Distribution of Analysis Confidence Scores')
                plt.grid(True, alpha=0.3)
                
                buffer = BytesIO()
                plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
                buffer.seek(0)
                charts['confidence_dist'] = base64.b64encode(buffer.read()).decode('utf-8')
                buffer.close()
                plt.close()
            
            # Chart 2: Protein length distribution
            lengths = [analysis.length for analysis in protein_analyses]
            if lengths:
                plt.figure(figsize=(10, 6))
                plt.hist(lengths, bins=20, color='lightgreen', alpha=0.7, edgecolor='black')
                plt.xlabel('Protein Length (amino acids)')
                plt.ylabel('Number of Proteins')
                plt.title('Distribution of Protein Lengths')
                plt.grid(True, alpha=0.3)
                
                buffer = BytesIO()
                plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
                buffer.seek(0)
                charts['length_dist'] = base64.b64encode(buffer.read()).decode('utf-8')
                buffer.close()
                plt.close()
            
            # Chart 3: Enzyme class distribution
            if 'enzyme_class_distribution' in comparative and comparative['enzyme_class_distribution']:
                ec_classes = list(comparative['enzyme_class_distribution'].keys())
                ec_counts = list(comparative['enzyme_class_distribution'].values())
                
                plt.figure(figsize=(10, 6))
                plt.bar(ec_classes, ec_counts, color='orange', alpha=0.7)
                plt.xlabel('Enzyme Class')
                plt.ylabel('Number of Proteins')
                plt.title('Distribution of Predicted Enzyme Classes')
                plt.xticks(rotation=45)
                
                buffer = BytesIO()
                plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
                buffer.seek(0)
                charts['enzyme_dist'] = base64.b64encode(buffer.read()).decode('utf-8')
                buffer.close()
                plt.close()
            
            # Chart 4: Biosynthetic pathway distribution
            if 'biosynthetic_pathway_distribution' in comparative and comparative['biosynthetic_pathway_distribution']:
                pathway_types = list(comparative['biosynthetic_pathway_distribution'].keys())
                pathway_counts = list(comparative['biosynthetic_pathway_distribution'].values())
                
                plt.figure(figsize=(12, 6))
                bars = plt.bar(pathway_types, pathway_counts, color=['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#ffeaa7', '#dda0dd'])
                plt.xlabel('Biosynthetic Pathway Type')
                plt.ylabel('Number of Proteins')
                plt.title('Distribution of Biosynthetic Pathway Types')
                plt.xticks(rotation=45, ha='right')
                
                # Add value labels on bars
                for bar, count in zip(bars, pathway_counts):
                    height = bar.get_height()
                    plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                            f'{count}', ha='center', va='bottom', fontsize=10, fontweight='bold')
                
                plt.tight_layout()
                
                buffer = BytesIO()
                plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
                buffer.seek(0)
                charts['pathway_dist'] = base64.b64encode(buffer.read()).decode('utf-8')
                buffer.close()
                plt.close()
            
            # Chart 5: Functional features summary
            feature_counts = {
                'Catalytic Sites': sum(len(p.catalytic_sites) for p in protein_analyses),
                'Binding Sites': sum(len(p.binding_sites) for p in protein_analyses),
                'Biosynthetic Functions': sum(len([f for f in p.enzyme_functions if f.is_biosynthetic]) for p in protein_analyses),
                'General Functions': sum(len([f for f in p.enzyme_functions if not f.is_biosynthetic]) for p in protein_analyses)
            }
            
            if any(feature_counts.values()):
                plt.figure(figsize=(8, 8))
                colors = ['#ff9999', '#66b3ff', '#ff6b6b', '#99ff99']
                plt.pie(feature_counts.values(), labels=feature_counts.keys(), autopct='%1.1f%%',
                       colors=colors)
                plt.title('Distribution of Predicted Functional Features')
                
                buffer = BytesIO()
                plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
                buffer.seek(0)
                charts['features_pie'] = base64.b64encode(buffer.read()).decode('utf-8')
                buffer.close()
                plt.close()
            
        except Exception as e:
            logger.warning(f"Failed to create some charts: {e}")
        
        return charts
    
    def _build_html_content(self, results: Dict, charts: Dict[str, str]) -> str:
        """Build complete HTML content"""
        protein_analyses = results['protein_analyses']
        comparative = results['comparative_analysis']
        
        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>PDB Enzymatic Function Analysis Report</title>
            <style>
                body {{
                    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                    line-height: 1.6;
                    margin: 0;
                    padding: 20px;
                    background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
                    color: #333;
                }}
                .container {{
                    max-width: 1400px;
                    margin: 0 auto;
                    background: white;
                    padding: 30px;
                    border-radius: 15px;
                    box-shadow: 0 10px 30px rgba(0,0,0,0.1);
                }}
                h1 {{
                    color: #2c3e50;
                    text-align: center;
                    border-bottom: 3px solid #3498db;
                    padding-bottom: 15px;
                    margin-bottom: 30px;
                }}
                h2 {{
                    color: #34495e;
                    border-left: 5px solid #3498db;
                    padding-left: 15px;
                    margin-top: 40px;
                }}
                h3 {{
                    color: #5d6d7e;
                    margin-top: 25px;
                }}
                .summary-grid {{
                    display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                    gap: 20px;
                    margin: 20px 0;
                }}
                .summary-card {{
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: white;
                    padding: 20px;
                    border-radius: 10px;
                    text-align: center;
                    box-shadow: 0 5px 15px rgba(0,0,0,0.1);
                }}
                .summary-card h3 {{
                    margin: 0 0 10px 0;
                    color: white;
                }}
                .summary-card .value {{
                    font-size: 2em;
                    font-weight: bold;
                    margin: 10px 0;
                }}
                .chart-container {{
                    text-align: center;
                    margin: 30px 0;
                    background: #f8f9fa;
                    padding: 20px;
                    border-radius: 10px;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.05);
                }}
                .chart-container img {{
                    max-width: 100%;
                    height: auto;
                    border-radius: 8px;
                }}
                table {{
                    width: 100%;
                    border-collapse: collapse;
                    margin: 20px 0;
                    box-shadow: 0 5px 15px rgba(0,0,0,0.08);
                    border-radius: 10px;
                    overflow: hidden;
                }}
                th {{
                    background: linear-gradient(135deg, #34495e 0%, #2c3e50 100%);
                    color: white;
                    padding: 15px;
                    text-align: left;
                    font-weight: bold;
                }}
                td {{
                    padding: 12px 15px;
                    border-bottom: 1px solid #ecf0f1;
                }}
                tr:nth-child(even) {{
                    background: #f8f9fa;
                }}
                tr:hover {{
                    background: #e3f2fd;
                    transform: scale(1.02);
                    transition: all 0.2s ease;
                }}
                .confidence-high {{ background: #27ae60; color: white; padding: 4px 8px; border-radius: 4px; }}
                .confidence-medium {{ background: #f39c12; color: white; padding: 4px 8px; border-radius: 4px; }}
                .confidence-low {{ background: #e74c3c; color: white; padding: 4px 8px; border-radius: 4px; }}
                .function-tag {{
                    background: #3498db;
                    color: white;
                    padding: 2px 6px;
                    border-radius: 12px;
                    font-size: 0.8em;
                    margin: 2px;
                    display: inline-block;
                }}
                .biosynthetic-tag {{
                    background: #e67e22;
                    color: white;
                    padding: 2px 6px;
                    border-radius: 12px;
                    font-size: 0.8em;
                    margin: 2px;
                    display: inline-block;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.2);
                }}
                .biosynthetic-indicator {{
                    background: linear-gradient(45deg, #ff6b6b, #feca57);
                    color: white;
                    padding: 3px 8px;
                    border-radius: 15px;
                    font-size: 0.7em;
                    font-weight: bold;
                    margin-right: 5px;
                    display: inline-block;
                    box-shadow: 0 2px 6px rgba(0,0,0,0.3);
                }}
                .site-tag {{
                    background: #e74c3c;
                    color: white;
                    padding: 2px 6px;
                    border-radius: 12px;
                    font-size: 0.8em;
                    margin: 2px;
                    display: inline-block;
                }}
                .key-findings {{
                    background: linear-gradient(135deg, #74b9ff 0%, #0984e3 100%);
                    color: white;
                    padding: 20px;
                    border-radius: 10px;
                    margin: 20px 0;
                }}
                .footer {{
                    text-align: center;
                    margin-top: 40px;
                    padding-top: 20px;
                    border-top: 2px solid #bdc3c7;
                    color: #7f8c8d;
                }}
                .expandable {{
                    cursor: pointer;
                    transition: all 0.3s ease;
                }}
                .expandable:hover {{
                    background: #e8f4f8;
                }}
                .details {{
                    display: none;
                    padding: 10px;
                    background: #f1f2f6;
                    border-radius: 5px;
                    margin-top: 5px;
                }}
            </style>
            <script>
                function toggleDetails(id) {{
                    var details = document.getElementById(id);
                    if (details.style.display === 'none' || details.style.display === '') {{
                        details.style.display = 'block';
                    }} else {{
                        details.style.display = 'none';
                    }}
                }}
            </script>
        </head>
        <body>
            <div class="container">
                <h1>🧬 PDB Enzymatic Function & Structural Analysis Report</h1>
                
                <div class="key-findings">
                    <h3>📊 Analysis Overview</h3>
                    <p><strong>Analysis Date:</strong> {results['analysis_timestamp']}</p>
                    <p><strong>Total PDB Files Processed:</strong> {results['total_files']}</p>
                    <p><strong>Successful Analyses:</strong> {results['successful_analyses']}</p>
                    <p><strong>Average Analysis Confidence:</strong> {comparative.get('average_confidence', 0):.3f}</p>
                </div>
                
                <div class="summary-grid">
                    <div class="summary-card">
                        <h3>Proteins Analyzed</h3>
                        <div class="value">{len(protein_analyses)}</div>
                    </div>
                    <div class="summary-card">
                        <h3>Catalytic Sites Found</h3>
                        <div class="value">{comparative.get('proteins_with_catalytic_sites', 0)}</div>
                    </div>
                    <div class="summary-card">
                        <h3>Biosynthetic Enzymes</h3>
                        <div class="value">{comparative.get('proteins_with_biosynthetic_functions', 0)}</div>
                    </div>
                    <div class="summary-card">
                        <h3>Functions Predicted</h3>
                        <div class="value">{comparative.get('proteins_with_predicted_functions', 0)}</div>
                    </div>
                </div>
        """
        
        # Add charts
        if 'confidence_dist' in charts:
            html += f"""
                <div class="chart-container">
                    <h3>📈 Analysis Confidence Distribution</h3>
                    <img src="data:image/png;base64,{charts['confidence_dist']}" alt="Confidence Distribution">
                    <p><em>Distribution of analysis confidence scores across all proteins</em></p>
                </div>
            """
        
        if 'length_dist' in charts:
            html += f"""
                <div class="chart-container">
                    <h3>📏 Protein Length Distribution</h3>
                    <img src="data:image/png;base64,{charts['length_dist']}" alt="Length Distribution">
                    <p><em>Distribution of protein lengths in the analyzed dataset</em></p>
                </div>
            """
        
        if 'enzyme_dist' in charts:
            html += f"""
                <div class="chart-container">
                    <h3>🧪 Enzyme Class Distribution</h3>
                    <img src="data:image/png;base64,{charts['enzyme_dist']}" alt="Enzyme Distribution">
                    <p><em>Distribution of predicted enzyme classes (EC numbers)</em></p>
                </div>
            """
        
        if 'pathway_dist' in charts:
            html += f"""
                <div class="chart-container">
                    <h3>🧬 Biosynthetic Pathway Distribution</h3>
                    <img src="data:image/png;base64,{charts['pathway_dist']}" alt="Biosynthetic Pathways">
                    <p><em>Distribution of biosynthetic pathway types in natural product enzymes</em></p>
                </div>
            """
        
        if 'features_pie' in charts:
            html += f"""
                <div class="chart-container">
                    <h3>⚙️ Functional Features Overview</h3>
                    <img src="data:image/png;base64,{charts['features_pie']}" alt="Features Pie Chart">
                    <p><em>Distribution of predicted functional features across all proteins</em></p>
                </div>
            """
        
        # Individual protein results table
        html += """
            <h2>🔬 Individual Protein Analysis Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>PDB ID</th>
                        <th>Length (aa)</th>
                        <th>Confidence</th>
                        <th>Catalytic Sites</th>
                        <th>Binding Sites</th>
                        <th>Predicted Functions</th>
                        <th>Key Features</th>
                    </tr>
                </thead>
                <tbody>
        """
        
        for i, analysis in enumerate(protein_analyses):
            # Determine confidence class
            conf = analysis.analysis_confidence
            conf_class = "confidence-high" if conf > 0.7 else "confidence-medium" if conf > 0.4 else "confidence-low"
            
            # Get catalytic sites summary
            catalytic_summary = f"{len(analysis.catalytic_sites)} sites" if analysis.catalytic_sites else "None"
            
            # Get binding sites summary  
            binding_summary = f"{len(analysis.binding_sites)} sites" if analysis.binding_sites else "None"
            
            # Get function summary
            if analysis.enzyme_functions:
                functions = []
                biosynthetic_count = 0
                for func in analysis.enzyme_functions[:3]:  # Limit to 3
                    if func.is_biosynthetic:
                        biosynthetic_count += 1
                        functions.append(f'<span class="biosynthetic-tag">{func.pathway_type}</span>')
                    else:
                        functions.append(f'<span class="function-tag">{func.enzyme_class}</span>')
                
                function_summary = ''.join(functions)
                if len(analysis.enzyme_functions) > 3:
                    function_summary += f'<span class="function-tag">+{len(analysis.enzyme_functions)-3} more</span>'
                
                # Add biosynthetic indicator
                if biosynthetic_count > 0:
                    function_summary = f'<span class="biosynthetic-indicator">🧬 BIOSYNTHETIC</span>' + function_summary
            else:
                function_summary = "None predicted"
            
            # Key features
            key_features = []
            if analysis.catalytic_sites:
                key_features.extend([f'<span class="site-tag">{site.site_type}</span>' 
                                   for site in analysis.catalytic_sites[:2]])
            if analysis.structural_features.unusual_features:
                key_features.append('<span class="site-tag">Unusual features</span>')
            
            key_features_html = ''.join(key_features) if key_features else "Standard structure"
            
            html += f"""
                <tr class="expandable" onclick="toggleDetails('details_{i}')">
                    <td><strong>{analysis.pdb_id}</strong></td>
                    <td>{analysis.length}</td>
                    <td><span class="{conf_class}">{conf:.3f}</span></td>
                    <td>{catalytic_summary}</td>
                    <td>{binding_summary}</td>
                    <td>{function_summary}</td>
                    <td>{key_features_html}</td>
                </tr>
                <tr>
                    <td colspan="7">
                        <div id="details_{i}" class="details">
                            <strong>Detailed Analysis for {analysis.pdb_id}:</strong><br>
            """
            
            # Add detailed information
            if analysis.catalytic_sites:
                html += "<strong>Catalytic Sites:</strong><br>"
                for site in analysis.catalytic_sites:
                    residues_str = ", ".join([f"{res[0]}{res[1]}" for res in site.residues])
                    html += f"• {site.site_type} (confidence: {site.confidence:.3f}): {residues_str}<br>"
            
            if analysis.enzyme_functions:
                html += "<strong>Predicted Functions:</strong><br>"
                for func in analysis.enzyme_functions:
                    ec_str = f" ({func.ec_number})" if func.ec_number else ""
                    biosynthetic_str = ""
                    if func.is_biosynthetic:
                        biosynthetic_str = f" 🧬 <em>Biosynthetic Pathway: {func.pathway_type}</em>"
                        if func.microbial_origin and func.microbial_origin != "Unknown":
                            biosynthetic_str += f" ({func.microbial_origin})"
                    
                    html += f"• {func.enzyme_class}{ec_str} (confidence: {func.confidence:.3f}){biosynthetic_str}<br>"
                    html += f"  Mechanism: {func.catalytic_mechanism}<br>"
                    
                    if func.substrate_prediction:
                        substrates = ", ".join(func.substrate_prediction[:3])  # Limit to 3
                        if len(func.substrate_prediction) > 3:
                            substrates += f" (+{len(func.substrate_prediction)-3} more)"
                        html += f"  Substrates: {substrates}<br>"
            
            if analysis.structural_features.unusual_features:
                html += f"<strong>Unusual Features:</strong> {', '.join(analysis.structural_features.unusual_features)}<br>"
            
            # Physicochemical properties
            if analysis.physicochemical_properties:
                mw = analysis.physicochemical_properties.get('molecular_weight', 0)
                pi = analysis.physicochemical_properties.get('isoelectric_point', 0)
                html += f"<strong>Properties:</strong> MW: {mw:.0f} Da, pI: {pi:.2f}<br>"
            
            html += """
                        </div>
                    </td>
                </tr>
            """
        
        html += """
                </tbody>
            </table>
        """
        
        # Comparative analysis section
        if comparative:
            html += f"""
                <h2>📊 Comparative Analysis</h2>
                <div class="key-findings">
                    <h3>🔍 Key Insights</h3>
                    <ul>
                        <li><strong>Analysis Success Rate:</strong> {(results['successful_analyses']/results['total_files']*100):.1f}% of PDB files successfully analyzed</li>
                        <li><strong>Biosynthetic Enzyme Detection:</strong> {(comparative.get('proteins_with_biosynthetic_functions', 0)/len(protein_analyses)*100):.1f}% of proteins identified as biosynthetic enzymes</li>
                        <li><strong>Functional Annotation Rate:</strong> {(comparative.get('proteins_with_predicted_functions', 0)/len(protein_analyses)*100):.1f}% of proteins have predicted enzymatic functions</li>
                        <li><strong>Catalytic Site Detection:</strong> {(comparative.get('proteins_with_catalytic_sites', 0)/len(protein_analyses)*100):.1f}% of proteins have identifiable catalytic sites</li>
                        <li><strong>Natural Product Pathways:</strong> {len(comparative.get('biosynthetic_pathway_distribution', {}))} different pathway types detected</li>
                        <li><strong>Average Protein Length:</strong> {comparative.get('length_distribution', {}).get('mean', 0):.0f} amino acids</li>
                    </ul>
                </div>
            """
            
            # Most common catalytic residues
            if 'catalytic_residue_distribution' in comparative:
                html += """
                    <h3>Most Common Catalytic Residues</h3>
                    <table>
                        <thead>
                            <tr><th>Residue</th><th>Frequency</th><th>Percentage</th></tr>
                        </thead>
                        <tbody>
                """
                
                total_catalytic = sum(comparative['catalytic_residue_distribution'].values())
                for residue, count in comparative['catalytic_residue_distribution'].items():
                    percentage = (count / total_catalytic * 100) if total_catalytic > 0 else 0
                    html += f"""
                        <tr>
                            <td><strong>{residue}</strong></td>
                            <td>{count}</td>
                            <td>{percentage:.1f}%</td>
                        </tr>
                    """
                
                html += """
                        </tbody>
                    </table>
                """
            
            # Biosynthetic pathway distribution
            if 'biosynthetic_pathway_distribution' in comparative and comparative['biosynthetic_pathway_distribution']:
                html += """
                    <h3>🧬 Biosynthetic Pathway Types Detected</h3>
                    <table>
                        <thead>
                            <tr><th>Pathway Type</th><th>Number of Proteins</th><th>Potential Products</th></tr>
                        </thead>
                        <tbody>
                """
                
                # Add pathway descriptions
                pathway_products = {
                    'Type I Polyketide': 'Macrolides, polyenes, complex polyketides',
                    'Type II Polyketide': 'Anthracyclines, tetracyclines, aromatic compounds',
                    'Type III Polyketide': 'Chalcones, pyrones, simple polyketides',
                    'NRPS Peptide': 'Cyclic peptides, lipopeptides, siderophores',
                    'Terpenoid': 'Monoterpenes, sesquiterpenes, diterpenes',
                    'P450 Oxidation': 'Hydroxylated metabolites, epoxides',
                    'Methylation': 'O-methylated, N-methylated compounds',
                    'Halogenation': 'Chlorinated, brominated metabolites',
                    'Glycosylation': 'Glycosides, nucleoside antibiotics',
                    'Bacterial metabolism': 'Metal-dependent oxidation products',
                    'Fungal metabolism': 'Phenolic oxidation products'
                }
                
                for pathway, count in comparative['biosynthetic_pathway_distribution'].items():
                    products = pathway_products.get(pathway, 'Various secondary metabolites')
                    html += f"""
                        <tr>
                            <td><strong>{pathway}</strong></td>
                            <td>{count}</td>
                            <td><em>{products}</em></td>
                        </tr>
                    """
                
                html += """
                        </tbody>
                    </table>
                """
        
        # Recommendations section
        html += f"""
            <div class="key-findings">
                <h3>💡 Recommendations & Next Steps</h3>
                <ul>
                    <li><strong>High Confidence Proteins:</strong> Focus experimental validation on proteins with confidence > 0.7</li>
                    <li><strong>Catalytic Site Validation:</strong> Use site-directed mutagenesis to validate predicted catalytic residues</li>
                    <li><strong>Functional Assays:</strong> Design enzymatic assays based on predicted EC classes and mechanisms</li>
                    <li><strong>Structural Validation:</strong> Compare predictions with experimental structures when available</li>
                    <li><strong>Comparative Studies:</strong> Analyze structurally similar proteins for functional relationships</li>
                    <li><strong>Literature Mining:</strong> Cross-reference predictions with published functional studies</li>
                </ul>
            </div>
            
            <div class="footer">
                <p>Generated by PDB Enzymatic Function Analyzer</p>
                <p>Report created on {results['analysis_timestamp']}</p>
                <p>Click on protein rows to expand detailed information</p>
            </div>
            
            </div>
        </body>
        </html>
        """
        
        return html
    
    def _generate_summary_report(self, results: Dict):
        """Generate a text summary report"""
        summary_file = self.output_dir / "analysis_summary.txt"
        
        lines = []
        lines.append("PDB ENZYMATIC FUNCTION ANALYSIS SUMMARY")
        lines.append("=" * 60)
        lines.append(f"Analysis Date: {results['analysis_timestamp']}")
        lines.append(f"Total PDB files: {results['total_files']}")
        lines.append(f"Successful analyses: {results['successful_analyses']}")
        lines.append(f"Failed analyses: {results['failed_analyses']}")
        lines.append("")
        
        if results['failed_files']:
            lines.append("Failed files:")
            for failed_file in results['failed_files']:
                lines.append(f"  - {failed_file}")
            lines.append("")
        
        comparative = results['comparative_analysis']
        if comparative:
            lines.append("COMPARATIVE ANALYSIS:")
            lines.append("-" * 30)
            lines.append(f"Proteins with catalytic sites: {comparative.get('proteins_with_catalytic_sites', 0)}")
            lines.append(f"Proteins with binding sites: {comparative.get('proteins_with_binding_sites', 0)}")
            lines.append(f"Proteins with predicted functions: {comparative.get('proteins_with_predicted_functions', 0)}")
            lines.append(f"Average analysis confidence: {comparative.get('average_confidence', 0):.3f}")
            lines.append("")
            
            if 'enzyme_class_distribution' in comparative:
                lines.append("Enzyme class distribution:")
                for ec_class, count in comparative['enzyme_class_distribution'].items():
                    lines.append(f"  {ec_class}: {count} proteins")
                lines.append("")
        
        # Top proteins by confidence
        protein_analyses = results['protein_analyses']
        if protein_analyses:
            sorted_proteins = sorted(protein_analyses, 
                                   key=lambda x: x.analysis_confidence, reverse=True)
            
            lines.append("TOP 10 PROTEINS BY ANALYSIS CONFIDENCE:")
            lines.append("-" * 40)
            for i, protein in enumerate(sorted_proteins[:10], 1):
                lines.append(f"{i:2d}. {protein.pdb_id} (confidence: {protein.analysis_confidence:.3f})")
                if protein.enzyme_functions:
                    func_names = [func.enzyme_class for func in protein.enzyme_functions]
                    lines.append(f"     Functions: {', '.join(func_names)}")
                lines.append("")
        
        with open(summary_file, 'w') as f:
            f.write('\n'.join(lines))
        
        logger.info(f"Saved summary report to {summary_file}")

def main():
    """Main command line interface"""
    parser = argparse.ArgumentParser(
        description='PDB Enzymatic Function and Structural Analysis Pipeline with Biosynthetic Enzyme Detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
ANALYSIS FEATURES:
  • Catalytic site prediction using structural patterns
  • Binding site detection and characterization  
  • General enzyme function prediction with EC classification
  • 🧬 BIOSYNTHETIC ENZYME ANALYSIS:
    - Polyketide synthases (Type I, II, III)
    - Nonribosomal peptide synthetases (NRPS)
    - Terpene synthases
    - Cytochrome P450s, methyltransferases, halogenases
    - Bacterial and fungal-specific enzyme features
  • Natural product pathway type prediction
  • Structural feature analysis and comparison
  • Physicochemical property calculation
  • Comprehensive HTML reporting with biosynthetic insights

BIOSYNTHETIC SPECIALIZATION:
  • Detects enzyme families involved in natural product biosynthesis
  • Predicts secondary metabolite pathway types
  • Identifies bacterial vs fungal enzyme characteristics
  • Suggests potential natural products and substrates
  • Specialized for antibiotic and secondary metabolite discovery

REQUIREMENTS:
  • biopython, numpy, scipy, matplotlib, scikit-learn
  • Optional: DSSP/mkdssp for enhanced secondary structure analysis

EXAMPLES:
  python pdb_analyzer.py /path/to/bacterial/pdb/folder
  python pdb_analyzer.py /path/to/fungal/structures -o biosynthetic_analysis --verbose
  python pdb_analyzer.py ./unknown_enzymes --confidence-threshold 0.5
        """
    )
    
    parser.add_argument('pdb_folder', help='Path to folder containing PDB files')
    parser.add_argument('-o', '--output_dir', default='pdb_analysis_results',
                       help='Output directory (default: pdb_analysis_results)')
    parser.add_argument('--confidence-threshold', type=float, default=0.3,
                       help='Minimum confidence threshold for predictions (default: 0.3)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        # Keep matplotlib font manager warnings suppressed even in debug mode
        logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)
    
    # Validate input
    pdb_folder = Path(args.pdb_folder)
    if not pdb_folder.exists() or not pdb_folder.is_dir():
        print(f"❌ Error: PDB folder '{args.pdb_folder}' not found or not a directory")
        sys.exit(1)
    
    pdb_files = list(pdb_folder.glob("*.pdb"))
    if not pdb_files:
        print(f"❌ Error: No PDB files found in '{args.pdb_folder}'")
        sys.exit(1)
    
    # Print header
    print("🧬 PDB ENZYMATIC FUNCTION & STRUCTURAL ANALYSIS PIPELINE")
    print("=" * 80)
    print(f"📂 PDB folder: {args.pdb_folder}")
    print(f"📁 Output directory: {args.output_dir}")
    print(f"🎯 Confidence threshold: {args.confidence_threshold}")
    print(f"📊 PDB files found: {len(pdb_files)}")
    print()
    
    try:
        # Initialize pipeline
        pipeline = PDBAnalysisPipeline(args.output_dir)
        
        # Show DSSP status
        dssp_status = "✅ DSSP available" if pipeline.structural_analyzer.dssp_available else "⚠️  DSSP not available (using phi/psi fallback)"
        print(f"🔬 Secondary structure method: {dssp_status}")
        print()

        # Run analysis
        start_time = time.time()
        results = pipeline.analyze_pdb_folder(args.pdb_folder)
        analysis_time = time.time() - start_time
        
        # Print results summary
        print("\n" + "=" * 80)
        print("✅ ANALYSIS COMPLETED!")
        print("=" * 80)
        print(f"⏱️  Total analysis time: {analysis_time:.1f} seconds")
        print(f"📊 Files processed: {results['successful_analyses']}/{results['total_files']}")
        print(f"⚡ Average time per protein: {analysis_time/results['successful_analyses']:.1f}s")
        
        if results['protein_analyses']:
            print(f"\n🔬 FUNCTIONAL ANALYSIS SUMMARY:")
            comparative = results['comparative_analysis']
            print(f"   Proteins with catalytic sites: {comparative.get('proteins_with_catalytic_sites', 0)}")
            print(f"   Proteins with binding sites: {comparative.get('proteins_with_binding_sites', 0)}")
            print(f"   🧬 Biosynthetic enzymes detected: {comparative.get('proteins_with_biosynthetic_functions', 0)}")
            print(f"   Proteins with predicted functions: {comparative.get('proteins_with_predicted_functions', 0)}")
            print(f"   Average analysis confidence: {comparative.get('average_confidence', 0):.3f}")
            
            # Show biosynthetic pathway types
            if comparative.get('biosynthetic_pathway_distribution'):
                print(f"\n🧬 DETECTED BIOSYNTHETIC PATHWAYS:")
                for pathway, count in comparative['biosynthetic_pathway_distribution'].items():
                    print(f"   {pathway}: {count} proteins")
            
            # Top predictions
            sorted_proteins = sorted(results['protein_analyses'], 
                                   key=lambda x: x.analysis_confidence, reverse=True)
            
            print(f"\n🎯 TOP 5 MOST CONFIDENT PREDICTIONS:")
            for i, protein in enumerate(sorted_proteins[:5], 1):
                functions = []
                biosynthetic_functions = []
                
                for func in protein.enzyme_functions:
                    if func.is_biosynthetic:
                        biosynthetic_functions.append(f"{func.pathway_type} ({func.enzyme_class})")
                    else:
                        functions.append(func.enzyme_class)
                
                all_functions = biosynthetic_functions + functions[:2]  # Prioritize biosynthetic
                
                print(f"   {i}. {protein.pdb_id} (confidence: {protein.analysis_confidence:.3f})")
                if biosynthetic_functions:
                    print(f"      🧬 Biosynthetic: {', '.join(biosynthetic_functions[:2])}")
                if functions and len(all_functions) < 3:
                    print(f"      Functions: {', '.join(functions[:2])}")
                if protein.catalytic_sites:
                    catalytic_types = [site.site_type for site in protein.catalytic_sites]
                    print(f"      Catalytic sites: {', '.join(catalytic_types[:2])}")
        
        if results['failed_analyses'] > 0:
            print(f"\n⚠️  Failed to analyze {results['failed_analyses']} files:")
            for failed_file in results['failed_files']:
                print(f"   - {Path(failed_file).name}")
        
        print(f"\n📁 GENERATED FILES:")
        print(f"   📊 pdb_analysis_results.json - Complete analysis data")
        print(f"   🌐 pdb_analysis_report.html - ⭐ MAIN INTERACTIVE REPORT")
        print(f"   📋 analysis_summary.txt - Quick summary")
        print(f"   📂 individual_reports/ - Detailed per-protein reports")
        
        print(f"\n💡 NEXT STEPS:")
        print(f"   1. Open pdb_analysis_report.html in your browser for interactive results")
        print(f"   2. Review high-confidence predictions for experimental validation")
        print(f"   3. 🧬 Prioritize biosynthetic enzymes for natural product discovery")
        print(f"   4. Check individual_reports/ for detailed protein analyses")
        print(f"   5. Use predicted catalytic sites for targeted mutagenesis studies")
        print(f"   6. Design functional assays based on predicted enzyme classes")
        if comparative.get('proteins_with_biosynthetic_functions', 0) > 0:
            print(f"   7. 🧬 Consider heterologous expression for biosynthetic enzymes")
            print(f"   8. 🧬 Analyze gene clusters for complete biosynthetic pathways")
            print(f"   9. 🧬 Test predicted substrates and cofactors for natural product formation")
       
    except KeyboardInterrupt:
        print("\n❌ Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        print(f"\n❌ Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
