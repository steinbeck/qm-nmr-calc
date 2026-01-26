"""Tests for canonical atom ordering functionality.

These tests ensure atom indices remain consistent across the conformer
lifecycle (generation -> optimization -> NMR parsing).
"""
import pytest
from qm_nmr_calc.atom_ordering import (
    canonical_atom_order,
    map_nwchem_to_canonical,
    get_atom_count,
)


class TestCanonicalAtomOrder:
    """Test canonical_atom_order() function."""

    def test_methane_ordering(self):
        """Methane (C) should have 5 atoms (1C + 4H) with deterministic ordering."""
        order = canonical_atom_order("C")

        # Should have 5 atoms total
        assert len(order) == 5

        # All indices should be unique (proper ranking)
        assert len(set(order)) == 5

        # Ranks should be 0-4 (canonical ranks are 0-based)
        assert set(order) == {0, 1, 2, 3, 4}

    def test_ethanol_ordering(self):
        """Ethanol (CCO) should have 9 atoms (2C + 1O + 6H)."""
        order = canonical_atom_order("CCO")

        # Should have 9 atoms total
        assert len(order) == 9

        # All indices should be unique
        assert len(set(order)) == 9

        # Ranks should be 0-8
        assert set(order) == set(range(9))

    def test_benzene_ordering(self):
        """Benzene (c1ccccc1) should have 12 atoms (6C + 6H)."""
        order = canonical_atom_order("c1ccccc1")

        # Should have 12 atoms total
        assert len(order) == 12

        # All indices should be unique
        assert len(set(order)) == 12

        # Ranks should be 0-11
        assert set(order) == set(range(12))

    def test_deterministic_ordering(self):
        """Same SMILES called multiple times should produce identical ordering."""
        smiles = "CCO"

        # Call 10 times
        orders = [canonical_atom_order(smiles) for _ in range(10)]

        # All should be identical
        first_order = orders[0]
        for order in orders[1:]:
            assert order == first_order, "Canonical ordering is not deterministic"

    def test_invalid_smiles(self):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid SMILES"):
            canonical_atom_order("not a valid smiles")


class TestMapNWChemToCanonical:
    """Test map_nwchem_to_canonical() function."""

    def test_methane_mapping(self):
        """Methane mapping should be 1-based NWChem index to canonical rank."""
        mapping = map_nwchem_to_canonical("C")

        # Should have 5 entries (1-based: 1 through 5)
        assert len(mapping) == 5
        assert set(mapping.keys()) == {1, 2, 3, 4, 5}

        # Values should be canonical ranks (0-4)
        assert set(mapping.values()) == {0, 1, 2, 3, 4}

    def test_ethanol_mapping(self):
        """Ethanol mapping should have 9 entries."""
        mapping = map_nwchem_to_canonical("CCO")

        # Should have 9 entries (1-based: 1 through 9)
        assert len(mapping) == 9
        assert set(mapping.keys()) == set(range(1, 10))

        # Values should be canonical ranks (0-8)
        assert set(mapping.values()) == set(range(9))

    def test_mapping_consistency(self):
        """Mapping should be consistent with canonical_atom_order."""
        smiles = "CCO"

        order = canonical_atom_order(smiles)
        mapping = map_nwchem_to_canonical(smiles)

        # For each NWChem 1-based index, the mapping should give
        # the canonical rank for that atom
        for nwchem_idx in mapping:
            rdkit_idx = nwchem_idx - 1  # Convert to 0-based
            canonical_rank = order[rdkit_idx]
            assert mapping[nwchem_idx] == canonical_rank


class TestGetAtomCount:
    """Test get_atom_count() function."""

    def test_methane_count(self):
        """Methane should have 5 total, 4 H, 1 heavy."""
        total, num_h, num_heavy = get_atom_count("C")
        assert total == 5
        assert num_h == 4
        assert num_heavy == 1

    def test_ethanol_count(self):
        """Ethanol should have 9 total, 6 H, 3 heavy."""
        total, num_h, num_heavy = get_atom_count("CCO")
        assert total == 9
        assert num_h == 6
        assert num_heavy == 3

    def test_benzene_count(self):
        """Benzene should have 12 total, 6 H, 6 heavy."""
        total, num_h, num_heavy = get_atom_count("c1ccccc1")
        assert total == 12
        assert num_h == 6
        assert num_heavy == 6


class TestSymmetricMolecules:
    """Test handling of symmetric molecules."""

    def test_benzene_symmetry(self):
        """Benzene canonical ordering should be stable despite symmetry."""
        # Run multiple times
        orders = [canonical_atom_order("c1ccccc1") for _ in range(5)]

        # All should be identical (deterministic despite symmetry)
        first = orders[0]
        for order in orders[1:]:
            assert order == first

    def test_different_smiles_same_molecule(self):
        """Different SMILES for same molecule should give equivalent canonical ordering.

        Note: This tests that the canonical ranking itself is consistent,
        not necessarily that atom order is identical (which depends on SMILES parsing).
        """
        # Ethanol written two ways
        smiles1 = "CCO"
        smiles2 = "OCC"

        order1 = canonical_atom_order(smiles1)
        order2 = canonical_atom_order(smiles2)

        # Both should have same length
        assert len(order1) == len(order2)

        # Both should use same set of ranks
        assert set(order1) == set(order2)

        # Note: The actual order list may differ because atom indexing
        # depends on SMILES parsing order, but the canonical ranks used
        # should be the same set
