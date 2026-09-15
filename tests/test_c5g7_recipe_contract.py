from __future__ import annotations

from pathlib import Path
import unittest

import h5py


EXPECTED_SCATTER_MGXS = "consistent scatter matrix"


class C5G7RecipeContractTests(unittest.TestCase):
    def test_recipe_declares_locked_scatter_mgxs_type(self) -> None:
        text = (_repo_root() / "scripts" / "c5g7_export_recipe.py").read_text(
            encoding="utf-8"
        )

        self.assertIn(f'"{EXPECTED_SCATTER_MGXS}"', text)
        self.assertNotIn('"reduced absorption"', text)
        self.assertIn("def scatter_mgxs_type", text)

    def test_direct_export_declares_ordinary_consistent_scatter(self) -> None:
        text = (_repo_root() / "scripts" / "export_c5g7_statepoint.py").read_text(
            encoding="utf-8"
        )

        self.assertIn(f'"{EXPECTED_SCATTER_MGXS}"', text)
        self.assertNotIn('"reduced absorption"', text)

    def test_accepted_hdf5_records_locked_scatter_mgxs_type(self) -> None:
        path = (
            _repo_root()
            / "examples"
            / "donjon_openmc2donjon"
            / "c5g7_assembly_p1_adf_production.h5"
        )

        with h5py.File(path, "r") as h5:
            self.assertEqual(h5.attrs["openmc_scatter_mgxs_type"], EXPECTED_SCATTER_MGXS)
            self.assertFalse(bool(h5.attrs["openmc_scatter_multiplicity_weighted"]))
            self.assertEqual(
                h5.attrs["openmc_scatter_balance_dataset"], "absorption"
            )
            self.assertEqual(
                h5.attrs["scatter_contract_source_statepoint_sha256"],
                "d71664009ddbb562c8fcb88639aba2282c7bb0036526be77ece7f9c640bc2b8a",
            )
            self.assertEqual(
                float(
                    h5.attrs[
                        "scatter_contract_scatter_nu_scatter_max_abs_rate_difference"
                    ]
                ),
                0.0,
            )
            for group in h5["mixtures"].values():
                self.assertEqual(
                    group.attrs["openmc_scatter_mgxs_type"],
                    EXPECTED_SCATTER_MGXS,
                )
                self.assertFalse(
                    bool(group.attrs["openmc_scatter_multiplicity_weighted"])
                )
                self.assertEqual(
                    group.attrs["openmc_scatter_balance_dataset"], "absorption"
                )

    def test_current_c5g7_docs_record_ordinary_contract_and_replay_boundary(
        self,
    ) -> None:
        for relative_path in (
            "docs/HANDOFF_NOTE.md",
            "docs/HANDOFF_SNAPSHOT.md",
            "docs/VALIDATION.md",
        ):
            text = (_repo_root() / relative_path).read_text(encoding="utf-8")
            with self.subTest(path=relative_path):
                self.assertIn("`consistent scatter matrix`", text)
                self.assertIn("transport-complete", text)
                self.assertIn("statepoint", text)
                self.assertIn("regenerat", text.lower())


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


if __name__ == "__main__":
    unittest.main()
