from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from openmc2donjon.handoff_summary import handoff_summary_payload
from openmc2donjon.scatter_contract import scatter_contract_from_preflight


class ScatterContractTests(unittest.TestCase):
    def test_extracts_complete_common_contract(self) -> None:
        contract = {
            "openmc_scatter_mgxs_type": "consistent nu-scatter matrix",
            "openmc_scatter_multiplicity_weighted": True,
            "openmc_scatter_balance_dataset": "reduced_absorption",
            "openmc_scatter_contract_declared": True,
            "openmc_scatter_contract_valid": True,
            "openmc_transport_mgxs_type": "nu-transport",
            "openmc_transport_contract_declared": True,
        }
        preflight = {"inputs": [{"ok": True, **contract}]}

        self.assertEqual(scatter_contract_from_preflight(preflight), contract)

    def test_refuses_partial_or_conflicting_contracts(self) -> None:
        partial = {
            "inputs": [
                {
                    "openmc_scatter_mgxs_type": "scatter matrix",
                    "openmc_scatter_multiplicity_weighted": False,
                }
            ]
        }
        conflicting = {
            "inputs": [
                {
                    "openmc_scatter_mgxs_type": "scatter matrix",
                    "openmc_scatter_multiplicity_weighted": False,
                    "openmc_scatter_balance_dataset": "absorption",
                    "openmc_scatter_contract_declared": True,
                    "openmc_scatter_contract_valid": True,
                    "openmc_transport_mgxs_type": "transport",
                    "openmc_transport_contract_declared": True,
                },
                {
                    "openmc_scatter_mgxs_type": "consistent nu-scatter matrix",
                    "openmc_scatter_multiplicity_weighted": True,
                    "openmc_scatter_balance_dataset": "reduced_absorption",
                    "openmc_scatter_contract_declared": True,
                    "openmc_scatter_contract_valid": True,
                    "openmc_transport_mgxs_type": "nu-transport",
                    "openmc_transport_contract_declared": True,
                },
            ]
        }

        self.assertIsNone(scatter_contract_from_preflight(partial))
        self.assertIsNone(scatter_contract_from_preflight(conflicting))

    def test_preserves_explicit_legacy_resolution(self) -> None:
        contract = {
            "openmc_scatter_mgxs_type": None,
            "openmc_scatter_multiplicity_weighted": False,
            "openmc_scatter_balance_dataset": "absorption",
            "openmc_scatter_contract_declared": False,
            "openmc_scatter_contract_valid": True,
            "openmc_transport_mgxs_type": "transport",
            "openmc_transport_contract_declared": False,
        }

        self.assertEqual(
            scatter_contract_from_preflight({"inputs": [contract]}),
            contract,
        )

    def test_handoff_summary_surfaces_hash_bundled_check_contract(self) -> None:
        contract = {
            "openmc_scatter_mgxs_type": "consistent nu-scatter matrix",
            "openmc_scatter_multiplicity_weighted": True,
            "openmc_scatter_balance_dataset": "reduced_absorption",
            "openmc_scatter_contract_declared": True,
            "openmc_scatter_contract_valid": True,
            "openmc_transport_mgxs_type": "nu-transport",
            "openmc_transport_contract_declared": True,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            check_summary = root / "check_summary.json"
            manifest = root / "manifest.json"
            check_summary.write_text(
                json.dumps({"inputs": [{"ok": True, **contract}]}),
                encoding="utf-8",
            )
            manifest.write_text(
                json.dumps({"artifact_count": 1, "artifacts": []}),
                encoding="utf-8",
            )

            payload = handoff_summary_payload(
                package_version="test",
                run_dir=root,
                recipe_path=root / "recipe.py",
                statepoint_path=None,
                hdf5_path=root / "mgxs.h5",
                output_path=root / "out.macrolib.txt",
                output_format="macrolib",
                summary={},
                run_summary_json=None,
                check_summary_json=check_summary,
                manifest_path=manifest,
                bundle_validation_summary_json=None,
                bundle_validation_passed=True,
                bundle_validation_decision=None,
                adf_enabled=False,
                sph_enabled=False,
            )

        self.assertEqual(payload["scatter_contract"], contract)

    def test_refuses_a_resolved_but_invalid_contract(self) -> None:
        preflight = {
            "inputs": [
                {
                    "openmc_scatter_mgxs_type": "consistent nu-scatter matrix",
                    "openmc_scatter_multiplicity_weighted": True,
                    "openmc_scatter_balance_dataset": "reduced_absorption",
                    "openmc_scatter_contract_declared": True,
                    "openmc_scatter_contract_valid": False,
                    "openmc_transport_mgxs_type": "transport",
                    "openmc_transport_contract_declared": True,
                }
            ]
        }

        self.assertIsNone(scatter_contract_from_preflight(preflight))


if __name__ == "__main__":
    unittest.main()
