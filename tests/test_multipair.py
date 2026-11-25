"""Offline unit tests. Do not trigger LightDock, PyMOL or network requests."""
import csv
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

BASE = Path(__file__).resolve().parent.parent


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, BASE / relative)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


prep = load("prepare_pair_protein", "workflow/scripts/prepare_pair_protein.py")
score = load("score_pair", "workflow/scripts/score_pair.py")
build = load("build_pair_complex", "workflow/scripts/build_pair_complex.py")
download = load("download_rcsb_pdb", "workflow/scripts/download_rcsb_pdb.py")


class TestMultiPair(unittest.TestCase):
    def test_original_ensemble_data_not_modified(self):
        # Spot-check original data and ensure it remains a 600-pose dataset.
        with (BASE / 'results/e2a_hpr/ensemble/docking_scores.csv').open() as f:
            original = list(csv.DictReader(f))
        self.assertEqual(len(original), 600)
        self.assertEqual(int(original[0]['model']), 9)
        self.assertAlmostEqual(float(original[0]['score']), -4.13045687)

    def test_pair_score_parser_on_actual_lightdock_files(self):
        source = BASE / 'results/e2a_hpr/ensemble/model_9'
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)
            for idx in (0, 1):
                swarm = dest / f'swarm_{idx}'
                swarm.mkdir()
                (swarm / 'gso_10.out').write_bytes(
                    (source / f'swarm_{idx}/gso_10.out').read_bytes()
                )
            csv_file = dest / 'docking_scores.csv'
            rows = score.collect(dest, 10, csv_file)
            self.assertEqual(len(rows), 20)
            self.assertEqual(rows[0]['swarm'], 0)
            self.assertEqual(rows[0]['pose_index'], 5)
            self.assertAlmostEqual(rows[0]['score'], -4.13045687)
            with csv_file.open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 20)

    def test_chain_split_on_real_best_complex(self):
        source = BASE / 'results/e2a_hpr/ensemble/model_9/swarm_0/lightdock_5.pdb'
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'best_complex_chainsAB.pdb'
            build.rewrite_two_chain_pdb(source, out)
            lines = [x for x in out.read_text().splitlines() if x.startswith('ATOM  ')]
            self.assertEqual({x[21] for x in lines}, {'A', 'B'})
            self.assertEqual(sum(x[21] == 'A' for x in lines), 1137)
            self.assertEqual(sum(x[21] == 'B' for x in lines), 639)

    def test_pdb_chain_selection_and_cleanup(self):
        src = BASE / 'data/raw/e2a_hpr/1F3G.pdb'
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / 'receptor.pdb'
            n = prep.extract_chain(src, dest, chain='A')
            self.assertEqual(n, 1115)
            self.assertEqual(sum(x.startswith('ATOM  ') for x in dest.read_text().splitlines()), n)
            self.assertTrue(dest.read_text().endswith('END\n'))
            self.assertTrue(all(x[21] == 'A' for x in dest.read_text().splitlines() if x.startswith('ATOM  ')))
            with self.assertRaises(ValueError):
                prep.extract_chain(src, Path(tmp) / 'bad.pdb', chain='Z')

    def test_nmr_model_selection(self):
        source = BASE / 'data/raw/e2a_hpr/1HDN.pdb'
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / 'hpr_model9.pdb'
            n = prep.extract_chain(source, destination, chain='A', model=9)
            self.assertGreater(n, 600)
            self.assertEqual(n, sum(x.startswith('ATOM  ') for x in destination.read_text().splitlines()))

    def test_generic_reconstruction_command_with_mocked_lightdock(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp) / 'automated'
            swarm = base / 'swarm_0'
            swarm.mkdir(parents=True)
            src = BASE / 'results/e2a_hpr/ensemble/model_9/swarm_0/lightdock_5.pdb'
            scores = base / 'docking_scores.csv'
            with scores.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=['swarm','pose_index','score'])
                writer.writeheader()
                writer.writerow({'swarm':0,'pose_index':5,'score':-4.13045687})
            def fake_run(args, cwd, check):
                self.assertIn('lgd_generate_conformations.py', args[0])
                self.assertIn('--setup', args)
                self.assertEqual(Path(cwd), swarm)
                self.assertTrue(check)
                (swarm / 'lightdock_5.pdb').write_bytes(src.read_bytes())
            final = base / 'best_complex_chainsAB.pdb'
            with patch.object(build.subprocess, 'run', side_effect=fake_run):
                build.build(scores, base, final, steps=10, glowworms=10)
            self.assertTrue(final.is_file())
            self.assertIn(' B ', final.read_text())

    def test_mocked_rcsb_download_and_invalid_id(self):
        class FakeResponse:
            def __enter__(self): return self
            def __exit__(self, exc_type, exc, tb): pass
            def read(self): return b'HEADER   TEST\nATOM      1  N   GLY A   1       0.0     0.0     0.0   1.00 0.00           N\nEND\n'
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / '1BRS.pdb'
            with patch.object(download, 'urlopen', return_value=FakeResponse()) as req:
                download.download('1BRS', dest)
            self.assertTrue(dest.is_file())
            self.assertEqual(req.call_args.args[0], 'https://files.rcsb.org/download/1BRS.pdb')
            with self.assertRaises(ValueError):
                download.download('../hijack', dest)


if __name__ == '__main__':
    unittest.main()
