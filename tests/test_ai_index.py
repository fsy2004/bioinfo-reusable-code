from __future__ import annotations
import json
import unittest
import tempfile
from unittest.mock import patch
from pathlib import Path
from scripts.build_ai_index import build, discover, ROOT
from scripts.query_modules import query


class AIIndexTests(unittest.TestCase):
    def test_nested_module_and_loose_template_discovery(self):
        files=['modules/17_population/01_cohorts/529_panel/529_panel.R',
               'modules/17_population/01_cohorts/529_panel/example_data/777_sample.csv',
               'modules/04_genetics/01_prep/028_old.R',
               'modules/_upstream_sources/999_ignore.py']
        groups=discover(Path('.'),files)
        self.assertEqual(len(groups),2)
        self.assertEqual(groups['modules/17_population/01_cohorts/529_panel'][0],'529')

    def test_all_indexed_paths_and_entrypoints_exist(self):
        records=json.loads(build()['modules.json'])['modules']
        self.assertGreater(len(records),180)
        self.assertEqual(len(records),len({r['id'] for r in records}))
        for row in records:
            self.assertTrue((ROOT/row['path']).exists(),row['path'])
            for entry in row['entrypoints']:
                self.assertTrue((ROOT/entry).is_file(),entry)
            self.assertTrue(row['metadata_source'])

    def test_query_finds_renamed_nhanes_and_retains_old_alias(self):
        records=json.loads(build()['modules.json'])['modules']
        matches=list(query(records,['NHANES']))
        row=next(r for r in matches if r['id']=='528')
        self.assertEqual(row['domain'],'17_population_epidemiology')
        self.assertTrue(row['aliases'])
        self.assertEqual(list(query(records,['528'],'04')),[])

    def test_build_is_deterministic(self):
        self.assertEqual(build(),build())

    def test_ignored_local_manifest_and_readme_are_not_consumed(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            framework=root/'modules/_framework'
            framework.mkdir(parents=True)
            (framework/'module_metadata.json').write_text('[{"id":"001","purpose":"public"}]',encoding='utf-8')
            (framework/'new_modules.json').write_text('[]',encoding='utf-8')
            folder=root/'modules/17_population_epidemiology/01_disease_burden/001_test'
            folder.mkdir(parents=True)
            (folder/'001_test.py').write_text('',encoding='utf-8')
            (folder/'module.json').write_text('{"summary":"private"}',encoding='utf-8')
            (folder/'README.md').write_text('private',encoding='utf-8')
            with patch('scripts.build_ai_index.visible_files',return_value=['modules/17_population_epidemiology/01_disease_burden/001_test/001_test.py']):
                row=json.loads(build(root)['modules.json'])['modules'][0]
            self.assertEqual(row['summary'],'public')
            self.assertIsNone(row['readme'])


if __name__=='__main__':
    unittest.main()
