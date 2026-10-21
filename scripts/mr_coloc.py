#!/usr/bin/env python3
"""Repository entry point for module611; reads only its registry and requested sources."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).resolve().parents[1] / 'modules/04_causal_inference_genetics/09_method_library/611_mr_coloc_library/method_library.py'), run_name='__main__')
