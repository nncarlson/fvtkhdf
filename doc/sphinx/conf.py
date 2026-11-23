# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from pathlib import Path
import re
import sys


source_root = Path(__file__).resolve().parents[2]


project = 'fVTKHDF'
copyright = '2026, Neil N. Carlson'
author = 'Neil N. Carlson'
# Keep CMakeLists.txt as the source of the project version. No compiler or
# configured CMake build is needed to work on the documentation.
version = re.search(
    r'project\(fVTKHDF\s+VERSION\s+(\d+\.\d+\.\d+)\b',
    (source_root / 'CMakeLists.txt').read_text(),
).group(1)
# This documentation branch describes the source at the release tag.
release = version
repo_revision = f'v{version}'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, str(Path(__file__).resolve().parent / '_ext'))
extensions = ['method_doc', 'sphinx.ext.extlinks']

extlinks = {
    'repo': (f'https://github.com/nncarlson/fvtkhdf/blob/{repo_revision}/%s', '%s'),
}

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'sphinx_rtd_theme'
#html_theme = 'sphinx_book_theme'
html_theme = 'furo'

html_sidebars = {
    "index": [] # This removes the sidebar ONLY from the index.html page
}
html_static_path = ['_static']
html_css_files = ['custom.css']

# get rid of navigation sidebar for single-page document
#html_sidebars = { '**': [] }

html_theme_options = {
    'nosidebar': False,
}
