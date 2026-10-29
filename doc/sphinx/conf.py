# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from pathlib import Path
import re
import subprocess


source_root = Path(__file__).resolve().parents[2]


def git_output(*args):
    try:
        result = subprocess.run(
            ['git', *args], cwd=source_root, capture_output=True, text=True,
            check=False,
        )
    except OSError:
        return ''
    return result.stdout.strip() if result.returncode == 0 else ''


project = 'fVTKHDF'
copyright = '2026, Neil N. Carlson'
author = 'Neil N. Carlson'
# Keep CMakeLists.txt as the source of the project version. No compiler or
# configured CMake build is needed to work on the documentation.
version = re.search(
    r'project\(fVTKHDF\s+VERSION\s+(\d+\.\d+\.\d+)\b',
    (source_root / 'CMakeLists.txt').read_text(),
).group(1)
release = f'{version}-dev'
if (source_root / '.git').exists():
    if git_output('describe', '--tags', '--exact-match') == f'v{version}':
        release = version
    else:
        revision = git_output('rev-parse', '--short', 'HEAD')
        if revision:
            release += f'+g{revision}'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

#extensions = ['sphinx_rtd_theme']

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
# html_static_path = ['_static']

# get rid of navigation sidebar for single-page document
#html_sidebars = { '**': [] }

html_theme_options = {
    'nosidebar': False,
}
