#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/reaktoro-pse/"
#################################################################################
# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import shutil
import sys
from importlib import metadata

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "_ext"))

# -- Project information -----------------------------------------------------

project = "Reaktoro-PSE"
copyright = "2024-2026, WaterTAP"
author = "WaterTAP"
release = version = metadata.version("reaktoro-pse")

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.napoleon",  # Google and NumPy-style docstrings
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.doctest",
    "nbsphinx",  # Jupyter notebooks as docs
    "rkt_docs",  # ReaktoroBlock options reference, docstring checks (docs/_ext)
]

mathjax3_config = {"chtml": {"displayAlign": "left", "displayIndent": "2em"}}

# labels are prefixed with document name, so headings in different pages do not clash
autosectionlabel_prefix_document = True
autodoc_warningiserror = False
autodoc_member_order = "bysource"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "apidoc/*tests*"]

# tutorials are rendered with outputs stored in notebooks, they are executed in CI
# with nbmake (see .github/workflows/checks.yml)
nbsphinx_execute = "never"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pyomo": ("https://pyomo.readthedocs.io/en/stable/", None),
    "idaes": ("https://idaes-pse.readthedocs.io/en/stable/", None),
    "watertap": ("https://watertap.readthedocs.io/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "images/reaktoro_pse.png"
html_theme_options = {"logo_only": False}

# -- Build hooks ---------------------------------------------------------------

DOCS_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(DOCS_DIR, "..", "src", "reaktoro_pse")


def run_apidoc(*args):
    """generate module reference (same approach as WaterTAP), set SKIP_APIDOC
    environment variable to skip it for faster local builds"""
    # NOTE the env var must be set before importing apidoc, or the options
    # will have no effect
    os.environ["SPHINX_APIDOC_OPTIONS"] = "members,show-inheritance"
    from sphinx.ext import apidoc

    apidoc.main(
        [
            SRC_DIR,
            os.path.join(SRC_DIR, "*tests"),
            os.path.join(SRC_DIR, "tutorials"),
            "-o",
            os.path.join(DOCS_DIR, "apidoc"),
            "--force",
            # parallel_tools, reaktoro_block_config, examples have no __init__.py
            "--implicit-namespaces",
        ]
    )


def copy_tutorials(*args):
    """notebooks live in the package (so they are tested with nbmake), Sphinx needs
    them inside docs source tree, so they are copied (with their images) on build"""
    src = os.path.join(SRC_DIR, "tutorials")
    dst = os.path.join(DOCS_DIR, "tutorials", "notebooks")
    shutil.copytree(src, dst, dirs_exist_ok=True)


def skip(app, what, name, obj, would_skip, options):
    """Do not skip constructors!"""
    if name == "__init__":
        return False
    return would_skip


def setup(app):
    copy_tutorials()
    if os.environ.get("SKIP_APIDOC", False):
        print("Skipping apidoc")
    else:
        app.connect("builder-inited", run_apidoc)
    app.connect("autodoc-skip-member", skip)
