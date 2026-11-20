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
"""Sphinx helpers for reaktoro-pse docs

- ``rkt-config`` directive renders options of a Pyomo ConfigDict (e.g. ReaktoroBlock
  CONFIG) from the code, so the options reference never goes out of date.
- docstrings that are not valid reStructuredText are shown as preformatted text, so
  docstring formatting can not break the docs build (it is built with warnings as
  errors).
"""
import importlib
import inspect

from docutils import nodes, utils
from docutils.frontend import get_default_settings
from docutils.parsers.rst import Directive, Parser, directives
from docutils.statemachine import StringList
from pyomo.common.config import ConfigDict

_settings = get_default_settings(Parser)
_settings.report_level = 5  # do not print anything while checking
_settings.halt_level = 5


def is_valid_rst(text):
    """True if text parses as reStructuredText without warnings or errors"""
    document = utils.new_document("<check>", _settings)
    messages = []
    document.reporter.attach_observer(
        lambda msg: messages.append(msg) if msg["level"] >= 2 else None
    )
    Parser().parse(text, document)
    # substitutions are resolved after parsing, |x| in plain text (e.g. norms) is an
    # undefined substitution
    substitutions = list(document.findall(nodes.substitution_reference))
    return not messages and not substitutions


def _is_bullet(line):
    return line.lstrip().startswith(("- ", "* "))


def _indent(line):
    return len(line) - len(line.lstrip())


def normalize_lists(text):
    """option docs are written as plain text with bullet lists that directly follow
    a paragraph, add the blank lines rst needs around lists"""
    lines = text.splitlines()
    out = []
    bullet_indent = None
    for line in lines:
        if not line.strip():
            bullet_indent = None
            out.append(line)
            continue
        if _is_bullet(line):
            if out and out[-1].strip() and bullet_indent is None:
                out.append("")
            elif bullet_indent is not None and _indent(line) != bullet_indent:
                # nested or dedented list
                out.append("")
            bullet_indent = _indent(line)
        elif bullet_indent is not None and _indent(line) <= bullet_indent:
            # paragraph after a list
            out.append("")
            bullet_indent = None
        out.append(line)
    return "\n".join(out)


def as_rst_lines(text):
    """dedented text as rst, or a literal block if it is not valid rst"""
    text = inspect.cleandoc(text or "")
    if not text:
        return []
    for candidate in (text, normalize_lists(text)):
        if is_valid_rst(candidate):
            return candidate.splitlines()
    return ["::", ""] + ["    " + line for line in text.splitlines()]


def _import(path):
    """import dotted path to an attribute, e.g. package.module.Class.CONFIG"""
    parts = path.split(".")
    for i in range(len(parts) - 1, 0, -1):
        try:
            obj = importlib.import_module(".".join(parts[:i]))
        except ImportError:
            continue
        for name in parts[i:]:
            obj = getattr(obj, name)
        return obj
    raise ImportError(f"Can not import {path}")


def _default(cv):
    try:
        value = cv.value()
    except Exception:
        return None
    if isinstance(value, ConfigDict):
        return None
    if isinstance(value, type):
        return value.__name__
    text = repr(value)
    return text if len(text) < 120 else text[:117] + "..."


def _domain(cv):
    domain = getattr(cv, "_domain", None)
    if domain is None:
        return None
    try:
        name = cv.domain_name()
    except Exception:
        name = getattr(domain, "__name__", None)
    return name


class ConfigReference(Directive):
    """render options of a ConfigDict

    .. rkt-config:: reaktoro_pse.reaktoro_block.ReaktoroBlockData.CONFIG
       :path: aqueous_phase
       :values-only:
       :exclude: bfgs_init_type, ...
    """

    required_arguments = 1
    option_spec = {
        "path": directives.unchanged,
        "values-only": directives.flag,
        "exclude": directives.unchanged,
    }

    def run(self):
        config = _import(self.arguments[0])
        if callable(config) and not isinstance(config, ConfigDict):
            config = config()
        path = self.options.get("path")
        if path:
            for key in path.split("."):
                config = config[key]
        exclude = {
            e.strip() for e in self.options.get("exclude", "").split(",") if e.strip()
        }
        values_only = "values-only" in self.options

        lines = []
        keys = list(config.keys())
        for key in keys:
            if key in exclude:
                continue
            cv = config.get(key)
            if (cv._description or "").startswith("Deprecated"):
                continue
            if key.endswith("_indexed") and key[: -len("_indexed")] in keys:
                # documented with the input it applies to
                continue
            if isinstance(cv, ConfigDict):
                if values_only:
                    continue
                lines += [
                    f"``{key}``",
                    f"   Group of options, see ``{key}`` section.",
                    "",
                ]
                continue
            lines.append(f"``{key}``")
            meta = []
            default = _default(cv)
            if default is not None:
                meta.append(f"*default:* ``{default}``")
            domain = _domain(cv)
            if domain:
                meta.append(f"*type:* ``{domain}``")
            if f"{key}_indexed" in keys:
                indexed_default = _default(config.get(f"{key}_indexed"))
                meta.append(f"*indexed option:* ``{key}_indexed={indexed_default}``")
            body = []
            if meta:
                body += [" | ".join(meta), ""]
            if cv._description:
                body += as_rst_lines(cv._description)
                body.append("")
            if cv._doc and cv._doc.strip() != (cv._description or "").strip():
                body += as_rst_lines(cv._doc)
                body.append("")
            lines += ["   " + line if line else "" for line in body]
            lines.append("")

        node = nodes.container(classes=["rkt-config"])
        self.state.nested_parse(StringList(lines), self.content_offset, node)
        return [node]


def safe_docstring(app, what, name, obj, options, lines):
    """IDAES process block classes get a generated docstring with their CONFIG,
    options are documented in the options reference instead. Other docstrings that
    are not valid rst are shown as preformatted text."""
    if not name.startswith("reaktoro_pse"):
        return
    if what == "class" and "Config args" in "\n".join(lines):
        lines[:] = [
            lines[0] if lines else "",
            "",
            "Options are documented in :doc:`/technical_reference/reaktoro_block_options`.",
        ]
        return
    text = "\n".join(lines)
    if text.strip() and not is_valid_rst(text):
        lines[:] = ["::", ""] + ["    " + line for line in lines]


def setup(app):
    app.add_directive("rkt-config", ConfigReference)
    # runs after napoleon (registered first), so google style sections are converted
    app.connect("autodoc-process-docstring", safe_docstring, priority=900)
    return {"parallel_read_safe": True}
