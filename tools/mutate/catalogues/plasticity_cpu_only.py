"""The same 68 substitutions as `plasticity.py`, against a build with no Vulkan.

§7 was written before §8 existed. The suite it was measured against had no GPU
solid-shell kernel in it and therefore no CPU-against-device cross-check, and that
check is one of the strongest instruments in the suite: it compares every
integration point of the element the plasticity feeds. Several mutants below die
on it *first*, which means today's suite would kill some mutants §7's could not,
and a rate re-derived against today's suite is not the same measurement.

Configuring the copy with `CMAKE_DISABLE_FIND_PACKAGE_Vulkan` removes exactly that
half of the suite, so the question "which of these kills does the suite owe to a
test that did not exist when the figure was published?" has an answer rather than
a caveat. Nothing else changes: same tree, same catalogue, same bound.
"""
import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location(
    "_plasticity_catalogue", pathlib.Path(__file__).with_name("plasticity.py"))
_source = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_source)

NAME = "plasticity-cpu-only"
MUTANTS = _source.MUTANTS
CONFIGURE = ["cmake", "-S", "{src}", "-B", "{build}", "-G", "Ninja",
             "-DCMAKE_BUILD_TYPE=RelWithDebInfo", "-DCMAKE_DISABLE_FIND_PACKAGE_Vulkan=ON"]
