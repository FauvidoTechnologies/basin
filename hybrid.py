"""Compatibility shim: `import hybrid` gives the `basin.core` module itself.

The runners in allen/ and neuromorpho/ set switches with `setattr(hybrid, ...)`; because
this module replaces itself with basin.core, those settings reach the code that reads them.
"""
import sys

from basin import core as _core

sys.modules[__name__] = _core
