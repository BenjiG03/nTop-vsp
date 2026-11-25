"""nTop to OpenVSP drivers: build OpenVSP models from nTop's section contours and run VSPAERO, Wave Drag and
Parasite Drag on them. The nTop custom blocks run the ``run_*`` modules with ``python -m vsp_ntop.<module>``."""
import os
import sys

__version__ = "0.4.0"


def conda_dll_dirs():
    """Put a conda Python's DLL folders on PATH. nTop starts the Python of NTOP_PYTHON_PATH without activating its
    environment, so numpy's BLAS/LAPACK DLLs are not found at its first matrix product and the process dies with no
    traceback."""
    if os.name != "nt":
        return
    dirs = [os.path.join(sys.prefix, *p) for p in (("Library", "mingw-w64", "bin"), ("Library", "usr", "bin"),
                                                   ("Library", "bin"), ("Scripts",))]
    base = getattr(sys, "base_prefix", sys.prefix)
    if base != sys.prefix:  # a venv made from a conda Python
        dirs += [os.path.join(base, *p) for p in (("Library", "mingw-w64", "bin"), ("Library", "bin"))]
    dirs = [d for d in dirs if os.path.isdir(d)]
    if not dirs:
        return
    os.environ["PATH"] = os.pathsep.join(dirs + [os.environ.get("PATH", "")])
    for d in dirs:
        try:
            os.add_dll_directory(d)
        except (AttributeError, OSError):
            pass


conda_dll_dirs()
