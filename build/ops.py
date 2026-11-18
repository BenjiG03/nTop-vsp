"""Thin typed wrappers around nTop block signatures, so graph code reads like the maths.

Every function adds one nTop block instance to a ``Graph`` and returns its id. Nothing here computes
a value: the arithmetic happens when nTop runs the notebook.
"""
from __future__ import annotations

from ntopfile import Graph, ModelInput, Ref  # noqa: F401

LENGTH = {"length": 1}
AREA = {"length": 2}
ANGLE = {"angle": 1}


class Ops:
    """Block builders bound to one graph. Helpers are hidden (top=False) unless named."""

    def __init__(self, g: Graph):
        self.g = g
        if not hasattr(g, "_constants"):
            g._constants = {}


    def _b(self, func, typ, inputs, name=None):
        return self.g.block(func, typ, list(inputs), name=name, top=name is not None)

    # --- constants ---------------------------------------------------------------------------
    def real(self, v, units=None):
        """A constant, shared with every other use of the same value and units in this graph."""
        key = (float(v), tuple(sorted((units or {}).items())))
        if key not in self.g._constants:
            self.g._constants[key] = self.g.real(v, units)
        return self.g._constants[key]

    def field(self, expr, units=None):
        """A field from nTop's expression language (x, y, z, + - * / ^, sin cos tan sqrt log min max)."""
        return self.g.leaf("real_field", {"expression": expr, "units": units or {}})

    def vector(self, x, y, z):
        return self.g.vector(x, y, z)

    # --- scalars -----------------------------------------------------------------------------
    def add(self, a, b, name=None):
        return self._b("add<real,real>", "real", [a, b], name)

    def sub(self, a, b, name=None):
        return self._b("subtract<real,real>", "real", [a, b], name)

    def mul(self, a, b, name=None):
        return self._b("multiply<real,real>", "real", [a, b], name)

    def div(self, a, b, name=None):
        return self._b("divide<real,real>", "real", [a, b], name)

    def pow(self, a, b, name=None):
        return self._b("pow<real,real>", "real", [a, b], name)

    def sq(self, a, name=None):
        return self._b("square<real>", "real", [a], name)

    def sqrt(self, a, name=None):
        return self._b("sqrt<real>", "real", [a], name)

    def sin(self, a, name=None):
        return self._b("sin<real>", "real", [a], name)

    def cos(self, a, name=None):
        return self._b("cos<real>", "real", [a], name)

    def max(self, a, b, name=None):
        return self._b("max<real,real>", "real", [a, b], name)

    def min(self, a, b, name=None):
        return self._b("min<real,real>", "real", [a, b], name)

    def select(self, choice, items, name=None):
        """The real in ``items`` picked by a dropdown input (constants or scalar inputs only)."""
        values = self._b("core.list<real>", "list<real>", items)
        return self._b("core.select_by_choice<choice,list_interface>", "real", [choice, values], name)

    def pick(self, choice, items, name=None):
        """As ``select``, for items that may be lists (one value per flight condition): the sum of each item times
        a 0/1 weight picked by the dropdown. A list cannot be an item of a list in nTop, so ``select`` cannot take
        them. The items must be finite."""
        out = None
        for i, item in enumerate(items):
            w = self.select(choice, [self.real(1.0 if j == i else 0.0) for j in range(len(items))])
            last = i == len(items) - 1
            term = self.mul(item, w, name=name if last and out is None else None)
            out = term if out is None else self.add(out, term, name=name if last else None)
        return out

    def supersonic(self, mach):
        """(weight, Mach for a supersonic formula, Mach for a subsonic one): the weight is 1 above Mach 1 and 0
        at or below it; each formula gets a Mach number on its own side of 1, so neither returns NaN."""
        one = self.real(1.0)
        w = self.min(self.max(self.mul(self.sub(mach, one), self.real(1e9)), self.real(0.0)), one)
        return w, self.max(mach, one), self.min(mach, one)

    def blend(self, w, a, b, name=None):
        """w a + (1 - w) b."""
        return self.add(self.mul(w, a), self.mul(self.sub(self.real(1.0), w), b), name=name)

    # --- vectors -----------------------------------------------------------------------------
    def vadd(self, a, b, name=None):
        return self._b("add<vector,vector>", "vector", [a, b], name)

    def vsub(self, a, b, name=None):
        return self._b("subtract<vector,vector>", "vector", [a, b], name)

    def vscale(self, a, s, name=None):
        return self._b("multiply<vector,real>", "vector", [a, s], name)

    def vdiv(self, a, s, name=None):
        return self._b("divide<vector,real>", "vector", [a, s], name)

    def dot(self, a, b, name=None):
        return self._b("dot_product<vector,vector>", "real", [a, b], name)

    def cross(self, a, b, name=None):
        return self._b("cross_product<vector,vector>", "vector", [a, b], name)

    def unit(self, a, name=None):
        return self.vdiv(a, self.sqrt(self.dot(a, a)), name)

    @staticmethod
    def comp(ref, c):
        """Component ``c`` ('x', 'y' or 'z') of a vector or point instance or model input."""
        if isinstance(ref, ModelInput):
            return ModelInput(ref.index, ref.prop + (c,))
        r = ref if isinstance(ref, Ref) else Ref(ref)
        return Ref(r.instance, r.prop + (c,))

    # --- scalar fields -----------------------------------------------------------------------
    def fadd(self, a, b, name=None):
        return self._b("add<real_field,real_field>", "real_field", [a, b], name)

    def fsub(self, a, b, name=None):
        return self._b("subtract<real_field,real_field>", "real_field", [a, b], name)

    def fmul(self, a, b, name=None):
        return self._b("multiply<real_field,real_field>", "real_field", [a, b], name)

    def fdiv(self, a, b, name=None):
        return self._b("divide<real_field,real_field>", "real_field", [a, b], name)

    def fmax(self, a, b, name=None):
        return self._b("max<real_field,real_field>", "real_field", [a, b], name)

    def fmin(self, a, b, name=None):
        return self._b("min<real_field,real_field>", "real_field", [a, b], name)

    def fsq(self, a, name=None):
        return self._b("square<real_field>", "real_field", [a], name)

    def fsqrt(self, a, name=None):
        return self._b("sqrt<real_field>", "real_field", [a], name)

    def fabs(self, a, name=None):
        return self._b("abs<real_field>", "real_field", [a], name)

    def fpow(self, a, b, name=None):
        return self._b("pow<real_field,real_field>", "real_field", [a, b], name)

    def flog(self, a, name=None):
        return self._b("log<real_field>", "real_field", [a], name)

    def fasin(self, a, name=None):
        return self._b("asin<real_field>", "real_field", [a], name)

    def fclamp(self, a, lo, hi, name=None):
        return self._b("clamp<real_field,real_field,real_field>", "real_field", [a, lo, hi], name)

    def remap(self, f, x, y, z, name=None):
        """f evaluated at (x(p), y(p), z(p)) instead of p."""
        return self._b("remap<real_field,real_field,real_field,real_field>", "real_field", [f, x, y, z], name)

    # --- vector fields -----------------------------------------------------------------------
    def vfield(self, x, y, z, name=None):
        return self._b("vector_field_from_components<real_field,real_field,real_field>", "vector_field", [x, y, z], name)

    def vconst(self, v, name=None):
        """A uniform vector field from a vector instance."""
        return self.vfield(self.comp(v, "x"), self.comp(v, "y"), self.comp(v, "z"), name)

    def grad(self, f, name=None):
        return self._b("vector_field_from_gradient<real_field>", "vector_field", [f], name)

    def fdot(self, a, b, name=None):
        return self._b("dot_product<vector_field,vector_field>", "real_field", [a, b], name)

    def fcross(self, a, b, name=None):
        return self._b("cross_product<vector_field,vector_field>", "vector_field", [a, b], name)

    def vfscale(self, s, v, name=None):
        return self._b("multiply<real_field,vector_field>", "vector_field", [s, v], name)

    def vfsub(self, a, b, name=None):
        return self._b("subtract<vector_field,vector_field>", "vector_field", [a, b], name)

    # --- integrals ---------------------------------------------------------------------------
    def flux(self, body, v, tol, name=None):
        """Integral over the body's surface of v . n dA (nTop's Surface Integral of a vector field)."""
        return self._b("surface_integral<implicit,vector_field,real>", "real", [body, v, tol], name)

    def surface_integral(self, body, f, tol, name=None):
        return self._b("surface_integral<implicit,real_field,real>", "real", [body, f, tol], name)


class More(Ops):
    """Further scalar and list blocks. Any scalar block also maps over lists (nTop's list processing)."""

    def atan(self, a, name=None):
        return self._b("atan<real>", "real", [a], name)

    def atan2(self, y, x, name=None):
        return self._b("atan2<real,real>", "real", [y, x], name)

    def asin(self, a, name=None):
        return self._b("asin<real>", "real", [a], name)

    def tan(self, a, name=None):
        return self._b("tan<real>", "real", [a], name)

    def exp(self, a, name=None):
        return self._b("exp<real>", "real", [a], name)

    def log(self, a, name=None):
        return self._b("log<real>", "real", [a], name)

    def abs(self, a, name=None):
        return self._b("abs<real>", "real", [a], name)

    def clamp(self, a, lo, hi, name=None):
        return self._b("clamp<real,real,real>", "real", [a, lo, hi], name)

    def neg(self, a, name=None):
        return self.mul(a, self.real(-1.0), name)

    def lst(self, items, name=None):
        return self._b("core.list<real>", "list<real>", items, name)

    def seq(self, start, end, count, name=None):
        return self._b("sequence_from_bounds<real,real,integer>", "list<real>", [start, end, self.g.integer(count)], name)

    def cat(self, a, b, name=None):
        return self._b("core.concatenate_lists<list_interface,list_interface>", "list<real>", [a, b], name)

    def lmax(self, a, name=None):
        return self._b("max<list<real>>", "real", [a], name)

    def lsum(self, a, count, name=None):
        """Sum of a list of ``count`` reals (count known when the graph is built)."""
        return self.mul(self._b("average<list<real>>", "real", [a]), self.real(float(count)), name)

    @staticmethod
    def at(ref, i):
        """Element i of a list instance."""
        r = ref if isinstance(ref, Ref) else Ref(ref)
        return Ref(r.instance, r.prop + (int(i),))

    def interp(self, x, xs, ys, clamp=True, name=None):
        """Piecewise-linear table lookup (nTop's Transfer Function); clamped or linear outside the table."""
        return self._b("transfer_function<real,list<real>,list<real>,interpolation_enum,extrapolation_enum>", "real",
                       [x, xs, ys, self.g.enum("interpolation_enum", 0), self.g.enum("extrapolation_enum", 1 if clamp else 0)], name)

    def table(self, values, name=None):
        return self.lst([self.real(float(v)) for v in values], name)

    def point(self, x, y, z, name=None):
        return self._b("point<real,real,real>", "point", [x, y, z], name)

    def var(self, typ, src, name=None):
        return self._b(f"core.var<{typ}>", typ, [src], name)
