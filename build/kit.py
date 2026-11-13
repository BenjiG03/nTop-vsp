"""Building pieces shared by the nTop -> OpenVSP custom blocks.

Every block here turns its inputs into lines of ``key = value`` text (SI units: metres, square metres, degrees,
m/s, kg/m^3, seconds), and the analysis blocks write those lines to a settings file, export the geometry, and
call a Python driver in ``vsp_ntop`` with Run Command. Readers bring the driver's CSV files back with Import
Table and Import Scalar Point Map.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass

from ntopfile import CustomBlock, Graph, InputSpec, ModelInput, Ref

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))).replace("\\", "/")
DEFAULT_VSP_DIR = ROOT + "/OpenVSP-3.56.0-win64"
DEFAULT_DRIVER_DIR = ROOT + "/vsp_ntop"
DEFAULT_WORK_DIR = ROOT + "/runs"
DEFAULT_PYTHON = "C:/Users/Benji/anaconda3/python.exe"

LENGTH = {"length": 1}
AREA = {"length": 2}
ANGLE = {"angle": 1}
TIME = {"time": 1}
VELOCITY = {"length": 1, "time": -1}
DENSITY = {"mass": 1, "length": -3}
PRESSURE = {"mass": 1, "length": -1, "time": -2}

RUN_COMMAND = "run_command<text,bool,bool,real,file_path>[5.44.0]"
MD = "$markdown-enabled$"
PRECISION = 12


# --- input defaults ---------------------------------------------------------------------------------------------

def d_real(v, units=None):
    return lambda g: g.real(v, units)


def d_int(v):
    return lambda g: g.integer(v)


def d_text(s):
    return lambda g: g.text(s)


def d_path(p):
    return lambda g: g.file_path(p)


def d_choice(labels, selected=0):
    return lambda g: g.choice(list(labels), selected)


def d_point_mm(x, y, z):
    return lambda g: g.point(x, y, z)


def d_angle_deg(v):
    return lambda g: g.real(math.radians(v), ANGLE)


# --- settings fields ----------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Field:
    """One block input and how it is written to the settings file.

    kind: "real" (value / reference unit), "int", "choice" (writes ``values[i]``), "point" (three keys _x _y _z in
    metres), "text", "points" (list<point>: three keys _x _y _z, values joined by ';'), "reals" (list<real>).
    """

    spec: InputSpec
    key: str
    kind: str = "real"
    reference: tuple = (1.0, None)  # (value, units) of one SI unit for "real"
    values: tuple = ()  # for "choice"
    default: object = None  # the default as written to the settings file (SI; choice: index)

    def default_lines(self):
        d = self.default
        if self.kind == "choice":
            return [f"{self.key} = {self.values[d]}"]
        if self.kind == "point":
            return [f"{self.key}_{c} = {_num(v)}" for c, v in zip("xyz", d)]
        if self.kind == "points":
            return [f"{self.key}_{c} = " + ";".join(_num(p[i]) for p in d) for i, c in enumerate("xyz")]
        if self.kind == "reals":
            return [f"{self.key} = " + ";".join(_num(v) for v in d)]
        return [f"{self.key} = {_num(d) if not isinstance(d, str) else d}"]


def _num(v):
    return f"{float(v):.12g}"


def default_text(fields):
    """The text a settings block outputs with every input at its default."""
    return "".join(line + chr(10) for f in fields for line in f.default_lines())


def F_real(name, key, default, units=None, description="", reference=None):
    ref = reference or (1.0, units)
    return Field(InputSpec(name, "real", description, units, d_real(default, units)), key, "real", ref, default=default)


def F_length(name, key, default_m, description=""):
    return Field(InputSpec(name, "real", description, LENGTH, d_real(default_m, LENGTH)), key, "real", (1.0, LENGTH),
                 default=default_m)


def F_angle(name, key, default_deg, description=""):
    return Field(InputSpec(name, "real", description, ANGLE, d_angle_deg(default_deg)), key, "real",
                 (math.radians(1.0), ANGLE), default=default_deg)


def F_int(name, key, default, description=""):
    return Field(InputSpec(name, "integer", description, None, d_int(default)), key, "int", default=default)


def F_choice(name, key, labels, values, default=0, description=""):
    return Field(InputSpec(name, "choice", description, None, d_choice(labels, default)), key, "choice",
                 values=tuple(str(v) for v in values), default=default)


def F_switch(name, key, default_on=False, description="", labels=("Off", "On")):
    return F_choice(name, key, labels, (0, 1), 1 if default_on else 0, description)


def F_text(name, key, default, description=""):
    return Field(InputSpec(name, "text", description, None, d_text(default)), key, "text", default=default)


def F_point(name, key, default_mm=(0, 0, 0), description=""):
    return Field(InputSpec(name, "point", description, LENGTH, d_point_mm(*default_mm)), key, "point",
                 default=tuple(v / 1000.0 for v in default_mm))


# --- graph helpers ----------------------------------------------------------------------------------------------

class Kit:
    """Text and file plumbing on one graph."""

    def __init__(self, g: Graph):
        self.g = g
        self._texts = {}

    def t(self, s):
        """A shared text constant."""
        if s not in self._texts:
            self._texts[s] = self.g.text(s)
        return self._texts[s]

    def cat(self, parts, name=None, sep=""):
        items = [self.t(p) if isinstance(p, str) else p for p in parts]
        return self.g.concat(items, separator=sep, top=name is not None, name=name)

    def num(self, value, ref_value=1.0, ref_units=None):
        """Text of ``value`` expressed in the unit ``ref_value ref_units`` (nTop converts units)."""
        ref = self.g.real(ref_value, ref_units)
        return self.g.block("text_from_scalar<real,real,integer>", "text", [value, ref, self.g.integer(PRECISION)],
                            top=False)

    def pick(self, choice, values):
        """The text in ``values`` picked by a dropdown."""
        lst = self.g.block("core.list<text>", "list<text>", [self.t(v) for v in values], top=False)
        return self.g.block("core.select_by_choice<choice,list_interface>", "text", [choice, lst], top=False)

    def join(self, list_of_text, sep=";"):
        return self.g.block("concatenate<list<text>,text>", "text", [list_of_text, self.t(sep)], top=False)

    def line(self, key, value_text):
        return self.cat([key + " = ", value_text, "\n"])

    def field_lines(self, fields, inp):
        """Settings lines for ``fields``; ``inp(field)`` gives the instance or model input of its value."""
        out = []
        for f in fields:
            v = inp(f)
            if f.kind == "real":
                out.append(self.line(f.key, self.num(v, *f.reference)))
            elif f.kind == "int":
                out.append(self.line(f.key, self.num(v)))
            elif f.kind == "choice":
                out.append(self.line(f.key, self.pick(v, f.values)))
            elif f.kind == "text":
                out.append(self.line(f.key, v))
            elif f.kind == "point":
                for c in "xyz":
                    out.append(self.line(f"{f.key}_{c}", self.num(prop(v, c), 1.0, LENGTH)))
            elif f.kind == "points":
                for c in "xyz":
                    out.append(self.line(f"{f.key}_{c}", self.join(self.num(prop(v, c), 1.0, LENGTH))))
            elif f.kind == "reals":
                out.append(self.line(f.key, self.join(self.num(v))))
            else:
                raise ValueError(f.kind)
        return out

    def run_command(self, command, workdir, timeout, name="Run Command"):
        g = self.g
        return g.block(RUN_COMMAND, "command_output", [command, g.boolean(True), g.boolean(True), timeout, workdir], name)

    def wait(self, wait_for, value, typ="text", name=None):
        return self.g.block("await<any,any>", typ, [wait_for, value], name=name, top=name is not None)

    def check(self, rc, what, name="Check"):
        """Stop with the driver's error text when its exit code is not zero."""
        g = self.g
        ok = g.block("equals<real,real>", "bool", [prop(rc, "exit code"), g.real(0.0)], top=False)
        msg = self.cat([what + " failed: ", prop(rc, "error")])
        return g.block("core.validate<bool,text,message_type_enum>", "text", [ok, msg, g.enum("message_type_enum", 1)], name)


def prop(ref, name):
    """Property ``name`` of an instance or model input (nTop maps it over lists)."""
    if isinstance(ref, ModelInput):
        return ModelInput(ref.index, ref.prop + (name,))
    r = ref if isinstance(ref, Ref) else Ref(ref)
    return Ref(r.instance, r.prop + (name,))


# --- inputs every analysis block has --------------------------------------------------------------------------

def run_inputs(default_case):
    return (
        InputSpec("Working Directory", "file_path",
                  "Folder for every file of the run (created if missing). Results stay there and the reader blocks "
                  "load them from it.", default=d_path(DEFAULT_WORK_DIR)),
        InputSpec("Case Name", "text", "Base name of the run's files. Give each concurrent run its own name.",
                  default=d_text(default_case)),
        InputSpec("OpenVSP Directory", "file_path", "Folder holding vsp.exe, vspscript.exe and vspaero.exe "
                  "(OpenVSP 3.56).", default=d_path(DEFAULT_VSP_DIR)),
        InputSpec("Python", "file_path", "Any Python 3.8+ interpreter; the driver scripts need only its standard library.",
                  default=d_path(DEFAULT_PYTHON)),
        InputSpec("Driver Directory", "file_path", "Folder holding the nTop-VSP driver scripts (vsp_ntop).",
                  default=d_path(DEFAULT_DRIVER_DIR)),
        InputSpec("Timeout", "real", "Longest time each external program may run.", TIME, d_real(7200.0, TIME)),
    )


RUN_INPUT_NAMES = ("Working Directory", "Case Name", "OpenVSP Directory", "Python", "Driver Directory", "Timeout")


def analysis(g: Graph, k: Kit, inp, script, exports, settings_lines, what):
    """The shared tail of every analysis block.

    ``exports(after_mkdir_base)`` adds the geometry export blocks and returns the instances to wait for; the
    settings file gets ``settings_lines`` plus the run's own keys. Returns the results base path (text),
    available once the driver has finished without error.
    """
    g.comment("**Run folder.** The Working Directory is created if it does not exist; every file of this run "
              "is named after the Case Name inside it.")
    work = k.cat([inp("Working Directory")], name="Working Directory Text")
    base = k.cat([work, "/", inp("Case Name")], name="Results Base")
    mkdir = k.run_command(k.cat(['"', inp("Python"), '" -c "import os,sys; os.makedirs(sys.argv[1], exist_ok=True)" "',
                                 work, '"']), inp("OpenVSP Directory"), inp("Timeout"), name="Create Working Directory")
    ready = k.wait(mkdir, base, name="Base After Folder Exists")
    waits = exports(ready)
    g.comment("**Settings file.** Every input becomes one `key = value` line in SI units (lengths in metres, angles "
              "in degrees). The driver script reads this file; it is kept in the Working Directory as a record of "
              "the run.")
    own = [
        k.line("work_dir", work), k.line("case", inp("Case Name")), k.line("vsp_dir", k.cat([inp("OpenVSP Directory")])),
        k.line("timeout_s", k.num(inp("Timeout"), 1.0, TIME)),
    ]
    text = g.block("core.list<text>", "list<text>", own + list(settings_lines), top=False)
    settings = g.block("concatenate<list<text>,text>", "text", [text, k.t("")], "Settings")
    settings_path = k.cat([ready, ".ntvsp.settings.txt"])
    written = g.block("export_text<file_path,text>", "file_path", [settings_path, settings], "Write Settings File")
    g.comment(f"**Run.** Run Command calls the Python driver `{script}` on the settings file once the geometry and "
              "settings are written (Await enforces the order). Its progress appears in nTop's log and in "
              "`<Case Name>.ntvsp.log`. A non-zero exit code stops the notebook with the driver's error message.")
    command = k.cat(['"', inp("Python"), '" "', inp("Driver Directory"), "/" + script + '" "', settings_path, '"'])
    after = command
    for w in [written] + list(waits):
        after = k.wait(w, after)
    rc = k.run_command(after, inp("Working Directory"), inp("Timeout"), name="Run Driver")
    ok = k.check(rc, what)
    return k.wait(ok, base, name="Results")


def make_block(key, display, description, specs, output_type, g, output):
    return CustomBlock(key, display, MD + description, tuple(specs), output_type, g, output, side_effects=True)


def model_inputs(specs):
    index = {s.name: i for i, s in enumerate(specs)}
    return lambda name: ModelInput(index[name])
