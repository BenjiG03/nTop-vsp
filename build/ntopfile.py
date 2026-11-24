"""Write nTop custom blocks and host notebooks as .ntop containers.

A notebook is a graph of block instances in JSON plus a table of constants. This module
keeps that graph as plain data: ``Graph`` accumulates instances and the functions at the
bottom turn it into container sections. Container reading and writing comes from
``implicit_io``, which also documents the format.
"""
from __future__ import annotations

import uuid
from dataclasses import dataclass, field, replace

from implicit_io.file_processing import nbgen
from implicit_io.file_processing import ntop_container as nc

TEMPLATE = nbgen.TEMPLATE


@dataclass(frozen=True)
class Ref:
    """A connection to another instance, optionally through one of its properties."""

    instance: int
    prop: tuple = ()


@dataclass(frozen=True)
class ModelInput:
    """A connection to the custom block's own input number ``index``."""

    index: int
    prop: tuple = ()


NONE = Ref(0)  # an optional block input left unconnected


def _connection(item) -> dict:
    if isinstance(item, ModelInput):
        target, model, chain = -1, item.index, list(item.prop)
    else:
        ref = item if isinstance(item, Ref) else Ref(item)
        target, model, chain = ref.instance, -1, list(ref.prop)
    return {
        "instanceId": target,
        "meta": {"expression": "", "is_default": False, "name": ""},
        "modelInputIdx": model,
        "propchain": chain,
    }


def _instance(iid, func, typ, inputs, name, units=None, auto=False, desc="") -> dict:
    return {
        "autoc": auto, "bldst": 0, "func": func, "id": iid,
        "inputs": [_connection(i) for i in inputs], "meta": {"desc": desc, "import": 0},
        "name": name, "optional": False, "paust": 0, "qual": 0, "type": typ,
        "units": units or {}, "unitsKnown": True,
    }


def _root(top) -> dict:
    root = _instance(100, "core.group<any>", "group", list(top), "root")
    root["qual"] = 2
    return root


@dataclass
class Graph:
    """Block instances and constants of one function, built in execution order."""

    code: list = field(default_factory=list)
    leaves: list = field(default_factory=list)
    top: list = field(default_factory=list)
    next_id: int = 101

    def _new(self) -> int:
        self.next_id += 1
        return self.next_id

    def leaf(self, typ, value, units=None, name=None) -> int:
        iid = self._new()
        self.code.append(_instance(iid, "", typ, [], name or f"{typ}_{iid}", units=units, auto=True))
        self.leaves.append({"id": str(iid), "type": typ, "value": value})
        return iid

    def real(self, v, units=None) -> int:
        return self.leaf("real", {"isFinite": True, "units": units or {}, "val": float(v)}, units=units)

    def integer(self, v) -> int:
        return self.leaf("integer", {"val": int(v)})

    def boolean(self, v) -> int:
        return self.leaf("bool", {"val": bool(v)})

    def text(self, s) -> int:
        return self.leaf("text", {"string": s})

    def enum(self, typ, v) -> int:
        return self.leaf(typ, {"enum": int(v)})

    def choice(self, labels, selected=0) -> int:
        return self.leaf("choice", {"choices": list(labels), "selected": int(selected)})

    def unit_length(self, unit="mm") -> int:
        return self.leaf("unit_length_enum", {"id": unit})

    def file_path(self, p) -> int:
        return self.leaf("file_path", {"val": p.replace("\\", "/")})

    def point(self, x, y, z) -> int:
        """Point with coordinates in millimetres."""
        return self.leaf("point", [{"isFinite": True, "val": c * nbgen.MM} for c in (x, y, z)], units={"length": 1})

    def vector(self, x, y, z) -> int:
        return self.leaf("vector", {"units": {}, "value": [{"isFinite": True, "val": float(c)} for c in (x, y, z)]})

    def block(self, func, typ, inputs, name=None, top=True, units=None) -> int:
        iid = self._new()
        # nTop renames duplicates on load, so unnamed helpers get unique names here.
        self.code.append(_instance(iid, func, typ, inputs, name or f"{func.split('<')[0]}_{iid}", units=units))
        if top:
            self.top.append(iid)
        return iid

    def comment(self, text) -> int:
        """A floating comment in the notebook body, shown in order with the top-level blocks (markdown)."""
        iid = self._new()
        inst = _instance(iid, "", "text", [], f"Floating Comment_{iid}", desc="$markdown-enabled$" + text)
        inst["isFloatingComment"] = True
        self.code.append(inst)
        self.top.append(iid)
        return iid

    def concat(self, parts, separator="", top=False, name=None) -> int:
        """Join text pieces (instance references) into one text."""
        pieces = self.block("core.list<text>", "list<text>", parts, top=False)
        return self.block("concatenate<list<text>,text>", "text", [pieces, self.text(separator)], name=name, top=top)


@dataclass(frozen=True)
class InputSpec:
    """One input of a custom block. ``default`` builds the input's default in its own graph."""

    name: str
    typ: str
    description: str
    units: dict | None = None
    default: object = None  # callable(Graph) -> instance id


_DEF_FLAGS = {
    "newSystemOutput": {"type": ""},
    "outputBaseUnits": {},
    "releasestage": 2,
    "runaspects.defaultname": "",
    "runaspects.disablelistprocessing": False,
    "runaspects.embeddable": False,
    "runaspects.implicitparam": 0,
    "runaspects.iterative": False,
    "runaspects.overloadpriority": 0,
    "runaspects.progstages": -1,
    "runaspects.rons": 1,
    "runaspects.searchable": True,
    "runaspects.startvisible": True,
    "unitsVariadicMultiplier": [0, 1],
}


def _function_name(key: str) -> str:
    """A stable function id derived from the block's name, so rebuilds stay compatible."""
    return "user_func_" + str(uuid.uuid5(uuid.NAMESPACE_URL, "ntop_vsp/" + key)).replace("-", "_")


def _definition(display, description, function, inputs, output_type, side_effects, output_units=None) -> dict:
    return {
        **_DEF_FLAGS,
        "description": description,
        "displayname": display,
        "function": function,
        "inputs": inputs,
        "output": {"default_type": output_type, "type": output_type},
        "outputBaseUnits": output_units or {},
        "runaspects.security": side_effects,
        "runaspects.sideeffects": side_effects,
    }


def _ntopfn(name, program, leaves) -> nc.Section:
    fn = nc.Section("json", "fn")
    fn.set_json(program)
    index = nc.Section("json", "index")
    index.set_json(leaves)
    return nc.Section("ntopfn", name, children=[fn, nc.Section("obj_container", "leaves", children=[index])])


def _input_section(name, spec: InputSpec, block_key, index) -> nc.Section:
    graph = Graph()
    value = spec.default(graph) if spec.default else None
    connection = Ref(value) if value is not None else Ref(-1)
    var_id = 101  # nTop's own input sections always use this id for the input variable
    var = _instance(var_id, f"core.var<{spec.typ}>", spec.typ, [connection], spec.name, desc=spec.description,
                    units=spec.units)
    if value is None:
        var["inputs"][0]["type"] = spec.typ
    code = [_root([var_id]), var] + graph.code
    function = _function_name(f"{block_key}/input/{index}")
    definition = _definition(spec.name, spec.description, function, [], spec.typ, False, spec.units)
    return _ntopfn(name, {"code": code, "def": definition, "dependencies": [], "fninf": None, "output": var_id},
                   graph.leaves)


# nTop's own output types of blocks the builders annotate loosely: a variable checks its block's type, a plain block
# does not (Export Text's output is a Text File, "text_file", not a path).
OUTPUT_TYPES = {"export_text": "text_file"}
LIST_PROPERTIES = {"points", "values"}  # properties that are lists (of a point map)


def _params(func):
    """Parameter types of a block function from its signature, ``name<a,b<c,d>,e>[version]``."""
    if "<" not in func:
        return []
    inner = func[func.index("<") + 1:func.rindex(">")]
    out, depth, cur = [], 0, ""
    for ch in inner:
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
            continue
        depth += ch == "<"
        depth -= ch == ">"
        cur += ch
    return out + ([cur] if cur else [])


def _takes_list(param):
    return param.startswith("list") or param in ("any", "list_interface")


def list_outputs(code, input_types):
    """Ids of the instances whose value is a list in a custom block's graph, by nTop's list processing: a block
    whose input takes one value but is given a list runs once per item. ``input_types``: the block's input types."""
    by_id = {inst["id"]: inst for inst in code}
    known = {}

    def is_list(conn):
        if conn["modelInputIdx"] >= 0:
            base = input_types[conn["modelInputIdx"]].startswith("list<")
            chain = conn["propchain"]
            return base and not (chain and chain[-1] == "size") or bool(set(chain) & LIST_PROPERTIES)
        target = conn["instanceId"]
        if target not in by_id:
            return False
        if set(conn["propchain"]) & LIST_PROPERTIES:
            return True
        if conn["propchain"] and conn["propchain"][-1] == "size":
            return False
        return value_is_list(target)

    def value_is_list(iid):
        if iid in known:
            return known[iid]
        known[iid] = False  # no cycles in a graph; guards recursion
        inst = by_id[iid]
        if not inst["func"]:
            result = inst["type"].startswith("list<")
        elif inst["func"].startswith(("core.var<", "await<")):
            result = is_list(inst["inputs"][-1])  # a variable or an Await passes its (last) input through
        elif inst["type"].startswith("list<") or inst["func"].startswith("core.group<"):
            result = inst["type"].startswith("list<")
        else:
            params = _params(inst["func"])
            result = any(is_list(c) and not _takes_list(p) for c, p in zip(inst["inputs"], params))
        known[iid] = result
        return result

    return {iid for iid in by_id if value_is_list(iid)}


def as_variables(top, code, output, input_types=()):
    """Every top-level block of a custom block's graph made a variable, as nTop shows them: each is wrapped in a
    ``core.var<type>`` that takes the block's name and its place at the top level, the block moves inside it, and
    every connection to the block (properties included) and the block's output go through the variable. A block run
    once per item of a list (``list_outputs``) gets a ``list<type>`` variable. Variables and floating comments stay
    as they are. Returns (top, code, output) without changing the graph."""
    code = [dict(inst, type=OUTPUT_TYPES.get(inst["func"].split("<")[0], inst["type"])) if inst["func"] else inst
            for inst in code]
    by_id = {inst["id"]: inst for inst in code}
    targets = [t for t in top if not by_id[t].get("isFloatingComment") and not by_id[t]["func"].startswith("core.var<")]
    if not targets:
        return list(top), list(code), output
    lists = list_outputs(code, list(input_types))
    next_id = max(by_id) + 1
    var_of, variables = {}, []
    for t in targets:
        block = by_id[t]
        typ = block["type"]
        if t in lists and not typ.startswith("list<"):
            typ = f"list<{typ}>"
        var_of[t] = next_id
        var = _instance(next_id, f"core.var<{typ}>", typ, [Ref(t)], block["name"], units=block.get("units") or None)
        variables.append(var)
        next_id += 1
    new_code = []
    for inst in code:
        inst = dict(inst)
        if inst["id"] in var_of:
            inst["name"] = f"{inst['name']} (block)"
        inputs = []
        for conn in inst["inputs"]:
            if conn["instanceId"] in var_of:
                conn = dict(conn, instanceId=var_of[conn["instanceId"]])
            inputs.append(conn)
        inst["inputs"] = inputs
        new_code.append(inst)
        if inst["id"] in var_of:  # the variable right after its block
            new_code.append(variables[targets.index(inst["id"])])
    return [var_of.get(t, t) for t in top], new_code, var_of.get(output, output)


@dataclass(frozen=True)
class CustomBlock:
    """A complete custom block, ready to save alone or embed in a notebook."""

    key: str
    display: str
    description: str
    inputs: tuple
    output_type: str
    graph: Graph
    output: int
    side_effects: bool = True
    function_id: str = ""  # keeps the id of a block saved from nTop; default derives one from the key
    dependencies: tuple = ()  # indices of the custom blocks this one calls, among the embedded sf sections
    requires: tuple = field(default=(), compare=False)  # the custom blocks this one calls (set dependencies on save)

    @property
    def signature(self) -> str:
        function = self.function_id or _function_name(self.key)
        return f"{function}<{','.join(s.typ for s in self.inputs)}>"

    def program(self) -> dict:
        inputs = [
            {
                "cardinality": "one", "description": s.description, "name": s.name, "parameterization": 1,
                "type": {"default_type": s.typ, "type": s.typ}, "unitsMult": [0, 1],
                **({"unitsReq": s.units} if s.units else {}),
            }
            for s in self.inputs
        ]
        definition = _definition(self.display, self.description, self.signature, inputs, self.output_type, self.side_effects)
        top, code, output = as_variables(self.graph.top, self.graph.code, self.output, [s.typ for s in self.inputs])
        return {
            "code": [_root(top)] + code,
            "def": definition,
            "dependencies": list(self.dependencies),
            "fninf": None,
            "output": output,
        }

    def sections(self, prefix="") -> list:
        """Sections of the function and its inputs, named for a file or embedded in a host."""
        main = _ntopfn(prefix or "main", self.program(), self.graph.leaves)
        stem = prefix or ""
        return [main] + [
            _input_section(f"{stem}inp{i}", s, self.key, i) for i, s in enumerate(self.inputs)
        ]


def _template_sections(top_count, path=TEMPLATE):
    """Container sections of a template notebook with its cache and view state cleared."""
    secs = nc.read(path)
    for s in secs:
        if s.name == "cache":
            s.children = [nc.Section("json", "index", b"[]")]
        elif s.name in ("open", "view"):
            s.set_json([])
        elif s.name == "sections":
            s.set_json({"decorations": [{"collapse": False, "description": "", "name": "Section 1"}], "poles": [0, top_count]})
    return secs


def embed_order(blocks) -> list:
    """``blocks`` and every block they require, each once, required blocks before the blocks that call them."""
    out = []

    def add(b):
        for r in b.requires:
            add(r)
        if all(x.key != b.key for x in out):
            out.append(b)

    for b in blocks:
        add(b)
    return out


def with_dependencies(blocks, first=0) -> list:
    """The blocks to embed as sf sections numbered from ``first``: ``blocks`` and what they require, each with
    ``dependencies`` set to the sf numbers of the blocks it calls."""
    order = embed_order(blocks)
    number = {b.key: first + i for i, b in enumerate(order)}
    return [replace(b, dependencies=tuple(number[r.key] for r in b.requires)) for b in order]


def save_block(block: CustomBlock, path, template=TEMPLATE, embedded=()) -> str:
    """Write the block as its own .ntop file, importable with File > Import.

    The custom blocks it calls (``block.requires``, or ``embedded`` in the order its ``dependencies`` index) are
    stored as sf sections, as nTop stores a custom block that uses other custom blocks.
    """
    if block.requires:
        embedded = with_dependencies(block.requires)
        number = {b.key: i for i, b in enumerate(embedded)}
        block = replace(block, dependencies=tuple(number[r.key] for r in block.requires))
    secs = _template_sections(len(block.graph.top), template)
    keep = [s for s in secs if s.name not in ("main",) and not s.name.startswith(("inp", "sf"))]
    head = [s for s in keep if s.name == "turhe"]
    tail = [s for s in keep if s.name != "turhe"]
    inner = []
    for i, child in enumerate(embedded):
        parts = child.sections(prefix=f"sf{i}")
        meta = nc.Section("json", f"sf{i}_metadata")
        meta.set_json({"local": False, "time": 0})
        inner += parts[:1] + [meta] + parts[1:]
    nc.write(path, head + inner + block.sections() + tail)
    return path


def save_host(graph: Graph, blocks, path, output=None, output_type="", name="host", template=TEMPLATE) -> str:
    """Write a notebook that embeds ``blocks`` and runs ``graph``."""
    secs = _template_sections(len(graph.top), template)
    keep = [s for s in secs if s.name != "main" and not s.name.startswith(("inp", "sf"))]
    head = [s for s in keep if s.name == "turhe"]
    tail = [s for s in keep if s.name != "turhe"]
    blocks = with_dependencies(blocks)
    embedded = []
    for i, block in enumerate(blocks):
        parts = block.sections(prefix=f"sf{i}")
        meta = nc.Section("json", f"sf{i}_metadata")
        meta.set_json({"local": False, "time": 0})
        embedded += parts[:1] + [meta] + parts[1:]
    definition = _definition(name, "", _function_name(f"host/{name}"), [], output_type, True)
    program = {
        "code": [_root(graph.top)] + graph.code,
        "def": definition,
        "dependencies": list(range(len(blocks))),
        "output": -1 if output is None else output,
    }
    nc.write(path, head + embedded + [_ntopfn("main", program, graph.leaves)] + tail)
    return path
