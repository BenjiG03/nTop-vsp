"""Derive notebooks from the X-59 hackathon starter: keep its parametric geometry and append blocks to it."""
from __future__ import annotations

import json
import os
import re

import ntopfile
from ntopfile import CustomBlock, Graph, InputSpec, ModelInput, nc

HERE = os.path.dirname(os.path.abspath(__file__))
STARTER = os.path.join(os.path.dirname(os.path.dirname(HERE)), "x59-hackathon-starter", "x59-hackathon-starter.ntop")
CREDIT = ("Geometry: X-59 Hackathon Starter, MIT License, Copyright (c) 2026 matthewmuellernTop "
          "(x59-hackathon-starter.ntop).")


class Derived:
    def __init__(self, path=STARTER):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"X-59 starter notebook not found at {path}")
        self.secs = nc.read(path)
        main = next(s for s in self.secs if s.name == "main")
        self.fn = next(c for c in main.children if c.name == "fn").json()
        leaves = next(c for c in main.children if c.name == "leaves")
        self.leaves = json.loads(next(c for c in leaves.children if c.name == "index").payload)
        self.code = self.fn["code"]
        self.root = self.code[0]
        self.by_id = {c["id"]: c for c in self.code}
        self.inputs = []  # (InputSpec) notebook inputs, in order

    def id_of(self, name) -> int:
        hits = [c["id"] for c in self.code if c.get("name") == name and not c.get("isFloatingComment")]
        if len(hits) != 1:
            raise KeyError(f"{len(hits)} blocks named {name!r} in the starter")
        return hits[0]

    def keep(self, names):
        """Drop every block that the named variables do not depend on."""
        want, stack = set(), [self.id_of(n) for n in names]
        while stack:
            i = stack.pop()
            if i in want or i not in self.by_id:
                continue
            want.add(i)
            stack += [c["instanceId"] for c in self.by_id[i]["inputs"] if c["instanceId"] > 0]
        self.root["inputs"] = [c for c in self.root["inputs"] if c["instanceId"] in want]
        self.code = [self.root] + [c for c in self.code[1:] if c["id"] in want]
        for c in self.code:
            c["bldst"] = 0
        self.by_id = {c["id"]: c for c in self.code}
        self.leaves = [e for e in self.leaves if int(e["id"]) in want and "binoffset" not in e]

    def expose(self, name, spec: InputSpec):
        """Make the variable ``name`` a notebook input (ntopcl sets it from its -j file)."""
        var = self.by_id[self.id_of(name)]
        var["inputs"][0] = ntopfile._connection(ModelInput(len(self.inputs)))
        self.inputs.append(spec)

    def graph(self) -> Graph:
        """A graph whose new blocks follow the starter's; existing blocks are referenced by id."""
        return Graph(next_id=max(self.by_id) + 1)

    def save(self, path, graph: Graph, blocks=(), output=None, output_type="", display="X-59 example", description=""):
        used = [int(m.group(1)) for s in self.secs for m in [re.fullmatch(r"sf(\d+)", s.name)] if m]
        first = max(used, default=-1) + 1
        blocks = ntopfile.with_dependencies(blocks, first)
        embedded = []
        for k, block in enumerate(blocks):
            parts = block.sections(prefix=f"sf{first + k}")
            meta = nc.Section("json", f"sf{first + k}_metadata")
            meta.set_json({"local": False, "time": 0})
            embedded += parts[:1] + [meta] + parts[1:]
        top = self.root["inputs"] + [ntopfile._connection(i) for i in graph.top]
        root = dict(self.root, inputs=top)
        host = CustomBlock("x59_example/" + display, display, description, tuple(self.inputs), output_type, Graph(), -1)
        made = host.program()["def"]
        definition = dict(self.fn["def"], description=description, displayname=display, inputs=made["inputs"],
                          output=made["output"])
        program = {
            "code": [root] + self.code[1:] + graph.code,
            "def": definition,
            "dependencies": list(self.fn["dependencies"]) + [first + k for k in range(len(blocks))],
            "output": -1 if output is None else output,
        }
        main = ntopfile._ntopfn("main", program, self.leaves + graph.leaves)
        inputs = [ntopfile._input_section(f"inp{i}", s, host.key, i) for i, s in enumerate(self.inputs)]
        out = []
        for s in self.secs:
            if s.name == "main":
                out += embedded + [main] + inputs
                continue
            if re.fullmatch(r"inp\d+", s.name):
                continue
            if s.name == "cache":
                s.children = [nc.Section("json", "index", b"[]")]
            elif s.name in ("open", "view"):
                s.set_json([])
            elif s.name == "sections":
                s.set_json({"decorations": [{"collapse": False, "description": "", "name": "Section 1"}],
                            "poles": [0, len(top)]})
            out.append(s)
        nc.write(path, out)
        return path
