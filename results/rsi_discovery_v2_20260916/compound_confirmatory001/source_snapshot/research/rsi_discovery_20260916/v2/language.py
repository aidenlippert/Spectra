"""Pure, resource-bounded scientific programs; arithmetic is not a parameter grid.

Run untrusted programs in the existing deadline-limited child, never the controller.
This allowlist is a language boundary, not an OS security sandbox.
"""
import ast
import builtins

from research.rsi_discovery_20260916.dialect import NODES as OLD_NODES

NODES = OLD_NODES | {"While", "SetComp", "Del", "Delete", "Is", "IsNot"}
ATTRS = {"append", "extend", "get", "items", "keys", "values", "pop", "copy",
         "setdefault", "update", "clear", "add", "discard", "sort", "index", "count",
         "numerator", "denominator"}
BUILTINS = {name: getattr(builtins, name) for name in
            ("len", "range", "enumerate", "zip", "sum", "min", "max", "abs", "round",
             "int", "float", "str", "bool", "list", "tuple", "dict", "set", "sorted",
             "reversed", "any", "all")}


def compile_program(source, capabilities):
    if not isinstance(source, str) or len(source) > 40000:
        raise ValueError("Program exceeds source envelope")
    tree = ast.parse(source)
    if not tree.body or any(not isinstance(n, ast.FunctionDef) for n in tree.body):
        raise ValueError("Only function definitions at module level")
    if len({n.name for n in tree.body}) != len(tree.body):
        raise ValueError("Duplicate function names")
    entry = [n for n in tree.body if n.name == "propose"]
    if len(entry) != 1 or [a.arg for a in entry[0].args.args] != ["payload"]:
        raise ValueError("A propose(payload) entry point is required")
    for node in ast.walk(tree):
        if type(node).__name__ not in NODES:
            raise ValueError("Unsupported syntax: " + type(node).__name__)
        if isinstance(node, ast.FunctionDef):
            if node.name.startswith("_") or node.decorator_list or node.returns:
                raise ValueError("Decorators, annotations and private names forbidden")
            if node.args.defaults or node.args.kw_defaults or node.args.vararg or node.args.kwarg:
                raise ValueError("Only explicit function arguments are allowed")
        if isinstance(node, ast.arg) and (node.arg.startswith("_") or node.annotation):
            raise ValueError("Private/annotated argument")
        if isinstance(node, ast.Name) and node.id.startswith("_"):
            raise ValueError("Private name")
        if isinstance(node, ast.Attribute) and node.attr not in ATTRS:
            raise ValueError("Unsupported attribute: " + node.attr)
        if isinstance(node, ast.Call) and not isinstance(node.func, (ast.Name, ast.Attribute)):
            raise ValueError("Indirect calls are forbidden")
        if isinstance(node, ast.Constant) and not isinstance(node.value, (str, int, float, bool, type(None))):
            raise ValueError("Unsupported constant")
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            if not isinstance(node.right, ast.Constant) or type(node.right.value) is not int or not 0 <= node.right.value <= 16:
                raise ValueError("Power must have a literal exponent in 0..16")
    scope = {"__builtins__": BUILTINS, **capabilities}
    exec(compile(tree, "<scientific-program-v2>", "exec"), scope)
    return scope["propose"]
