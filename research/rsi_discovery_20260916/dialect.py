"""A small executable Python dialect with explicit numerical capabilities.

No imports, attributes beyond an allowlist, I/O, introspection or dynamic code.
This is a restricted language, not a sandbox for arbitrary Python packages.
"""
import ast

NODES = {"Module", "FunctionDef", "arguments", "arg", "Return", "Assign", "AugAssign", "Expr",
         "For", "If", "IfExp", "Break", "Continue", "Pass", "Name", "Load", "Store", "Constant",
         "List", "Tuple", "Dict", "Set", "ListComp", "DictComp", "GeneratorExp", "comprehension",
         "Subscript", "Slice", "Attribute", "Call", "keyword", "BinOp", "UnaryOp", "BoolOp", "Compare",
         "Add", "Sub", "Mult", "Div", "FloorDiv", "Mod", "Pow", "USub", "UAdd", "Not", "And", "Or",
         "Eq", "NotEq", "Lt", "LtE", "Gt", "GtE", "In", "NotIn"}
ATTRIBUTES = {"append", "extend", "get", "items", "keys", "values", "numerator", "denominator"}
BUILTINS = {k: v for k, v in vars(__import__("builtins")).items()
            if k in {"len", "range", "enumerate", "zip", "sum", "min", "max", "abs", "round", "int",
                     "float", "str", "list", "tuple", "sorted", "reversed", "any", "all"}}


def compile_candidate(source, capabilities):
    if not isinstance(source, str) or len(source) > 24000:
        raise ValueError("Candidate source exceeds the declared dialect envelope")
    tree = ast.parse(source)
    if len(tree.body) != 1 or not isinstance(tree.body[0], ast.FunctionDef) or tree.body[0].name != "propose":
        raise ValueError("One propose(payload) function required")
    fn = tree.body[0]
    if fn.decorator_list or fn.returns or fn.args.defaults or fn.args.kw_defaults or fn.args.vararg or fn.args.kwarg:
        raise ValueError("Decorators, annotations and default arguments are forbidden")
    if len(fn.args.args) != 1 or fn.args.args[0].arg != "payload" or fn.args.kwonlyargs or fn.args.posonlyargs:
        raise ValueError("Signature must be propose(payload)")
    for node in ast.walk(tree):
        if type(node).__name__ not in NODES:
            raise ValueError("Unsupported syntax: " + type(node).__name__)
        if isinstance(node, ast.FunctionDef) and node is not fn:
            raise ValueError("Nested functions are not supported")
        if isinstance(node, ast.Name) and node.id.startswith("_"):
            raise ValueError("Private names are forbidden")
        if isinstance(node, ast.arg) and node.annotation:
            raise ValueError("Annotations are forbidden")
        if isinstance(node, ast.Attribute) and node.attr not in ATTRIBUTES:
            raise ValueError("Unsupported attribute: " + node.attr)
        if isinstance(node, ast.Call) and not isinstance(node.func, (ast.Name, ast.Attribute)):
            raise ValueError("Indirect calls are forbidden")
        if isinstance(node, ast.Constant) and not isinstance(node.value, (str, int, float, bool, type(None))):
            raise ValueError("Unsupported constant")
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            if not isinstance(node.right, ast.Constant) or type(node.right.value) is not int or not 0 <= node.right.value <= 16:
                raise ValueError("Only small constant powers are supported")
    scope = {"__builtins__": BUILTINS, **capabilities}
    exec(compile(tree, "<retained-candidate>", "exec"), scope)
    return scope["propose"]
