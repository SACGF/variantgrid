"""
pylint plugin: drop the `if TYPE_CHECKING:` method stubs django-model-utils declares on `InheritanceManagerMixin`.

Their `...` bodies are all pylint sees, so every `objects.filter()` / `.all()` / `.none()` on a model with an
`InheritanceManager` was reported as assignment-from-no-return. Without the stubs those calls infer from Django's
`Manager`, as on any other model. Loaded by `scripts/linting/run_pylint.sh`.
"""
import astroid
from astroid import nodes


def _is_type_checking_block(node: nodes.NodeNG) -> bool:
    return isinstance(node, nodes.If) and node.test.as_string() == "TYPE_CHECKING"


def _strip_type_checking_stubs(cls: nodes.ClassDef) -> nodes.ClassDef:
    for block in [n for n in cls.body if _is_type_checking_block(n)]:
        cls.body.remove(block)
        for stub in block.body:
            if isinstance(stub, nodes.FunctionDef) and stub in cls.locals.get(stub.name, []):
                cls.locals[stub.name].remove(stub)
                if not cls.locals[stub.name]:
                    del cls.locals[stub.name]
    return cls


def _is_inheritance_manager_mixin(cls: nodes.ClassDef) -> bool:
    return cls.name == "InheritanceManagerMixin" and cls.root().name == "model_utils.managers"


def register(_linter):
    """ pylint plugin entry point - the work is the astroid transform registered at import """


astroid.MANAGER.register_transform(nodes.ClassDef, _strip_type_checking_stubs, _is_inheritance_manager_mixin)
