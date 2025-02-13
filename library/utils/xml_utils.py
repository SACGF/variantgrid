import inspect
from typing import Callable, Optional, Any, Union

from lxml import etree
from lxml.etree import Element


class PP:
    # Path Predicate
    def __init__(self, tag: str, text: Optional[str] = None, **kwargs):
        self.tag = tag
        self.text = text
        self.atts = kwargs

    def matches(self, elem):
        if self.tag != elem.tag:
            return False
        if self.text and self.text != elem.text:
            return False
        if self.atts:
            for key, value in self.atts.items():
                if elem.get(key) != value:
                    return False
        return True

    @staticmethod
    def convert(test):
        if isinstance(test, str):
            return PP(tag=test)
        else:
            return test

    def __repr__(self):
        output = self.tag
        if self.atts:
            pairs = [f"{key}={value}" for key, value in self.atts.items()]
            output += "[" + ",".join(pairs) + "]"
        if self.text:
            output += f"(\"{self.text}\")"
        return output


class PathPredicates:

    def __init__(self, path: list[Union[PP, str]]):
        self._path = [PP.convert(elem) for elem in path]

    def with_prefix(self, prefix: list[Union[PP, str]]) -> 'PathPredicates':
        return PathPredicates(path=prefix + self._path)

    def queue_test(self, elem: Element) -> Optional['PathPredicates']:
        if self._path[0].matches(elem):
            return PathPredicates(path=self._path[1:])
        else:
            return None

    @property
    def is_satisfied(self):
        return not self._path

    def __repr__(self):
        return "/".join(repr(pp) for pp in self._path)


class ParserMethod:

    def __init__(self, method: Callable, path: PathPredicates, on_start: bool, debug: Optional[str] = None):
        self.method = method
        self.path = path
        self.on_start = on_start
        if not debug:
            debug = f"{path} {method.__name__}"
        self.debug = debug

    def test_and_remainder(self, elem) -> Optional['ParserMethod']:
        if sub_test := self.path.queue_test(elem):
            return ParserMethod(method=self.method, path=sub_test, on_start=self.on_start, debug=self.debug)

    def __repr__(self):
        return self.debug

    @property
    def is_satisfied(self):
        return self.path.is_satisfied


def parser_path(*args, on_start=False):
    """
    Extend XmlParser and annotate methods with parser_path.
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            return method(*args, **kwargs)
        # have to cache the line number of the source method, otherwise we just get the line number of this wrapper
        wrapper.is_parser = True
        wrapper.path = PathPredicates(args)
        wrapper.on_start = on_start
        wrapper.__name__ = method.__name__
        return wrapper
    return decorator


class XmlParser:
    """
    A utility for helping us parse streaming XML data (as in open tag/close tag events) into classes.
    Can be run on things like clinvar's multi GB XML output, or just a few kb of XML.

    Annotate methods with a parser_path (a series of cascading path predicates), if the predicate is matched
    at the stage of processing the XML, then the method is called.

    There's an expectation that we're only populating one kind of object. When the XML indicates that there should
    be a new instance, call `set_yieldable` with the new object. Then other methods should parse it.
    """

    @classmethod
    def get_parser_methods(cls, prefix: list) -> list[ParserMethod]:
        if not hasattr(cls, 'parser_methods'):
            parser_methods = [func for _, func in inspect.getmembers(cls, lambda x: getattr(x, 'is_parser', False))]
            cls.parser_methods = [ParserMethod(method=pm, path=pm.path.with_prefix(prefix), on_start=pm.on_start) for pm in parser_methods]

            if not cls.parser_methods:
                raise ValueError(f"XmlParser class {cls} has no @parser_path")

        return cls.parser_methods

    # determine if we only save tagName
    def __init__(self, prefix: Optional[list[Union[PP, str]]] = None):
        self._prefix = prefix or []
        self._stack: list = []
        self._candidates: list[ParserMethod] = []
        self._execute: list[ParserMethod] = []
        self._yieldable: Optional[Any] = None

    def set_yieldable(self, obj):
        self._clear()
        self.post_record_parse(obj)
        self._yieldable = obj

    def post_record_parse(self, obj):
        pass

    def yield_me(self):
        if temp := self._yieldable:
            self._yieldable = None
            return temp

    def finish(self):
        pass

    def parse(self, source):
        self._stack = []
        self._candidates = [self.__class__.get_parser_methods(self._prefix)]
        self._execute = []

        context = etree.iterparse(source, events=('start', 'end'), huge_tree=True, recover=True, encoding="utf-8")
        for event, elem in context:
            if event == "start":
                self._push(elem)
            elif event == "end":
                self._pop()
            if result := self.yield_me():
                yield result

        self.finish()
        if result := self.yield_me():
            yield result

    @property
    def _peek_candidates(self) -> list[ParserMethod]:
        return self._candidates[-1]

    def _push(self, elem: Element):
        self._stack.append(elem)
        execute_later = []
        remaining_candidates = []
        evaluated = [cand for cand in [cand.test_and_remainder(elem) for cand in self._peek_candidates]
                     if cand]
        for candidate in evaluated:
            if candidate.is_satisfied:
                if candidate.on_start:
                    # execute now
                    # print(f"Executing ({candidate}) {candidate.method.__name__}")
                    candidate.method(self, elem)
                else:
                    execute_later.append(candidate)
            else:
                remaining_candidates.append(candidate)

        self._candidates.append(remaining_candidates)
        self._execute.append(execute_later)

    def _pop(self) -> list[ParserMethod]:
        if last_execute := self._execute.pop():
            for method in last_execute:
                candidate = self._peek
                # print(f"Executing ({method}) {method.method.__name__}")
                method.method(self, candidate)
        self._stack.pop()
        self._candidates.pop()

    @property
    def tag_names(self) -> list[str]:
        return [item.tag for item in self._stack]

    @property
    def _peek(self) -> Element:
        return self._stack[-1]

    def _clear(self):
        # really not sure when the right time to call this is
        if self._stack:
            self._stack[0].clear()
