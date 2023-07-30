import inspect
from typing import List, Callable, Optional

from lxml import etree


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

    def __init__(self, path: List):
        self._path = [PP.convert(elem) for elem in path]

    def queue_test(self, elem) -> Optional['PathPredicates']:
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

    def __init__(self, method, path: PathPredicates, on_start):
        self.method = method
        self.path = path
        self.on_start = on_start

    def test_and_remainder(self, elem) -> Optional['ParserMethod']:
        if sub_test := self.path.queue_test(elem):
            return ParserMethod(method=self.method, path=sub_test, on_start=self.on_start)

    def __repr__(self):
        return f"PM:{self.path}"

    @property
    def is_satisfied(self):
        return self.path.is_satisfied


def parser_path(*args, on_start=True):
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
        return wrapper
    return decorator


class XmlParser:

    @staticmethod
    def get_parser_methods(klass):
        if not hasattr(klass, 'parser_methods'):
            parser_methods = [func for _, func in inspect.getmembers(klass, lambda x: getattr(x, 'is_parser', False))]
            klass.parser_methods = [ParserMethod(method=pm, path=pm.path, on_start=pm.on_start) for pm in parser_methods]

            if not klass.parser_methods:
                raise ValueError(f"XmlParser class {klass} has no @parser_path")

        return klass.parser_methods

    # determine if we only save tagName
    def __init__(self):
        self._stack: List = list()
        self._candidates: List[List[Callable]] = list()
        self._execute: List[List[Callable]] = list()
        self._yieldable: Optional = None

    def set_yieldable(self, obj):
        self._clear()
        self._yieldable = obj

    def yield_me(self):
        if temp := self._yieldable:
            self._yieldable = None
            return temp

    def finish(self):
        pass

    def parse(self, source):
        self._stack = list()
        self._candidates = [XmlParser.get_parser_methods(self.__class__)]
        self._execute = list()

        context = etree.iterparse(source, events=('start', 'end'))
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
    def _peek_candidates(self):
        return self._candidates[-1]

    def _push(self, elem):
        self._stack.append(elem)
        execute_later = list()
        remaining_candidates = list()
        evaluated = [cand for cand in [cand.test_and_remainder(elem) for cand in self._peek_candidates]
                     if cand]
        for candidate in evaluated:
            if candidate.is_satisfied:
                if candidate.on_start:
                    # execute now
                    candidate.method(self, elem)
                else:
                    execute_later.append(candidate)
            else:
                remaining_candidates.append(candidate)

        self._candidates.append(remaining_candidates)
        self._execute.append(execute_later)

    def _pop(self):
        if last_execute := self._execute.pop():
            for method in last_execute:
                method.method(self, self._peek)
        self._stack.pop()
        self._candidates.pop()

    @property
    def tag_names(self):
        return [item.tag for item in self._stack]

    @property
    def _peek(self):
        return self._stack[-1]

    def _clear(self):
        # really not sure when the right time to call this is
        self._stack[0].clear()