import importlib


def nice_class_name(obj_or_klass) -> str:
    if isinstance(obj_or_klass, type):
        klass = obj_or_klass
    else:
        klass = obj_or_klass.__class__
    return klass.__name__


def full_class_name(klass):
    return klass.__module__ + '.' + klass.__name__


def import_class(full_class_path):
    modulename, classname = full_class_path.rsplit('.', 1)
    # logging.debug("import_class: modulename = %s", modulename)
    mod = importlib.import_module(modulename)
    return getattr(mod, classname)


def get_subclasses(cls):
    """returns all subclasses of argument, cls"""
    if issubclass(cls, type):
        subclasses = cls.__subclasses__(cls)
    else:
        subclasses = cls.__subclasses__()
    for subclass in subclasses:
        subclasses.extend(get_subclasses(subclass))
    return subclasses


def get_all_subclasses(cls):
    """ From https://stackoverflow.com/a/17246726 """
    all_subclasses = set()

    for subclass in cls.__subclasses__():
        all_subclasses.add(subclass)
        all_subclasses.update(get_all_subclasses(subclass))

    return all_subclasses
