from matplotlib.patches import Rectangle
from matplotlib.ticker import Formatter


def format_strip_zeros(val, num_decimals=2):
    if num_decimals < 1:
        return str(val)
    format_str = f".{num_decimals}f"
    return format(val, format_str).rstrip('0').rstrip('.')


def format_m_and_k(value, num_decimals=2):
    """ Format decimal to (M)illions and (K)ilos """
    if value >= 1000000:
        label = format_strip_zeros(value / 1000000, num_decimals) + "M"
    elif value >= 1000:
        label = format_strip_zeros(value / 1000, num_decimals) + "K"
    else:
        label = format_strip_zeros(value, num_decimals)
    return label


class MandKIntFormatter(Formatter):
    """ Formats 10000 -> 10K, 2000000 -> 2M """

    def __call__(self, value, pos=None):
        label = str(value)
        if isinstance(value, int) or float(value).is_integer():
            if value % 1000000 == 0:
                label = "%dM" % (value // 1000000)
            elif value % 1000 == 0:
                label = "%dK" % (value // 1000)

        return label


class ForceMandKIntFormatter(Formatter):
    """ Formats 10000 -> 10K, 2000000 -> 2M
        Converts non-round numbers, using 2 decimal places
    """

    def __call__(self, value, pos=None):
        return format_m_and_k(value)


def legend_patches_and_labels(legend):
    """ legend: [('label', 'color'), ('label2', 'color2')] """
    patches = []
    labels = []
    for (name, color) in legend:
        labels.append(name)
        shape = Rectangle((0, 0), 1, 1, fc=color)
        patches.append(shape)

    return patches, labels
