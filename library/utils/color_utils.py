def rgb_hex_to_tuples(rgb: str) -> bytes:
    rgb = rgb.replace('#', '')
    return bytes.fromhex(rgb)


def rgb_to_hex(red: int, green: int, blue: int) -> str:
    return f"#{red:02x}{green:02x}{blue:02x}"


def rgb_invert(rgb: str) -> str:
    red, green, blue = rgb_hex_to_tuples(rgb)
    inverted = (255 - red, 255 - green, 255 - blue)
    return rgb_to_hex(*inverted)


# Relative luminance where "black text on this" and "white text on this" have equal WCAG contrast
_CONTRAST_MIDPOINT = 0.1791


def rgb_relative_luminance(rgb: str) -> float:
    """ WCAG 2.x relative luminance (0=black, 1=white) """
    linear = []
    for channel in rgb_hex_to_tuples(rgb):
        c = channel / 255
        linear.append(c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4)
    red, green, blue = linear
    return 0.2126 * red + 0.7152 * green + 0.0722 * blue


def rgb_contrasting_text(rgb: str) -> str:
    """ Black or white text - whichever is more readable on an rgb background.
        Inverting the background instead leaves mid-tone colors invisible (#808080 inverts to #7f7f7f) """
    return "#000000" if rgb_relative_luminance(rgb) > _CONTRAST_MIDPOINT else "#ffffff"
