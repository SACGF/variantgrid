def rgb_hex_to_tuples(rgb: str):
    rgb = rgb.replace('#', '')
    return bytes.fromhex(rgb)


def rgb_to_hex(red, green, blue):
    return "#%02x%02x%02x" % (red, green, blue)


def rgb_invert(rgb):
    red, green, blue = rgb_hex_to_tuples(rgb)
    inverted = (255 - red, 255 - green, 255 - blue)
    return rgb_to_hex(*inverted)

