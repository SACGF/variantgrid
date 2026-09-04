import re

from django.template.loader import render_to_string
from django.test import TestCase

from analysis.models.nodes.node_display import NodeChip
from analysis.models.nodes.node_types import get_node_types_hash

SPRITE_TEMPLATE = "uicore/tags/svg_icon_sprite.html"


def get_sprite_symbol_ids() -> set[str]:
    return set(re.findall(r'<symbol id="([^"]+)"', render_to_string(SPRITE_TEMPLATE)))


class NodeDisplayTestCase(TestCase):
    """ A new node class can't ship without the metadata the canvas card needs """

    def test_every_node_class_has_card_metadata(self):
        sprite_symbols = get_sprite_symbol_ids()
        for label, node_class in get_node_types_hash().items():
            with self.subTest(node_class=label):
                icon = node_class.get_node_class_icon()
                self.assertNotEqual(bool(icon.fa), bool(icon.symbol),
                                    f"{label} icon must set exactly one of fa/symbol: {icon}")
                if icon.symbol:
                    self.assertIn(icon.symbol, sprite_symbols)

                self.assertTrue(node_class.get_node_class_label_short())

                node = node_class()
                self.assertEqual(icon, node.get_node_icon())
                self.assertTrue(node.get_node_strip_label())
                chips = node.get_node_chips()
                self.assertIsInstance(chips, list)
                for chip in chips:
                    self.assertIsInstance(chip, NodeChip)

    def test_sample_node_patient_shape_symbols_exist(self):
        """ SampleNode picks its badge by string, so the sprite has to carry all 4 combinations """
        sprite_symbols = get_sprite_symbol_ids()
        for sex in ["male", "female"]:
            for deceased in ["", "-deceased"]:
                self.assertIn(f"node-icon-sample-{sex}{deceased}", sprite_symbols)
