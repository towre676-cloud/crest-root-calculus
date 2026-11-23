"""Check vector manufacturing geometry at its final file-format boundary."""

import json
import math
from pathlib import Path
import unittest
from xml.etree import ElementTree

import ezdxf
import fitz

HERE = Path(__file__).resolve().parent
SPEC = json.loads((HERE / "design.json").read_text())
MM_PER_PT = 25.4 / 72
SVG = "{http://www.w3.org/2000/svg}"


class FinalDrawings(unittest.TestCase):
    def match(self, actual, expected, tolerance=0.001):
        self.assertAlmostEqual(actual, expected, delta=tolerance)

    def check_centres(self, centres, keys):
        self.assertEqual(len(centres), 2 * len(keys))
        for index, key in enumerate(keys):
            left, right = centres[2 * index : 2 * index + 2]
            self.match(left[1], right[1])
            self.match(math.dist(left, right), SPEC["spacings"][key])

    def test_pdf_circles_and_scales(self):
        for name, dims in (("a4", (210, 297)), ("letter", (215.9, 279.4))):
            with self.subTest(name=name), fitz.open(HERE / f"build-{name}.pdf") as doc:
                self.assertEqual(len(doc), 6)
                self.match(doc[0].rect.width * MM_PER_PT, dims[0])
                self.match(doc[0].rect.height * MM_PER_PT, dims[1])
                for page, diameter, keys in (
                    (doc[2], 3.025, ["ground"]),
                    (doc[3], 3.15, ["crank", "coupler", "rocker"]),
                ):
                    boxes = [p["rect"] * MM_PER_PT for p in page.get_drawings()]
                    circles = [
                        b
                        for b in boxes
                        if abs(b.width - diameter) < 0.001 and abs(b.height - diameter) < 0.001
                    ]
                    centres = [((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2) for b in circles]
                    self.check_centres(centres, keys)
                    self.assertTrue(any(abs(b.width - 50) < 0.001 and b.height < 0.001 for b in boxes))
                    self.assertTrue(any(abs(b.height - 50) < 0.001 and b.width < 0.001 for b in boxes))
                    if keys == ["ground"]:
                        [outside] = [
                            b for b in boxes if abs(b.width - 174) < 0.001 and abs(b.height - 207) < 0.001
                        ]
                        self.match(centres[0][0] - outside.x0, 62)
                        self.match(outside.y1 - centres[0][1], 62)

    def test_actual_size_inspection_gauge(self):
        gauge_file = HERE.parent / "inspection" / "angle-gauge.pdf"
        with fitz.open(gauge_file) as document:
            self.assertEqual(len(document), 1)
            page = document[0]
            self.match(page.rect.width * MM_PER_PT, 195)
            self.match(page.rect.height * MM_PER_PT, 250)
            boxes = [item["rect"] * MM_PER_PT for item in page.get_drawings()]
            for width, height in ((50, 0), (0, 50), (174, 207)):
                self.assertTrue(
                    any(abs(box.width - width) < 0.001 and abs(box.height - height) < 0.001 for box in boxes)
                )
            pivots = [
                (box.x0 + box.width / 2, box.y0 + box.height / 2)
                for box in boxes
                if abs(box.width - 16) < 0.001 and abs(box.height - 16) < 0.001
            ]
            self.check_centres(pivots, ["ground"])
            self.assertIn("backwards along A-O", page.get_text())

    def test_native_svg_and_dxf(self):
        for name, width, height, keys in (
            ("base", 174, 207, ["ground"]),
            ("links", 168, 104, ["crank", "coupler", "rocker"]),
        ):
            with self.subTest(name=name):
                doc = ezdxf.readfile(HERE / f"cut-{name}.dxf")
                self.assertEqual(doc.units, ezdxf.units.MM)
                model = doc.modelspace()
                dxf_bores = list(model.query('CIRCLE[layer=="BORE"]'))
                xy = [(c.dxf.center.x, c.dxf.center.y) for c in dxf_bores]
                poly = list(model.query('LWPOLYLINE[layer=="PROFILE"]'))
                self.assertEqual(len(poly), len(keys))
                self.assertTrue(all(path.closed for path in poly))
                self.check_centres(xy, keys)
                root = ElementTree.parse(HERE / f"cut-{name}.svg").getroot()
                self.match(float(root.get("width").removesuffix("mm")), width + 2)
                self.match(float(root.get("height").removesuffix("mm")), height + 2)
                svg_bores = root.findall(SVG + "circle")
                self.assertEqual(len(svg_bores), len(dxf_bores))
                for dxf, svg in zip(dxf_bores, svg_bores):
                    self.match(dxf.dxf.radius, float(svg.get("r")))
                    self.match(dxf.dxf.center.x, float(svg.get("cx")) - 1)
                    self.match(dxf.dxf.center.y, height + 1 - float(svg.get("cy")))
                if name == "base":
                    self.assertEqual(xy, [(62, 62), (155.8, 62)])
                else:
                    for outline, key in zip(poly, keys):
                        pts = list(outline.get_points("xyb"))
                        self.assertEqual([p[2] for p in pts], [0, 1, 0, 1])
                        expected_length = SPEC["spacings"][key]
                        if key == "rocker":
                            expected_length += SPEC["links"]["pointer_extension"]
                        self.match(pts[1][0] - pts[0][0], expected_length)
                        self.match(pts[2][1] - pts[1][1], SPEC["links"]["width"])


if __name__ == "__main__":
    unittest.main()
