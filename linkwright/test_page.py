"""Check the actual offline UI. Use --chromium to select a local Chromium executable."""
import argparse
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromium")
    parser.add_argument("--receipt", type=Path)
    args = parser.parse_args()
    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=args.chromium, args=["--no-sandbox"])
        page = browser.new_page(viewport={"width": 1280, "height": 1000}, offline=True)
        release_design = json.loads((HERE / "kit" / "design.json").read_text(encoding="utf-8"))
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto((HERE / "index.html").as_uri())
        animated = page.evaluate("""() => ({
            ground:kitModel.B0.x,crank:kitModel.a,coupler:kitModel.b,rocker:kitModel.c,
            extension:kitModel.pointerExtension, width:Number(document.getElementById('kitRocker').getAttribute('stroke-width'))
        })""")
        for name in ["ground","crank","coupler","rocker"]:
            assert animated[name] == release_design["spacings"][name], name
        assert animated["extension"] == release_design["links"]["pointer_extension"]
        assert animated["width"] == release_design["links"]["width"]
        receipt = page.evaluate((HERE / "browser_scenario.js").read_text(encoding="utf-8"))
        widths = [1280, 390, 320]
        for width in widths:
            page.set_viewport_size({"width": width, "height": 900})
            assert page.evaluate("document.documentElement.scrollWidth <= innerWidth"), width
        for button, suffix in [("saveSvg", ".svg"), ("saveJson", ".json")]:
            with page.expect_download() as event:
                page.locator("#" + button).click()
            assert event.value.suggested_filename.endswith(suffix)
            downloaded = Path(event.value.path())
            if suffix == ".svg":
                from xml.etree import ElementTree
                svg = ElementTree.parse(downloaded).getroot()
                assert svg.get("width").endswith("mm") and svg.get("height").endswith("mm")
                assert "50 mm vertical" in "".join(svg.itertext())
            else:
                exported = json.loads(downloaded.read_text(encoding="utf-8"))
                assert exported["version"] == 3 and exported["units"] == "mm"
                assert exported["pathGap"]["dimensionUnits"] == "mm at current cut scale"
        assert not errors, errors
        receipt["viewportWidths"] = widths
        receipt["referenceAnimationSpec"] = release_design["revision"]
        receipt["downloaded"] = ["parsed 1:1 SVG", "parsed v3 JSON"]
        browser.close()
    if args.receipt:
        args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
