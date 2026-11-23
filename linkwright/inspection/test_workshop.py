"""Exercise the released 3D viewer, controls, mobile sizing and physical-record round trip offline."""
import argparse
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]

def fixture():
    values=dict(buildId='synthetic-browser-test',inspectionDate='2026-11-23',instruments='test gauges',
                printError=.04,lengthError=.02,gapError=.02,angleError=.5,
                ground=93.8,crank=46.6,coupler=84.1,rocker=111.3,goBlade=.55,limitingPair='test shaft / deck',
                limitingPhase=57,contact='none',cyclesCW=10,cyclesCCW=10,
                fit='yes',layers='yes',operation='yes',retention='yes')
    for sheet in ('Base','Links','Gauge'):
        for axis in 'XY':values['scale'+sheet+axis]=50
    for joint in 'ODAB':values['endplay'+joint]=.25
    for i in range(3):
        for key,v in [('rockerLow',101.33),('rockerHigh',161.32),('crankAtLow',56.61),('crankAtHigh',288.08)]:
            values[key+str(i)]=v
    return values

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chromium', required=True)
    parser.add_argument('--screenshots',type=Path)
    args=parser.parse_args()
    with sync_playwright() as pw:
        browser=pw.chromium.launch(executable_path=args.chromium,
                 args=['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
        page=browser.new_page(offline=True,viewport={'width':1080,'height':900})
        errors=[]
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto((ROOT/'index.html').as_uri())
        page.locator('a[href="assembly-view.html"]').first.click()
        assert page.locator('#model').is_visible() and not page.locator('#fallback').is_visible()
        assert page.evaluate("document.querySelector('canvas').getContext('webgl').getError()") == 0
        p0=page.locator('#pins b').nth(2).get_attribute('style')
        for id,val in [('theta','120'),('spread','75')]:
            page.locator('#'+id).fill(val)
        assert page.locator('#thetaValue').inner_text()=='120°'
        assert '75%' in page.locator('#spreadValue').inner_text()
        assert p0!=page.locator('#pins b').nth(2).get_attribute('style')
        page.locator('button[data-joint="A"]').click()
        assert page.locator('button[data-joint="A"]').get_attribute('aria-pressed')=='true'
        assert '10 mm spacer' in page.locator('#jointInfo').inner_text()
        page.locator('button[data-joint="all"]').click()
        page.locator('#model').focus();page.keyboard.press('ArrowRight')
        page.locator('#top').click();page.locator('#home').click()
        page.locator('#play').click();page.wait_for_function("Number(document.querySelector('#theta').value)>0")
        page.locator('#play').click();page.locator('#home').click()
        if args.screenshots:
            args.screenshots.mkdir(parents=True,exist_ok=True)
            page.locator('main').screenshot(path=str(args.screenshots/'assembly.png'))
        page.locator('a[href="commissioning.html"]').click()
        page.locator('button[type=submit]').click()
        assert page.locator('#result').get_attribute('data-state')=='incomplete'
        for name,val in fixture().items():
            control=page.locator('[name="'+name+'"]')
            if control.evaluate("e=>e.tagName")== 'SELECT': control.select_option(str(val))
            else: control.fill(str(val))
        page.locator('button[type=submit]').click()
        assert page.locator('#result').get_attribute('data-state')=='pass'
        with page.expect_download() as dl:page.locator('#save').click()
        result=json.loads(Path(dl.value.path()).read_text())
        assert result['format']=='linkwright.physical-inspection/2' and result['assessment']['state']=='pass'
        page.locator('[name=goBlade]').fill('.5')
        assert page.locator('#result').get_attribute('data-state')=='review'
        page.locator('[name=contact]').select_option('observed')
        assert page.locator('#result').get_attribute('data-state')=='fail'
        # Reload prior evidence through the actual upload widget and recompute its classification.
        page.once('dialog',lambda dialog:dialog.accept())
        page.locator('#open').set_input_files({'name':'test.json','mimeType':'application/json','buffer':json.dumps(result).encode()})
        page.wait_for_function("document.querySelector('#result').dataset.state==='pass'")
        for width in (390,320):
            page.set_viewport_size({'width':width,'height':860})
            for path in ('commissioning.html','assembly-view.html','kit/assembly-guide.html','custom-guide.html','verification/engineering.html'):
                page.goto((ROOT/path).as_uri())
                assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth'),(path,width)
        assert not errors,errors
        browser.close()
        print('Passed live WebGL controls, offline sign-off JSON, and all five documents at 320 and 390 px.')
if __name__=='__main__':main()
